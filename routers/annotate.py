from fastapi import APIRouter, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from typing import Optional, Literal
from pathlib import Path
import asyncio
import pandas as pd
import duckdb
import re
import shlex

router = APIRouter()

# ====== Config ======
PLANTDB_PARQUET = Path("data/moleculesdb.parquet").resolve()
# Set absolute path on Windows if needed, e.g.:
# RSCRIPT_BIN = r"C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
RSCRIPT_BIN = "Rscript"

H_ATOM = 1.007825  # as requested

# ====== Minimal ASCII-only R script for formula prediction ======
R_MIN_SCRIPT = r"""
suppressPackageStartupMessages({
  load_or_install <- function(pkg){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "http://cran.us.r-project.org")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  pkgs <- c("Rdisop")
  invisible(lapply(pkgs, load_or_install))
})

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(name, default=NULL) {
  hit <- grep(paste0("^--", name, "="), args, value=TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit)
}

in_csv  <- get_arg("in")
out_csv <- get_arg("out")
if (is.null(in_csv) || !file.exists(in_csv)) stop(paste("Input CSV not found:", in_csv))
if (is.null(out_csv)) stop("Missing --out")

inp <- tryCatch(read.csv(in_csv, stringsAsFactors=FALSE, check.names=FALSE),
                error=function(e) stop(paste("Cannot read input CSV:", e$message)))
if (!all(c("id","mass") %in% names(inp))) stop("Input CSV must contain 'id' and 'mass'")

extract_tab <- function(cand) {
  if (is.null(cand)) return(NULL)
  if (isS4(cand)) {
    sn <- try(slotNames(cand), silent=TRUE)
    if (!inherits(sn, "try-error") && "results" %in% sn) return(cand@results)
  }
  if (is.list(cand) && !is.null(cand$results)) return(cand$results)
  if (is.data.frame(cand)) return(cand)
  return(NULL)
}

predict_formula_one <- function(mass) {
  if (is.na(mass)) return(NA_character_)
  cand <- tryCatch(Rdisop::decomposeMass(mass, mzabs=0.005),
                   error=function(e) NULL)
  tab <- extract_tab(cand)
  if (is.null(tab) || nrow(tab)==0) return(NA_character_)
  if (!("exactmass" %in% names(tab)) && "mass" %in% names(tab)) tab$exactmass <- tab$mass
  if ("dbd" %in% names(tab)) {
    tab <- tab[!is.na(tab$dbd) & tab$dbd > 0, , drop=FALSE]
  } else if ("RDBE" %in% names(tab)) {
    tab <- tab[!is.na(tab$RDBE) & tab$RDBE > 0, , drop=FALSE]
  }
  if (nrow(tab)==0) return(NA_character_)
  tab$err <- abs(tab$exactmass - mass)
  tab <- tab[order(tab$err), , drop=FALSE]
  tab$formula[1]
}

res <- data.frame(id = inp$id, Formula = NA_character_, stringsAsFactors=FALSE)
n <- nrow(inp)
for (i in seq_len(n)) {
  res$Formula[i] <- predict_formula_one(as.numeric(inp$mass[i]))
  if (i %% 1000 == 0) cat(".. Rdisop processed", i, "rows\n")
}
write.csv(res, out_csv, row.names=FALSE)
"""

# ====== Request model ======
class AnnotParams(BaseModel):
    directory: str = Field(min_length=1, description="Folder containing Features_Matrix.csv")
    table: str = Field(min_length=1, description="Base name (no .csv), e.g., Features_Matrix")
    polarity: Literal["positive", "negative"]
    analyte_csv: Optional[str] = Field(
        None, description="Optional analyte CSV with columns: MolecularFormula, ExactMass, CompoundName"
    )

# ====== Helpers ======
def _safe(p: str) -> Path:
    return Path(p).expanduser().resolve()

def _ensure_exists(path: Path, what: str):
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"{what} not found: {path}")

def is_molecular_row(tag) -> bool:
    """
    Keep rows where isotopes is:
      - empty/NA, OR
      - a molecular-ion pattern "[number][M]" optionally followed by +/- (e.g., "[10][M]+", "[1][M]-", "[2][M]")
    """
    if tag is None:
        return True
    s = str(tag).strip()
    if s == "":
        return True
    return re.match(r"^\[\d+\]\[M\][+-]?$", s) is not None

# ====== Core streaming worker ======
def _make_stream(params: AnnotParams):
    async def stream():
        directory = _safe(params.directory)
        _ensure_exists(directory, "Directory")

        features_csv = (directory / f"{params.table}.csv").resolve()
        _ensure_exists(features_csv, "Features CSV")
        _ensure_exists(PLANTDB_PARQUET, "DB Parquet")

        analyte_path: Optional[Path] = None
        if params.analyte_csv:
            analyte_path = _safe(params.analyte_csv)
            _ensure_exists(analyte_path, "Analyte CSV")

        ft_ready_csv    = (directory / "FT_Ready.csv").resolve()
        r_in_csv        = (directory / f"{params.table}__r_in.csv").resolve()
        r_out_csv       = (directory / f"{params.table}__r_out.csv").resolve()
        r_src           = (directory / f"{params.table}__r.R").resolve()
        fm_ready_csv    = (directory / "Features_Matrix_Ready.csv").resolve()
        fm_annot_csv    = (directory / "Features_Matrix_Ready_Annotated.csv").resolve()

        try:
            # ---------- Load & filter ----------
            yield "=== Step 1: Load & filter to molecular ions ===\n"
            yield f"Loading: {features_csv}\n"
            df = pd.read_csv(features_csv)

            # Ensure mz exists
            if "mz" not in df.columns:
                yield "❌ Features CSV must contain 'mz' column.\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            # RT → minutes column
            if "RT_min" not in df.columns:
                if "rt" in df.columns:
                    df["RT_min"] = pd.to_numeric(df["rt"], errors="coerce") / 60.0
                elif "RT" in df.columns:
                    df["RT_min"] = pd.to_numeric(df["RT"], errors="coerce")
                else:
                    df["RT_min"] = pd.NA

            # Filter isotopes to empty or [number][M](+/-)
            if "isotopes" in df.columns:
                keep = df["isotopes"].apply(is_molecular_row)
                kept, dropped = int(keep.sum()), int((~keep).sum())
                yield f"Isotope filter: kept {kept}; dropped {dropped}\n"
                df = df[keep].copy()
            else:
                yield "No 'isotopes' column; keeping all rows as molecular ions.\n"

            if len(df) == 0:
                yield "No rows remain after filtering. Writing empty outputs.\n"
                pd.DataFrame(columns=["RT_min","mz","ExactMass"]).to_csv(ft_ready_csv, index=False)
                pd.DataFrame(columns=["RT_min","mz","ExactMass","Formula"]).to_csv(fm_ready_csv, index=False)
                pd.DataFrame(columns=["RT_min","mz","ExactMass","Formula","Candidates","AnalyteList"]).to_csv(fm_annot_csv, index=False)
                yield f"✅ Wrote: {ft_ready_csv}\n"
                yield f"✅ Wrote: {fm_ready_csv}\n"
                yield f"✅ Wrote: {fm_annot_csv}\n"
                yield "__ANNOT_DONE__ EXIT_CODE=0\n"
                return

            # ---------- ExactMass by polarity (using H_ATOM as requested) ----------
            yield "=== Step 2: Compute ExactMass from mz and polarity ===\n"
            if params.polarity == "positive":
                df["ExactMass"] = pd.to_numeric(df["mz"], errors="coerce") - H_ATOM
            else:
                df["ExactMass"] = pd.to_numeric(df["mz"], errors="coerce") + H_ATOM

            before = len(df)
            df = df[df["ExactMass"].notna()].copy()
            after = len(df)
            yield f"ExactMass non-null: {after} / {before}\n"
            try:
                prev = df[["mz","RT_min","ExactMass"]].head(5).to_string(index=False)
                yield "Preview (mz, RT_min, ExactMass):\n" + prev + "\n"
            except Exception:
                pass

            # Save FT_Ready.csv
            df_ready = df.copy()
            df_ready[["RT_min","mz","ExactMass"]].to_csv(ft_ready_csv, index=False)
            yield f"✅ Wrote: {ft_ready_csv}\n"

            # ---------- Rdisop formulas ----------
            yield "=== Step 3: Rdisop formula prediction ===\n"
            df_ready = df_ready.reset_index(drop=True)
            df_ready["_rid"] = df_ready.index + 1
            r_in = df_ready[["_rid","ExactMass"]].rename(columns={"_rid":"id","ExactMass":"mass"})
            r_in.to_csv(r_in_csv, index=False)
            r_src.write_text(R_MIN_SCRIPT, encoding="utf-8")
            yield f"Wrote R input: {r_in_csv}\n"
            yield f"Wrote R script: {r_src}\n"

            cmd = [RSCRIPT_BIN, "--encoding=utf8", str(r_src), f"--in={r_in_csv}", f"--out={r_out_csv}"]
            yield "Command: " + " ".join(shlex.quote(p) for p in cmd) + "\n"

            try:
                proc = await asyncio.create_subprocess_exec(
                    *cmd,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.STDOUT,
                )
            except FileNotFoundError:
                yield f"❌ Rscript not found: {RSCRIPT_BIN}\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            async for chunk in proc.stdout:
                yield chunk.decode(errors="ignore")
            r_code = await proc.wait()
            if r_code != 0:
                yield f"❌ Rdisop failed (exit {r_code}).\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            if not r_out_csv.exists():
                yield f"❌ Missing R output: {r_out_csv}\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            # Merge Formula -> Features_Matrix_Ready.csv
            yield "=== Step 4: Merge formulas and write Features_Matrix_Ready.csv ===\n"
            rform = pd.read_csv(r_out_csv)  # columns: id, Formula
            if not {"id","Formula"}.issubset(rform.columns):
                yield "❌ R output must contain 'id' and 'Formula'.\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return
            rform = rform.rename(columns={"id":"_rid"})
            df_ready = df_ready.merge(rform, on="_rid", how="left")

            # Keep RT_min, mz, ExactMass, Formula + intensity columns
            intensity_cols = [c for c in df.columns if re.search(r"(?i)(intensity|into|maxo|area)", c)]
            cols_out = ["RT_min","mz","ExactMass","Formula"] + intensity_cols
            df_ready[cols_out].to_csv(fm_ready_csv, index=False)
            yield f"✅ Wrote: {fm_ready_csv}\n"

            # ---------- Join with DB + optional analyte list ----------
            yield "=== Step 5: Join with DB and Analyte list ===\n"
            con = duckdb.connect()
            con.register("fm", df_ready[cols_out])  # RT_min, mz, ExactMass, Formula, intensities...

            db_path_sql = str(PLANTDB_PARQUET).replace("'", "''")
            con.execute(f"""
                CREATE OR REPLACE VIEW mols AS
                SELECT * FROM read_parquet('{db_path_sql}')
            """)

            con.execute("""
                CREATE OR REPLACE VIEW dbagg AS
                SELECT MolecularFormula AS formula,
                       string_agg(CompoundName, ' / ') AS Candidates
                FROM mols
                WHERE MolecularFormula IS NOT NULL
                GROUP BY MolecularFormula
            """)

            analyte_join_sql = ""
            analyte_sel_sql  = ""
            if analyte_path:
                analyte_path_sql = str(analyte_path).replace("'", "''")
                con.execute(f"""
                    CREATE OR REPLACE VIEW analyte AS
                    SELECT * FROM read_csv_auto('{analyte_path_sql}', HEADER=TRUE)
                """)
                con.execute("""
                    CREATE OR REPLACE VIEW aagg AS
                    SELECT MolecularFormula AS formula,
                           string_agg(CompoundName, ' / ') AS AnalyteList
                    FROM analyte
                    WHERE MolecularFormula IS NOT NULL
                    GROUP BY MolecularFormula
                """)
                analyte_join_sql = "LEFT JOIN aagg A ON A.formula = fm.Formula"
                analyte_sel_sql  = ", A.AnalyteList"

            sql = f"""
            SELECT
              fm.RT_min,
              fm.mz,
              fm.ExactMass,
              fm.Formula,
              {"fm.* EXCLUDE (RT_min, mz, ExactMass, Formula)," if len(intensity_cols) > 0 else ""}
              D.Candidates
              {analyte_sel_sql}
            FROM fm
            LEFT JOIN dbagg D ON D.formula = fm.Formula
            {analyte_join_sql}
            """
            out_df = con.execute(sql).fetch_df()
            out_df.to_csv(fm_annot_csv, index=False)

            yield f"✅ Wrote: {fm_annot_csv}\n"
            yield "__ANNOT_DONE__ EXIT_CODE=0\n"
        except Exception as e:
            yield f"❌ Error: {e}\n"
            yield "__ANNOT_DONE__ EXIT_CODE=1\n"
    return stream

# ====== Endpoints ======
@router.post("/annotation/run")
async def annotation_run(params: AnnotParams):
    return StreamingResponse(_make_stream(params)(), media_type="text/plain; charset=utf-8")

# Back-compat alias if your page still posts to /annotate/run
@router.post("/annotate/run")
async def annotate_run_alias(params: AnnotParams):
    return StreamingResponse(_make_stream(params)(), media_type="text/plain; charset=utf-8")
