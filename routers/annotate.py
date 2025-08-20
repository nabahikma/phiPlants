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
import subprocess

router = APIRouter()

# ====== Config ======
PLANTDB_PARQUET = Path("data/moleculesdb.parquet").resolve()
# Set absolute path on Windows if Rscript isn't in PATH, e.g.:
# RSCRIPT_BIN = r"C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
RSCRIPT_BIN = "Rscript"

# Adduct mass shifts (Da) and charge
ADDUCT_INFO = {
    "[M+H]+":   (1.007276, +1),
    "[M+Na]+":  (22.989218, +1),
    "[M+K]+":   (38.963158, +1),
    "[M+NH4]+": (18.033823, +1),
    "[M+2H]2+": (2 * 1.007276, +2),

    "[M-H]-":    (-1.007276, -1),
    "[M+Cl]-":   (34.969402, -1),
    "[M+FA-H]-": (46.00548 - 1.007276, -1),
    "[M+HCOO]-": (46.00548 - 1.007276, -1),
}

# ====== Minimal R script (ASCII only) for formula prediction ======
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

cat("Rdisop: reading masses from", in_csv, "\n")
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
  cand <- tryCatch(Rdisop::decomposeMass(mass, mzabs=0.005, elements=c(C=0,H=0,N=0,O=0,P=0,S=0)),
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

res <- data.frame(id = inp$id, predicted_formula = NA_character_, stringsAsFactors=FALSE)
n <- nrow(inp)
for (i in seq_len(n)) {
  res$predicted_formula[i] <- predict_formula_one(as.numeric(inp$mass[i]))
  if (i %% 1000 == 0) cat(".. Rdisop processed", i, "rows\n")
}
write.csv(res, out_csv, row.names=FALSE)
cat("Rdisop: wrote formulas to", out_csv, "\n")
"""

# ====== Request model ======
class AnnotParams(BaseModel):
    directory: str = Field(min_length=1, description="Folder containing Features CSV")
    table: str = Field(min_length=1, description="Base name of features CSV (no .csv). E.g., Features_Matrix")
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

def choose_adduct(raw_adduct: Optional[str], polarity: str) -> str:
    if raw_adduct and isinstance(raw_adduct, str) and raw_adduct.strip():
        first = re.split(r"[;,]", raw_adduct)[0].strip()
        if first:
            return first
    return "[M+H]+" if polarity == "positive" else "[M-H]-"

def neutral_mass_from_mz(mz_val, adduct: str) -> Optional[float]:
    info = ADDUCT_INFO.get(adduct)
    try:
        mz = float(mz_val)
    except Exception:
        return None
    if not info:
        return None
    shift, z = info
    return (abs(z) * mz - shift) / abs(z)

def is_molecular_isotope_tag(tag: Optional[str]) -> bool:
    """Keep empty/NA or tags like [M]+ / [M]- / [10][M]+ ; drop others."""
    if tag is None:
        return True
    s = str(tag).strip()
    if s == "":
        return True
    return re.match(r"^\[\d*\]\[M\][+-]?$", s) is not None

# ====== Endpoint (streams logs) ======
@router.post("/annotate/run")
async def annotate_run(params: AnnotParams):
    directory = _safe(params.directory)
    _ensure_exists(directory, "Directory")

    features_csv = (directory / f"{params.table}.csv").resolve()
    _ensure_exists(features_csv, "Features CSV")

    _ensure_exists(PLANTDB_PARQUET, "DB Parquet")

    analyte_path: Optional[Path] = None
    if params.analyte_csv:
        analyte_path = _safe(params.analyte_csv)
        _ensure_exists(analyte_path, "Analyte CSV")

    in_r_csv  = (directory / f"{params.table}__rdisop_in.csv").resolve()
    out_r_csv = (directory / f"{params.table}__rdisop_out.csv").resolve()
    r_src     = (directory / f"{params.table}__rdisop.R").resolve()
    final_out = (directory / "Features_Matrix_Annotated.csv").resolve()

    async def stream():
        try:
            # ---------- Python preprocessing ----------
            yield "=== Annotation: Python preprocessing ===\n"
            yield f"Loading: {features_csv}\n"
            df = pd.read_csv(features_csv)

            if "mz" not in df.columns:
                yield "❌ Features CSV must contain 'mz' column.\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            # Normalize RT → minutes
            if "RT_min" not in df.columns:
                if "rt" in df.columns:
                    df["RT_min"] = pd.to_numeric(df["rt"], errors="coerce") / 60.0
                elif "RT" in df.columns:
                    df["RT_min"] = pd.to_numeric(df["RT"], errors="coerce")
                else:
                    df["RT_min"] = pd.NA

            # Isotope filter
            if "isotopes" in df.columns:
                keep_mask = df["isotopes"].apply(is_molecular_isotope_tag)
                kept = int(keep_mask.sum())
                dropped = int((~keep_mask).sum())
                yield f"Isotope filter: kept {kept}; dropped {dropped}\n"
                df = df[keep_mask].copy()
            else:
                yield "No 'isotopes' column; treating all rows as molecular ions.\n"

            if len(df) == 0:
                yield "No rows remain after isotope filter; writing empty annotated CSV.\n"
                pd.DataFrame(columns=["RT_min","mz","adduct","ExactMass","PredictedFormula","isotopes"]).to_csv(final_out, index=False)
                yield f"✅ Wrote: {final_out}\n"
                yield "__ANNOT_DONE__ EXIT_CODE=0\n"
                return

            # Choose adduct + compute neutral mass
            adduct_exists = "adduct" in df.columns
            df["adduct_used"] = [
                choose_adduct(df.at[i, "adduct"] if adduct_exists else None, params.polarity)
                for i in df.index
            ]
            df["ExactMass"] = [neutral_mass_from_mz(df.at[i, "mz"], df.at[i, "adduct_used"]) for i in df.index]

            # Assign row ids for round-trip to R
            df["_rid"] = range(1, len(df) + 1)
            masses_for_r = df[["_rid", "ExactMass"]].rename(columns={"_rid": "id", "ExactMass": "mass"})

            # ---------- Write R inputs ----------
            masses_for_r.to_csv(in_r_csv, index=False)
            r_src.write_text(R_MIN_SCRIPT, encoding="utf-8")
            yield f"Wrote R input: {in_r_csv}\n"
            yield f"Wrote R script: {r_src}\n"

            # ---------- Run Rdisop ----------
            yield "=== Rdisop: predicting formulas ===\n"
            cmd = [RSCRIPT_BIN, "--encoding=utf8", str(r_src), f"--in={in_r_csv}", f"--out={out_r_csv}"]
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

            # stream R output live
            async for chunk in proc.stdout:
                yield chunk.decode(errors="ignore")
            r_code = await proc.wait()
            if r_code != 0:
                yield f"❌ Rdisop failed (exit {r_code}).\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            if not out_r_csv.exists():
                yield f"❌ Missing R output: {out_r_csv}\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            # ---------- Merge formulas back ----------
            yield "=== Python postprocessing: merging & joining ===\n"
            rform = pd.read_csv(out_r_csv)
            if not {"id","predicted_formula"}.issubset(rform.columns):
                yield "❌ R output must contain 'id' and 'predicted_formula'.\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return
            rform = rform.rename(columns={"id": "_rid", "predicted_formula": "PredictedFormula"})
            df = df.merge(rform, on="_rid", how="left")

            # Column order: key + intensities + rest
            intensity_cols = [c for c in df.columns if re.search(r"(?i)(intensity|into|maxo|area)", c)]
            key_cols = ["RT_min","mz","adduct_used","ExactMass","PredictedFormula"]
            maybe_iso = ["isotopes"] if "isotopes" in df.columns else []
            rest = [c for c in df.columns if c not in set(key_cols + maybe_iso + intensity_cols + ["_rid"])]
            final_cols = key_cols + maybe_iso + intensity_cols + rest
            df_final = df[final_cols].rename(columns={"adduct_used":"adduct"})

            # ---------- DuckDB joins ----------
            con = duckdb.connect()
            con.execute("CREATE OR REPLACE VIEW feats AS SELECT * FROM df_final")
            con.execute("CREATE OR REPLACE VIEW mols  AS SELECT * FROM read_parquet(?)", [str(PLANTDB_PARQUET)])

            con.execute("""
                CREATE OR REPLACE VIEW dbagg AS
                SELECT MolecularFormula AS formula, string_agg(CompoundName, ' / ') AS Candidates
                FROM mols
                WHERE MolecularFormula IS NOT NULL
                GROUP BY MolecularFormula
            """)

            analyte_join_sql = ""
            analyte_sel_sql  = ""
            if analyte_path:
                analyte_name = analyte_path.stem
                con.execute("CREATE OR REPLACE VIEW analyte AS SELECT * FROM read_csv_auto(?, HEADER=TRUE)", [str(analyte_path)])
                con.execute(f"""
                    CREATE OR REPLACE VIEW aagg AS
                    SELECT MolecularFormula AS formula, string_agg(CompoundName, ' / ') AS "{analyte_name}"
                    FROM analyte
                    WHERE MolecularFormula IS NOT NULL
                    GROUP BY MolecularFormula
                """)
                analyte_join_sql = "LEFT JOIN aagg A ON A.formula = feats.PredictedFormula"
                analyte_sel_sql  = f', A."{analyte_name}"'

            sql = f"""
            SELECT
              feats.RT_min,
              feats.mz,
              feats.adduct,
              feats.ExactMass,
              feats.PredictedFormula,
              {"feats.isotopes," if "isotopes" in df_final.columns else ""}
              feats.* EXCLUDE (RT_min, mz, adduct, ExactMass, PredictedFormula{", isotopes" if "isotopes" in df_final.columns else ""}),
              D.Candidates
              {analyte_sel_sql}
            FROM feats
            LEFT JOIN dbagg D ON D.formula = feats.PredictedFormula
            {analyte_join_sql}
            """
            out_df = con.execute(sql).fetch_df()
            out_df.to_csv(final_out, index=False)

            yield f"✅ Wrote: {final_out}\n"
            yield "__ANNOT_DONE__ EXIT_CODE=0\n"
        except Exception as e:
            yield f"❌ Error: {e}\n"
            yield "__ANNOT_DONE__ EXIT_CODE=1\n"

    return StreamingResponse(stream(), media_type="text/plain; charset=utf-8")
