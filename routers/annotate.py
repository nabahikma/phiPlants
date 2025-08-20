# routers/annotation.py
from fastapi import APIRouter, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from typing import Optional, Literal, Dict, Tuple
from pathlib import Path
import asyncio
import pandas as pd
import duckdb
import re
import shlex

router = APIRouter()

# ========= Config =========
PLANTDB_PARQUET = Path("data/moleculesdb.parquet").resolve()
# If Rscript isn't in PATH on Windows, set the full path, e.g.:
# RSCRIPT_BIN = r"C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
RSCRIPT_BIN = "Rscript"

# Monoisotopic constants (Da)
PROTON = 1.007276
NA_ION = 22.989218
K_ION  = 38.963158
NH4_ION = 18.033823
CL_ION = 34.969402
BR_ION = 78.918885
H2O = 18.010565
ACN = 41.026549      # CH3CN
MEOH = 32.026215     # CH3OH
FORMIC_ACID = 46.005480
ACETIC_ACID = 60.021129
FORMATE = 45.017444  # HCOO- anion mass? Often use FA-H net ~59? We will use FA pathway via acid - H+
ACETATE = 59.013851  # CH3COO-

# ========= Expanded adduct list (shift, z, n) =========
# shift is net adduct mass in Da (what you add/subtract to n*M before dividing by |z|)
# z is the integer charge (sign matters), n is stoichiometry (e.g., 2 for [2M+H]+)
ADDUCTS: Dict[str, Tuple[float, int, int]] = {
    # --- Positive, monomer (n=1) ---
    "[M+H]+":        (PROTON, +1, 1),
    "[M+Na]+":       (NA_ION, +1, 1),
    "[M+K]+":        (K_ION, +1, 1),
    "[M+NH4]+":      (NH4_ION, +1, 1),
    "[M+2H]2+":      (2*PROTON, +2, 1),
    "[M+3H]3+":      (3*PROTON, +3, 1),

    "[M+ACN+H]+":    (ACN + PROTON, +1, 1),
    "[M+CH3OH+H]+":  (MEOH + PROTON, +1, 1),
    "[M+H-H2O]+":    (PROTON - H2O, +1, 1),
    "[M+Na+ACN]+":   (NA_ION + ACN, +1, 1),

    # --- Positive, dimers (n=2) ---
    "[2M+H]+":       (PROTON, +1, 2),
    "[2M+Na]+":      (NA_ION, +1, 2),
    "[2M+K]+":       (K_ION, +1, 2),
    "[2M+NH4]+":     (NH4_ION, +1, 2),
    "[2M+ACN+H]+":   (ACN + PROTON, +1, 2),

    # --- Positive, trimers (n=3) ---
    "[3M+H]+":       (PROTON, +1, 3),
    "[3M+Na]+":      (NA_ION, +1, 3),

    # --- Negative, monomer (n=1) ---
    "[M-H]-":        (-PROTON, -1, 1),
    "[M-2H]2-":      (-2*PROTON, -2, 1),
    "[M+Cl]-":       (CL_ION, -1, 1),
    "[M+Br]-":       (BR_ION, -1, 1),

    # Formate / formic acid paths
    "[M+FA-H]-":     (FORMIC_ACID - PROTON, -1, 1),  # common notation
    "[M+HCOO]-":     (FORMATE, -1, 1),

    # Acetate / acetic acid paths
    "[M+Ac-H]-":     (ACETIC_ACID - PROTON, -1, 1),
    "[M+CH3COO]-":   (ACETATE, -1, 1),

    # --- Negative, dimers (n=2) ---
    "[2M-H]-":       (-PROTON, -1, 2),
    "[2M+Cl]-":      (CL_ION, -1, 2),
    "[2M+FA-H]-":    (FORMIC_ACID - PROTON, -1, 2),
    "[2M+HCOO]-":    (FORMATE, -1, 2),
    "[2M+Ac-H]-":    (ACETIC_ACID - PROTON, -1, 2),

    # --- Negative, trimers (n=3) ---
    "[3M-H]-":       (-PROTON, -1, 3),
}

# ========= Minimal ASCII-only R script (Rdisop only) =========
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
  # No 'elements=' constraint: let Rdisop search fully, then filter DBE>0
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

res <- data.frame(id = inp$id, predicted_formula = NA_character_, stringsAsFactors=FALSE)
n <- nrow(inp)
for (i in seq_len(n)) {
  res$predicted_formula[i] <- predict_formula_one(as.numeric(inp$mass[i]))
  if (i %% 1000 == 0) cat(".. Rdisop processed", i, "rows\n")
}
write.csv(res, out_csv, row.names=FALSE)
"""

# ========= Request model =========
class AnnotParams(BaseModel):
    directory: str = Field(min_length=1, description="Folder containing Features CSV")
    table: str = Field(min_length=1, description="Base name (no .csv), e.g., Features_Matrix")
    polarity: Literal["positive", "negative"]
    analyte_csv: Optional[str] = Field(
        None, description="Optional analyte CSV with columns: MolecularFormula, ExactMass, CompoundName"
    )

# ========= Helpers =========
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
    """
    General: m/z = (n*M + shift) / |z|  =>  M = (|z|*m/z - shift) / n
    """
    if adduct not in ADDUCTS:
        return None
    shift, z, n = ADDUCTS[adduct]
    try:
        mz = float(mz_val)
    except Exception:
        return None
    return ((abs(z) * mz) - shift) / max(n, 1)

def is_molecular_isotope_tag(tag: Optional[str]) -> bool:
    """Keep empty/NA or tags like [M]+ / [M]- / [10][M]+ ; drop others (e.g., [M+1]+)."""
    if tag is None:
        return True
    s = str(tag).strip()
    if s == "":
        return True
    return re.match(r"^\[\d*\]\[M\][+-]?$", s) is not None

# ========= Core streaming worker =========
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

        in_r_csv  = (directory / f"{params.table}__r_in.csv").resolve()
        out_r_csv = (directory / f"{params.table}__r_out.csv").resolve()
        r_src     = (directory / f"{params.table}__r.R").resolve()
        with_formula = (directory / f"{params.table}_with_Formula.csv").resolve()
        final_out = (directory / "Features_Matrix_Annotated.csv").resolve()

        try:
            # ---------- Python preprocessing ----------
            yield "=== Annotation: Python preprocessing ===\n"
            yield f"Loading: {features_csv}\n"
            df = pd.read_csv(features_csv)

            if "mz" not in df.columns:
                yield "❌ Features CSV must contain 'mz' column.\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            # RT → minutes
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
                kept, dropped = int(keep_mask.sum()), int((~keep_mask).sum())
                yield f"Isotope filter: kept {kept}; dropped {dropped}\n"
                df = df[keep_mask].copy()
            else:
                yield "No 'isotopes' column; treating all rows as molecular ions.\n"

            if len(df) == 0:
                yield "No rows remain after isotope filter.\n"
                # still write empty outputs for UX consistency
                pd.DataFrame(columns=["RT_min","mz","adduct","ExactMass","PredictedFormula","isotopes"]).to_csv(with_formula, index=False)
                pd.DataFrame(columns=["RT_min","mz","adduct","ExactMass","PredictedFormula","isotopes","Candidates"]).to_csv(final_out, index=False)
                yield f"✅ Wrote: {with_formula}\n"
                yield f"✅ Wrote: {final_out}\n"
                yield "__ANNOT_DONE__ EXIT_CODE=0\n"
                return

            # Choose adduct + calculate neutral ExactMass
            adduct_exists = "adduct" in df.columns
            df["adduct_used"] = [
                choose_adduct(df.at[i, "adduct"] if adduct_exists else None, params.polarity)
                for i in df.index
            ]
            df["ExactMass"] = [neutral_mass_from_mz(df.at[i, "mz"], df.at[i, "adduct_used"]) for i in df.index]

            # Drop NA ExactMass, log status and preview
            before = len(df)
            df = df[pd.to_numeric(df["ExactMass"], errors="coerce").notna()].copy()
            after = len(df)
            yield f"ExactMass non-null: {after} / {before}\n"
            try:
                prev = df[["mz","adduct_used","ExactMass"]].head(5).to_string(index=False)
                yield "Preview (mz, adduct_used, ExactMass):\n" + prev + "\n"
            except Exception:
                pass

            if after == 0:
                yield "No valid ExactMass values to submit to R. Aborting.\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            # Prepare compact CSV for Rdisop (id,mass)
            df["_rid"] = range(1, len(df) + 1)
            masses_for_r = df[["_rid", "ExactMass"]].rename(columns={"_rid": "id", "ExactMass": "mass"})
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

            # ---------- Merge back & save with formula ----------
            yield "=== Python: merging formulas & saving with_Formula CSV ===\n"
            rform = pd.read_csv(out_r_csv)
            if not {"id","predicted_formula"}.issubset(rform.columns):
                yield "❌ R output must contain 'id' and 'predicted_formula'.\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return
            rform = rform.rename(columns={"id": "_rid", "predicted_formula": "PredictedFormula"})
            df = df.merge(rform, on="_rid", how="left")

            cols_with = ["RT_min","mz","adduct_used","ExactMass","PredictedFormula"] + (["isotopes"] if "isotopes" in df.columns else [])
            df_with = df[cols_with + [c for c in df.columns if c not in cols_with + ["_rid"]]].rename(columns={"adduct_used":"adduct"})
            df_with.to_csv(with_formula, index=False)
            yield f"✅ Wrote: {with_formula}\n"

            # ---------- DuckDB joins (DB + optional analyte list) ----------
            yield "=== Python: joining with DB & analyte list ===\n"
            con = duckdb.connect()
            con.register("feats", df_with)  # register pandas DF explicitly

            # inline file paths safely (escape single quotes)
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
                analyte_name = analyte_path.stem
                analyte_path_sql = str(analyte_path).replace("'", "''")
                con.execute(f"""
                    CREATE OR REPLACE VIEW analyte AS
                    SELECT * FROM read_csv_auto('{analyte_path_sql}', HEADER=TRUE)
                """)
                con.execute(f"""
                    CREATE OR REPLACE VIEW aagg AS
                    SELECT MolecularFormula AS formula,
                           string_agg(CompoundName, ' / ') AS "{analyte_name}"
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
              {"feats.isotopes," if "isotopes" in df_with.columns else ""}
              feats.* EXCLUDE (RT_min, mz, adduct, ExactMass, PredictedFormula{", isotopes" if "isotopes" in df_with.columns else ""}),
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
    return stream

# ========= Endpoints =========
@router.post("/annotation/run")
async def annotation_run(params: AnnotParams):
    return StreamingResponse(_make_stream(params)(), media_type="text/plain; charset=utf-8")

# Back-compat alias if your page still posts to /annotate/run
@router.post("/annotate/run")
async def annotate_run_alias(params: AnnotParams):
    return StreamingResponse(_make_stream(params)(), media_type="text/plain; charset=utf-8")
