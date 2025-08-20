# routers/annotate.py
from fastapi import APIRouter, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from typing import Optional, Literal
from pathlib import Path
import asyncio
import tempfile
import json
import duckdb

router = APIRouter()

# Your bundled molecules DB (Parquet)
PLANTDB_PARQUET = Path("data/moleculesdb.parquet").resolve()

# Set to absolute path to Rscript.exe on Windows if needed, e.g.:
# RSCRIPT_BIN = r"C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
RSCRIPT_BIN = "Rscript"


class AnnotParams(BaseModel):
    directory: str = Field(min_length=1, description="Folder where Features_Matrix.csv lives")
    table: str = Field(min_length=1, description="Base name of features CSV (no .csv). E.g., Features_Matrix")
    polarity: Literal["positive", "negative"]
    # Default = DB. Optional = Analyte List *in addition* to DB.
    analyte_csv: Optional[str] = Field(
        None, description="Optional analyte CSV with columns: MolecularFormula, ExactMass, CompoundName"
    )


def _safe(p: str) -> Path:
    return Path(p).expanduser().resolve()


def _ensure_exists(path: Path, what: str):
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"{what} not found: {path}")


# ===== Embedded R stage (isotope filter + adduct + Rdisop formulas) =====
R_SCRIPT = r'''
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

features_csv <- get_arg("features")
polarity     <- get_arg("polarity", "positive")
out_csv      <- get_arg("out")

if (is.null(features_csv) || !file.exists(features_csv)) {
  stop(paste("Features CSV not found:", features_csv))
}
if (is.null(out_csv)) {
  stop("Missing --out argument")
}

cat("Rdisop annotation starting...\n")
cat("Features:", features_csv, "\nPolarity:", polarity, "\nOut:", out_csv, "\n")

# Coalesce helper (must be defined early)
`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else ifelse(is.na(a), b, a)

# Read features
df <- tryCatch({
  read.csv(features_csv, stringsAsFactors = FALSE, check.names = FALSE)
}, error=function(e) stop(paste("Failed to read features CSV:", e$message)))

# Check essential columns
if (!("mz" %in% names(df))) stop("Features CSV must contain column 'mz'.")

# Normalize RT to minutes if possible
if ("rt" %in% names(df)) {
  df$RT_min <- suppressWarnings(as.numeric(df$rt)) / 60
} else if ("RT" %in% names(df)) {
  df$RT_min <- suppressWarnings(as.numeric(df$RT))
} else if (!("RT_min" %in% names(df))) {
  df$RT_min <- NA_real_
}

# --- Isotope filter: keep only molecular ions; empty/NA is molecular ion ---
if ("isotopes" %in% names(df)) {
  iso <- trimws(as.character(df$isotopes))   # normalize
  empty_iso <- is.na(iso) | iso == ""        # keep empty as molecular ion

  # Accept patterns: [M]+, [M]- and with optional leading index like [10][M]+
  # Regex: optional [digits] then [M] then optional + or -
  molecular_pattern <- "^\\[\\d*\\]\\[M\\][+-]?$"

  is_molecular <- grepl(molecular_pattern, iso)
  keep <- empty_iso | is_molecular

  kept <- sum(keep, na.rm=TRUE)
  dropped <- sum(!keep, na.rm=TRUE)
  cat("Isotope filter: kept", kept, "rows; dropped", dropped, "non-molecular isotopes\n")

  df <- df[keep, , drop = FALSE]
}

# If nothing remains, write empty CSV (with headers) and exit 0
if (nrow(df) == 0) {
  cat("No rows remain after isotope filter; writing empty CSV and exiting.\n")
  res <- df
  if (!("RT_min" %in% names(res))) res$RT_min <- numeric(0)
  if (!("mz" %in% names(res)))     res$mz     <- numeric(0)
  if (!("isotopes" %in% names(res))) res$isotopes <- character(0)
  res$adduct_used      <- character(0)
  res$ExactMass        <- numeric(0)
  res$PredictedFormula <- character(0)
  # keep likely intensity columns too
  int_cols <- grep("(?i)intensity|into|maxo|area", names(res), perl=TRUE, value=TRUE)
  keep_names <- unique(c(
    "RT_min","mz","adduct_used","ExactMass","PredictedFormula","isotopes",
    int_cols, setdiff(names(res), c("rt"))
  ))
  res <- res[, intersect(keep_names, names(res)), drop=FALSE]
  write.csv(res, out_csv, row.names = FALSE)
  cat("Rdisop stage complete (empty) ->", out_csv, "\n")
  quit(status = 0)
}

# Choose adduct: first in 'adduct' column if present, otherwise default by polarity
choose_adduct <- function(adduct, polarity) {
  if (is.na(adduct) || adduct == "") {
    return(ifelse(polarity == "positive", "[M+H]+", "[M-H]-"))
  }
  parts <- strsplit(adduct, ";|,")[[1]]
  trimws(parts[1])
}

# Adduct mass shifts (Da) and charge; extend as needed
adduct_info <- list(
  "[M+H]+"   = c(1.007276,  1),
  "[M+Na]+"  = c(22.989218, 1),
  "[M+K]+"   = c(38.963158, 1),
  "[M+NH4]+" = c(18.033823, 1),
  "[M+2H]2+" = c(2*1.007276, 2),

  "[M-H]-"    = c(-1.007276, -1),
  "[M+Cl]-"   = c(34.969402, -1),
  "[M+FA-H]-" = c(46.00548 - 1.007276, -1),
  "[M+HCOO]-" = c(46.00548 - 1.007276, -1)
)

neutral_mass <- function(mz, adduct){
  if (is.na(mz) || is.na(adduct)) return(NA_real_)
  info <- adduct_info[[adduct]]
  if (is.null(info)) return(NA_real_)
  shift <- info[1]; z <- info[2]
  # neutral mass: M = (|z|*mz) - shift ; then / |z|
  M <- (abs(z) * as.numeric(mz)) - shift
  M / abs(z)
}

# Predict formula with Rdisop within ±0.005 Da, require DBE > 0
predict_formula <- function(neutral_da){
  if (is.na(neutral_da)) return(NA_character_)
  window <- 0.005
  cand <- tryCatch({
    Rdisop::decomposeMass(neutral_da, mzabs=window, elements=c(C=0,H=0,N=0,O=0,P=0,S=0))
  }, error=function(e) NULL)
  if (is.null(cand) || nrow(cand@results) == 0) return(NA_character_)
  tab <- cand@results
  tab$err <- abs(tab$exactmass - neutral_da)
  tab <- tab[order(tab$err), , drop=FALSE]
  # DBE (double bond equivalents) > 0
  tab <- tab[tab$dbd > 0, , drop=FALSE]
  if (nrow(tab) == 0) return(NA_character_)
  tab$formula[1]
}

# Build result
res <- df
has_adduct_col <- "adduct" %in% names(df)
res$adduct_used      <- NA_character_
res$ExactMass        <- NA_real_
res$PredictedFormula <- NA_character_

n <- nrow(df)
for (i in seq_len(n)) {
  mz <- suppressWarnings(as.numeric(df$mz[i]))
  ad <- if (has_adduct_col) choose_adduct(df$adduct[i] %||% NA_character_, polarity) else choose_adduct(NA_character_, polarity)
  nm <- neutral_mass(mz, ad)
  pf <- predict_formula(nm)
  res$adduct_used[i]      <- ad
  res$ExactMass[i]        <- nm
  res$PredictedFormula[i] <- pf
  if (i %% 1000 == 0) cat(".. processed", i, "rows\n")
}

# Keep key + intensity columns first
int_cols <- grep("(?i)intensity|into|maxo|area", names(res), perl=TRUE, value=TRUE)
keep_names <- unique(c(
  "RT_min","mz","adduct_used","ExactMass","PredictedFormula",
  if ("isotopes" %in% names(res)) "isotopes" else character(0),
  int_cols,
  setdiff(names(res), c("rt"))
))
res <- res[, intersect(keep_names, names(res)), drop=FALSE]

# Write out
write.csv(res, out_csv, row.names = FALSE)
cat("Rdisop stage complete ->", out_csv, "\n")
'''

@router.post("/annotate/run")
async def annotate_run(params: AnnotParams):
    directory = _safe(params.directory)
    if not directory.exists():
        raise HTTPException(status_code=400, detail=f"Directory not found: {directory}")

    features_csv = (directory / f"{params.table}.csv").resolve()
    _ensure_exists(features_csv, "Features CSV")

    _ensure_exists(PLANTDB_PARQUET, "DB Parquet")

    analyte_csv_path: Optional[Path] = None
    if params.analyte_csv:
        analyte_csv_path = _safe(params.analyte_csv)
        _ensure_exists(analyte_csv_path, "Analyte CSV")

    # Outputs
    r_out_csv = (directory / f"{params.table}__rdisop_stage.csv").resolve()
    final_out = (directory / "Features_Matrix_Annotated.csv").resolve()

    # Write embedded R to a temp file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as tf:
        tf.write(R_SCRIPT)
        r_path = Path(tf.name)

    cmd = [
        RSCRIPT_BIN, str(r_path),
        f"--features={str(features_csv)}",
        f"--polarity={params.polarity}",
        f"--out={str(r_out_csv)}",
    ]

    async def stream():
        # 1) Run Rdisop stage
        try:
            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.STDOUT,
            )
        except FileNotFoundError:
            yield f"❌ Rscript not found at: {RSCRIPT_BIN}\n"
            yield "__ANNOT_DONE__ EXIT_CODE=1\n"
            return

        yield "=== Annotation (Rdisop) started ===\n"
        async for chunk in proc.stdout:
            yield chunk.decode(errors="ignore")
        r_code = await proc.wait()
        if r_code != 0:
            yield f"❌ R stage failed (exit {r_code}).\n"
            yield "__ANNOT_DONE__ EXIT_CODE=1\n"
            return

        # 2) Python/DuckDB stage: match by PredictedFormula to DB + optional Analyte CSV
        try:
            yield "=== Python stage: matching by PredictedFormula ===\n"
            if not r_out_csv.exists():
                yield f"❌ Intermediate CSV missing: {r_out_csv}\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            con = duckdb.connect()

            # Views
            con.execute("CREATE OR REPLACE VIEW feats AS SELECT * FROM read_csv_auto(?, HEADER=TRUE)", [str(r_out_csv)])
            con.execute("CREATE OR REPLACE VIEW mols  AS SELECT * FROM read_parquet(?)", [str(PLANTDB_PARQUET)])

            # Aggregate candidates by formula from DB
            con.execute("""
                CREATE OR REPLACE VIEW dbagg AS
                SELECT MolecularFormula AS formula, string_agg(CompoundName, ' / ') AS Candidates
                FROM mols
                WHERE MolecularFormula IS NOT NULL
                GROUP BY MolecularFormula
            """)

            analyte_join_sql = ""
            analyte_sel_sql  = ""
            if analyte_csv_path:
                analyte_name = analyte_csv_path.stem
                con.execute("CREATE OR REPLACE VIEW analyte AS SELECT * FROM read_csv_auto(?, HEADER=TRUE)", [str(analyte_csv_path)])
                con.execute(f"""
                    CREATE OR REPLACE VIEW aagg AS
                    SELECT MolecularFormula AS formula, string_agg(CompoundName, ' / ') AS "{analyte_name}"
                    FROM analyte
                    WHERE MolecularFormula IS NOT NULL
                    GROUP BY MolecularFormula
                """)
                analyte_join_sql = "LEFT JOIN aagg A ON A.formula = feats.PredictedFormula"
                analyte_sel_sql  = f', A."{analyte_name}"'

            # Final table
            sql = f"""
            SELECT
              feats.RT_min         AS RT_min,
              feats.mz             AS mz,
              feats.adduct_used    AS adduct,
              feats.ExactMass      AS ExactMass,
              feats.PredictedFormula,
              feats.isotopes,
              feats.* EXCLUDE (RT_min, mz, adduct_used, ExactMass, PredictedFormula, isotopes),
              D.Candidates
              {analyte_sel_sql}
            FROM feats
            LEFT JOIN dbagg D ON D.formula = feats.PredictedFormula
            {analyte_join_sql}
            """
            out_df = con.execute(sql).fetch_df()
            out_df.to_csv(final_out, index=False)

            yield f"✅ Wrote: {final_out}\n"
            yield json.dumps({"ok": True, "output_csv": str(final_out)}) + "\n"
            yield "__ANNOT_DONE__ EXIT_CODE=0\n"
        except Exception as e:
            yield f"❌ Error: {e}\n"
            yield "__ANNOT_DONE__ EXIT_CODE=1\n"
        finally:
            try:
                r_path.unlink(missing_ok=True)
            except Exception:
                pass

    return StreamingResponse(stream(), media_type="text/plain; charset=utf-8")
