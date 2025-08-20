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

# Your bundled DB (Parquet)
PLANTDB_PARQUET = Path("data/moleculesdb.parquet").resolve()

# Set to absolute path to Rscript.exe on Windows if needed
RSCRIPT_BIN = "Rscript"

class AnnotParams(BaseModel):
    directory: str = Field(min_length=1, description="Folder where Features_Matrix.csv lives")
    table: str = Field(min_length=1, description="Base name of features CSV (no .csv). E.g., Features_Matrix")
    polarity: Literal["positive", "negative"]
    # DB = default; CSV is optional and *in addition* to DB (as you asked)
    analyte_csv: Optional[str] = Field(None, description="Optional analyte list CSV (MolecularFormula, ExactMass, CompoundName)")

def _safe(p: str) -> Path:
    return Path(p).expanduser().resolve()

def _exist_file(p: Path, kind: str):
    if not p.exists():
        raise HTTPException(status_code=404, detail=f"{kind} not found: {p}")

R_SCRIPT = r'''
suppressPackageStartupMessages({
  load_or_install <- function(pkg){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "http://cran.us.r-project.org")
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  pkgs <- c("Rdisop", "jsonlite")
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

# Read features
df <- tryCatch({
  read.csv(features_csv, stringsAsFactors = FALSE, check.names = FALSE)
}, error=function(e) stop(paste("Failed to read features CSV:", e$message)))

# Columns we expect / handle
has_mz <- "mz" %in% names(df)
if (!has_mz) stop("Features CSV must contain column 'mz'.")

# Normalize RT_min
if ("rt" %in% names(df)) {
  # many pipelines store rt in seconds; export in minutes
  df$RT_min <- as.numeric(df$rt) / 60
} else if ("RT" %in% names(df)) {
  df$RT_min <- as.numeric(df$RT)
} else if (!"RT_min" %in% names(df)) {
  df$RT_min <- NA_real_
}

# Helper: keep only isotopes that are [M]+ or [M]- ; drop others
if ("isotopes" %in% names(df)) {
  ok_iso <- grepl("^\\s*\\[M[\\+\\-]?\\]\\s*$", df$isotopes %||% "")
  # Keep if NA (no isotope info) or matches [M]+ / [M]-
  keep <- is.na(df$isotopes) | ok_iso
  dropped <- sum(!keep, na.rm=TRUE)
  if (dropped > 0) cat("Dropping", dropped, "non-[M]+/- isotope rows\n")
  df <- df[keep, , drop=FALSE]
}

# Helper: choose adduct
choose_adduct <- function(adduct, polarity) {
  if (is.na(adduct) || adduct == "") {
    return(ifelse(polarity == "positive", "[M+H]+", "[M-H]-"))
  }
  # If multiple (e.g. "[M+H]+;[M+Na]+"), take first
  parts <- strsplit(adduct, ";|,")[[1]]
  trimws(parts[1])
}

# Map adduct to (mass_shift, charge). Add more as needed.
adduct_info <- list(
  "[M+H]+"   = c(1.007276,  1),
  "[M+Na]+"  = c(22.989218, 1),
  "[M+K]+"   = c(38.963158, 1),
  "[M+NH4]+" = c(18.033823, 1),
  "[M+2H]2+" = c(2*1.007276, 2),

  "[M-H]-"   = c(-1.007276, -1),
  "[M+Cl]-"  = c(34.969402, -1),
  "[M+FA-H]-"= c(46.00548 - 1.007276, -1), # ~ formate adduct
  "[M+HCOO]-"= c(46.00548 - 1.007276, -1)
)

neutral_mass <- function(mz, adduct){
  if (is.na(mz) || is.na(adduct)) return(NA_real_)
  info <- adduct_info[[adduct]]
  if (is.null(info)) return(NA_real_)
  shift <- info[1]; z <- info[2]
  M <- (abs(z) * mz) - shift
  M / abs(z)
}

# Coalesce helper
`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else ifelse(is.na(a), b, a)

# Predict formula using Rdisop within 5 mDa window; require DBE>0 and valid
predict_formula <- function(neutral_mass_da, ppm_fallback=5){
  if (is.na(neutral_mass_da)) return(NA_character_)
  # search window 0.005 Da
  window <- 0.005
  # Rdisop::decomposeMass returns candidate formulas; we'll pick the closest with DBE>0
  # We restrict elements to CHNOPS by default (adjust if needed)
  cand <- tryCatch({
    Rdisop::decomposeMass(neutral_mass_da, mzabs=window, elements=c(C=0, H=0, N=0, O=0, P=0, S=0))
  }, error=function(e) NULL)

  if (is.null(cand) || nrow(cand@results) == 0) return(NA_character_)
  tab <- cand@results
  # tab has: formula, exactmass, score, dbd (double bond equivalents)
  tab$err <- abs(tab$exactmass - neutral_mass_da)
  tab <- tab[order(tab$err), , drop=FALSE]
  # Filter: DBE > 0
  tab <- tab[tab$dbd > 0, , drop=FALSE]
  if (nrow(tab) == 0) return(NA_character_)
  # pick best
  tab$formula[1]
}

# Process rows
n <- nrow(df)
res <- df

# Determine adduct source column, if present
has_adduct <- "adduct" %in% names(df)
res$adduct_used <- NA_character_
res$ExactMass <- NA_real_
res$PredictedFormula <- NA_character_

for (i in seq_len(n)) {
  mz <- suppressWarnings(as.numeric(df$mz[i]))
  ad <- NA_character_
  if (has_adduct) ad <- choose_adduct(df$adduct[i] %||% NA_character_, polarity) else ad <- choose_adduct(NA_character_, polarity)
  nm <- neutral_mass(mz, ad)
  pf <- predict_formula(nm)
  res$adduct_used[i] <- ad
  res$ExactMass[i] <- nm
  res$PredictedFormula[i] <- pf
  if (i %% 1000 == 0) cat(".. processed", i, "rows\n")
}

# Keep key columns first, preserve any intensity columns
int_cols <- grep("(?i)intensity|into|maxo|area", names(df), perl=TRUE, value=TRUE)

keep_names <- unique(c(
  "RT_min", "mz", "adduct_used", "ExactMass", "PredictedFormula",
  if ("isotopes" %in% names(df)) "isotopes" else character(0),
  int_cols,
  setdiff(names(df), c("rt")) # keep the rest after key columns
))

res <- res[ , keep_names, drop=FALSE]

# Write out intermediate (R stage) CSV to out_csv
write.csv(res, out_csv, row.names = FALSE)
cat("Rdisop stage complete ->", out_csv, "\n")
'''

@router.post("/annotate/run")
async def annotate_run(params: AnnotParams):
    directory = _safe(params.directory)
    _exist_file(directory, "Directory") if directory.exists() else (_ for _ in ()).throw(HTTPException(400, f"Directory not found: {directory}"))

    features_csv = (directory / f"{params.table}.csv").resolve()
    _exist_file(features_csv, "Features CSV")

    if not PLANTDB_PARQUET.exists():
        raise HTTPException(404, f"DB Parquet not found: {PLANTDB_PARQUET}")

    analyte_csv_path: Optional[Path] = None
    if params.analyte_csv:
        analyte_csv_path = _safe(params.analyte_csv)
        _exist_file(analyte_csv_path, "Analyte CSV")

    # Output names
    r_out_csv = (directory / f"{params.table}__rdisop_stage.csv").resolve()
    final_out = (directory / "Features_Matrix_Annotated.csv").resolve()  # fixed name as requested

    # Write R to temp file
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
        # 1) Run the R stage (Rdisop, adduct/isotope handling)
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

        # 2) Python stage: join on PredictedFormula to DB and optional analyte list
        try:
            yield "=== Python stage: matching by formula ===\n"
            if not r_out_csv.exists():
                yield f"❌ Intermediate CSV missing: {r_out_csv}\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            con = duckdb.connect()

            # R output
            con.execute("CREATE OR REPLACE VIEW feats AS SELECT * FROM read_csv_auto(?, HEADER=TRUE)", [str(r_out_csv)])
            # DB parquet
            con.execute("CREATE OR REPLACE VIEW mols AS SELECT * FROM read_parquet(?)", [str(PLANTDB_PARQUET)])

            # Aggregate DB candidates by PredictedFormula
            db_cand = con.execute("""
                SELECT
                  PredictedFormula,
                  string_agg(CompoundName, ' / ') AS Candidates
                FROM mols
                WHERE MolecularFormula IS NOT NULL
                GROUP BY PredictedFormula := MolecularFormula
            """).fetch_df()

            # Left join with feats
            con.execute("CREATE OR REPLACE VIEW db_cand AS SELECT * FROM db_cand_temp",)
        except duckdb.CatalogException:
            # duckdb can't create view from a pandas df directly; load via values
            pass
        finally:
            # Use a single SQL to produce final table with optional analyte matches
            con = duckdb.connect()
            con.execute("CREATE OR REPLACE VIEW feats AS SELECT * FROM read_csv_auto(?, HEADER=TRUE)", [str(r_out_csv)])
            con.execute("CREATE OR REPLACE VIEW mols AS SELECT * FROM read_parquet(?)", [str(PLANTDB_PARQUET)])

            # Optional analyte list
            analyte_join_sql = ""
            analyte_sel_sql = ""
            if analyte_csv_path:
                analyte_name = analyte_csv_path.stem  # column name requested
                con.execute("CREATE OR REPLACE VIEW analyte AS SELECT * FROM read_csv_auto(?, HEADER=TRUE)", [str(analyte_csv_path)])
                analyte_join_sql = f"""
                    LEFT JOIN (
                      SELECT MolecularFormula AS a_formula,
                             string_agg(CompoundName, ' / ') AS "{analyte_name}"
                      FROM analyte
                      WHERE MolecularFormula IS NOT NULL
                      GROUP BY MolecularFormula
                    ) A ON A.a_formula = feats.PredictedFormula
                """
                analyte_sel_sql = f', A."{analyte_name}"'

            # Final select (DB candidates by formula, plus optional analyte list)
            sql = f"""
            WITH DBAGG AS (
              SELECT MolecularFormula AS f, string_agg(CompoundName, ' / ') AS Candidates
              FROM mols
              WHERE MolecularFormula IS NOT NULL
              GROUP BY MolecularFormula
            )
            SELECT
              feats.RT_min AS RT_min,
              feats.mz,
              feats.adduct_used AS adduct,
              feats.ExactMass,
              feats.PredictedFormula,
              feats.isotopes,
              feats.* EXCLUDE (RT_min, mz, adduct_used, ExactMass, PredictedFormula, isotopes),
              D.Candidates
              {analyte_sel_sql}
            FROM feats
            LEFT JOIN DBAGG D ON D.f = feats.PredictedFormula
            {analyte_join_sql}
            """
            out_df = con.execute(sql).fetch_df()
            out_df.to_csv(final_out, index=False)

            yield f"✅ Wrote: {final_out}\n"
            yield json.dumps({"ok": True, "output_csv": str(final_out)}) + "\n"
            yield "__ANNOT_DONE__ EXIT_CODE=0\n"

        # Cleanup temp R file
        try:
            r_path.unlink(missing_ok=True)
        except Exception:
            pass

    return StreamingResponse(stream(), media_type="text/plain; charset=utf-8")
