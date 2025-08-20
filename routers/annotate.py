from fastapi import APIRouter
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from typing import Literal, Optional
from pathlib import Path
import pandas as pd
import duckdb
import re
import os

router = APIRouter()

# === Config ===
# Point this to your parquet DB (or set env var PLANTDB_PARQUET)
PLANTDB_PARQUET = Path(
    os.getenv("PLANTDB_PARQUET", r"C:\Users\Administrator\PycharmProjects\PlantResearch\data\moleculesdb.parquet")
).resolve()

H_ATOM = 1.007825  # proton mass per your workflow
MASS_TOL = 0.002   # ± 0.002 Da window for DB join


class AnnotSemiParams(BaseModel):
    directory: str = Field(..., description="Folder that contains the features CSV")
    table: str = Field(..., description="CSV base name without .csv (e.g., 'Features_Matrix')")
    polarity: Literal["positive", "negative"]
    # optional: if your features file name differs, you can pass it in full:
    features_csv_override: Optional[str] = None


ADDUCT_RE = re.compile(r"(\[[^\]]+\][+-]?\d*)\s+([\d\.]+)")

def split_adducts(adduct_str: str):
    """
    From CAMERA 'adduct' column: extract pairs like "[M+H]+ 123.456"
    Returns list of (adduct_token, mass_value_string)
    """
    if adduct_str is None:
        return []
    s = str(adduct_str)
    if not s or s.strip() == "" or s.strip().lower() == "na":
        return []
    return ADDUCT_RE.findall(s)


def is_molecular_ion_isotope_tag(tag: Optional[str]) -> bool:
    """
    True if isotopes is empty/NA, or a molecular-ion tag like "[10][M]" or "[2][M]+".
    We treat these rows as OK to default to [M+H]+ / [M-H]- when adduct is missing.
    """
    if tag is None:
        return True
    s = str(tag).strip()
    if s == "" or s.lower() == "na":
        return True
    # [number][M] optionally followed by +/-
    return re.match(r"^\[\d+\]\[M\][+-]?$", s) is not None


def _safe_path(p: str) -> Path:
    return Path(p).expanduser().resolve()


def _fmt(x) -> str:
    try:
        return str(x)
    except Exception:
        return repr(x)


def _build_stream(params: AnnotSemiParams):
    async def stream():
        directory = _safe_path(params.directory)
        features_csv = (
            _safe_path(params.features_csv_override)
            if params.features_csv_override
            else (directory / f"{params.table}.csv").resolve()
        )
        out_csv = (directory / f"{params.table}_Annotated.csv").resolve()

        # Pre-flight checks BEFORE yielding (safe to stop with messages)
        if not directory.exists():
            yield f"❌ Directory not found: {directory}\n"
            yield "__ANNOT_DONE__ EXIT_CODE=1\n"
            return
        if not features_csv.exists():
            yield f"❌ Features CSV not found: {features_csv}\n"
            yield "__ANNOT_DONE__ EXIT_CODE=1\n"
            return

        db_exists = PLANTDB_PARQUET.exists()
        if not db_exists:
            yield f"⚠ DB parquet not found: {PLANTDB_PARQUET}\n"
            yield "   → Will skip Candidates join and still produce output.\n"

        # Begin processing
        yield "=== Semi-automated annotation ===\n"
        yield f"Loading features: {features_csv}\n"
        try:
            df = pd.read_csv(features_csv)
        except Exception as e:
            yield f"❌ Could not read CSV: {e}\n"
            yield "__ANNOT_DONE__ EXIT_CODE=1\n"
            return

        # Ensure mz exists
        if "mz" not in df.columns:
            yield "❌ CSV must include 'mz' column.\n"
            yield "__ANNOT_DONE__ EXIT_CODE=1\n"
            return

        # Ensure isotopes column exists (optional)
        if "isotopes" not in df.columns:
            df["isotopes"] = ""

        # RT in minutes
        if "RT_min" not in df.columns:
            if "rt" in df.columns:
                df["RT_min"] = pd.to_numeric(df["rt"], errors="coerce") / 60.0
            elif "RT" in df.columns:
                df["RT_min"] = pd.to_numeric(df["RT"], errors="coerce")
            else:
                df["RT_min"] = pd.NA

        # Split CAMERA adducts and keep FIRST pair only
        yield "Splitting CAMERA adducts and picking first...\n"
        adduct_1 = []
        mass_1 = []
        for s in df.get("adduct", "").astype(str).tolist():
            pairs = split_adducts(s)
            if pairs:
                a, m = pairs[0]
                adduct_1.append(a)
                try:
                    mass_1.append(float(m))
                except Exception:
                    mass_1.append(float("nan"))
            else:
                adduct_1.append("")
                mass_1.append(float("nan"))

        df["adduct_1"] = adduct_1
        df["adduct_mass_1"] = mass_1

        # Default adduct & ExactMass
        yield "Applying default adducts where appropriate and computing ExactMass...\n"
        pol = params.polarity
        mz = pd.to_numeric(df["mz"], errors="coerce")

        # Start ExactMass with NaN; fill from adduct_mass_1 if present
        df["ExactMass"] = pd.Series([pd.NA] * len(df), dtype="float64")

        # Case 1: first adduct mass available -> trust CAMERA value
        has_mass = pd.notna(df["adduct_mass_1"])
        df.loc[has_mass, "ExactMass"] = df.loc[has_mass, "adduct_mass_1"]

        # Case 2: no adduct mass AND row is eligible for default (not an isotope peak)
        need_default = (~has_mass) & df["isotopes"].apply(is_molecular_ion_isotope_tag)
        if pol == "positive":
            df.loc[need_default, "adduct_1"] = df.loc[need_default, "adduct_1"].replace("", "[M+H]+")
            df.loc[need_default, "ExactMass"] = mz[need_default] - H_ATOM
        else:
            df.loc[need_default, "adduct_1"] = df.loc[need_default, "adduct_1"].replace("", "[M-H]-")
            df.loc[need_default, "ExactMass"] = mz[need_default] + H_ATOM

        # Report stats
        n_total = len(df)
        n_from_cam = int(has_mass.sum())
        n_defaulted = int(need_default.sum())
        n_blank = int(df["ExactMass"].isna().sum())
        yield f"Rows: total={n_total}, from CAMERA mass={n_from_cam}, defaulted={n_defaulted}, still blank={n_blank}\n"

        # Build output frame: keep intensities too
        intensity_cols = [c for c in df.columns if re.search(r"(?i)(intensity|into|maxo|area)", c)]
        cols_out = ["RT_min", "mz", "adduct_1", "ExactMass"] + intensity_cols

        out = df.copy()
        out = out[cols_out].copy()

        # Annotate with DB Candidates (±0.002 Da on ExactMass), if DB exists
        if db_exists:
            yield f"Joining DB parquet (±{MASS_TOL} Da)...\n"
            con = duckdb.connect()
            con.register("ft", out[out["ExactMass"].notna()].copy())

            db_path_sql = str(PLANTDB_PARQUET).replace("'", "''")
            con.execute(f"""
                CREATE OR REPLACE VIEW mols AS
                SELECT CompoundName, MolecularFormula, ExactMass
                FROM read_parquet('{db_path_sql}')
                WHERE ExactMass IS NOT NULL
            """)

            # Find candidates within mass tolerance, aggregate names
            cand_df = con.execute(f"""
                WITH hits AS (
                    SELECT
                        ft.ROW_NUMBER() OVER () AS rid,
                        ft.ExactMass AS q_mass,
                        mols.CompoundName,
                        mols.ExactMass  AS db_mass,
                        abs(mols.ExactMass - ft.ExactMass) AS dm
                    FROM ft
                    JOIN mols
                      ON abs(mols.ExactMass - ft.ExactMass) <= {MASS_TOL}
                )
                SELECT q_mass,
                       string_agg(CompoundName, ' / ' ORDER BY dm, CompoundName) AS Candidates
                FROM hits
                GROUP BY q_mass
            """).fetch_df()

            out = out.merge(cand_df, left_on="ExactMass", right_on="q_mass", how="left").drop(columns=["q_mass"])
        else:
            out["Candidates"] = pd.NA

        # Save final
        out.to_csv(out_csv, index=False)
        yield f"✅ Wrote: {out_csv}\n"
        yield "__ANNOT_DONE__ EXIT_CODE=0\n"

    return stream


@router.post("/annotation/semi-run")
async def annotation_semi_run(params: AnnotSemiParams):
    """
    Semi-automated annotation:
      - Split CAMERA adducts; pick first adduct+mass
      - Default to [M+H]+ / [M-H]- when no adduct and not an isotope peak
      - ExactMass from adduct mass or neutral mz ± 1.007825
      - Join to Parquet DB within ±0.002 Da → Candidates
      - Writes <table>_Annotated.csv
      - Streams logs with __ANNOT_DONE__ EXIT_CODE=*
    """
    return StreamingResponse(_build_stream(params)(), media_type="text/plain; charset=utf-8")
