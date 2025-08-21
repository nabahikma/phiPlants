# routers/annotate.py
from __future__ import annotations

import os
import re
import bisect
import subprocess
from pathlib import Path
from typing import Optional, Tuple, List, Literal, AsyncGenerator

import numpy as np
import pandas as pd
from fastapi import APIRouter, Form
from fastapi.responses import StreamingResponse, FileResponse, JSONResponse
from pydantic import BaseModel, Field

router = APIRouter()

# =========================
# Mass tables & helpers
# =========================

# Neutral atom monoisotopic masses (u)
ELEMENT_MASS = {
    'H': 1.00782503223, 'C': 12.0, 'N': 14.00307400443, 'O': 15.99491461957,
    'F': 18.998403163, 'Na': 22.9897692820, 'Mg': 23.985041699, 'Al': 26.98153853,
    'Si': 27.97692653465, 'P': 30.97376199842, 'S': 31.9720711744, 'Cl': 34.968852682,
    'K': 38.9637064864, 'Ca': 39.962590863, 'Br': 78.9183376, 'I': 126.9044719,
}

# Charged adduct ion masses (for m/z arithmetic)
ION_MASS = {
    'H': 1.007276466812,   # H+
    'Na': 22.989218,       # Na+
    'K': 38.963158,        # K+
    'NH4': 18.033823,      # NH4+
    'Cl': 34.969402,       # Cl-
}

# Neutral losses (uncharged)
NEUTRAL_LOSSES = {
    'H2O': 18.010564684, 'NH3': 17.026549101, 'CO2': 43.989829, 'C2H4': 28.031300,
    'H': ELEMENT_MASS['H'], 'H2': 2 * ELEMENT_MASS['H']
}

MASS_TOL_DA_DEFAULT = 0.002  # ±0.002 Da

def formula_mass(formula: str) -> Optional[float]:
    """Neutral mass from a formula like C2H4 or CH3OH. Returns None if unparsable."""
    if not isinstance(formula, str) or not formula.strip():
        return None
    if formula in ION_MASS:  # ion keyword, not a neutral formula
        return None
    total = 0.0
    pos = 0
    for m in re.finditer(r'([A-Z][a-z]?)(\d*)', formula):
        el, count = m.group(1), m.group(2)
        if el not in ELEMENT_MASS:
            return None
        n = int(count) if count else 1
        total += ELEMENT_MASS[el] * n
        pos = m.end()
    if pos != len(formula):
        return None
    return total

def parse_first_adduct(adduct_str: str) -> Tuple[Optional[str], Optional[dict]]:
    """
    Parse ONLY the FIRST adduct token from a possibly messy 'adduct' cell, e.g.:
      [M+H]+, [M-H]-, [2M+Na+3K-H]3+, [M+H-C2H4]+
    Returns (normalized_pattern, info) where info={n, z, plus[], minus[]}.
    """
    if not isinstance(adduct_str, str) or not adduct_str.strip():
        return None, None
    s = adduct_str.strip()
    m = re.search(r'\[([^\]]+)\]\s*([0-9]*)([+-])', s)
    if not m:
        return None, None
    inside = m.group(1)
    zmag = int(m.group(2)) if m.group(2) else 1
    sign = m.group(3)
    z = zmag if sign == '+' else -zmag
    pattern_str = f'[{inside}]{zmag if zmag != 1 else ""}{sign}'

    # split tokens inside [] like +H, +Na, -H2O, +3K, +C2H4, and "2M"
    tokens = re.findall(r'([+-]?)([^+-]+)', inside)
    n = 1
    plus, minus = [], []
    for sgn, tok in tokens:
        tok = tok.strip()
        if not tok:
            continue
        if tok.endswith('M'):  # 2M, 3M
            num = tok[:-1]
            n = int(num) if num else 1
            continue
        m2 = re.match(r'(\d+)([A-Za-z0-9]+)', tok)  # 3K, 2Na, 2H2O, etc.
        count, base = (int(m2.group(1)), m2.group(2)) if m2 else (1, tok)
        (minus if sgn == '-' else plus).append((base, count))
    return pattern_str, {'n': n, 'z': z, 'plus': plus, 'minus': minus}

def compute_neutral_mass_from_adduct(mz: float, adduct_pattern: str) -> Optional[float]:
    """
    Compute neutral ExactMass (M) from m/z and first adduct:
      M = (|z|*mz - sum(plus) + sum(minus)) / n
    Supports multimers, charges, common ions & neutral losses.
    """
    pattern_str, info = parse_first_adduct(adduct_pattern)
    if info is None:
        return None
    n, z = info['n'], info['z']

    plus_mass = 0.0
    for base, count in info['plus']:
        if base in ION_MASS:
            plus_mass += ION_MASS[base] * count
        else:
            plus_mass += (formula_mass(base) or 0.0) * count

    minus_mass = 0.0
    for base, count in info['minus']:
        # Special case: deprotonation '-H' with negative charge → use H+ mass
        if base == 'H' and z < 0:
            minus_mass += ION_MASS['H'] * count
        elif base in NEUTRAL_LOSSES:
            minus_mass += NEUTRAL_LOSSES[base] * count
        else:
            minus_mass += (formula_mass(base) or 0.0) * count

    M = ((abs(z) * mz) - plus_mass + minus_mass) / max(n, 1)
    return M

def is_non_molecular_isotope(val) -> bool:
    """
    True → skip (isotope like M+1, [M+2], etc.)
    False → keep (empty/NA, 'M', '[M]', or CAMERA style '[10][M]' / '[2][M]+').
    """
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return False
    s = str(val).strip()
    if not s or s.lower() == 'na':
        return False
    if s in {'M', '[M]'}:
        return False
    # CAMERA-like molecular ion tags (keep)
    if re.fullmatch(r'\[\d+\]\[M\][+-]?', s):
        return False
    # M+1, [M+2], etc. → skip
    if re.search(r'\[?M\+\d+\]?', s):
        return True
    return False

# =========================
# Pydantic payload
# =========================

class AnnotateParams(BaseModel):
    directory: str = Field(..., description="Folder containing the features CSV")
    table: str = Field(..., description="Base name without .csv (e.g., 'Features_Matrix')")
    polarity: Literal['positive', 'negative']
    analyte_csv: Optional[str] = Field(None, description="Optional PlantDB CSV path")
    tol: float = Field(MASS_TOL_DA_DEFAULT, description="Mass tolerance in Da (±)")

# =========================
# Core annotation
# =========================

def match_candidates_by_mass(
    exact_masses: np.ndarray,
    plant_df: pd.DataFrame,
    tol_da: float
) -> Tuple[List[Optional[str]], List[Optional[str]]]:
    """
    For each ExactMass, find PlantDB rows with |db.ExactMass - ExactMass| <= tol_da.
    Returns (CandidateFormula, Candidates).
    """
    plant = plant_df.copy()

    # Normalize mass column: prefer ExactMass, fall back to MonoMass
    if 'ExactMass' not in plant.columns and 'MonoMass' in plant.columns:
        plant = plant.rename(columns={'MonoMass': 'ExactMass'})

    plant['ExactMass'] = pd.to_numeric(plant['ExactMass'], errors='coerce')
    plant = plant.dropna(subset=['ExactMass']).sort_values('ExactMass').reset_index(drop=True)

    pmasses = plant['ExactMass'].to_numpy()
    have_formula = 'MolecularFormula' in plant.columns
    have_name = 'CompoundName' in plant.columns

    out_formula: List[Optional[str]] = []
    out_names: List[Optional[str]] = []

    for M in exact_masses:
        if M is None or (isinstance(M, float) and np.isnan(M)):
            out_formula.append(None)
            out_names.append(None)
            continue
        lo, hi = M - tol_da, M + tol_da
        i = bisect.bisect_left(pmasses, lo)
        j = bisect.bisect_right(pmasses, hi)
        if i >= j:
            out_formula.append(None)
            out_names.append(None)
            continue

        sub = plant.iloc[i:j]

        # Names (dedup, joined with '/')
        if have_name:
            seen = set()
            names = []
            for n in sub['CompoundName']:
                if pd.isna(n):
                    continue
                n = str(n)
                if n not in seen:
                    seen.add(n)
                    names.append(n)
            out_names.append('/'.join(names) if names else None)
        else:
            out_names.append(None)

        # One formula only (first non-empty)
        if have_formula:
            f = next((x for x in sub['MolecularFormula'].tolist() if pd.notna(x) and str(x).strip()), None)
            out_formula.append(f)
        else:
            out_formula.append(None)

    return out_formula, out_names

def annotate_frame_with_plantdb(
    features_df: pd.DataFrame,
    plant_df: Optional[pd.DataFrame],
    polarity: str,
    tol_da: float
) -> pd.DataFrame:
    df = features_df.copy()

    # 1) Filter isotopes (keep only molecular ion / empty)
    if 'isotopes' in df.columns:
        df = df[~df['isotopes'].map(is_non_molecular_isotope)].copy()

    # 2) Normalize/assign adduct; 3) Compute ExactMass from m/z + adduct
    if 'mz' not in df.columns:
        raise ValueError("Features CSV must include an 'mz' column.")

    pol = polarity.lower()
    default_adduct = '[M+H]+' if pol.startswith('pos') else '[M-H]-'

    used_adducts: List[str] = []
    exact_masses: List[float] = []

    # Iterate row-wise to pair each mz with its own adduct cell
    for row in df.itertuples(index=False):
        mz_val = float(getattr(row, 'mz'))
        raw_adduct = getattr(row, 'adduct') if hasattr(row, 'adduct') else None
        norm, _ = parse_first_adduct(raw_adduct) if isinstance(raw_adduct, str) else (None, None)
        if norm is None:
            norm = default_adduct

        M = compute_neutral_mass_from_adduct(mz_val, norm)
        if M is None or np.isnan(M):
            # conservative fallback if parsing fails
            M = mz_val - ION_MASS['H'] if pol.startswith('pos') else mz_val + ION_MASS['H']

        used_adducts.append(norm)
        exact_masses.append(M)

    df['Adduct'] = used_adducts
    df['ExactMass'] = exact_masses

    # 4) Candidates from PlantDB (if provided)
    if plant_df is not None:
        cand_formula, cand_names = match_candidates_by_mass(df['ExactMass'].to_numpy(), plant_df, tol_da)
        df['CandidateFormula'] = cand_formula
        df['Candidates'] = cand_names
    else:
        df['CandidateFormula'] = None
        df['Candidates'] = None

    # 5) Output columns & RT minutes
    if 'rt' in df.columns:
        df['RT in min'] = pd.to_numeric(df['rt'], errors='coerce') / 60.0
    elif 'RT' in df.columns:
        df['RT in min'] = pd.to_numeric(df['RT'], errors='coerce')
    else:
        df['RT in min'] = np.nan

    meta = {
        'mz','mzmin','mzmax','rt','rtmin','rtmax','npeaks','blank','sample','ms_level',
        'isotopes','adduct','pcgroup','feature_id','has_ms2',
        'RT in min','ExactMass','CandidateFormula','Candidates','Adduct'
    }

    intensity_cols = [c for c in df.columns
                      if c.endswith('.mzML') or (c not in meta and np.issubdtype(df[c].dtype, np.number))]

    for c in ['pcgroup', 'feature_id', 'has_ms2']:
        if c not in df.columns:
            df[c] = None

    out_cols = ['RT in min','mz','ExactMass','CandidateFormula','Candidates','Adduct'] \
               + intensity_cols + ['pcgroup','feature_id','has_ms2']
    out_cols = [c for c in out_cols if c in df.columns]
    return df[out_cols].copy()

# =========================
# I/O helpers
# =========================

def _safe_path(p: str) -> Path:
    return Path(p).expanduser().resolve()

def _project_root() -> Path:
    # routers/ -> project root
    return Path(__file__).resolve().parents[1]

def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)

def _load_plantdb_csv(directory: Path, analyte_csv: Optional[str]) -> Optional[pd.DataFrame]:
    """
    Search order:
      1) analyte_csv param (if provided)
      2) env var PLANTDB_CSV
      3) <directory>/PlantDB_ESI.csv
      4) <directory>/PlantDB.csv
      5) <project_root>/data/PlantDB_ESI.csv
      6) <project_root>/data/PlantDB.csv
    Returns a DataFrame or None. Accepts MonoMass as fallback to ExactMass.
    """
    cand_paths: List[Path] = []

    if analyte_csv:
        cand_paths.append(_safe_path(analyte_csv))

    env_p = os.getenv("PLANTDB_CSV")
    if env_p:
        cand_paths.append(_safe_path(env_p))

    cand_paths.append(directory / "PlantDB_ESI.csv")
    cand_paths.append(directory / "PlantDB.csv")

    root = _project_root()
    cand_paths.append(root / "data" / "PlantDB_ESI.csv")
    cand_paths.append(root / "data" / "PlantDB.csv")

    for p in cand_paths:
        try:
            if p.exists():
                df = pd.read_csv(p)
                # normalize mass column
                if "ExactMass" not in df.columns and "MonoMass" in df.columns:
                    df = df.rename(columns={"MonoMass": "ExactMass"})
                if "ExactMass" not in df.columns:
                    continue
                df["__source_path__"] = str(p)
                return df
        except Exception:
            continue
    return None

def _csv_path(directory: str, file: str) -> Path:
    return Path(directory).expanduser().resolve() / file

def _load_df_for_update(csv_path: Path) -> pd.DataFrame:
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    if 'feature_id' in df.columns:
        df['feature_id'] = df['feature_id'].astype(str)
    return df

def _find_r_script(script_name: str | None = None) -> Optional[Path]:
    """
    Locate an R script under <project_root>/Scripts.
    Priority:
      1) Scripts/<script_name> (if provided)
      2) Scripts/finalReport.R
      3) Scripts/generate_final_report.R
      4) Scripts/GenerateReport.R
      5) first *.R in Scripts/
    """
    scripts_dir = _project_root() / "Scripts"
    if script_name:
        p = scripts_dir / script_name
        if p.exists():
            return p
    for cand in ["finalReport.R", "generate_final_report.R", "GenerateReport.R", "report.R"]:
        p = scripts_dir / cand
        if p.exists():
            return p
    if scripts_dir.exists():
        rs = sorted(scripts_dir.glob("*.R"))
        if rs:
            return rs[0]
    return None

# =========================
# Streaming annotation endpoint used by UI
# =========================

@router.post("/annotate/run")
def annotate_run(params: AnnotateParams):
    """
    POST JSON:
    {
      "directory": "C:\\path\\to\\folder",
      "table": "Features_Matrix",
      "polarity": "positive" | "negative",
      "analyte_csv": "C:\\path\\to\\PlantDB_ESI.csv",  # optional
      "tol": 0.002
    }
    Streams logs, writes <directory>/Features_Matrix_Annotated.csv,
    and prints '__ANNOT_DONE__ EXIT_CODE=*' at the end.
    """

    async def stream() -> AsyncGenerator[str, None]:
        try:
            directory = _safe_path(params.directory)
            features_csv = (directory / f"{params.table}.csv").resolve()
            out_csv = (directory / "Features_Matrix_Annotated.csv").resolve()

            if not directory.exists():
                yield f"❌ Directory not found: {directory}\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return
            if not features_csv.exists():
                yield f"❌ Features CSV not found: {features_csv}\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            yield "=== Annotation (PlantDB CSV) ===\n"
            yield f"Loading features: {features_csv}\n"
            try:
                feat = _read_csv(features_csv)
            except Exception as e:
                yield f"❌ Could not read features CSV: {e}\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            plant_df = _load_plantdb_csv(directory, params.analyte_csv)
            if plant_df is None:
                yield "⚠ PlantDB CSV not found (analyte_csv / PLANTDB_CSV / Directory / project data).\n"
            else:
                src = plant_df.get("__source_path__", "<unknown>")
                n_rows = len(plant_df)
                cols_preview = ", ".join(list(plant_df.columns)[:6])
                yield f"Loaded PlantDB CSV: {src}\n"
                yield f"Rows: {n_rows} | Columns (head): {cols_preview}\n"
                # keep a clean copy without helper col
                plant_df = plant_df.drop(columns=["__source_path__"], errors="ignore")

            yield "Parsing adducts, filtering isotopes, computing ExactMass...\n"
            try:
                out = annotate_frame_with_plantdb(
                    features_df=feat,
                    plant_df=plant_df,
                    polarity=params.polarity,
                    tol_da=float(params.tol or MASS_TOL_DA_DEFAULT)
                )
            except Exception as e:
                yield f"❌ Annotation failed: {e}\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            try:
                out.to_csv(out_csv, index=False)
            except Exception as e:
                yield f"❌ Failed to write output CSV: {e}\n"
                yield "__ANNOT_DONE__ EXIT_CODE=1\n"
                return

            n_rows_out = len(out)
            n_with_cand = int(out['Candidates'].notna().sum()) if 'Candidates' in out.columns else 0
            yield f"✅ Wrote: {out_csv}\n"
            yield f"Rows annotated: {n_rows_out}, with candidates: {n_with_cand}\n"
            yield "__ANNOT_DONE__ EXIT_CODE=0\n"
        except Exception as e:
            yield f"❌ Unexpected error: {e}\n"
            yield "__ANNOT_DONE__ EXIT_CODE=1\n"

    return StreamingResponse(stream(), media_type="text/plain; charset=utf-8")

# =========================
# Annotation editing + report endpoints
# =========================

@router.post("/annotate/update_annotation")
def update_annotation(
    directory: str = Form(...),
    file: str = Form(...),
    feature_id: str = Form(...),
    selected_annotation: str = Form(...)
):
    """
    Persist the selected candidate into the CSV under 'SelectedAnnotation' column.
    """
    try:
        csv_path = _csv_path(directory, file)
        df = _load_df_for_update(csv_path)
        if 'feature_id' not in df.columns:
            return JSONResponse({"ok": False, "error": "feature_id column not found"}, status_code=400)

        if 'SelectedAnnotation' not in df.columns:
            df['SelectedAnnotation'] = pd.NA

        mask = df['feature_id'].astype(str) == str(feature_id)
        if not mask.any():
            return JSONResponse({"ok": False, "error": f"feature_id {feature_id} not found"}, status_code=404)

        df.loc[mask, 'SelectedAnnotation'] = selected_annotation
        tmp = csv_path.with_suffix('.tmp.csv')
        df.to_csv(tmp, index=False)
        tmp.replace(csv_path)

        return JSONResponse({"ok": True})
    except Exception as e:
        return JSONResponse({"ok": False, "error": str(e)}, status_code=500)


def _sanitize_header_val(s: str, maxlen: int = 1024) -> str:
    """Make a safe, single-line ASCII-ish header value."""
    if not s:
        return ""
    s = re.sub(r'[\r\n]+', ' | ', s)            # no newlines in headers
    s = ''.join(ch if 32 <= ord(ch) <= 126 else ' ' for ch in s)
    return s[:maxlen].strip()

@router.post("/annotate/generate_final_report")
def generate_final_report(payload: dict):
    """
    JSON:
      { "directory": "...", "file": "Features_Matrix_Annotated.csv", "feature_ids": ["..."] }
    Writes <directory>/finalReport.csv, then calls R if available:
      save_final_plots_and_report.R <mzml_folder> <finalReport.csv> <Final_Report.xlsx>
    Returns finalReport.csv and writes Final_Report_R.log with stdout/stderr.
    """
    try:
        directory = payload.get("directory")
        file = payload.get("file")
        feature_ids = payload.get("feature_ids") or []
        script_name = payload.get("script_name")  # optional

        if not directory or not file:
            return JSONResponse({"error": "directory and file are required"}, status_code=400)
        if not feature_ids:
            return JSONResponse({"error": "feature_ids is empty"}, status_code=400)

        src_csv = _csv_path(directory, file)
        df = _load_df_for_update(src_csv)
        sel = df[df['feature_id'].astype(str).isin([str(x) for x in feature_ids])].copy()
        if sel.empty:
            return JSONResponse({"error": "No matching features found in CSV"}, status_code=404)

        out_dir = Path(directory).expanduser().resolve()
        out_csv = out_dir / "finalReport.csv"
        sel.to_csv(out_csv, index=False)

        r_exec = os.getenv("RSCRIPT_EXE", "Rscript")
        r_script = _find_r_script(script_name)

        headers = {}
        if r_script is not None:
            sb = r_script.name.lower()
            if "save_final_plots_and_report" in sb:
                # IMPORTANT: pass the CSV with feature_id column
                output_xlsx_name = "Final_Report.xlsx"
                cmd = [r_exec, str(r_script), str(out_dir), str(out_csv), output_xlsx_name]
            else:
                # fallback signature
                cmd = [r_exec, str(r_script), str(out_csv), str(out_dir)]

            proc = subprocess.run(cmd, capture_output=True, text=True)

            # write log next to data
            log_path = out_dir / "Final_Report_R.log"
            try:
                with open(log_path, "w", encoding="utf-8") as f:
                    f.write("Command:\n" + " ".join(cmd) + "\n\n")
                    f.write(f"Exit code: {proc.returncode}\n\n")
                    f.write("=== STDOUT ===\n" + (proc.stdout or "") + "\n\n")
                    f.write("=== STDERR ===\n" + (proc.stderr or "") + "\n")
            except Exception:
                pass

            headers = {
                "X-R-Status": "0" if proc.returncode == 0 else str(proc.returncode),
                "X-R-Stderr": _sanitize_header_val(proc.stderr or ""),
                "X-R-Log": log_path.name,
            }
        else:
            headers = {"X-R-Status": "Script not found in Scripts/"}

        return FileResponse(
            path=str(out_csv),
            filename="finalReport.csv",
            media_type="text/csv",
            headers=headers
        )
    except Exception as e:
        return JSONResponse({"error": str(e)}, status_code=500)


