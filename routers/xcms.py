import uvicorn
from fastapi import FastAPI, Request, HTTPException
from starlette.responses import HTMLResponse
from starlette.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
import asyncio
import shlex
from fastapi import APIRouter, Form
import subprocess



router = APIRouter()



# ===========================
# /xcms/start (R embedded)
# ===========================
import asyncio
import tempfile
from pathlib import Path
from typing import Optional
from fastapi import HTTPException
from pydantic import BaseModel, Field, validator

RSCRIPT_BIN = "Rscript"  # e.g., "/usr/local/bin/Rscript"

class XcmsParams(BaseModel):
    output_name: str = Field(..., min_length=1, max_length=128)
    directory: str = Field(..., min_length=1)
    polarity: str = Field(..., regex="^(positive|negative)$")
    ppm: float = Field(..., gt=0)
    snr: float = Field(..., gt=0)
    noise: float = Field(..., ge=0)
    min_peakwidth: float = Field(2.0, gt=0)
    max_peakwidth: float = Field(30.0, gt=0)
    mzdiff: float = Field(0.01)
    timeout_sec: Optional[int] = Field(0, ge=0)

    @validator("output_name")
    def safe_output(cls, v):
        import re
        if not re.match(r"^[A-Za-z0-9._-]+$", v):
            raise ValueError("output_name allows letters, numbers, ., _, - only")
        return v

    @validator("max_peakwidth")
    def peakwidth_order(cls, v, values):
        if "min_peakwidth" in values and v <= values["min_peakwidth"]:
            raise ValueError("max_peakwidth must be greater than min_peakwidth")
        return v

def _ensure_ready(params: XcmsParams):
    d = Path(params.directory).expanduser().resolve()
    if not d.exists() or not d.is_dir():
        raise HTTPException(status_code=400, detail=f"Directory not found: {d}")
    if not any(p.suffix.lower() == ".mzml" for p in d.glob("*.mzml")):
        # also check *.mzML (mixed case)
        if not any(p.suffix == ".mzML" for p in d.glob("*.mzML")):
            raise HTTPException(status_code=400, detail=f"No mzML files found in: {d}")
    return d

R_SCRIPT_EMBEDDED = r'''
# ==== Safe Library Loader ====
load_or_install <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

pkgs <- c("xcms","CAMERA","tibble","ggplot2","purrr","MSnbase","dplyr","gridExtra")
invisible(lapply(pkgs, load_or_install))

# ==== Parse CLI Args ====
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default=NULL) {
  hit <- grep(paste0("^--", name, "="), args, value=TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit)
}

data_dir      <- get_arg("dir")
polarity_mode <- get_arg("polarity", "negative")
output_name   <- get_arg("output", "processed_data")
ppm           <- as.numeric(get_arg("ppm", 5))
snthr         <- as.numeric(get_arg("snr", 5))
noise         <- as.numeric(get_arg("noise", 30))
min_pw        <- as.numeric(get_arg("min_peakwidth", 4))
max_pw        <- as.numeric(get_arg("max_peakwidth", 50))
mzdiff        <- as.numeric(get_arg("mzdiff", -0.001))

if (is.null(data_dir) || !dir.exists(data_dir)) {
  stop(paste("Invalid data_dir:", data_dir))
}

# ==== Load mzML (MS1 + MS2) ====
files <- list.files(data_dir, pattern="mzML$", full.names=TRUE, ignore.case=TRUE)
if (length(files) == 0) stop("No mzML files found.")

raw_data <- readMSData(files, mode = "onDisk", msLevel. = 1:2)
sample_classes <- ifelse(grepl("blank", tolower(basename(files))), "blank", "sample")
pData(raw_data)$sample_group <- factor(sample_classes)

# Optional: save raw_data
saveRDS(raw_data, file.path(data_dir, "raw_data.rds"))

# ==== XCMS Processing ====
cwp <- CentWaveParam(
  ppm = ppm,
  peakwidth = c(min_pw, max_pw),
  snthr = snthr,
  prefilter = c(3, 1000),
  mzCenterFun = "wMeanApex3",
  integrate = 1,
  mzdiff = mzdiff,
  noise = noise
)

xdata <- findChromPeaks(raw_data, param = cwp)
xdata <- groupChromPeaks(xdata, param = PeakDensityParam(sampleGroups = sample_classes))
xdata <- adjustRtime(xdata, param = PeakGroupsParam(minFraction = 0.5, smooth = "loess"))
xdata <- groupChromPeaks(xdata, param = PeakDensityParam(sampleGroups = sample_classes))
xdata <- fillChromPeaks(xdata)

saveRDS(xdata, file.path(data_dir, "xdata.rds"))

# ==== CAMERA Annotation ====
xset <- as(xdata, "xcmsSet")
an <- xsAnnotate(xset)
an <- groupFWHM(an)
an <- findIsotopes(an)
an <- groupCorr(an)
an <- findAdducts(an, polarity = polarity_mode)

annotated_peaks <- getPeaklist(an)
annotated_peaks$feature_id <- paste0("ID", seq_len(nrow(annotated_peaks)))

# ==== MS2 Presence Check (±10s, ±0.01 m/z) ====
ms2_idx <- which(msLevel(raw_data) == 2)
ms2_rt  <- rtime(raw_data)[ms2_idx]
ms2_pmz <- precursorMz(raw_data)[ms2_idx]

has_ms2_vec <- vapply(seq_len(nrow(annotated_peaks)), function(i){
  mz <- annotated_peaks$mz[i]
  rt <- annotated_peaks$rt[i]
  rt_ok <- which(ms2_rt >= (rt - 10) & ms2_rt <= (rt + 10))
  if (length(rt_ok) == 0) return(FALSE)
  any(ms2_pmz[rt_ok] >= (mz - 0.01) & ms2_pmz[rt_ok] <= (mz + 0.01))
}, logical(1))

annotated_peaks$has_ms2 <- has_ms2_vec

# ==== Save outputs ====
out_csv <- file.path(data_dir, paste0(output_name, ".csv"))
write.csv(annotated_peaks, file = out_csv, row.names = FALSE)
saveRDS(an, file.path(data_dir, "camera_annotated.rds"))

cat("SUCCESS: Wrote CSV to", out_csv, "\n")
'''

@router.post("/xcms/start")
async def xcms_start(params: XcmsParams):
    data_dir = _ensure_ready(params)

    # Write embedded R to a temp file for this run
    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False) as tf:
        tf.write(R_SCRIPT_EMBEDDED)
        r_path = Path(tf.name)

    cmd = [
        RSCRIPT_BIN, str(r_path),
        f"--dir={str(data_dir)}",
        f"--polarity={params.polarity}",
        f"--output={params.output_name}",
        f"--ppm={params.ppm}",
        f"--snr={params.snr}",
        f"--noise={params.noise}",
        f"--min_peakwidth={params.min_peakwidth}",
        f"--max_peakwidth={params.max_peakwidth}",
        f"--mzdiff={params.mzdiff}",
    ]

    try:
        proc = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT,
        )
    except FileNotFoundError:
        # cleanup temp file
        try: r_path.unlink(missing_ok=True)
        except Exception: pass
        raise HTTPException(status_code=500, detail=f"Rscript binary not found: {RSCRIPT_BIN}")

    try:
        if params.timeout_sec and params.timeout_sec > 0:
            try:
                stdout, _ = await asyncio.wait_for(proc.communicate(), timeout=params.timeout_sec)
            except asyncio.TimeoutError:
                proc.kill()
                try:
                    await proc.wait()
                finally:
                    raise HTTPException(status_code=504, detail=f"xcms run exceeded {params.timeout_sec}s and was terminated.")
        else:
            stdout, _ = await proc.communicate()
    finally:
        # best-effort cleanup of the temp R file
        try: r_path.unlink(missing_ok=True)
        except Exception: pass

    log_text = (stdout or b"").decode(errors="ignore")

    if proc.returncode != 0:
        tail = log_text[-2000:] if log_text else "No logs captured."
        raise HTTPException(status_code=500, detail=f"xcms failed (exit {proc.returncode}).\n\n{tail}")

    return {"ok": True, "output_name": params.output_name}
