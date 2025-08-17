# routers/xcms.py
from fastapi import APIRouter, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field, field_validator, model_validator
from typing import Optional, Literal
from pathlib import Path
import asyncio
import tempfile

router = APIRouter()

# Set to absolute Rscript.exe if needed on Windows
RSCRIPT_BIN = "Rscript"


class XcmsParams(BaseModel):
    output_name: str = Field(min_length=1, max_length=128)
    directory: str = Field(min_length=1)
    polarity: Literal["positive", "negative"]
    ppm: float = Field(gt=0)
    snr: float = Field(gt=0)
    noise: float = Field(ge=0)
    min_peakwidth: float = Field(2.0, gt=0)
    max_peakwidth: float = Field(30.0, gt=0)
    mzdiff: float = Field(0.01)
    timeout_sec: Optional[int] = Field(0, ge=0, description="0 = no timeout")

    @field_validator("output_name")
    @classmethod
    def safe_output(cls, v: str) -> str:
        import re
        if not re.match(r"^[A-Za-z0-9._-]+$", v):
            raise ValueError("output_name allows letters, numbers, ., _, - only")
        return v

    @model_validator(mode="after")
    def check_peakwidths(self):
        if self.max_peakwidth <= self.min_peakwidth:
            raise ValueError("max_peakwidth must be greater than min_peakwidth")
        return self


def _ensure_ready(params: XcmsParams) -> Path:
    d = Path(params.directory).expanduser().resolve()
    if not d.exists() or not d.is_dir():
        raise HTTPException(status_code=400, detail=f"Directory not found: {d}")
    # quick presence check (case-insensitive)
    if not any(p.suffix.lower() == ".mzml" for p in d.glob("*")):
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
ppm           <- as.numeric(get_arg("ppm", 15))
snthr         <- as.numeric(get_arg("snr", 6))
noise         <- as.numeric(get_arg("noise", 0))
min_pw        <- as.numeric(get_arg("min_peakwidth", 2))
max_pw        <- as.numeric(get_arg("max_peakwidth", 30))
mzdiff        <- as.numeric(get_arg("mzdiff", 0.01))

cat("PARAMS:\n",
    " dir=", data_dir, "\n",
    " polarity=", polarity_mode, "\n",
    " ppm=", ppm, "\n",
    " snthr=", snthr, "\n",
    " noise=", noise, "\n",
    " peakwidth=[", min_pw, ",", max_pw, "]\n",
    " mzdiff=", mzdiff, "\n", sep="")

if (is.null(data_dir) || !dir.exists(data_dir)) {
  stop(paste("Invalid data_dir:", data_dir))
}

# ==== Load mzML (MS1 + MS2) ====
files <- list.files(data_dir, pattern="mzML$", full.names=TRUE, ignore.case=TRUE)
cat("Found", length(files), "mzML file(s)\n")
if (length(files) == 0) stop("No mzML files found.")

raw_data <- MSnbase::readMSData(files, mode = "onDisk", msLevel. = 1:2)

ms1_count <- sum(MSnbase::msLevel(raw_data) == 1)
ms2_count <- sum(MSnbase::msLevel(raw_data) == 2)
cat("Spectra counts: MS1=", ms1_count, " MS2=", ms2_count, "\n", sep="")
if (ms1_count == 0) {
  stop("No MS1 spectra found. Check msconvert filters (remove activation/analyzer/mzWindow/scanTime).")
}

sample_classes <- ifelse(grepl("blank", tolower(basename(files))), "blank", "sample")
pData(raw_data)$sample_group <- factor(sample_classes)

# Optional: save raw_data
# saveRDS(raw_data, file.path(data_dir, "raw_data.rds"))

# ==== XCMS Processing ====
cwp <- xcms::CentWaveParam(
  ppm = ppm,
  peakwidth = c(min_pw, max_pw),
  snthr = snthr,
  prefilter = c(3, 300),            # slightly relaxed for robustness
  mzCenterFun = "wMeanApex3",
  integrate = 1,
  mzdiff = -0.001,
  noise = noise
)

cat("Detecting mass traces / peaks ...\n")
xdata <- xcms::findChromPeaks(raw_data, param = cwp)
cat("Chrom peaks found per sample:", paste(xcms::chromPeaks(xdata, bySample=TRUE) |> lengths(), collapse=", "), "\n")

xdata <- xcms::groupChromPeaks(xdata, param = xcms::PeakDensityParam(sampleGroups = sample_classes))
xdata <- xcms::adjustRtime(xdata, param = xcms::PeakGroupsParam(minFraction = 0.5, smooth = "loess"))
xdata <- xcms::groupChromPeaks(xdata, param = xcms::PeakDensityParam(sampleGroups = sample_classes))
xdata <- xcms::fillChromPeaks(xdata)

# saveRDS(xdata, file.path(data_dir, "xdata.rds"))

# ==== CAMERA Annotation ====
xset <- as(xdata, "xcmsSet")
an <- CAMERA::xsAnnotate(xset)
an <- CAMERA::groupFWHM(an)
an <- CAMERA::findIsotopes(an)
an <- CAMERA::groupCorr(an)
an <- CAMERA::findAdducts(an, polarity = polarity_mode)

annotated_peaks <- CAMERA::getPeaklist(an)
annotated_peaks$feature_id <- paste0("ID", seq_len(nrow(annotated_peaks)))

# ==== MS2 Presence Check (±10s, ±0.01 m/z) ====
ms2_idx <- which(MSnbase::msLevel(raw_data) == 2)
ms2_rt  <- MSnbase::rtime(raw_data)[ms2_idx]
ms2_pmz <- MSnbase::precursorMz(raw_data)[ms2_idx]

has_ms2_vec <- vapply(seq_len(nrow(annotated_peaks)), function(i){
  mz <- annotated_peaks$mz[i]
  rt <- annotated_peaks$rt[i]
  rt_ok <- which(ms2_rt >= (rt - 10) & ms2_rt <= (rt + 10))
  if (length(rt_ok) == 0) return(FALSE)
  any(ms2_pmz[rt_ok] >= (mz - 0.01) & ms2_pmz[rt_ok] <= (mz + 0.01))
}, logical(1))

annotated_peaks$has_ms2 <- has_ms2_vec

# ==== Save outputs ====
out_csv <- file.path(data_dir, paste0(Features_Matrix, ".csv"))
utils::write.csv(annotated_peaks, file = out_csv, row.names = FALSE)

cat("SUCCESS: Wrote CSV to ", out_csv, "\n", sep="")
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
        f"--ppm={params.ppm}",
        f"--snr={params.snr}",
        f"--noise={params.noise}",
        f"--min_peakwidth={params.min_peakwidth}",
        f"--max_peakwidth={params.max_peakwidth}",
    ]

    async def stream():
        try:
            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.STDOUT,
            )
        except FileNotFoundError:
            yield f"❌ Rscript not found at {RSCRIPT_BIN}\n"
            # cleanup
            try: r_path.unlink(missing_ok=True)
            except Exception: pass
            return

        # Optional timeout handling while streaming
        done = False
        async def reader_task():
            async for chunk in proc.stdout:
                yield chunk.decode(errors="ignore")

        # We stream directly without wait_for so logs arrive as they are produced
        async for line in proc.stdout:
            yield line.decode(errors="ignore")

        code = await proc.wait()
        yield f"\n__XCMS_DONE__ EXIT_CODE={code}\n"

        # cleanup the temp R file
        try: r_path.unlink(missing_ok=True)
        except Exception: pass

    return StreamingResponse(stream(), media_type="text/plain; charset=utf-8")
