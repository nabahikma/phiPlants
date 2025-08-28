from fastapi import APIRouter, Form
from starlette.responses import StreamingResponse
from pathlib import Path
from typing import List
import subprocess
import os
import shutil

router = APIRouter()


def collect_d_folders(parent: Path) -> List[Path]:
    if not parent.exists() or not parent.is_dir():
        return []
    return [p for p in parent.iterdir() if p.is_dir() and p.suffix.lower() == ".d"]


def _find_msconvert_exe() -> str:
    # 1) Env override (portable/pwiz/msconvert.exe in run_app.py)
    exe = os.getenv("MSCONVERT_EXE")
    if exe and Path(exe).exists():
        return exe
    # 2) PATH fallback
    for name in ("msconvert.exe", "msconvert"):
        hit = shutil.which(name)
        if hit:
            return hit
    return ""


@router.post("/run-msconvert/")
async def run_msconvert(
    polarity: str = Form(...),
    scan_time: str = Form(...),   # e.g. "60,600" OR "[60,600]" (we normalize below)
    parent_folder: str = Form(...),
    mz_window: str = Form(...),   # e.g. "100,500" OR "[100,500]"
):
    parent_path = Path(parent_folder).expanduser().resolve()
    d_folders = collect_d_folders(parent_path)

    if not d_folders:
        def err_stream():
            yield f"❌ No .d folders found under: {parent_path}\n"
            yield "__MSCONVERT_DONE__ EXIT_CODE=1\n"
        return StreamingResponse(err_stream(), media_type="text/plain; charset=utf-8")

    msconvert_exe = "portable/pwiz/msconvert.exe"
    if not msconvert_exe:
        def err_stream2():
            yield "❌ msconvert.exe not found.\n"
            yield "   Set MSCONVERT_EXE env var or ensure ProteoWizard is on PATH.\n"
            yield "__MSCONVERT_DONE__ EXIT_CODE=1\n"
        return StreamingResponse(err_stream2(), media_type="text/plain; charset=utf-8")

    # Normalize bracketed ranges for filters
    st = scan_time.strip()
    if not (st.startswith("[") and st.endswith("]")):
        st = f"[{st}]"
    mw = mz_window.strip()
    if not (mw.startswith("[") and mw.endswith("]")):
        mw = f"[{mw}]"

    out_dir = str(parent_path)

    # Build arg list (no shell, robust quoting)
    args = [
        msconvert_exe,
        "--filter", "peakPicking vendor snr=10 msLevel=1-2",
        "--ext", ".mzML",
        "--zlib",
        "-v",
        "-o", out_dir,
        "--filter", f"polarity {polarity}",
        "--filter", "msLevel 1-2",
        "--filter", f"scanTime {st}",
        "--filter", "activation CID",
        "--filter", "analyzer TOF",
        "--filter", f"mzWindow {mw}",
        "--filter", "MS2Denoise 6 0 true",
        "--filter", "titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState>",
    ] + [str(p) for p in d_folders]

    def stream_logs():
        # Echo command for debugging
        yield "Command: " + " ".join(f'"{a}"' if " " in a else a for a in args) + "\n"
        try:
            process = subprocess.Popen(
                args,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                shell=False,
            )
        except Exception as e:
            yield f"❌ Failed to start msconvert: {e}\n"
            yield "__MSCONVERT_DONE__ EXIT_CODE=1\n"
            return

        for line in iter(process.stdout.readline, ''):
            yield line
        process.stdout.close()
        process.wait()
        code = process.returncode
        yield f"\n__MSCONVERT_DONE__ EXIT_CODE={code}\n"

    return StreamingResponse(stream_logs(), media_type="text/plain; charset=utf-8")
