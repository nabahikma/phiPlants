from fastapi import APIRouter, Form
import subprocess
from starlette.responses import RedirectResponse, StreamingResponse
from pathlib import Path
from typing import List

router = APIRouter()


def collect_d_folders(parent: Path) -> List[Path]:

    if not parent.exists() or not parent.is_dir():
        return []
    return [p for p in parent.iterdir() if p.is_dir() and p.suffix.lower() == ".d"]



@router.post("/run-msconvert/")
async def run_msconvert(
        polarity: str = Form(...),
        scan_time: str = Form(...),
        parent_folder: str = Form(...),
        mz_window: str = Form(...)
):
    parent_path = Path(parent_folder).expanduser().resolve()
    d_folders = collect_d_folders(parent_path)
    if not d_folders:
        return RedirectResponse(url="/processing", status_code=303)
    out_dir = parent_folder

    exe = r'"C:\\pwiz\\msconvert.exe"'
    filters = " ".join([
        r'--filter "peakPicking vendor snr=10 msLevel=1-2"',
        r'--ext .mzML',
        r'--zlib',
        r'-v',
        f'-o "{str(out_dir)}"',
        f'--filter "polarity {polarity}"',
        r'--filter "msLevel 1-2"',
        f'--filter "scanTime [{scan_time}]"',
        r'--filter "activation CID"',
        r'--filter "analyzer TOF"',
        f'--filter "mzWindow [{mz_window}]"',
        r'--filter "MS2Denoise 6 0 true"',
        r'--filter "titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState>"',
    ])
    # Quote every .d path and append
    inputs = " ".join([f'"{str(p)}"' for p in d_folders])
    command = f"{exe} {filters} {inputs}"
    # Debug prints (server logs)
    print("Parent path:", str(parent_path))
    print("Found .d folders:", [str(p) for p in d_folders])
    print("Output folder:", str(out_dir))
    print("Command:", command)

    # Stream msconvert output
    def stream_logs():
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            shell=True
        )
        for line in iter(process.stdout.readline, ''):
            yield line
        process.stdout.close()
        process.wait()
        code = process.returncode  # ✅ Get actual exit code
        yield f"\n__MSCONVERT_DONE__ EXIT_CODE={code}\n"

    return StreamingResponse(stream_logs(), media_type="text/plain")