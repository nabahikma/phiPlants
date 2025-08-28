# run_app.py
import os, sys, threading, time, webbrowser, shutil
from pathlib import Path
import uvicorn

# Resolve base folder even when running as onefile EXE
BASE = Path(getattr(sys, "_MEIPASS", Path(__file__).parent)).resolve()
os.chdir(BASE)

# Portable tool paths
PORTABLE_R = BASE / "portable" / "R"
PORTABLE_RSCRIPT = PORTABLE_R / "bin" / "Rscript.exe"
PORTABLE_R_LIB = PORTABLE_R / "library"

PORTABLE_PWIZ = BASE / "portable" / "pwiz"
PORTABLE_MSCONVERT = PORTABLE_PWIZ / "msconvert.exe"

# Add portable bins to PATH so subprocess can find them (fallbacks if env not set elsewhere)
os.environ["PATH"] = str(PORTABLE_PWIZ) + os.pathsep + str(PORTABLE_R / "bin") + os.pathsep + os.environ.get("PATH", "")

# Environment defaults used by your routers
os.environ.setdefault("PLANTDB_CSV", str((BASE / "data" / "PlantDB_ESI.csv").resolve()))
os.environ.setdefault("RSCRIPT_EXE", str(PORTABLE_RSCRIPT))       # FastAPI will call this
os.environ.setdefault("R_HOME", str(PORTABLE_R))
os.environ.setdefault("R_USER_LIBS", str(PORTABLE_R_LIB))         # where packages live
os.environ.setdefault("MSCONVERT_EXE", str(PORTABLE_MSCONVERT))   # msconvert path for routers/msconvert.py
os.environ.setdefault("PORT", "8000")

# Sanity logs (optional)
def _log_env():
    print("== Portable tool check ==")
    print("Rscript:", os.environ["RSCRIPT_EXE"], "- exists?", (Path(os.environ["RSCRIPT_EXE"]).exists()))
    print("msconvert:", os.environ["MSCONVERT_EXE"], "- exists?", (Path(os.environ["MSCONVERT_EXE"]).exists()))
    print("R libs:", os.environ["R_USER_LIBS"], "- exists?", (Path(os.environ["R_USER_LIBS"]).exists()))
_log_env()

# Import your FastAPI app
from main import app  # noqa

def _open_browser(url: str):
    def run():
        time.sleep(1.5)
        try: webbrowser.open(url)
        except: pass
    threading.Thread(target=run, daemon=True).start()

if __name__ == "__main__":
    url = f"http://127.0.0.1:{os.environ['PORT']}"
    _open_browser(url)
    uvicorn.run(app, host="127.0.0.1", port=int(os.environ["PORT"]))
import uvicorn
from fastapi import FastAPI, Request, Query
from starlette.responses import HTMLResponse
from starlette.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles
from typing import Literal
from fastapi.responses import PlainTextResponse
import csv
from pathlib import Path
# Routers
from routers import msconvert, xcms, annotate

app = FastAPI()

# Include routers
app.include_router(msconvert.router)
app.include_router(xcms.router)
app.include_router(annotate.router)

# Templates & static
templates = Jinja2Templates(directory="templates")
app.mount("/static", StaticFiles(directory="static"), name="static")


@app.get("/", response_class=HTMLResponse)
async def mainpage(request: Request):
    return templates.TemplateResponse("mainpage.html", {"request": request})


@app.get("/firstStep_msConvert", response_class=HTMLResponse)
async def firststep_msconvert(request: Request):
    return templates.TemplateResponse("firstStep_msConvert.html", {"request": request})


@app.get("/processing", response_class=HTMLResponse)
async def processing_page(request: Request):
    return templates.TemplateResponse("processing.html", {"request": request})


# Require directory & polarity to be passed in the URL
@app.get("/Processing_xcms", response_class=HTMLResponse)
async def processing_xcms(
    request: Request,
    directory: str = Query(..., description="Absolute path to mzML folder"),
    polarity: Literal["positive", "negative"] = Query(...)
):
    return templates.TemplateResponse(
        "Processing_xcms.html",
        {"request": request, "directory": directory, "polarity": polarity},
    )


@app.get("/annotation_part", response_class=HTMLResponse)
async def annotation_part(request: Request, directory: str, polarity: str, table: str):
    ctx = {"request": request, "directory": directory, "polarity": polarity, "table": table}
    return templates.TemplateResponse("annotation_part.html", ctx)


@app.get("/Learn", response_class=HTMLResponse)
async def Learn_Page(request: Request):
    return templates.TemplateResponse("LearnPage.html", {"request": request})


@app.get("/static-proxy", response_class=PlainTextResponse)
async def static_proxy(directory: str, file: str):
    p = Path(directory).expanduser().resolve() / file
    if not p.exists():
        return PlainTextResponse("Not found", status_code=404)
    return PlainTextResponse(p.read_text(encoding="utf-8"), media_type="text/csv; charset=utf-8")


@app.get("/ShowTable", response_class=HTMLResponse)
async def show_table(
    request: Request,
    directory: str = Query(..., description="Folder containing the CSV"),
    file: str = Query(..., description="CSV filename inside that folder")
):
    return templates.TemplateResponse(
        "ShowTable.html",
        {"request": request, "directory": directory, "file": file}
    )


@app.get("/OpenAnalyzedProject")
async def OpenAnalyzedProject(request: Request):
    return templates.TemplateResponse("OpenAnalyzedProject.html", {"request": request})


if __name__ == '__main__':
    uvicorn.run(app, host="127.0.0.1", port=8000)
