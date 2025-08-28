import os
import sys
import socket
import threading
import webbrowser
from pathlib import Path
import uvicorn
from typing import Literal

from fastapi import FastAPI, Request, Query, HTTPException
from starlette.responses import HTMLResponse, PlainTextResponse
from starlette.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles

# Routers
from routers import msconvert, xcms, annotate


# ---------- Paths that work in source & bundled EXE ----------

def resource_path(rel: str) -> str:
    """
    Absolute path to a bundled resource (PyInstaller-friendly).
    Falls back to the project folder when running from source.
    """
    base = getattr(sys, "_MEIPASS", Path(__file__).parent.resolve())
    return str(Path(base) / rel)

def app_root() -> Path:
    """
    A writable folder next to the EXE (or the project folder in dev).
    Use this for outputs/logs, not resource_path() which may be read-only.
    """
    return Path(sys.executable).parent if getattr(sys, "frozen", False) else Path(__file__).parent.resolve()


# ---------- Important directories ----------

TEMPLATES_DIR = resource_path("templates")
STATIC_DIR    = resource_path("static")
DATA_DIR      = resource_path("data")   # bundle your read-only CSVs here
BIN_DIR       = resource_path("bin")    # bundle helper EXEs here
RUNTIME_DIR   = str(app_root() / "runtime_data")  # writable outputs here
Path(RUNTIME_DIR).mkdir(exist_ok=True)

# Make paths discoverable by routers / subprocess code
os.environ.setdefault("PHIAPP_DATA_DIR", DATA_DIR)
os.environ.setdefault("PHIAPP_BIN_DIR", BIN_DIR)
os.environ.setdefault("PHIAPP_RUNTIME_DIR", RUNTIME_DIR)
# Ensure bundled executables are on PATH for subprocess calls
os.environ["PATH"] = BIN_DIR + os.pathsep + os.environ.get("PATH", "")


# ---------- FastAPI app ----------

app = FastAPI()

# Include routers
app.include_router(msconvert.router)
app.include_router(xcms.router)
app.include_router(annotate.router)

# Templates & static
templates = Jinja2Templates(directory=TEMPLATES_DIR)
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


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


# Serve a CSV from ANY absolute path (use carefully)
@app.get("/static-proxy", response_class=PlainTextResponse)
async def static_proxy(directory: str, file: str):
    p = Path(directory).expanduser().resolve() / file
    if not p.exists():
        return PlainTextResponse("Not found", status_code=404)
    return PlainTextResponse(p.read_text(encoding="utf-8"), media_type="text/csv; charset=utf-8")


# Serve a CSV that you BUNDLE inside ./data (safe: prevents path traversal)
@app.get("/static-data", response_class=PlainTextResponse)
async def static_data(file: str = Query(..., description="Relative path inside bundled data/")):
    base = Path(DATA_DIR).resolve()
    candidate = (base / file).resolve()
    if base not in candidate.parents and candidate != base:
        raise HTTPException(status_code=400, detail="Invalid path")
    if not candidate.exists():
        raise HTTPException(status_code=404, detail="Not found")
    return PlainTextResponse(candidate.read_text(encoding="utf-8"), media_type="text/csv; charset=utf-8")


# Simple health check
@app.get("/health", response_class=PlainTextResponse)
async def health():
    return PlainTextResponse("ok")


# ---------- Nice UX: open browser + pick free port ----------

def _port_is_free(port: int) -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(("127.0.0.1", port)) != 0

def _open_browser_later(url: str, delay: float = 1.5) -> None:
    def _open():
        try:
            webbrowser.open(url)
        except Exception:
            pass
    threading.Timer(delay, _open).start()


if __name__ == '__main__':
    port = int(os.environ.get("PORT", "8000"))
    if not _port_is_free(port):
        for candidate in range(port + 1, port + 11):
            if _port_is_free(candidate):
                port = candidate
                break
    url = f"http://127.0.0.1:{port}"
    _open_browser_later(url, delay=2.0)
    uvicorn.run(app, host="127.0.0.1", port=port)
