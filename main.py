import uvicorn
from fastapi import FastAPI, Form
from starlette.requests import Request
from starlette.responses import HTMLResponse, RedirectResponse, StreamingResponse
from starlette.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles
import subprocess
from pathlib import Path
from routers import msconvert


app = FastAPI()

app.include_router(msconvert.router)
templates = Jinja2Templates(directory="templates")
app.mount("/static", StaticFiles(directory="static"), name="static")


@app.get("/", response_class=HTMLResponse)
async def mainpage(request: Request):
    return templates.TemplateResponse("mainpage.html", {"request": request})


@app.get("/firstStep_msConvert", response_class=HTMLResponse)
async def firststep_msconvert(request: Request):
    return templates.TemplateResponse("firstStep_msConvert.html", {"request": request})


@app.get("/processing", response_class=HTMLResponse)
async def processing_page():
    # Serve the processing page HTML
    html_path = Path("templates/processing.html")
    return HTMLResponse(content=html_path.read_text(), status_code=200)


@app.get("/Learn", response_class=HTMLResponse)
async def learn(request: Request):
    return templates.TemplateResponse("LearnPage.html", {"request": request})


if __name__ == '__main__':
    uvicorn.run(app, host="127.0.0.1", port=8000)