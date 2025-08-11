import uvicorn
from fastapi import FastAPI, Request, HTTPException
from starlette.responses import HTMLResponse
from starlette.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
import asyncio
import shlex
from routers import msconvert, xcms
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
async def processing_page(request: Request):
    return templates.TemplateResponse("processing.html", {"request": request})


# ⬇️ Processing_xcms takes directory & polarity as query params and shows them
@app.get("/Processing_xcms", response_class=HTMLResponse)
async def processing_xcms(request: Request, directory: str, polarity: str):
    # you can sanitize/validate here if you need to
    ctx = {
        "request": request,
        "directory": directory,
        "polarity": polarity
    }
    return templates.TemplateResponse("Processing_xcms.html", ctx)






@app.get("/annotation_part", response_class=HTMLResponse)
async def annotation_part(request: Request, directory: str, polarity: str, table: str):
    # Show a simple landing page for annotation step (customize as you like)
    ctx = {"request": request, "directory": directory, "polarity": polarity, "table": table}
    return templates.TemplateResponse("annotation_part.html", ctx)


if __name__ == '__main__':
    uvicorn.run(app, host="127.0.0.1", port=8000)
