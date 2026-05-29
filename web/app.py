import asyncio
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor
from contextlib import asynccontextmanager
from io import BytesIO

import pandas as pd
from fastapi import FastAPI, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse

from core import SST_DESCRIPTIONS, classify_sst, read_mcif

_MSG_CSV = os.path.join(os.path.dirname(__file__), "msg_list.cvs")
_df = None

_executor = ThreadPoolExecutor(max_workers=1)


@asynccontextmanager
async def lifespan(app: FastAPI):
    global _df
    _df = pd.read_csv(_MSG_CSV, dtype=str).set_index("UNI_NUM")
    yield
    _executor.shutdown(wait=False)


app = FastAPI(lifespan=lifespan)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/analyze", response_class=PlainTextResponse)
async def analyze(file: UploadFile):
    suffix = os.path.splitext(file.filename or "")[1].lower()
    if suffix not in (".mcif", ".cif"):
        raise HTTPException(status_code=422, detail="File must be .mcif or .cif")

    content = await file.read()

    loop = asyncio.get_running_loop()
    try:
        result = await asyncio.wait_for(
            loop.run_in_executor(_executor, _run, content, file.filename),
            timeout=180.0,
        )
    except asyncio.TimeoutError:
        raise HTTPException(status_code=500, detail="Analysis timed out (>3 min)")
    except (Exception, SystemExit) as exc:
        msg = str(exc)
        if not msg or msg == "1":
            msg = "Analysis failed (check that the file contains magnetic atoms and valid magnetic moments)"
        raise HTTPException(status_code=500, detail=msg)

    return result


def _bns_symbol(msg) -> str:
    try:
        row = _df.loc[_df["BNS_NUM"] == msg.bns_number].iloc[0]
        return row["BNS_SYM"]
    except (KeyError, IndexError):
        return msg.bns_number


def _run(content: bytes, filename: str) -> str:
    mcif = BytesIO(content)

    sst_key, sog, msg_wo_soc, msg_w_soc, centrosym, theta_i, compensated = classify_sst(mcif)

    bns_wo = _bns_symbol(msg_wo_soc)
    bns_w = _bns_symbol(msg_w_soc)
    sst_desc = SST_DESCRIPTIONS.get(sst_key, "unknown")

    lines = [
        "FINDMAGSYM Analysis Report",
        "=" * 42,
        f"Input file: {filename}",
        "",
        "MSG without SOC",
        f"  BNS symbol : {bns_wo}",
        f"  BNS number : {msg_wo_soc.bns_number}",
        f"  MSG type   : {msg_wo_soc.type}",
        "",
        "MSG with SOC",
        f"  BNS symbol : {bns_w}",
        f"  BNS number : {msg_w_soc.bns_number}",
        f"  MSG type   : {msg_w_soc.type}",
        "",
        "Spin-only group",
        f"  Type       : {sog.spin_only_group_type}",
        "",
        "Classification",
        f"  SST class    : {sst_key}",
        f"  Description  : {sst_desc}",
        f"  Centrosymmetric : {'yes' if centrosym else 'no'}",
        f"  Has ΘI         : {'yes' if theta_i else 'no'}",
        f"  Compensated     : {'yes' if compensated else 'no'}",
    ]
    return "\n".join(lines) + "\n"
