"""HTTP API for Protein Atlas.

Run locally with:
    uvicorn api:app --reload

This file only wires HTTP routes to two framework-free layers:

  * atlas_engine.py   - the one validated sequence-analysis engine
  * atlas_services.py - UniProt paging, structures, PROSITE, dataset statistics

Nothing scientific is calculated here.  Errors are translated into explicit,
human-readable states; a raw traceback never reaches the browser.
"""
from __future__ import annotations

import logging
import os
import threading
from collections import OrderedDict
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any, Optional

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

import atlas_aminoacids as aminoacids
import atlas_engine as engine
import atlas_services as services
import protein_analyzer as pa

API_VERSION = "v1"
APP_VERSION = "3.0.0"
log = logging.getLogger("atlas")

http_client = services.HttpClient()
uniprot = services.UniProtClient(http_client)
structures = services.StructureService(http_client)
references = services.ReferenceStructures(http_client)
prosite = services.PrositeManager(http_client)
dataset_stats = services.DatasetStatsManager(http_client)

_analysis_cache: "OrderedDict[tuple, dict[str, Any]]" = OrderedDict()
_analysis_lock = threading.Lock()


# --------------------------------------------------------------------------
# request models
# --------------------------------------------------------------------------
class SequenceBody(BaseModel):
    sequence: str = Field(..., max_length=2_000_000)


class AnalyzeBody(SequenceBody):
    hydropathy_window: int = pa.DEFAULT_HYDROPATHY_WINDOW


class ScanBody(SequenceBody):
    exclude_frequent: bool = False


class PatternBody(SequenceBody):
    pattern: str = Field(default="", max_length=500)


class CompareBody(BaseModel):
    first: str = Field(default="", max_length=40)
    second: str = Field(default="", max_length=40)
    first_sequence: Optional[str] = Field(default=None, max_length=2_000_000)
    second_sequence: Optional[str] = Field(default=None, max_length=2_000_000)
    kmer_size: int = pa.KMER_SIZE


# --------------------------------------------------------------------------
# application
# --------------------------------------------------------------------------
@asynccontextmanager
async def lifespan(application: FastAPI):
    try:
        prosite.start()                               # background: the server is usable at once
        if os.getenv("ATLAS_PREBUILD_STATS", "1") != "0":
            dataset_stats.status("human_reviewed")    # starts the one-off background build if needed
    except Exception:                                 # noqa: BLE001 - startup must never fail on this
        log.exception("background start-up tasks could not be started")
    yield


app = FastAPI(
    title="Protein Atlas API", version=APP_VERSION, lifespan=lifespan,
    description="Backend for Protein Atlas: UniProt, PROSITE and structure services "
                "plus one validated sequence-analysis engine.",
)
app.add_middleware(
    CORSMiddleware, allow_origins=os.getenv("CORS_ALLOW_ORIGINS", "*").split(","),
    allow_methods=["*"], allow_headers=["*"],
)


@app.middleware("http")
async def keep_frontend_fresh(request: Request, call_next):
    response = await call_next(request)
    if not request.url.path.startswith("/api/"):
        response.headers["Cache-Control"] = "no-cache"   # never serve a stale page/script after an update
    return response


def _problem(status: int, detail: str, kind: str, **extra: Any) -> JSONResponse:
    return JSONResponse(status_code=status, content={"detail": detail, "kind": kind, **extra})


@app.exception_handler(engine.InvalidSequence)
async def _invalid_input(request: Request, error: engine.InvalidSequence):
    return _problem(422, str(error), "invalid_input")


@app.exception_handler(services.UpstreamError)
async def _upstream(request: Request, error: services.UpstreamError):
    if error.kind == "not_found":
        return _problem(404, str(error), "not_found")
    if error.kind == "bad_request":
        return _problem(422, str(error), "invalid_input")
    return _problem(503, "An external database is unavailable right now. " + str(error),
                    "database_unavailable")


@app.exception_handler(services.PageOutOfRange)
async def _page_range(request: Request, error: services.PageOutOfRange):
    return _problem(404, str(error), "page_out_of_range", total_pages=error.total_pages)


@app.exception_handler(services.PageUnreachable)
async def _page_unreachable(request: Request, error: services.PageUnreachable):
    return _problem(409, str(error), "page_unreachable")


@app.exception_handler(RequestValidationError)
async def _bad_request(request: Request, error: RequestValidationError):
    parts = []
    for item in error.errors():
        where = ".".join(str(p) for p in item.get("loc", []) if p not in ("body", "query", "path"))
        parts.append(f"{where}: {item.get('msg', 'invalid value')}" if where else str(item.get("msg", "invalid")))
    return _problem(422, "; ".join(parts) or "The request could not be understood.", "invalid_input")


@app.exception_handler(Exception)
async def _unexpected(request: Request, error: Exception):
    log.exception("unhandled error on %s", request.url.path)
    return _problem(500, engine_calc_error(), "calculation_error")


def engine_calc_error() -> str:
    return "Analysis unavailable \u2014 calculation error"


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def protein_summary(record: dict[str, Any]) -> dict[str, Any]:
    return {k: record.get(k) for k in
            ("id", "entry_name", "name", "gene", "organism", "length", "reviewed")}


def analysis_for(record: dict[str, Any], window: int) -> dict[str, Any]:
    key = (record["id"], window, len(record["sequence"]))
    with _analysis_lock:
        if key in _analysis_cache:
            _analysis_cache.move_to_end(key)
            return _analysis_cache[key]
    result = engine.analyze_sequence(record["sequence"], window)
    with _analysis_lock:
        _analysis_cache[key] = result
        while len(_analysis_cache) > 24:
            _analysis_cache.popitem(last=False)
    return result


def reference_row(code: str) -> dict[str, Any]:
    code = code.upper()
    base = engine.AMINO_ACID_BY_CODE.get(code)
    if base is None:
        raise services.UpstreamError("not_found", f"'{code}' is not one of the 20 standard amino acids.")
    return {**base, **aminoacids.extras(code)}


# --------------------------------------------------------------------------
# status
# --------------------------------------------------------------------------
@app.get("/health")
def health() -> dict[str, Any]:
    info = prosite.info()
    return {"status": "ok", "version": APP_VERSION, "api_version": API_VERSION,
            "prosite": {"state": info.get("state"), "patterns": info.get("patterns_usable")}}


@app.get(f"/api/{API_VERSION}/status")
def status_report() -> dict[str, Any]:
    return {
        "version": APP_VERSION, "engine_version": engine.ENGINE_VERSION,
        "prosite": prosite.info(),
        "default_query": services.DEFAULT_QUERY,
        "page_sizes": list(services.PAGE_SIZES),
        "kmer_limit": engine.KMER_LIMIT,
        "dataset_scopes": {k: v["label"] for k, v in services.SCOPES.items()},
    }


@app.get(f"/api/{API_VERSION}/prosite/status")
def prosite_status() -> dict[str, Any]:
    return prosite.info()


@app.post(f"/api/{API_VERSION}/prosite/retry")
def prosite_retry() -> dict[str, Any]:
    prosite.retry()
    return prosite.info()


# --------------------------------------------------------------------------
# proteins
# --------------------------------------------------------------------------
@app.get(f"/api/{API_VERSION}/proteins")
def list_proteins(
    query: Optional[str] = Query(default=None, max_length=300),
    page: int = Query(default=1, ge=1),
    page_size: int = Query(default=24),
) -> dict[str, Any]:
    if page_size not in services.PAGE_SIZES:
        raise HTTPException(status_code=422, detail="page_size must be one of "
                            + ", ".join(str(n) for n in services.PAGE_SIZES))
    return uniprot.page(query or "", page, page_size)


@app.get(f"/api/{API_VERSION}/proteins/{{protein_id}}")
def get_protein(protein_id: str) -> dict[str, Any]:
    return uniprot.get_protein(protein_id)


@app.get(f"/api/{API_VERSION}/proteins/{{protein_id}}/analysis")
def protein_analysis(protein_id: str,
                     hydropathy_window: int = Query(default=pa.DEFAULT_HYDROPATHY_WINDOW)) -> dict[str, Any]:
    record = uniprot.get_protein(protein_id)
    return {"protein": protein_summary(record), "analysis": analysis_for(record, hydropathy_window)}


@app.get(f"/api/{API_VERSION}/proteins/{{protein_id}}/prosite")
def protein_prosite(protein_id: str, exclude_frequent: bool = False) -> dict[str, Any]:
    record = uniprot.get_protein(protein_id)
    return prosite.scan(record["sequence"], exclude_frequent)


@app.get(f"/api/{API_VERSION}/proteins/{{protein_id}}/structures")
def protein_structures(protein_id: str) -> dict[str, Any]:
    return structures.list_for(uniprot.get_protein(protein_id))


@app.get(f"/api/{API_VERSION}/structures/{{source}}/{{ident}}/file")
def structure_file(source: str, ident: str) -> dict[str, Any]:
    text, fmt, meta = structures.fetch_file(source, ident, protein_loader=uniprot.get_protein)
    return {"format": fmt, "data": text, "meta": meta}


# --------------------------------------------------------------------------
# raw-sequence tools (Sequence Lab / Pattern Lab)
# --------------------------------------------------------------------------
@app.post(f"/api/{API_VERSION}/analyze")
def analyze_sequence(body: AnalyzeBody) -> dict[str, Any]:
    return engine.analyze_sequence(body.sequence, body.hydropathy_window)


@app.post(f"/api/{API_VERSION}/prosite/scan")
def prosite_scan(body: ScanBody) -> dict[str, Any]:
    return prosite.scan(body.sequence, body.exclude_frequent)


@app.post(f"/api/{API_VERSION}/pattern")
def pattern_search(body: PatternBody) -> dict[str, Any]:
    return engine.search_pattern(body.sequence, body.pattern)


@app.post(f"/api/{API_VERSION}/compare")
def compare(body: CompareBody) -> dict[str, Any]:
    sides = []
    for accession, raw in ((body.first, body.first_sequence), (body.second, body.second_sequence)):
        if raw and raw.strip():
            sides.append(({"id": "pasted", "name": "Pasted sequence", "length": len(engine.normalize_sequence(raw))},
                          raw))
        elif accession.strip():
            record = uniprot.get_protein(accession)
            sides.append((protein_summary(record), record["sequence"]))
        else:
            raise engine.InvalidSequence("Choose two proteins (accessions) or paste two sequences to compare.")
    (first, first_seq), (second, second_seq) = sides
    result = engine.compare_proteins(first_seq, second_seq, body.kmer_size)
    return {"first": first, "second": second, **result}


# --------------------------------------------------------------------------
# reference + dataset
# --------------------------------------------------------------------------
@app.get(f"/api/{API_VERSION}/reference/amino-acids")
def amino_acid_reference() -> list[dict[str, Any]]:
    return [reference_row(code) for code in pa.AMINO_ACIDS]


@app.get(f"/api/{API_VERSION}/reference/amino-acids/{{code}}")
def amino_acid_detail(code: str) -> dict[str, Any]:
    row = reference_row(code)
    return {"residue": row, "structure": references.get(row["code"])}


@app.get(f"/api/{API_VERSION}/dataset")
def dataset_overview() -> dict[str, Any]:
    return {"default_query": services.DEFAULT_QUERY,
            "scopes": [{"key": k, **v} for k, v in services.SCOPES.items()],
            "prosite": prosite.info()}


@app.get(f"/api/{API_VERSION}/dataset/stats")
def dataset_statistics(scope: str = "human_reviewed", refresh: bool = False) -> dict[str, Any]:
    return dataset_stats.status(scope, refresh)


@app.get(f"/api/{API_VERSION}/dataset/universe")
def dataset_universe() -> dict[str, Any]:
    return dataset_stats.universe()


# --------------------------------------------------------------------------
# front end
# --------------------------------------------------------------------------
FRONTEND_DIR = Path(__file__).resolve().parent / "frontend"
if not (FRONTEND_DIR / "index.html").is_file():
    raise RuntimeError(
        f"The 'frontend' folder with index.html was not found next to api.py ({FRONTEND_DIR}). "
        "Extract the whole project zip and keep the folder structure.")
app.mount("/static", StaticFiles(directory=FRONTEND_DIR), name="static")
app.mount("/", StaticFiles(directory=FRONTEND_DIR, html=True), name="frontend")
