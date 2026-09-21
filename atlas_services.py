"""Protein Atlas - service layer (no web-framework imports, fully testable).

* HttpClient            - urllib wrapper with retries and explicit error kinds
* UniProtClient         - server-side pagination (cursor based), protein records
* StructureService      - experimental (PDB) -> AlphaFold -> SWISS-MODEL -> ESMFold pipeline
* PrositeManager        - loads / caches the PROSITE database in the background
* DatasetStatsManager   - streams a UniProt dataset once, caches complete statistics
* ReferenceStructures   - 3D coordinates for the 20 amino acids (PubChem / RCSB)
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from collections import OrderedDict
from math import ceil
from pathlib import Path
from typing import Any, Callable, NamedTuple, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, quote, urlencode, urlparse
from urllib.request import Request as UrlRequest, urlopen

import atlas_aminoacids as local_aa
import atlas_engine as engine

USER_AGENT = "ProteinAtlas/3.0 (+https://github.com/protein-sequence-analyzer)"

UNIPROT_BASE = os.getenv("UNIPROT_API_BASE_URL", "https://rest.uniprot.org").rstrip("/")
PROSITE_URLS = [u for u in os.getenv(
    "PROSITE_URL",
    "https://ftp.expasy.org/databases/prosite/prosite.dat,"
    "https://ftp.ebi.ac.uk/pub/databases/prosite/prosite.dat").split(",") if u.strip()]
DEFAULT_QUERY = os.getenv("UNIPROT_QUERY", "reviewed:true AND organism_id:9606")
CACHE_DIR = Path(os.getenv("ATLAS_CACHE_DIR", str(Path(__file__).resolve().parent / ".atlas_cache")))


# ====================================================================
# HTTP
# ====================================================================

class UpstreamError(Exception):
    """A remote service failed. ``kind`` is not_found | bad_request | unavailable."""

    def __init__(self, kind: str, message: str, status: int = 0) -> None:
        super().__init__(message)
        self.kind = kind
        self.status = status


class Response(NamedTuple):
    status: int
    headers: dict
    body: bytes

    def text(self) -> str:
        data = self.body
        if data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
        return data.decode("utf-8", errors="replace")

    def json(self) -> Any:
        try:
            return json.loads(self.text())
        except json.JSONDecodeError as error:
            raise UpstreamError("unavailable", f"Malformed response: {error}") from error


class HttpClient:
    def __init__(self, timeout: float = 25.0, retries: int = 2) -> None:
        self.timeout = timeout
        self.retries = retries

    def request(self, url: str, params: Optional[dict] = None, headers: Optional[dict] = None,
                data: Optional[bytes] = None, method: Optional[str] = None,
                timeout: Optional[float] = None, retries: Optional[int] = None) -> Response:
        if params:
            url += ("&" if "?" in url else "?") + urlencode(params)
        head = {"User-Agent": USER_AGENT, "Accept": "*/*"}
        head.update(headers or {})
        request = UrlRequest(url, data=data, headers=head,
                             method=method or ("POST" if data is not None else "GET"))
        attempts = (self.retries if retries is None else retries) + 1
        last: Optional[Exception] = None
        for attempt in range(attempts):
            try:
                with urlopen(request, timeout=timeout or self.timeout) as response:
                    body = response.read()
                    return Response(response.status,
                                    {k.lower(): v for k, v in response.headers.items()}, body)
            except HTTPError as error:
                if error.code == 404:
                    raise UpstreamError("not_found", f"Not found (HTTP 404): {url}", 404) from error
                if error.code in (400, 422):
                    raise UpstreamError("bad_request", f"Rejected by the remote service (HTTP {error.code})",
                                        error.code) from error
                last = error
                if error.code not in (429, 500, 502, 503, 504):
                    break
            except (URLError, TimeoutError, ConnectionError, OSError) as error:
                last = error
            if attempt < attempts - 1:
                time.sleep(0.6 * (attempt + 1))
        raise UpstreamError("unavailable", f"Remote service unavailable: {last}") from last


# ====================================================================
# UniProt records + server-side pagination
# ====================================================================

PAGE_SIZES = (12, 24, 48, 96)
CHUNK = 480                       # divisible by every supported page size
SEARCH_FIELDS = "accession,id,protein_name,gene_names,organism_name,length"
WALK_BUDGET = 12                  # max sequential chunk requests for one page jump


class PageOutOfRange(Exception):
    def __init__(self, total_pages: int) -> None:
        super().__init__(f"Page is out of range (1 to {total_pages}).")
        self.total_pages = total_pages


class PageUnreachable(Exception):
    """UniProt only pages sequentially; this jump would need too many requests."""


def _entry_name(entry: dict) -> str:
    desc = entry.get("proteinDescription") or {}
    candidates = []
    rec = desc.get("recommendedName") or {}
    candidates.append((rec.get("fullName") or {}).get("value"))
    for group in ("submissionNames", "alternativeNames"):
        items = desc.get(group) or []
        if items:
            candidates.append((items[0].get("fullName") or {}).get("value"))
    candidates.append(entry.get("uniProtkbId"))
    for c in candidates:
        if c:
            return str(c)
    return ""


def _entry_gene(entry: dict) -> str:
    genes = entry.get("genes") or []
    if not genes:
        return ""
    g = genes[0]
    for key in ("geneName",):
        if (g.get(key) or {}).get("value"):
            return g[key]["value"]
    for key in ("orderedLocusNames", "orfNames"):
        items = g.get(key) or []
        if items and items[0].get("value"):
            return items[0]["value"]
    return ""


def summary_from_entry(entry: dict) -> dict[str, Any]:
    seq = entry.get("sequence") or {}
    length = seq.get("length")
    if length is None and seq.get("value"):
        length = len(seq["value"])
    entry_type = str(entry.get("entryType", "")).lower()
    return {
        "id": entry.get("primaryAccession", ""),
        "entry_name": entry.get("uniProtkbId", ""),
        "name": _entry_name(entry),
        "gene": _entry_gene(entry),
        "organism": (entry.get("organism") or {}).get("scientificName", ""),
        "length": int(length) if length is not None else None,
        "reviewed": ("reviewed" in entry_type and "unreviewed" not in entry_type)
        if entry_type else None,
    }


def _parse_pdb_refs(entry: dict) -> list[dict[str, Any]]:
    refs = []
    for x in entry.get("uniProtKBCrossReferences") or []:
        if x.get("database") != "PDB":
            continue
        props = {p.get("key"): p.get("value") for p in x.get("properties") or []}
        res = None
        m = re.match(r"([\d.]+)", props.get("Resolution", "") or "")
        if m:
            res = float(m.group(1))
        chains = props.get("Chains", "") or ""
        span = re.search(r"=(\d+)-(\d+)", chains)
        start, end = (int(span.group(1)), int(span.group(2))) if span else (None, None)
        refs.append({"id": x.get("id"), "method": props.get("Method", ""),
                     "resolution": res, "chains": chains.split("=")[0] if chains else "",
                     "start": start, "end": end})
    return refs


class UniProtClient:
    def __init__(self, http: HttpClient, base_url: str = UNIPROT_BASE) -> None:
        self.http = http
        self.base = base_url.rstrip("/")
        self._lock = threading.RLock()
        self._plans: dict[str, dict] = {}
        self._cursors: dict[tuple, dict[int, Optional[str]]] = {}
        self._chunks: "OrderedDict[tuple, list]" = OrderedDict()
        self._proteins: "OrderedDict[str, dict]" = OrderedDict()

    # ---- single records -------------------------------------------------
    def get_protein(self, protein_id: str) -> dict[str, Any]:
        key = protein_id.strip().upper()
        if not re.fullmatch(r"[A-Z0-9_\-]{2,30}", key):
            raise UpstreamError("bad_request", "That does not look like a UniProt accession or entry name.")
        with self._lock:
            if key in self._proteins:
                self._proteins.move_to_end(key)
                return self._proteins[key]
        try:
            entry = self.http.request(f"{self.base}/uniprotkb/{quote(key)}.json",
                                      headers={"Accept": "application/json"}).json()
        except UpstreamError as error:
            if error.kind in ("bad_request", "not_found"):
                raise UpstreamError("not_found", f"'{protein_id}' was not found in UniProt.") from error
            raise
        record = self.record_from_entry(entry)
        with self._lock:
            self._proteins[key] = record
            self._proteins[record["id"].upper()] = record
            while len(self._proteins) > 400:
                self._proteins.popitem(last=False)
        return record

    @staticmethod
    def record_from_entry(entry: dict) -> dict[str, Any]:
        seq = entry.get("sequence") or {}
        sequence = seq.get("value", "")
        info = summary_from_entry(entry)
        function_text = ""
        for comment in entry.get("comments") or []:
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts") or []
                if texts:
                    function_text = texts[0].get("value", "")
                    break
        keywords = [k.get("name") for k in entry.get("keywords") or [] if k.get("name")]
        af_ids = [x.get("id") for x in entry.get("uniProtKBCrossReferences") or []
                  if x.get("database") == "AlphaFoldDB"]
        header = (f"sp|{info['id']}|{info['entry_name']} {info['name']} "
                  f"OS={info['organism']} GN={info['gene']}")
        return {**info, "length": len(sequence), "sequence": sequence, "header": header,
                "function": function_text, "keywords": keywords[:12],
                "pdb_refs": _parse_pdb_refs(entry), "alphafold_ids": af_ids,
                "molecular_weight_uniprot": seq.get("molWeight")}

    # ---- listing / paging --------------------------------------------------
    def _fetch_chunk(self, query: str, sort: Optional[str], cursor: Optional[str]):
        params = {"query": query, "format": "json", "size": CHUNK, "fields": SEARCH_FIELDS}
        if sort:
            params["sort"] = sort
        if cursor:
            params["cursor"] = cursor
        response = self.http.request(f"{self.base}/uniprotkb/search", params=params,
                                     headers={"Accept": "application/json"})
        payload = response.json()
        total = int(response.headers.get("x-total-results", payload.get("totalResults", 0)) or 0)
        link = response.headers.get("link", "")
        match = re.search(r"[?&]cursor=([^&>]+)", link) if 'rel="next"' in link else None
        next_cursor = match.group(1) if match else None
        items = [summary_from_entry(e) for e in payload.get("results", [])]
        return items, next_cursor, total

    def _get_chunk(self, query: str, direction: str, sort_mode: bool, idx: int) -> list:
        key = (query, direction, idx)
        with self._lock:
            if key in self._chunks:
                self._chunks.move_to_end(key)
                return self._chunks[key]
            cursors = self._cursors.setdefault((query, direction), {0: None})
            known = max(i for i in cursors if i <= idx)
        if idx - known > WALK_BUDGET:
            raise PageUnreachable(
                "UniProt can only be paged sequentially, and that page is too far away to reach "
                "quickly. Use Next / Previous, or narrow your search.")
        sort = None
        if sort_mode:
            sort = "accession asc" if direction == "asc" else "accession desc"
        for i in range(known, idx + 1):
            with self._lock:
                cached = self._chunks.get((query, direction, i))
                cursor = cursors.get(i)
            if cached is not None and (i + 1) in cursors:
                continue
            items, next_cursor, _ = self._fetch_chunk(query, sort, cursor)
            with self._lock:
                self._chunks[(query, direction, i)] = items
                while len(self._chunks) > 40:
                    self._chunks.popitem(last=False)
                if next_cursor:
                    cursors[i + 1] = next_cursor
        with self._lock:
            if key not in self._chunks:
                raise UpstreamError("unavailable", "UniProt did not return that page of results.")
            return self._chunks[key]

    def _plan(self, query: str) -> dict:
        with self._lock:
            plan = self._plans.get(query)
        if plan:
            return plan
        items, next_cursor, total = self._fetch_chunk(query, None, None)
        sort_mode = total > CHUNK * WALK_BUDGET
        if sort_mode:
            try:
                items, next_cursor, total = self._fetch_chunk(query, "accession asc", None)
            except UpstreamError as error:
                if error.kind != "bad_request":
                    raise
                sort_mode = False          # sorting refused: keep default order, sequential paging only
        with self._lock:
            self._chunks[(query, "asc", 0)] = items
            cursors = self._cursors.setdefault((query, "asc"), {0: None})
            if next_cursor:
                cursors[1] = next_cursor
            plan = {"total": total, "sorted": sort_mode}
            self._plans[query] = plan
            if len(self._plans) > 200:
                self._plans.pop(next(iter(self._plans)))
        return plan

    def _range(self, query: str, plan: dict, start: int, end: int) -> list:
        total, sort_mode = plan["total"], plan["sorted"]

        def collect(direction: str, a: int, b: int) -> list:
            out: list = []
            for idx in range(a // CHUNK, (b - 1) // CHUNK + 1):
                chunk = self._get_chunk(query, direction, sort_mode, idx)
                lo = max(a, idx * CHUNK) - idx * CHUNK
                hi = min(b, (idx + 1) * CHUNK) - idx * CHUNK
                out.extend(chunk[lo:hi])
            return out

        if not sort_mode:
            return collect("asc", start, end)

        def cost(direction: str, a: int) -> int:
            cursors = self._cursors.get((query, direction), {0: None})
            idx = a // CHUNK
            return idx - max(i for i in cursors if i <= idx)

        if cost("asc", start) <= cost("desc", total - end):
            return collect("asc", start, end)
        return list(reversed(collect("desc", total - end, total - start)))

    def page(self, query: str, page: int, page_size: int) -> dict[str, Any]:
        if page_size not in PAGE_SIZES:
            raise ValueError(f"page_size must be one of {', '.join(map(str, PAGE_SIZES))}")
        query = (query or "").strip() or DEFAULT_QUERY
        plan = self._plan(query)
        total = plan["total"]
        total_pages = max(1, ceil(total / page_size)) if total else 1
        if total == 0:
            return {"items": [], "page": 1, "page_size": page_size, "total_items": 0,
                    "total_pages": 1, "query": query, "ordering": "relevance"}
        if page < 1 or page > total_pages:
            raise PageOutOfRange(total_pages)
        start = (page - 1) * page_size
        end = min(total, start + page_size)
        items = self._range(query, plan, start, end)
        return {"items": items, "page": page, "page_size": page_size, "total_items": total,
                "total_pages": total_pages, "query": query,
                "ordering": "accession" if plan["sorted"] else "relevance",
                "source": "UniProtKB"}

    def total_for(self, query: str) -> Optional[int]:
        try:
            response = self.http.request(f"{self.base}/uniprotkb/search",
                                         params={"query": query, "format": "json", "size": 1,
                                                 "fields": "accession"})
            return int(response.headers.get("x-total-results", 0))
        except (UpstreamError, ValueError):
            return None


# ====================================================================
# Structures
# ====================================================================

ALLOWED_MODEL_HOSTS = ("alphafold.ebi.ac.uk", "swissmodel.expasy.org")
ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMFOLD_MAX = 400


def _safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", text)


class StructureService:
    def __init__(self, http: HttpClient, cache_dir: Path = CACHE_DIR) -> None:
        self.http = http
        self.dir = Path(cache_dir) / "structures"
        self._lock = threading.Lock()
        self._listing_cache: dict[str, tuple[float, dict]] = {}

    # ---- disk cache -----------------------------------------------------
    def _cache_get(self, name: str) -> Optional[tuple[str, str]]:
        for fmt in ("pdb", "cif", "sdf"):
            path = self.dir / f"{_safe_name(name)}.{fmt}"
            if path.exists() and path.stat().st_size > 0:
                return path.read_text(encoding="utf-8", errors="replace"), fmt
        return None

    def _cache_put(self, name: str, fmt: str, text: str) -> None:
        try:
            self.dir.mkdir(parents=True, exist_ok=True)
            (self.dir / f"{_safe_name(name)}.{fmt}").write_text(text, encoding="utf-8")
        except OSError:
            pass                                    # caching is best-effort

    # ---- listing -----------------------------------------------------------
    def list_for(self, protein: dict) -> dict[str, Any]:
        acc = protein["id"]
        with self._lock:
            cached = self._listing_cache.get(acc)
        if cached and time.time() - cached[0] < 3600:
            return cached[1]

        with ThreadPoolExecutor(max_workers=2) as pool:
            af_future = pool.submit(self._alphafold_source, acc)
            sm_future = pool.submit(self._swissmodel_source, acc)
            af_source, sm_source = af_future.result(), sm_future.result()
        sources = [self._pdb_source(protein), af_source, sm_source, self._esm_source(protein)]
        recommended = None
        for src in sources:
            if src["status"] == "found" and src["entries"]:
                recommended = {"source": src["key"], "id": src["entries"][0]["id"]}
                break
        if recommended is None and sources[-1]["status"] == "available":
            recommended = None      # prediction is offered but never auto-run
        found = [s["label"] for s in sources if s["status"] == "found"]
        if found:
            summary = "Structures located in: " + "; ".join(found) + "."
        else:
            checked = [s["label"] for s in sources if s["status"] in ("none", "found")]
            unavailable = [s["label"] for s in sources if s["status"] == "unavailable"]
            summary = ("No experimental structure or model was located for this accession in "
                       + (", ".join(checked) or "the structure databases") + ".")
            if unavailable:
                summary += " Could not be checked: " + ", ".join(unavailable) + "."
        listing = {"accession": acc, "sources": sources, "recommended": recommended,
                   "summary": summary, "sequence_length": len(protein.get("sequence", ""))}
        if not any(s["status"] == "unavailable" for s in sources[:3]):
            with self._lock:
                self._listing_cache[acc] = (time.time(), listing)
        return listing

    def _pdb_source(self, protein: dict) -> dict:
        refs = protein.get("pdb_refs") or []
        length = max(1, protein.get("length") or len(protein.get("sequence", "")) or 1)
        entries = []
        for r in refs:
            span = (r["end"] - r["start"] + 1) if r.get("start") and r.get("end") else None
            entries.append({**r, "coverage_percent": (span / length * 100.0) if span else None,
                            "label": f"{r['id']} · {r.get('method') or 'experimental'}"
                                     + (f" · {r['resolution']:.2f} Å" if r.get("resolution") else "")})

        def rank(e: dict) -> tuple:
            method = (e.get("method") or "").lower()
            nmr = 1 if "nmr" in method else 0
            return (-(e.get("coverage_percent") or 0), nmr, e.get("resolution") or 99.0)

        entries.sort(key=rank)
        return {"key": "pdb", "label": "Experimental structures (PDB)",
                "kind": "experimental", "status": "found" if entries else "none",
                "count": len(entries), "entries": entries[:25],
                "message": (f"{len(entries)} experimental structure(s) cross-referenced by UniProt."
                            if entries else "UniProt lists no PDB entries for this accession.")}

    def _alphafold_source(self, acc: str) -> dict:
        base = {"key": "alphafold", "label": "AlphaFold DB", "kind": "predicted", "entries": [], "count": 0}
        try:
            data = self.http.request(f"https://alphafold.ebi.ac.uk/api/prediction/{quote(acc)}",
                                     headers={"Accept": "application/json"}, retries=1, timeout=15).json()
        except UpstreamError as error:
            if error.kind in ("not_found", "bad_request"):
                return {**base, "status": "none",
                        "message": "AlphaFold DB has no prediction for this accession."}
            return {**base, "status": "unavailable", "message": "Structure database unavailable"}
        entries = []
        for item in data if isinstance(data, list) else []:
            if not item.get("pdbUrl"):
                continue
            entries.append({
                "id": item.get("entryId") or f"AF-{acc}-F1", "label":
                    f"{item.get('entryId', 'AlphaFold model')} · mean pLDDT "
                    f"{item['globalMetricValue']:.1f}" if item.get("globalMetricValue") is not None
                    else item.get("entryId", "AlphaFold model"),
                "mean_plddt": item.get("globalMetricValue"),
                "start": item.get("uniprotStart"), "end": item.get("uniprotEnd"),
                "version": item.get("latestVersion")})
            entries[-1]["_pdb_url"] = item["pdbUrl"]
            entries[-1]["_cif_url"] = item.get("cifUrl")
        with self._lock:
            self._listing_cache[f"_af:{acc}"] = (time.time(), {"entries": entries})
        public = [{k: v for k, v in e.items() if not k.startswith("_")} for e in entries]
        return {**base, "status": "found" if public else "none", "count": len(public),
                "entries": public,
                "message": ("AlphaFold DB predicted model (per-residue confidence is stored as pLDDT)."
                            if public else "AlphaFold DB has no prediction for this accession.")}

    def _swissmodel_source(self, acc: str) -> dict:
        base = {"key": "swissmodel", "label": "SWISS-MODEL Repository", "kind": "homology model",
                "entries": [], "count": 0}
        try:
            data = self.http.request(
                f"https://swissmodel.expasy.org/repository/uniprot/{quote(acc)}.json",
                headers={"Accept": "application/json"}, retries=1, timeout=15).json()
        except UpstreamError as error:
            if error.kind in ("not_found", "bad_request"):
                return {**base, "status": "none", "message": "SWISS-MODEL holds no model for this accession."}
            return {**base, "status": "unavailable", "message": "Structure database unavailable"}
        structures = (data.get("result") or {}).get("structures") or []
        entries = []
        for i, s in enumerate(structures[:12]):
            url = s.get("coordinates")
            if not url:
                continue
            ident = s.get("identity")
            entries.append({"id": f"{acc}:{i}", "start": s.get("from"), "end": s.get("to"),
                            "identity": ident, "coverage": s.get("coverage"),
                            "template": s.get("template"),
                            "label": f"Model {i + 1} · template {s.get('template', '?')}"
                                     + (f" · {ident:.0f}% identity" if isinstance(ident, (int, float)) else ""),
                            "_url": url})
        with self._lock:
            self._listing_cache[f"_sm:{acc}"] = (time.time(), {"entries": entries})
        public = [{k: v for k, v in e.items() if not k.startswith("_")} for e in entries]
        return {**base, "status": "found" if public else "none", "count": len(public),
                "entries": public,
                "message": ("Template-based homology models." if public
                            else "SWISS-MODEL holds no model for this accession.")}

    def _esm_source(self, protein: dict) -> dict:
        n = len(protein.get("sequence", ""))
        base = {"key": "esmfold", "label": "ESMFold (on-demand prediction)", "kind": "predicted",
                "count": 0, "entries": []}
        if n == 0:
            return {**base, "status": "none", "message": "No sequence available."}
        if n > ESMFOLD_MAX:
            return {**base, "status": "too_long",
                    "message": f"ESMFold predicts sequences up to {ESMFOLD_MAX} aa; this one has {n}."}
        return {**base, "status": "available",
                "entries": [{"id": protein["id"], "label": f"Predict with ESMFold ({n} aa)"}],
                "message": "Runs a fresh prediction from the sequence. Low-confidence for very short peptides."}

    # ---- coordinate files -----------------------------------------------------
    def fetch_file(self, source: str, ident: str, protein_loader: Optional[Callable[[str], dict]] = None
                   ) -> tuple[str, str, dict]:
        """Return (text, format, meta). Raises UpstreamError."""
        if source == "pdb":
            if not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", ident):
                raise UpstreamError("bad_request", "Invalid PDB identifier.")
            name = f"pdb_{ident.upper()}"
            hit = self._cache_get(name)
            if hit:
                return hit[0], hit[1], {"source": "PDB", "id": ident.upper(), "cached": True}
            last: Optional[UpstreamError] = None
            for url, fmt in ((f"https://files.rcsb.org/download/{ident.upper()}.pdb", "pdb"),
                             (f"https://files.rcsb.org/download/{ident.upper()}.cif", "cif")):
                try:
                    text = self.http.request(url, retries=1, timeout=60).text()
                    self._cache_put(name, fmt, text)
                    return text, fmt, {"source": "PDB", "id": ident.upper()}
                except UpstreamError as error:
                    last = error
            raise last or UpstreamError("unavailable", "PDB download failed.")

        if source == "alphafold":
            acc = ident
            if not re.fullmatch(r"[A-Za-z0-9\-_]{2,40}", acc):
                raise UpstreamError("bad_request", "Invalid AlphaFold identifier.")
            acc = acc.upper()
            base_acc = re.sub(r"^AF-|-F\d+$", "", acc, flags=re.I)
            name = f"af_{base_acc}"
            hit = self._cache_get(name)
            if hit:
                return hit[0], hit[1], {"source": "AlphaFold DB", "id": acc, "cached": True}
            self._alphafold_source(base_acc)
            listing = self._listing_cache.get(f"_af:{base_acc}", (0, {"entries": []}))[1]["entries"]
            if not listing:
                raise UpstreamError("not_found", "AlphaFold DB has no prediction for this accession.")
            entry = listing[0]
            url = entry["_pdb_url"]
            self._check_host(url)
            text = self.http.request(url, retries=1, timeout=60).text()
            self._cache_put(name, "pdb", text)
            return text, "pdb", {"source": "AlphaFold DB", "id": entry["id"],
                                 "mean_plddt": entry.get("mean_plddt")}

        if source == "swissmodel":
            m = re.fullmatch(r"([A-Za-z0-9_]{2,20}):(\d{1,2})", ident)
            if not m:
                raise UpstreamError("bad_request", "Invalid SWISS-MODEL identifier.")
            acc, index = m.group(1).upper(), int(m.group(2))
            name = f"sm_{acc}_{index}"
            hit = self._cache_get(name)
            if hit:
                return hit[0], hit[1], {"source": "SWISS-MODEL", "id": ident, "cached": True}
            self._swissmodel_source(acc)
            listing = self._listing_cache.get(f"_sm:{acc}", (0, {"entries": []}))[1]["entries"]
            match = next((e for e in listing if e["id"] == f"{acc}:{index}"), None)
            if not match:
                raise UpstreamError("not_found", "That SWISS-MODEL entry no longer exists.")
            self._check_host(match["_url"])
            text = self.http.request(match["_url"], retries=1, timeout=60).text()
            self._cache_put(name, "pdb", text)
            return text, "pdb", {"source": "SWISS-MODEL", "id": ident,
                                 "template": match.get("template")}

        if source == "esmfold":
            if protein_loader is None:
                raise UpstreamError("bad_request", "ESMFold needs a sequence.")
            protein = protein_loader(ident)
            text = self.predict(protein["sequence"])
            return text, "pdb", {"source": "ESMFold", "id": ident, "predicted": True}
        raise UpstreamError("bad_request", "Unknown structure source.")

    @staticmethod
    def _check_host(url: str) -> None:
        parsed = urlparse(url)
        if parsed.scheme != "https" or parsed.hostname not in ALLOWED_MODEL_HOSTS:
            raise UpstreamError("bad_request", "Structure URL is not on an allowed host.")

    def predict(self, sequence: str) -> str:
        sequence = engine.normalize_sequence(sequence)
        if len(sequence) > ESMFOLD_MAX:
            raise UpstreamError("bad_request", f"ESMFold predicts sequences up to {ESMFOLD_MAX} aa.")
        if re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", sequence):
            raise UpstreamError("bad_request", "ESMFold accepts the 20 standard amino acids only.")
        name = "esm_" + hashlib.sha1(sequence.encode()).hexdigest()[:20]
        hit = self._cache_get(name)
        if hit:
            return hit[0]
        text = self.http.request(ESMFOLD_URL, data=sequence.encode(), method="POST",
                                 headers={"Content-Type": "text/plain"}, timeout=120, retries=0).text()
        if "ATOM" not in text:
            raise UpstreamError("unavailable", "ESMFold returned no coordinates.")
        self._cache_put(name, "pdb", text)
        return text


# ====================================================================
# Reference structures for the 20 amino acids
# ====================================================================

class ReferenceStructures:
    def __init__(self, http: HttpClient, cache_dir: Path = CACHE_DIR) -> None:
        self.http = http
        self.dir = Path(cache_dir) / "aminoacids"

    def get(self, code: str) -> dict[str, Any]:
        row = engine.AMINO_ACID_BY_CODE.get(code.upper())
        if not row:
            raise UpstreamError("not_found", "Unknown amino-acid code.")
        try:                                   # built-in verified model: always available, offline
            model = local_aa.structure(row["code"])
            return {"code": row["code"], "format": "sdf", "data": model["data"],
                    "source": model["source"], "note": model["note"]}
        except (KeyError, OSError, ValueError):
            pass
        path = self.dir / f"{row['code']}.sdf"
        if path.exists() and path.stat().st_size > 0:
            return {"code": row["code"], "format": "sdf", "data": path.read_text(encoding="utf-8"),
                    "source": "cache"}
        attempts = [
            (f"PubChem CID {row['pubchem_cid']} (3D conformer)",
             f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{row['pubchem_cid']}/record/SDF/",
             {"record_type": "3d", "response_type": "display"}),
            (f"RCSB chemical component {row['three'].upper()} (ideal coordinates)",
             f"https://files.rcsb.org/ligands/download/{row['three'].upper()}_ideal.sdf", None),
        ]
        errors = []
        for label, url, params in attempts:
            try:
                text = self.http.request(url, params=params, retries=1, timeout=30).text()
                if text.count("\n") > 8 and ("V2000" in text or "V3000" in text):
                    try:
                        self.dir.mkdir(parents=True, exist_ok=True)
                        path.write_text(text, encoding="utf-8")
                    except OSError:
                        pass
                    return {"code": row["code"], "format": "sdf", "data": text, "source": label}
                errors.append(f"{label}: unexpected response")
            except UpstreamError as error:
                errors.append(f"{label}: {error}")
        raise UpstreamError("unavailable", "Structure database unavailable. " + " | ".join(errors))


# ====================================================================
# PROSITE
# ====================================================================

class PrositeManager:
    """Downloads/caches/parses PROSITE in the background so the server starts instantly."""

    def __init__(self, http: HttpClient, cache_dir: Path = CACHE_DIR,
                 urls: Optional[list[str]] = None, max_age_days: float = 30) -> None:
        self.http = http
        self.path = Path(cache_dir) / "prosite.dat"
        self.urls = urls or PROSITE_URLS
        self.max_age = max_age_days * 86400
        self.state = "idle"            # idle | loading | ready | failed
        self.error: Optional[str] = None
        self.db: Optional[engine.PrositeDatabase] = None
        self._thread: Optional[threading.Thread] = None
        self._lock = threading.Lock()

    def start(self) -> None:
        with self._lock:
            if self.state == "loading" or (self.state == "ready" and self.db):
                return
            self.state = "loading"
            self.error = None
            self._thread = threading.Thread(target=self._load, name="prosite-loader", daemon=True)
            self._thread.start()

    def retry(self) -> None:
        with self._lock:
            if self.state == "failed":
                self.state = "idle"
        self.start()

    def wait(self, timeout: float = 0) -> None:
        if self._thread:
            self._thread.join(timeout)

    def _load(self) -> None:
        try:
            text = None
            source = ""
            if self.path.exists() and time.time() - self.path.stat().st_mtime < self.max_age:
                text = self.path.read_text(encoding="utf-8", errors="replace")
                source = f"local cache ({self.path.name})"
            if not text:
                last: Optional[Exception] = None
                for url in self.urls:
                    try:
                        text = self.http.request(url.strip(), timeout=120, retries=1).text()
                        source = url.strip()
                        break
                    except UpstreamError as error:
                        last = error
                if not text:
                    if self.path.exists():         # stale cache beats nothing
                        text = self.path.read_text(encoding="utf-8", errors="replace")
                        source = f"stale local cache ({self.path.name})"
                    else:
                        raise RuntimeError(str(last) if last else "PROSITE download failed")
                else:
                    try:
                        self.path.parent.mkdir(parents=True, exist_ok=True)
                        self.path.write_text(text, encoding="utf-8")
                    except OSError:
                        pass
            db = engine.PrositeDatabase.from_text(text, source)
            if db.patterns_available == 0:
                raise RuntimeError("The PROSITE file contained no usable pattern records.")
            self.db = db
            self.state = "ready"
        except Exception as error:      # noqa: BLE001 - report every failure as state
            self.state = "failed"
            self.error = str(error)

    def scan(self, sequence: str, exclude_frequent: bool = False) -> dict[str, Any]:
        if self.state == "idle":
            self.start()
        if self.state == "loading":
            return engine.motif_status(
                "loading", "The PROSITE database is still loading. This happens once; try again in a moment.")
        if self.state == "failed" or not self.db:
            return engine.motif_status(
                "database_unavailable",
                "PROSITE database unavailable — no scan was performed. " + (self.error or ""))
        try:
            return self.db.scan(sequence, exclude_frequent)
        except engine.InvalidSequence:
            raise
        except Exception as error:      # noqa: BLE001
            return engine.motif_status("failed", f"The PROSITE scan failed: {error}")

    def info(self) -> dict[str, Any]:
        return {"state": self.state, "error": self.error, **(self.db.info() if self.db else {})}


# ====================================================================
# Dataset statistics (complete, cached, built once in the background)
# ====================================================================

SCOPES = {
    "human_reviewed": {"label": "Human · reviewed (Swiss-Prot)",
                       "query": "organism_id:9606 AND reviewed:true"},
    "human_all": {"label": "Human · all UniProtKB entries", "query": "organism_id:9606"},
    "swissprot": {"label": "All reviewed proteins (Swiss-Prot)", "query": "reviewed:true"},
}


class DatasetStatsManager:
    def __init__(self, http: HttpClient, cache_dir: Path = CACHE_DIR, max_age_days: float = 14,
                 base_url: str = UNIPROT_BASE) -> None:
        self.http = http
        self.dir = Path(cache_dir) / "dataset_stats"
        self.max_age = max_age_days * 86400
        self.base = base_url.rstrip("/")
        self._jobs: dict[str, dict] = {}
        self._lock = threading.Lock()
        self._universe: Optional[tuple[float, dict]] = None

    def _cache_path(self, scope: str) -> Path:
        return self.dir / f"{scope}.json"

    def _read_cache(self, scope: str) -> Optional[dict]:
        path = self._cache_path(scope)
        try:
            if path.exists() and time.time() - path.stat().st_mtime < self.max_age:
                return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            pass
        return None

    def status(self, scope: str, refresh: bool = False) -> dict[str, Any]:
        if scope not in SCOPES:
            raise UpstreamError("bad_request", "Unknown dataset scope.")
        meta = {"scope": scope, **SCOPES[scope]}
        with self._lock:
            job = self._jobs.get(scope)
        if job and job["state"] == "building":
            return {**meta, "state": "building", "progress": dict(job["progress"])}
        if refresh:
            self.start(scope)
            return {**meta, "state": "building", "progress": dict(self._jobs[scope]["progress"])}
        if job and job["state"] == "ready":
            return {**meta, "state": "ready", "built_at": job["built_at"], "stats": job["stats"]}
        cached = self._read_cache(scope)
        if cached:
            with self._lock:
                self._jobs[scope] = {"state": "ready", "built_at": cached["built_at"],
                                     "stats": cached["stats"], "progress": {}}
            return {**meta, "state": "ready", "built_at": cached["built_at"], "stats": cached["stats"]}
        if job and job["state"] == "failed":
            return {**meta, "state": "failed", "error": job["error"]}
        self.start(scope)
        return {**meta, "state": "building", "progress": dict(self._jobs[scope]["progress"])}

    def start(self, scope: str) -> None:
        with self._lock:
            job = self._jobs.get(scope)
            if job and job["state"] == "building":
                return
            self._jobs[scope] = {"state": "building",
                                 "progress": {"processed": 0, "total": None, "percent": 0.0,
                                              "phase": "starting"}}
        threading.Thread(target=self._build, args=(scope,), name=f"stats-{scope}", daemon=True).start()

    def start_sync(self, scope: str) -> None:
        with self._lock:
            self._jobs[scope] = {"state": "building",
                                 "progress": {"processed": 0, "total": None, "percent": 0.0,
                                              "phase": "starting"}}
        self._build(scope)

    def _open_stream(self, query: str):
        params = {"query": query, "format": "fasta", "compressed": "true"}
        url = f"{self.base}/uniprotkb/stream?{urlencode(params)}"
        request = UrlRequest(url, headers={"User-Agent": USER_AGENT, "Accept": "text/plain"})
        return urlopen(request, timeout=90)

    def _build(self, scope: str) -> None:
        job = self._jobs[scope]
        try:
            query = SCOPES[scope]["query"]
            try:
                total = int(self.http.request(f"{self.base}/uniprotkb/search",
                                              params={"query": query, "format": "json", "size": 1,
                                                      "fields": "accession"}).headers.get("x-total-results", 0))
            except (UpstreamError, ValueError):
                total = None
            job["progress"].update({"total": total, "phase": "downloading and analysing sequences"})
            aggregator = engine.DatasetAggregator()
            response = self._open_stream(query)
            try:
                head = response.read(2)
                stream: Any = io.BufferedReader(_Prefixed(head, response))
                if head == b"\x1f\x8b":
                    stream = gzip.GzipFile(fileobj=stream)
                lines = io.TextIOWrapper(stream, encoding="utf-8", errors="replace")
                for _header, sequence in engine.iter_fasta(lines):
                    aggregator.add(sequence)
                    n = aggregator.n
                    if n % 500 == 0:
                        job["progress"].update({
                            "processed": n,
                            "percent": (n / total * 100.0) if total else 0.0})
            finally:
                response.close()
            job["progress"].update({"phase": "computing statistics"})
            stats = aggregator.finish()
            if stats.get("status") != "success":
                raise RuntimeError(stats.get("message", "No sequences were returned."))
            stats["scope"] = scope
            stats["query"] = query
            stats["reported_total"] = total
            built_at = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())
            try:
                self.dir.mkdir(parents=True, exist_ok=True)
                self._cache_path(scope).write_text(
                    json.dumps({"built_at": built_at, "stats": stats}), encoding="utf-8")
            except OSError:
                pass
            job.update({"state": "ready", "stats": stats, "built_at": built_at})
        except Exception as error:          # noqa: BLE001
            job.update({"state": "failed", "error": f"Dataset statistics unavailable: {error}"})

    def universe(self) -> dict[str, Any]:
        if self._universe and time.time() - self._universe[0] < 86400:
            return self._universe[1]
        out = {}
        for key, query in (("uniprotkb_total", "*"), ("reviewed_total", "reviewed:true"),
                           ("human_total", "organism_id:9606"),
                           ("human_reviewed_total", "organism_id:9606 AND reviewed:true")):
            try:
                out[key] = int(self.http.request(
                    f"{self.base}/uniprotkb/search",
                    params={"query": query, "format": "json", "size": 1, "fields": "accession"},
                    retries=1).headers.get("x-total-results", 0))
            except (UpstreamError, ValueError):
                out[key] = None
        if any(v is not None for v in out.values()):
            self._universe = (time.time(), out)
        return out


class _Prefixed(io.RawIOBase):
    """Raw stream that replays already-read bytes before the underlying response."""

    def __init__(self, head: bytes, rest) -> None:
        self._head = head
        self._rest = rest

    def readable(self) -> bool:
        return True

    def readinto(self, buffer) -> int:
        if self._head:
            n = min(len(buffer), len(self._head))
            buffer[:n] = self._head[:n]
            self._head = self._head[n:]
            return n
        data = self._rest.read(len(buffer))
        buffer[:len(data)] = data
        return len(data)
