"""Tests for the service layer using a fake HTTP client (no network needed)."""
import gzip
import io
import json
import sys
import time
import types

try:                                    # the CLI needs `rich`; the tests do not
    import rich  # noqa: F401
except ImportError:                     # pragma: no cover
    class _Any:
        def __init__(self, *a, **k): pass
        def __getattr__(self, n): return _Any()
        def __call__(self, *a, **k): return _Any()
    for _name in ("rich", "rich.box", "rich.console", "rich.panel", "rich.rule", "rich.table", "rich.text"):
        _mod = types.ModuleType(_name)
        _mod.__getattr__ = lambda n: _Any
        sys.modules[_name] = _mod
    sys.modules["rich"].box = sys.modules["rich.box"]
    for _n, _m in (("Console", "rich.console"), ("Panel", "rich.panel"), ("Rule", "rich.rule"),
                   ("Table", "rich.table"), ("Text", "rich.text")):
        setattr(sys.modules[_m], _n, _Any)

import pytest

import atlas_engine as E
import atlas_services as S

MINI_PROSITE = """CC   Release 2099_01 of 01-Jan-2099
//
ID   ASN_GLYCOSYLATION; PATTERN.
AC   PS00001;
DE   N-glycosylation site.
PA   N-{P}-[ST]-{P}.
//
ID   PKC_PHOSPHO_SITE; PATTERN.
AC   PS00005;
DE   Protein kinase C phosphorylation site.
PA   [ST]-x-[RK].
CC   /SKIP-FLAG=TRUE;
//
"""


def entry(i, sequence=None, reviewed=True):
    acc = f"P{i:05d}"
    out = {"primaryAccession": acc, "uniProtkbId": f"PROT{i}_HUMAN",
           "entryType": "UniProtKB reviewed (Swiss-Prot)" if reviewed else "UniProtKB unreviewed (TrEMBL)",
           "proteinDescription": {"recommendedName": {"fullName": {"value": f"Protein number {i}"}}},
           "genes": [{"geneName": {"value": f"G{i}"}}],
           "organism": {"scientificName": "Homo sapiens"},
           "sequence": {"length": 100 + i}}
    if sequence:
        out["sequence"] = {"value": sequence, "length": len(sequence)}
    return out


class FakeHttp:
    """Routes are (predicate(url, params), handler(url, params) -> Response)."""

    def __init__(self, routes=None):
        self.routes = routes or []
        self.calls = []

    def request(self, url, params=None, headers=None, data=None, method=None, timeout=None, retries=None):
        self.calls.append((url, dict(params or {})))
        for match, handler in self.routes:
            if match(url, params or {}):
                result = handler(url, params or {})
                if isinstance(result, Exception):
                    raise result
                return result
        raise S.UpstreamError("not_found", "no route", 404)


def ok_json(payload, headers=None):
    return S.Response(200, {k.lower(): v for k, v in (headers or {}).items()},
                      json.dumps(payload).encode())


def uniprot_search_route(total):
    """A fake UniProt search API with cursor paging, X-Total-Results and Link headers."""
    def handler(url, params):
        size = int(params.get("size", 25))
        start = int(params["cursor"]) if params.get("cursor") else 0
        order = list(range(total))
        if params.get("sort") == "accession desc":
            order.reverse()
        chunk = order[start:start + size]
        headers = {"X-Total-Results": str(total)}
        if start + size < total:
            headers["Link"] = (f'<https://rest.uniprot.org/uniprotkb/search?cursor={start + size}'
                               f'&query=x&size={size}>; rel="next"')
        return ok_json({"results": [entry(i) for i in chunk]}, headers)
    return (lambda url, params: url.endswith("/uniprotkb/search")), handler


# ---------------------------------------------------------------- pagination
def test_pagination_pages_are_distinct_contiguous_and_totals_are_real():
    client = S.UniProtClient(FakeHttp([uniprot_search_route(1000)]))
    p1 = client.page("reviewed:true", 1, 24)
    p2 = client.page("reviewed:true", 2, 24)
    assert p1["total_items"] == 1000 and p1["total_pages"] == 42 and p1["page_size"] == 24
    assert [i["id"] for i in p1["items"]] == [f"P{n:05d}" for n in range(24)]
    assert [i["id"] for i in p2["items"]] == [f"P{n:05d}" for n in range(24, 48)]
    last = client.page("reviewed:true", 42, 24)
    assert [i["id"] for i in last["items"]][-1] == "P00999" and len(last["items"]) == 1000 - 41 * 24


def test_pagination_going_back_and_jumping_and_out_of_range():
    http = FakeHttp([uniprot_search_route(2000)])
    client = S.UniProtClient(http)
    third = client.page("q", 3, 24)
    assert third["items"][0]["id"] == "P00048"
    back = client.page("q", 1, 24)
    assert back["items"][0]["id"] == "P00000"
    jump = client.page("q", 60, 24)             # inside chunk 2, reached by walking chunk cursors
    assert jump["items"][0]["id"] == "P01416"
    with pytest.raises(S.PageOutOfRange):
        client.page("q", 500, 24)
    with pytest.raises(ValueError):
        client.page("q", 1, 25)


def test_large_result_sets_use_sorted_paging_and_can_reach_the_last_page():
    http = FakeHttp([uniprot_search_route(20000)])
    client = S.UniProtClient(http)
    first = client.page("big", 1, 24)
    assert first["ordering"] == "accession" and first["total_pages"] == 834
    last = client.page("big", 834, 24)          # far away: served by paging backwards from the end
    assert last["items"][-1]["id"] == "P19999" and last["items"][0]["id"] == "P19992"[:0] + f"P{20000 - (20000 - 833 * 24):05d}"
    assert len(last["items"]) == 20000 - 833 * 24
    assert sum(1 for _, p in http.calls if p.get("sort") == "accession desc") >= 1
    assert len(http.calls) < 10


def test_unreachable_page_is_reported_not_faked():
    # sorting refused -> only sequential paging is possible, far jumps must say so
    def refuse_sort(url, params):
        if params.get("sort"):
            return S.UpstreamError("bad_request", "sort refused", 400)
        return uniprot_search_route(20000)[1](url, params)
    client = S.UniProtClient(FakeHttp([(lambda u, p: u.endswith("/uniprotkb/search"), refuse_sort)]))
    assert client.page("q", 1, 24)["ordering"] == "relevance"
    with pytest.raises(S.PageUnreachable):
        client.page("q", 800, 24)


def test_empty_result_set_is_a_valid_page():
    def empty(url, params):
        return ok_json({"results": []}, {"X-Total-Results": "0"})
    client = S.UniProtClient(FakeHttp([(lambda u, p: True, empty)]))
    page = client.page("nothing", 1, 24)
    assert page["items"] == [] and page["total_items"] == 0 and page["total_pages"] == 1


def test_get_protein_names_for_reviewed_and_unreviewed_entries():
    trembl = {"primaryAccession": "A0A000", "uniProtkbId": "A0A000_HUMAN",
              "entryType": "UniProtKB unreviewed (TrEMBL)",
              "proteinDescription": {"submissionNames": [{"fullName": {"value": "Submitted name"}}]},
              "organism": {"scientificName": "Homo sapiens"}, "sequence": {"value": "MKTIIALSYI", "length": 10},
              "uniProtKBCrossReferences": [{"database": "PDB", "id": "1ABC", "properties": [
                  {"key": "Method", "value": "X-ray"}, {"key": "Resolution", "value": "2.0 A"},
                  {"key": "Chains", "value": "A=1-10"}]}]}
    http = FakeHttp([(lambda u, p: u.endswith("A0A000.json"), lambda u, p: ok_json(trembl)),
                     (lambda u, p: u.endswith("P00001.json"), lambda u, p: ok_json(entry(1, "MKTIIALSYIFC")))])
    client = S.UniProtClient(http)
    a = client.get_protein("a0a000")
    assert a["name"] == "Submitted name" and a["reviewed"] is False and a["pdb_refs"][0]["id"] == "1ABC"
    b = client.get_protein("P00001")
    assert b["name"] == "Protein number 1" and b["gene"] == "G1" and b["length"] == 12
    client.get_protein("P00001")
    assert sum(1 for u, _ in http.calls if u.endswith("P00001.json")) == 1      # cached
    with pytest.raises(S.UpstreamError) as info:
        client.get_protein("ZZZ999")
    assert info.value.kind == "not_found"


# ---------------------------------------------------------------- structures
def _protein(pdb=True):
    refs = ([{"id": "1TUP", "method": "X-ray", "resolution": 2.2, "chains": "A/B", "start": 94, "end": 312},
             {"id": "2MNB", "method": "NMR", "resolution": None, "chains": "A", "start": 1, "end": 60}]
            if pdb else [])
    return {"id": "P04637", "sequence": "M" * 393, "length": 393, "pdb_refs": refs}


def test_structure_listing_distinguishes_found_none_and_unavailable(tmp_path):
    af = [{"entryId": "AF-P04637-F1", "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v6.pdb",
           "globalMetricValue": 76.4, "uniprotStart": 1, "uniprotEnd": 393, "latestVersion": 6}]
    http = FakeHttp([(lambda u, p: "alphafold.ebi.ac.uk/api" in u, lambda u, p: ok_json(af)),
                     (lambda u, p: "swissmodel" in u, lambda u, p: S.UpstreamError("not_found", "x", 404))])
    listing = S.StructureService(http, tmp_path).list_for(_protein())
    by_key = {s["key"]: s for s in listing["sources"]}
    assert by_key["pdb"]["status"] == "found" and by_key["pdb"]["entries"][0]["id"] == "1TUP"
    assert by_key["alphafold"]["status"] == "found" and by_key["swissmodel"]["status"] == "none"
    assert by_key["esmfold"]["status"] == "available"
    assert listing["recommended"] == {"source": "pdb", "id": "1TUP"}
    assert "Structures located in" in listing["summary"]

    down = FakeHttp([(lambda u, p: True, lambda u, p: S.UpstreamError("unavailable", "timeout"))])
    bad = S.StructureService(down, tmp_path).list_for(_protein(pdb=False))
    assert all(s["status"] == "unavailable" for s in bad["sources"] if s["key"] in ("alphafold", "swissmodel"))
    assert bad["recommended"] is None and "Could not be checked" in bad["summary"]

    nothing = FakeHttp([])
    none = S.StructureService(nothing, tmp_path).list_for(_protein(pdb=False))
    assert none["recommended"] is None and "No experimental structure or model was located" in none["summary"]


def test_no_alphafold_model_does_not_hide_an_experimental_structure(tmp_path):
    http = FakeHttp([])                             # every model database answers "not found"
    listing = S.StructureService(http, tmp_path).list_for(_protein())
    assert listing["recommended"] == {"source": "pdb", "id": "1TUP"}
    assert {s["key"]: s["status"] for s in listing["sources"]}["alphafold"] == "none"


def test_structure_files_are_downloaded_once_and_cif_is_the_fallback(tmp_path):
    calls = []

    def route(url, params):
        calls.append(url)
        if url.endswith(".pdb"):
            return S.UpstreamError("not_found", "no pdb", 404)
        return S.Response(200, {}, b"data_1ABC\n#\n")
    service = S.StructureService(FakeHttp([(lambda u, p: "files.rcsb.org" in u, route)]), tmp_path)
    text, fmt, meta = service.fetch_file("pdb", "1abc")
    assert fmt == "cif" and text.startswith("data_1ABC") and meta["id"] == "1ABC"
    n = len(calls)
    assert service.fetch_file("pdb", "1ABC")[1] == "cif" and len(calls) == n       # cached on disk
    with pytest.raises(S.UpstreamError):
        service.fetch_file("pdb", "not-an-id")


def test_alphafold_file_uses_the_url_from_the_api_not_a_hard_coded_version(tmp_path):
    af = [{"entryId": "AF-P04637-F1", "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v9.pdb",
           "globalMetricValue": 80.0}]
    seen = []

    def file_route(url, params):
        seen.append(url)
        return S.Response(200, {}, b"ATOM      1  N   MET A   1\n")
    http = FakeHttp([(lambda u, p: "/api/prediction/" in u, lambda u, p: ok_json(af)),
                     (lambda u, p: u.endswith("model_v9.pdb"), file_route)])
    text, fmt, meta = S.StructureService(http, tmp_path).fetch_file("alphafold", "P04637")
    assert fmt == "pdb" and "ATOM" in text and seen == ["https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v9.pdb"]
    assert meta["mean_plddt"] == 80.0


def test_structure_downloads_refuse_unknown_hosts(tmp_path):
    with pytest.raises(S.UpstreamError):
        S.StructureService._check_host("https://evil.example.com/x.pdb")
    with pytest.raises(S.UpstreamError):
        S.StructureService._check_host("http://alphafold.ebi.ac.uk/x.pdb")


# ---------------------------------------------------------------- reference structures
def test_amino_acid_reference_structures_work_offline_for_all_twenty(tmp_path):
    refs = S.ReferenceStructures(FakeHttp([]), tmp_path)      # any network call would raise
    for code in "ACDEFGHIKLMNPQRSTVWY":
        model = refs.get(code)
        assert model["format"] == "sdf" and "V2000" in model["data"] and model["data"].rstrip().endswith("$$$$")
    with pytest.raises(S.UpstreamError):
        refs.get("Z")


# ---------------------------------------------------------------- PROSITE
def test_prosite_manager_states_and_honest_scan_results(tmp_path):
    http = FakeHttp([(lambda u, p: "prosite.dat" in u, lambda u, p: S.Response(200, {}, MINI_PROSITE.encode()))])
    manager = S.PrositeManager(http, tmp_path, urls=["https://ftp.expasy.org/databases/prosite/prosite.dat"])
    assert manager.state == "idle"
    manager.start()
    manager.wait(5)
    assert manager.state == "ready"
    result = manager.scan("AANASAASRKAAAT")
    assert result["status"] == "completed" and result["patterns_tested"] == 2
    assert {m["id"] for m in result["matches"]} == {"ASN_GLYCOSYLATION", "PKC_PHOSPHO_SITE"}
    zero = manager.scan("GGGGGGGG")
    assert zero["status"] == "completed" and zero["patterns_tested"] == 2 and zero["matches"] == []
    assert manager.info()["release"] == "2099_01"
    assert (tmp_path / "prosite.dat").exists()
    again = S.PrositeManager(FakeHttp([]), tmp_path, urls=["https://x"])     # cache: no network needed
    again.start()
    again.wait(5)
    assert again.state == "ready"


def test_prosite_unavailable_is_not_reported_as_no_match(tmp_path):
    http = FakeHttp([(lambda u, p: True, lambda u, p: S.UpstreamError("unavailable", "timeout"))])
    manager = S.PrositeManager(http, tmp_path, urls=["https://a", "https://b"])
    manager.start()
    manager.wait(5)
    assert manager.state == "failed"
    result = manager.scan("MKTIIALSYIFC")
    assert result["status"] == "database_unavailable" and result["patterns_tested"] == 0
    assert "no scan was performed" in result["message"]


# ---------------------------------------------------------------- dataset statistics
def _seqs(n):
    return [("MKTIIALSYIFCLVFADYKDDDDK" * (1 + i % 5))[: 20 + 7 * (i % 15)] for i in range(n)]


def _fasta(n):
    lines = []
    for i, seq in enumerate(_seqs(n)):
        lines += [f">sp|P{i:05d}|X_HUMAN Protein {i}", seq[:30], seq[30:]]
    return "\n".join(lines) + "\n"


def test_dataset_stats_cover_every_streamed_sequence_and_are_cached(tmp_path):
    http = FakeHttp([(lambda u, p: u.endswith("/uniprotkb/search"),
                      lambda u, p: ok_json({"results": []}, {"X-Total-Results": "120"}))])
    manager = S.DatasetStatsManager(http, tmp_path)
    manager._open_stream = lambda query: io.BytesIO(gzip.compress(_fasta(120).encode()))
    manager.start_sync("human_reviewed")
    status = manager.status("human_reviewed")
    assert status["state"] == "ready"
    stats = status["stats"]
    assert stats["total_proteins"] == 120 and stats["reported_total"] == 120
    assert sum(stats["length"]["histogram"]["counts"]) == 120
    lengths = [len(x) for x in _seqs(120)]
    assert stats["length"]["min"] == min(lengths) and stats["length"]["max"] == max(lengths)
    assert stats["length"]["mean"] == pytest.approx(sum(lengths) / 120)
    assert stats["total_residues"] == sum(lengths)
    assert abs(sum(r["percent"] for r in stats["amino_acids"]) - 100.0) < 1e-6
    fresh = S.DatasetStatsManager(FakeHttp([]), tmp_path)      # a new process reads the cache
    cached = fresh.status("human_reviewed")
    assert cached["state"] == "ready" and cached["stats"]["total_proteins"] == 120


def test_dataset_stats_failure_is_reported_as_a_state(tmp_path):
    manager = S.DatasetStatsManager(FakeHttp([]), tmp_path)

    def boom(query):
        raise OSError("network down")
    manager._open_stream = boom
    manager.start_sync("human_reviewed")
    status = manager.status("human_reviewed")
    assert status["state"] == "failed" and "unavailable" in status["error"]
    with pytest.raises(S.UpstreamError):
        manager.status("no_such_scope")
