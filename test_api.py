"""HTTP-level tests (need fastapi + httpx).  Every upstream service is replaced by a fake, so no
internet access is used.  Run:  python -m pytest -q"""
import gzip
import io

import pytest
from fastapi.testclient import TestClient

import api
import atlas_services as S
from test_services import FakeHttp, MINI_PROSITE, _fasta, entry, ok_json, uniprot_search_route

P53 = ("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQ"
       "KTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVE"
       "YLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNT"
       "SSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD")


@pytest.fixture()
def client(tmp_path, monkeypatch):
    p53 = {**entry(4637, P53), "primaryAccession": "P04637",
           "uniProtKBCrossReferences": [{"database": "PDB", "id": "1TUP", "properties": [
               {"key": "Method", "value": "X-ray"}, {"key": "Resolution", "value": "2.20 A"},
               {"key": "Chains", "value": "A/B/C=94-312"}]}]}
    http = FakeHttp([
        (lambda u, p: u.endswith("/uniprotkb/search"), uniprot_search_route(300)[1]),
        (lambda u, p: u.endswith("P04637.json"), lambda u, p: ok_json(p53)),
        (lambda u, p: u.endswith("P01308.json"), lambda u, p: ok_json({**entry(1308, "GIVEQCCTSICSLYQLENYCN"), "primaryAccession": "P01308"})),
        (lambda u, p: "prosite.dat" in u, lambda u, p: S.Response(200, {}, MINI_PROSITE.encode())),
    ])
    monkeypatch.setattr(api, "http_client", http)
    monkeypatch.setattr(api, "uniprot", S.UniProtClient(http))
    monkeypatch.setattr(api, "structures", S.StructureService(http, tmp_path))
    monkeypatch.setattr(api, "references", S.ReferenceStructures(http, tmp_path))
    prosite = S.PrositeManager(http, tmp_path, urls=["https://x/prosite.dat"])
    prosite.start()
    prosite.wait(5)
    monkeypatch.setattr(api, "prosite", prosite)
    stats = S.DatasetStatsManager(http, tmp_path)
    stats._open_stream = lambda q: io.BytesIO(gzip.compress(_fasta(60).encode()))
    stats.start_sync("human_reviewed")
    monkeypatch.setattr(api, "dataset_stats", stats)
    return TestClient(api.app)      # no context manager: start-up downloads are not triggered


def test_health_and_status(client):
    assert client.get("/health").json()["status"] == "ok"
    body = client.get("/api/v1/status").json()
    assert body["prosite"]["state"] == "ready" and 24 in body["page_sizes"]


def test_pagination_has_real_totals_and_distinct_pages(client):
    one = client.get("/api/v1/proteins?page=1&page_size=24").json()
    two = client.get("/api/v1/proteins?page=2&page_size=24").json()
    assert one["total_items"] == 300 and one["total_pages"] == 13
    assert {i["id"] for i in one["items"]}.isdisjoint({i["id"] for i in two["items"]})
    missing = client.get("/api/v1/proteins?page=999&page_size=24")
    assert missing.status_code == 404 and missing.json()["kind"] == "page_out_of_range" and missing.json()["total_pages"] == 13
    assert client.get("/api/v1/proteins?page_size=25").status_code == 422


def test_protein_analysis_and_prosite(client):
    body = client.get("/api/v1/proteins/P04637/analysis").json()
    assert body["protein"]["id"] == "P04637"
    comp = body["analysis"]["composition"]
    assert comp["most_abundant"]["count"] > 0 and comp["least_abundant"]["count"] > 0
    scan = client.get("/api/v1/proteins/P04637/prosite").json()
    assert scan["status"] == "completed" and scan["patterns_tested"] == 2


def test_sequence_lab_pattern_and_validation_messages(client):
    lab = client.post("/api/v1/analyze", json={"sequence": "GGGGG"}).json()
    assert lab["composition"]["most_abundant"]["residues"][0]["name"] == "Glycine"
    hits = client.post("/api/v1/pattern", json={"sequence": "MTKSGGGGSRYLL", "pattern": "G"}).json()
    assert [(m["match"], m["start"], m["end"], m["sequence"]) for m in hits["matches"]][0] == (1, 5, 5, "G")
    bad = client.post("/api/v1/analyze", json={"sequence": "MKT1@"})
    assert bad.status_code == 422 and bad.json()["kind"] == "invalid_input"


def test_compare_explains_invalid_k(client):
    ok = client.post("/api/v1/compare", json={"first": "P04637", "second": "P01308", "kmer_size": 3})
    assert ok.status_code == 200 and ok.json()["similarity"]["k"] == 3
    for k in (0, 25, 99):
        r = client.post("/api/v1/compare", json={"first": "P04637", "second": "P01308", "kmer_size": k})
        assert r.status_code == 422 and len(r.json()["detail"]) > 20 and "failed" not in r.json()["detail"].lower()


def test_unknown_protein_and_upstream_outage(client, monkeypatch):
    assert client.get("/api/v1/proteins/ZZZ999").status_code == 404
    down = FakeHttp([(lambda u, p: True, lambda u, p: S.UpstreamError("unavailable", "timeout"))])
    monkeypatch.setattr(api, "uniprot", S.UniProtClient(down))
    r = client.get("/api/v1/proteins/P04637")
    assert r.status_code == 503 and r.json()["kind"] == "database_unavailable"


def test_reference_dataset_and_frontend(client):
    rows = client.get("/api/v1/reference/amino-acids").json()
    assert len(rows) == 20 and all("formula" in r and "pk1" in r for r in rows)
    trp = client.get("/api/v1/reference/amino-acids/W").json()
    assert trp["structure"]["format"] == "sdf" and "V2000" in trp["structure"]["data"]
    stats = client.get("/api/v1/dataset/stats?scope=human_reviewed").json()
    assert stats["state"] == "ready" and stats["stats"]["total_proteins"] == 60
    page = client.get("/")
    assert page.status_code == 200 and "Protein Atlas" in page.text and "Pattern Lab" in page.text and "structure-modal" in page.text
    assert client.get("/static/js/core.js").status_code == 200
