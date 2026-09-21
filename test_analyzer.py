import numpy as np

import protein_analyzer as analyzer


class FakeResponse:
    def __init__(self, payload):
        self.payload = payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False

    def read(self):
        import json
        return json.dumps(self.payload).encode("utf-8")


def test_empty_analysis_helpers_return_safe_results():
    assert analyzer.aromaticity("") == 0.0
    assert analyzer.composition_table("").empty
    assert analyzer.class_summary_table("").empty
    assert analyzer.dataset_composition([])["total_residues"] == 0


def test_log_histogram_ignores_non_positive_values(capsys):
    analyzer.display_histogram(
        np.array([0.0, -1.0]), "Empty", "Value", "#000000", log_x=True)

    assert "No positive data available" in capsys.readouterr().out


def test_remote_search_returns_paginated_records(monkeypatch):
    payload = {
        "totalResults": 2_000_000,
        "results": [{
            "primaryAccession": "P04637",
            "uniProtkbId": "P53_HUMAN",
            "proteinDescription": {
                "recommendedName": {
                    "fullName": {"value": "Cellular tumor antigen p53"}
                }
            },
            "genes": [{"geneName": {"value": "TP53"}}],
            "organism": {"scientificName": "Homo sapiens"},
            "sequence": {"value": "MKTIIALSYIFCLVFADYKDDDDK"},
        }],
    }
    monkeypatch.setattr(analyzer, "urlopen", lambda request, timeout: FakeResponse(payload))

    records, total = analyzer.fetch_remote_proteins("*", 0, 100)

    assert total == 2_000_000
    assert records[0]["id"] == "P04637"
    assert records[0]["gene"] == "TP53"


def test_load_more_remote_proteins_appends_next_page(monkeypatch):
    analyzer.REMOTE_QUERY = "*"
    analyzer.REMOTE_PAGE_SIZE = 2
    analyzer.REMOTE_OFFSET = 1
    analyzer.REMOTE_TOTAL = 3
    proteins = [analyzer.build_protein_record(
        "sp|P04637|P53_HUMAN TP53 OS=Homo sapiens GN=TP53",
        "MKT")]

    monkeypatch.setattr(
        analyzer,
        "fetch_remote_proteins",
        lambda query, offset, size: (
            [analyzer.build_protein_record(
                "sp|Q99999|TEST_HUMAN Test protein OS=Homo sapiens GN=TEST",
                "MKT")],
            3,
        ),
    )

    analyzer.load_more_remote_proteins(proteins)

    assert len(proteins) == 2
    assert analyzer.REMOTE_OFFSET == 2