"""Tests for the scientific engine (no network, no web framework needed)."""
import sys
import types

try:                                    # the CLI uses `rich`; tests do not need it installed
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

import json
import math
import random

import pytest

import atlas_engine as E


def contains_bad(value):
    """True if any string/number in a nested structure is NaN/None-like text."""
    if isinstance(value, float):
        return not math.isfinite(value)
    if isinstance(value, dict):
        return any(contains_bad(v) for v in value.values())
    if isinstance(value, list):
        return any(contains_bad(v) for v in value)
    if isinstance(value, str):
        return value in ("undefined", "NaN", "null")
    return False


MOTSC = "MRWQEMGYIFYPRKLR"


def test_overview_values_are_real_for_a_short_peptide():
    a = E.analyze_sequence(MOTSC)
    comp = a["composition"]
    assert comp["most_abundant"]["count"] == 3
    assert comp["most_abundant"]["residues"][0]["name"] == "Arginine"
    assert comp["least_abundant"]["count"] == 1
    assert {r["code"] for r in comp["least_abundant"]["residues"]} == set("EFGIKLPQW")
    assert comp["families_total"] == len(MOTSC)
    assert abs(sum(f["percent"] for f in comp["families"]) - 100.0) < 1e-9


def test_all_twenty_residues_have_count_and_percent():
    a = E.analyze_sequence(MOTSC)
    rows = [r for r in a["composition"]["rows"] if r["family"] != "Non-standard"]
    assert [r["code"] for r in rows] == list("ACDEFGHIKLMNPQRSTVWY")
    assert all(isinstance(r["count"], int) and r["percent"] is not None for r in rows)
    assert sum(r["count"] for r in rows) == len(MOTSC)


def test_mass_breakdown_reconciles_with_molecular_weight():
    a = E.analyze_sequence(MOTSC)
    mw = a["physicochemical"]["molecular_weight"]
    brk = a["physicochemical"]["mass_breakdown"]
    assert brk["reconciles"] is True
    assert abs(brk["total_da"] - mw["value_da"]) < 0.005
    assert abs(sum(r["mass_da"] for r in brk["rows"]) + brk["water_da"] - mw["value_da"]) < 0.005
    assert all(r["count"] > 0 and r["mass_da"] > 0 for r in brk["rows"])
    assert 2174.0 < mw["value_da"] < 2175.5          # published MOTS-c mass is about 2174.6 Da


def test_sequence_lab_and_protein_pages_share_one_engine():
    assert E.analyze_sequence("GGGGG")["composition"]["most_abundant"]["count"] == 5
    g = E.analyze_sequence("GGGGG")
    assert g["composition"]["least_abundant"]["only_one_type_present"] is True
    assert g["composition"]["most_abundant"]["residues"][0]["code"] == "G"
    assert E.analyze_sequence(MOTSC) == E.analyze_sequence(MOTSC.lower())


def test_results_are_json_safe_and_have_no_bad_values():
    for seq in (MOTSC, "GGGGG", "MKTIIALSYIFCLVFADYKDDDDK", "ACDEFGHIKLMNPQRSTVWYXUB"):
        a = E.analyze_sequence(seq)
        json.dumps(a, allow_nan=False)
        assert not contains_bad(a)


def test_hydropathy_shorter_than_window_is_explained_not_empty():
    a = E.analyze_sequence("GGGGG")
    h = a["hydropathy"]
    assert h["status"] == "no_data" and "shorter than the window" in h["message"]
    assert len(h["per_residue"]) == 5
    long_a = E.analyze_sequence("A" * 60 + "L" * 25 + "K" * 30)
    assert long_a["hydropathy"]["status"] == "success"
    assert long_a["hydropathy"]["membrane_segments"]


def test_invalid_input_is_reported_clearly():
    with pytest.raises(E.InvalidSequence):
        E.analyze_sequence("MKT123!!")
    with pytest.raises(E.InvalidSequence):
        E.analyze_sequence("   ")
    with pytest.raises(E.InvalidSequence):
        E.analyze_sequence("MKTIIALSYIF", window=2)


def test_fasta_and_whitespace_input_is_accepted():
    assert E.normalize_sequence(">sp|P1|X test\nMKT\nIIA\n") == "MKTIIA"
    assert E.normalize_sequence("mkt iia 10\n20") == "MKTIIA"


def test_custom_pattern_table_fields_line_up():
    r = E.search_pattern(MOTSC, "G")
    assert r["count"] == 1
    m = r["matches"][0]
    assert (m["match"], m["start"], m["end"], m["sequence"]) == (1, 7, 7, "G")


def test_prosite_pattern_semantics():
    def hits(pattern, seq):
        return [(h["start"], h["end"]) for h in E._find_hits(seq, E.compile_prosite(pattern))]
    assert hits("N-{P}-[ST]-{P}.", "AANASAANGTAAANPSAA") == [(3, 6), (8, 11)]
    assert hits("[ST]-x-[RK].", "ASAKTTR") == [(2, 4), (5, 7)]              # overlapping allowed
    assert hits("[STAGCN]-[RKH]-[LIVMAFY]>.", "AAAAGKL") == [(5, 7)]        # C-terminal anchor
    assert hits("[STAGCN]-[RKH]-[LIVMAFY]>.", "AAAAGKLA") == []
    assert hits("<M-x(2)-[LK]", "MAALK") == [(1, 4)]
    assert hits("<M-x(2)-[LK]", "AMAALK") == []                             # N-terminal anchor
    assert hits("A-[GS>]", "CCAGCA") == [(3, 4), (6, 6)]                    # '>' inside a class
    assert E.compile_prosite("bad-[") is None


PROSITE_SAMPLE = """CC   Release 2026_03 of 01-Jan-26
ID   ASN_GLYCOSYLATION; PATTERN.
AC   PS00001;
DE   N-glycosylation site.
PA   N-{P}-[ST]-{P}.
CC   /SKIP-FLAG=TRUE;
//
ID   SOME_PROFILE; MATRIX.
AC   PS50000;
DE   profile record (no pattern).
//
ID   PKC_PHOSPHO_SITE; PATTERN.
AC   PS00005;
DE   Protein kinase C phosphorylation site.
PA   [ST]-x-[RK].
//
ID   BROKEN; PATTERN.
AC   PS99999;
DE   broken.
PA   Z-Z-(
//
"""


def test_prosite_scan_reports_what_it_really_did():
    db = E.PrositeDatabase.from_text(PROSITE_SAMPLE, "test")
    info = db.info()
    assert info["release"] == "2026_03"
    assert info["patterns_usable"] == 2 and info["patterns_unsupported"] == 1
    full = db.scan("AANASAAAKSAR")
    assert full["status"] == "completed" and full["patterns_tested"] == 2
    assert full["patterns_matched"] == 2 and full["total_hits"] == 2
    none = db.scan("GGGGG")
    assert none["status"] == "completed" and none["patterns_matched"] == 0
    assert "2 PROSITE patterns tested" in none["message"]
    filtered = db.scan("AANASAAAKSAR", exclude_frequent=True)
    assert filtered["patterns_tested"] == 1


def test_kmer_comparison_validation_and_meaning():
    other = "MKTIIALSYIFCLVFADYKDDDDK"
    with pytest.raises(E.InvalidSequence, match="limited to 20"):
        E.compare_sequences(MOTSC, other, 21)
    with pytest.raises(E.InvalidSequence, match="shorter sequence"):
        E.compare_sequences(MOTSC, other, 17)
    with pytest.raises(E.InvalidSequence, match="at least 1"):
        E.compare_sequences(MOTSC, other, 0)
    r = E.compare_sequences(MOTSC, other, 5)
    assert r["shared_kmers"] == 0 and "genuine result" in r["interpretation"]
    same = E.compare_sequences(other, other, 5)
    assert same["word_similarity_percent"] == pytest.approx(100.0)
    assert same["shared_kmers"] == same["union_kmers"]
    partial = E.compare_sequences("AAAAKLMNPQ", "KLMNPQRSTV", 3)
    assert partial["shared_kmers"] > 0
    assert partial["profile"][0]["k"] == 1


def test_reference_data_is_complete_and_consistent():
    assert len(E.AMINO_ACID_REFERENCE) == 20
    for row in E.AMINO_ACID_REFERENCE:
        assert row["formula_average_mass_da"] == pytest.approx(row["free_mass_da"], abs=0.03)
        assert row["monoisotopic_free_mass_da"] < row["free_mass_da"] + 0.5
        assert row["smiles"] and row["pubchem_cid"] and row["three"] and row["polarity"]
        assert row["residue_mass_da"] == pytest.approx(row["free_mass_da"] - 18.015, abs=0.01)


def test_dataset_aggregator_covers_every_sequence_and_matches_engine():
    random.seed(7)
    seqs = ["".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(random.randint(20, 800)))
            for _ in range(300)]
    agg = E.DatasetAggregator()
    for s in seqs:
        agg.add(s)
    stats = agg.finish()
    assert stats["total_proteins"] == 300
    assert stats["length"]["min"] == min(map(len, seqs)) and stats["length"]["max"] == max(map(len, seqs))
    assert sum(stats["length"]["histogram"]["counts"]) == 300
    assert stats["length"]["median"] == sorted(map(len, seqs))[149] * 0.5 + sorted(map(len, seqs))[150] * 0.5
    mean_pi = sum(E.analyze_sequence(s)["physicochemical"]["isoelectric_point"]["value"] for s in seqs[:40]) / 40
    small = E.DatasetAggregator()
    for s in seqs[:40]:
        small.add(s)
    assert small.finish()["isoelectric_point"]["mean"] == pytest.approx(mean_pi, abs=1e-6)
    assert not contains_bad(stats)
    assert sum(r["percent"] for r in stats["amino_acids"]) == pytest.approx(100.0)


def test_compare_bundle_uses_the_analysis_engine_and_a_dot_plot():
    a = "MKTIIALSYIFCLVFADYKDDDDKMKTIIALSYIFC"
    b = "GGGGSRYLLMKTIIALSYIFCGGGG"
    bundle = E.compare_proteins(a, b, 4)
    props = {p["key"]: p for p in bundle["properties"]}
    assert props["length"]["first"] == len(a) and props["length"]["second"] == len(b)
    solo = E.analyze_sequence(a)
    assert props["molecular_weight"]["first"] == solo["physicochemical"]["molecular_weight"]["value_da"]
    assert props["gravy"]["first"] == solo["physicochemical"]["gravy"]
    assert len(bundle["composition"]) == 20 and all(isinstance(r["first_percent"], float) for r in bundle["composition"])
    assert bundle["similarity"]["shared_kmers"] > 0
    dots = bundle["dotplot"]
    assert dots["points"] and all(a[p["y"] - 1:p["y"] - 1 + dots["k"]] == b[p["x"] - 1:p["x"] - 1 + dots["k"]]
                                  for p in dots["points"])


def test_compare_bundle_rejects_bad_k_with_specific_message():
    import pytest
    with pytest.raises(E.InvalidSequence, match="shorter sequence"):
        E.compare_proteins("MKTIIALSYIFC", "MKTIIAL", 9)
    with pytest.raises(E.InvalidSequence, match="limited to 20"):
        E.compare_proteins("ACDEFGHIKLMNPQRSTVWY" * 3, "ACDEFGHIKLMNPQRSTVWY" * 3, 25)


def test_dotplot_densifies_low_complexity_input_and_reports_it():
    dots = E.dotplot_points("A" * 200, "A" * 200, 2)
    assert dots["k"] > 2 and dots["truncated"] in (True, False)
    assert len(dots["points"]) <= 6000
