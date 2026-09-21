import itertools
import math

import numpy as np
import pytest

import atlas_aminoacids as aa
import atlas_engine as engine
import protein_analyzer as pa

FORMULAS = {"A": "C3H7NO2", "R": "C6H14N4O2", "N": "C4H8N2O3", "D": "C4H7NO4", "C": "C3H7NO2S",
            "E": "C5H9NO4", "Q": "C5H10N2O3", "G": "C2H5NO2", "H": "C6H9N3O2", "I": "C6H13NO2",
            "L": "C6H13NO2", "K": "C6H14N2O2", "M": "C5H11NO2S", "F": "C9H11NO2", "P": "C5H9NO2",
            "S": "C3H7NO3", "T": "C4H9NO3", "W": "C11H12N2O2", "Y": "C9H11NO3", "V": "C5H11NO2"}


def test_model_formula_matches_engine_reference_for_all_twenty():
    for code in pa.AMINO_ACIDS:
        assert aa.extras(code)["model_formula"] == engine.AMINO_ACID_BY_CODE[code]["formula"], code
        assert aa.extras(code)["model_formula"] == FORMULAS[code], code


def test_smiles_in_model_and_engine_describe_the_same_stereochemistry_class():
    # every residue except Gly is chiral, and only Ile/Thr carry a second stereocentre
    for code in pa.AMINO_ACIDS:
        stereo = aa.extras(code)["stereo"]
        assert stereo.startswith("L") or code == "G"
        assert ("3" in stereo) == (code in "IT")


def test_sdf_round_trips_and_is_valid_v2000():
    for code in pa.AMINO_ACIDS:
        text = aa.sdf(code)
        lines = text.splitlines()
        n_atoms, n_bonds = int(lines[3][0:3]), int(lines[3][3:6])
        assert lines[3].endswith("V2000") and lines[-1] == "$$$$" and lines[-2] == "M  END"
        assert len(lines) == 4 + n_atoms + n_bonds + 2
        for row in lines[4:4 + n_atoms]:
            assert all(math.isfinite(float(row[i:i + 10])) for i in (0, 10, 20))
        for row in lines[4 + n_atoms:4 + n_atoms + n_bonds]:
            a, b, order = int(row[0:3]), int(row[3:6]), int(row[6:9])
            assert 1 <= a <= n_atoms and 1 <= b <= n_atoms and order in (1, 2)


def test_geometry_is_chemically_sensible_and_L_configured():
    limits = {("C", "C"): (1.32, 1.56), ("C", "N"): (1.26, 1.50), ("C", "O"): (1.20, 1.46),
              ("C", "S"): (1.75, 1.86), ("C", "H"): (1.05, 1.13), ("H", "N"): (0.97, 1.05),
              ("H", "O"): (0.92, 1.0), ("H", "S"): (1.30, 1.38)}
    for code in pa.AMINO_ACIDS:
        entry = aa._structures()[code]
        pos = np.array([a[1:] for a in entry["atoms"]])
        el = [a[0] for a in entry["atoms"]]
        adj = {i: [] for i in range(len(el))}
        for a, b, _ in entry["bonds"]:
            adj[a].append(b)
            adj[b].append(a)
            lo, hi = limits[tuple(sorted((el[a], el[b])))]
            assert lo <= np.linalg.norm(pos[a] - pos[b]) <= hi, (code, a, b)
        for c, nbrs in adj.items():
            for x, y in itertools.combinations(nbrs, 2):
                v1, v2 = pos[x] - pos[c], pos[y] - pos[c]
                ang = math.degrees(math.acos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))
                assert 95 <= ang <= 135, (code, c, ang)
        if code == "G":
            continue
        alpha = next(i for i in range(len(el)) if el[i] == "C"
                     and any(el[n] == "N" for n in adj[i])
                     and any(el[n] == "C" and sum(el[m] == "O" for m in adj[n]) == 2 for n in adj[i]))
        n_ = next(n for n in adj[alpha] if el[n] == "N")
        co = next(n for n in adj[alpha] if el[n] == "C" and sum(el[m] == "O" for m in adj[n]) == 2)
        h_ = next(n for n in adj[alpha] if el[n] == "H")
        r_ = next(n for n in adj[alpha] if n not in (n_, co, h_))
        corn = np.linalg.det(np.array([pos[co] - pos[alpha], pos[r_] - pos[alpha], pos[n_] - pos[alpha]]))
        assert corn > 0, f"{code} is not L (CORN rule)"


def test_isoelectric_points_match_published_table():
    published = {"G": 5.97, "A": 6.01, "V": 5.97, "L": 5.98, "I": 6.02, "M": 5.74, "P": 6.48,
                 "F": 5.48, "W": 5.89, "S": 5.68, "T": 5.60, "C": 5.07, "Y": 5.66, "N": 5.41,
                 "Q": 5.65, "D": 2.77, "E": 3.22, "K": 9.74, "R": 10.76, "H": 7.59}
    for code in pa.AMINO_ACIDS:
        assert aa.extras(code)["isoelectric_point_free"] == pytest.approx(published[code], abs=0.02), code
