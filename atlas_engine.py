"""Protein Atlas - the single, validated sequence-analysis engine.

Every screen in the web app (protein pages, Sequence Lab, Compare, Export)
gets its numbers from this one module, so the same sequence can never
produce contradictory results in two places.

Design rules
------------
* Every section of a result carries an explicit ``status``:
  ``success`` | ``no_data`` | ``invalid_input`` | ``error``.
* Values are always JSON-safe: NaN / infinity become ``None`` and the
  front end turns ``None`` into an honest "Data unavailable" message.
* Constants and the core formulas (pI, charge, GRAVY, ...) come from
  ``protein_analyzer`` so the CLI and the web app cannot drift apart.

This module has no web-framework dependency, so it is easy to test.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from typing import Any, Iterable, Iterator, Optional

import numpy as np

import protein_analyzer as pa

ENGINE_VERSION = "3.0"
KMER_LIMIT = 20          # documented upper limit for k-mer comparison
MAX_SEQUENCE_LENGTH = 1_000_000

# --------------------------------------------------------------------
# Reference chemistry for the 20 canonical amino acids
# --------------------------------------------------------------------

_ATOM_AVG = {"C": 12.011, "H": 1.008, "N": 14.007, "O": 15.999, "S": 32.06}
_ATOM_MONO = {"C": 12.0, "H": 1.00782503207, "N": 14.0030740048,
              "O": 15.99491461956, "S": 31.97207100}

# code: (three-letter, free-amino-acid formula, SMILES (L-form),
#        PubChem CID, side-chain H-bond role, donor atoms, acceptor atoms,
#        side-chain pKa or None, ring type, one-line note)
_AA_TABLE: dict[str, tuple] = {
    "A": ("Ala", "C3H7NO2", "C[C@H](N)C(=O)O", 5950,
          "none", [], [], None, "Aliphatic",
          "Smallest chiral residue; a methyl side chain."),
    "R": ("Arg", "C6H14N4O2", "N[C@@H](CCCNC(=N)N)C(=O)O", 6322,
          "donor", ["NE", "NH1", "NH2"], [], 12.5, "Aliphatic chain + guanidinium",
          "Guanidinium group is protonated across the physiological pH range."),
    "N": ("Asn", "C4H8N2O3", "N[C@@H](CC(N)=O)C(=O)O", 6267,
          "donor and acceptor", ["ND2"], ["OD1"], None, "Non-aromatic",
          "Amide side chain; the usual N-glycosylation attachment site."),
    "D": ("Asp", "C4H7NO4", "N[C@@H](CC(=O)O)C(=O)O", 5960,
          "acceptor", [], ["OD1", "OD2"], 3.9, "Non-aromatic",
          "Carboxylate side chain; often binds metal ions and active-site protons."),
    "C": ("Cys", "C3H7NO2S", "N[C@@H](CS)C(=O)O", 5862,
          "weak donor / weak acceptor", ["SG"], ["SG"], 8.5, "Non-aromatic",
          "Thiol side chain; pairs of Cys can form disulfide bonds."),
    "Q": ("Gln", "C5H10N2O3", "N[C@@H](CCC(N)=O)C(=O)O", 5961,
          "donor and acceptor", ["NE2"], ["OE1"], None, "Non-aromatic",
          "Amide side chain, one methylene longer than Asn."),
    "E": ("Glu", "C5H9NO4", "N[C@@H](CCC(=O)O)C(=O)O", 33032,
          "acceptor", [], ["OE1", "OE2"], 4.1, "Non-aromatic",
          "Carboxylate side chain, one methylene longer than Asp."),
    "G": ("Gly", "C2H5NO2", "NCC(=O)O", 750,
          "none", [], [], None, "No side chain (H)",
          "Achiral; gives backbones exceptional flexibility."),
    "H": ("His", "C6H9N3O2", "N[C@@H](Cc1c[nH]cn1)C(=O)O", 6274,
          "donor and acceptor", ["ND1", "NE2"], ["ND1", "NE2"], 6.5,
          "Aromatic (imidazole)",
          "Imidazole ring has a pKa near physiological pH, so it is a common catalytic residue."),
    "I": ("Ile", "C6H13NO2", "CC[C@H](C)[C@H](N)C(=O)O", 6306,
          "none", [], [], None, "Aliphatic",
          "Beta-branched hydrophobic residue with a second stereocentre."),
    "L": ("Leu", "C6H13NO2", "CC(C)C[C@H](N)C(=O)O", 6106,
          "none", [], [], None, "Aliphatic",
          "Very common hydrophobic residue in helices and membrane segments."),
    "K": ("Lys", "C6H14N2O2", "NCCCC[C@H](N)C(=O)O", 5962,
          "donor", ["NZ"], [], 10.8, "Aliphatic chain + amine",
          "Long flexible side chain ending in a positively charged amine."),
    "M": ("Met", "C5H11NO2S", "CSCC[C@H](N)C(=O)O", 6137,
          "weak acceptor", [], ["SD"], None, "Aliphatic (thioether)",
          "Thioether side chain; the usual start residue of translation."),
    "F": ("Phe", "C9H11NO2", "N[C@@H](Cc1ccccc1)C(=O)O", 6140,
          "none", [], [], None, "Aromatic",
          "Benzyl side chain; buried in hydrophobic cores."),
    "P": ("Pro", "C5H9NO2", "OC(=O)[C@@H]1CCCN1", 145742,
          "none", [], [], None, "Aliphatic (cyclic)",
          "Side chain loops back onto the backbone nitrogen, restricting rotation."),
    "S": ("Ser", "C3H7NO3", "N[C@@H](CO)C(=O)O", 5951,
          "donor and acceptor", ["OG"], ["OG"], None, "Non-aromatic",
          "Small hydroxyl side chain; frequent phosphorylation site."),
    "T": ("Thr", "C4H9NO3", "C[C@@H](O)[C@H](N)C(=O)O", 6288,
          "donor and acceptor", ["OG1"], ["OG1"], None, "Non-aromatic",
          "Beta-branched hydroxyl residue with a second stereocentre."),
    "W": ("Trp", "C11H12N2O2", "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", 6305,
          "donor", ["NE1"], [], None, "Aromatic (indole)",
          "Largest residue; dominates absorbance at 280 nm."),
    "Y": ("Tyr", "C9H11NO3", "N[C@@H](Cc1ccc(O)cc1)C(=O)O", 6057,
          "donor and acceptor", ["OH"], ["OH"], 10.1, "Aromatic",
          "Phenol side chain; absorbs at 280 nm and is a phosphorylation site."),
    "V": ("Val", "C5H11NO2", "CC(C)[C@H](N)C(=O)O", 6287,
          "none", [], [], None, "Aliphatic",
          "Beta-branched hydrophobic residue."),
}

_POLARITY = {
    "Nonpolar": "Nonpolar (hydrophobic)",
    "Aromatic": "Nonpolar / aromatic",
    "Polar": "Polar, uncharged",
    "Acidic": "Charged, negative (acidic)",
    "Basic": "Charged, positive (basic)",
}


def _parse_formula(formula: str) -> dict[str, int]:
    return {el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula)}


def _formula_string(atoms: dict[str, int]) -> str:
    parts = []
    for el in ("C", "H", "N", "O", "S"):
        if atoms.get(el):
            parts.append(el if atoms[el] == 1 else f"{el}{atoms[el]}")
    return "".join(parts)


def _mass(atoms: dict[str, int], table: dict[str, float]) -> float:
    return sum(table[el] * n for el, n in atoms.items())


def _build_reference() -> list[dict[str, Any]]:
    rows = []
    for code in pa.AMINO_ACIDS:
        (three, formula, smiles, cid, hb_role, donors, acceptors, pka,
         ring, note) = _AA_TABLE[code]
        atoms = _parse_formula(formula)
        residue_atoms = dict(atoms)
        residue_atoms["H"] -= 2
        residue_atoms["O"] -= 1
        family = pa.AMINO_ACID_CLASSES[code]
        charge = {"Acidic": -1, "Basic": +1}.get(family, 0)
        if code == "H":
            charge_text = "About +0.1 at pH 7 (partially protonated)"
        elif charge == -1:
            charge_text = "-1 at pH 7"
        elif charge == 1:
            charge_text = "+1 at pH 7"
        else:
            charge_text = "0 at pH 7 (neutral)"
        rows.append({
            "code": code,
            "three": three,
            "name": pa.AMINO_ACID_NAMES[code],
            "family": family,
            "side_chain_class": family,
            "polarity": _POLARITY[family],
            "charge_ph7": charge,
            "charge_text": charge_text,
            "side_chain_pka": pka,
            "formula": formula,
            "residue_formula": _formula_string(residue_atoms),
            "free_mass_da": pa.AMINO_ACID_MASSES[code],
            "residue_mass_da": round(pa.AMINO_ACID_MASSES[code] - pa.WATER_MASS, 3),
            "monoisotopic_free_mass_da": round(_mass(atoms, _ATOM_MONO), 4),
            "monoisotopic_residue_mass_da": round(_mass(residue_atoms, _ATOM_MONO), 4),
            "formula_average_mass_da": round(_mass(atoms, _ATOM_AVG), 2),
            "hydropathy": pa.HYDROPATHY_SCALE[code],
            "ring_type": ring,
            "aromatic": code in pa.AROMATIC_RESIDUES or code == "H",
            "hbond_role": hb_role,
            "hbond_donor_atoms": donors,
            "hbond_acceptor_atoms": acceptors,
            "smiles": smiles,
            "pubchem_cid": cid,
            "pubchem_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
            "note": note,
        })
    return rows


AMINO_ACID_REFERENCE: list[dict[str, Any]] = _build_reference()
AMINO_ACID_BY_CODE: dict[str, dict[str, Any]] = {r["code"]: r for r in AMINO_ACID_REFERENCE}
THREE_LETTER = {code: row["three"] for code, row in AMINO_ACID_BY_CODE.items()}
THREE_LETTER.update({"U": "Sec", "O": "Pyl", "X": "Xaa", "B": "Asx", "Z": "Glx", "J": "Xle"})

CLASSIFICATION_NOTE = (
    "Families: Nonpolar = G A V L I M P; Aromatic = F W Y; Polar (uncharged) = S T C N Q; "
    "Acidic = D E; Basic = K R H. Non-standard codes (U O X B Z J) are reported separately."
)
MASS_CONVENTION = (
    "Average (not monoisotopic) masses. Protein mass = sum of residue masses "
    "(free amino acid mass minus one water, 18.015 Da) plus one water for the termini."
)


class InvalidSequence(ValueError):
    """Raised when input cannot be analysed as a protein sequence."""


# --------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------

def clean_json(value: Any) -> Any:
    """Recursively replace NaN / infinity with None so JSON is always valid."""
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, (np.floating,)):
        v = float(value)
        return v if math.isfinite(v) else None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, np.ndarray):
        return clean_json(value.tolist())
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(v) for v in value]
    return value


def normalize_sequence(raw: str) -> str:
    """Accept plain sequence, FASTA text, or sequence with spaces/digits."""
    if raw is None:
        raise InvalidSequence("No sequence was provided.")
    lines = [ln.strip() for ln in str(raw).splitlines()]
    lines = [ln for ln in lines if ln and not ln.startswith(">")]
    text = "".join(lines)
    text = re.sub(r"[\s\d]+", "", text).upper()
    if text.endswith("*"):
        text = text[:-1]
    if not text:
        raise InvalidSequence("The sequence is empty.")
    if len(text) > MAX_SEQUENCE_LENGTH:
        raise InvalidSequence(f"Sequences are limited to {MAX_SEQUENCE_LENGTH:,} residues.")
    return text


def _pct(part: float, whole: float) -> Optional[float]:
    return (part / whole * 100.0) if whole else None


def _residue_ref(code: str) -> dict[str, str]:
    return {"code": code, "name": pa.AMINO_ACID_NAMES.get(code, code),
            "three": THREE_LETTER.get(code, code)}


# --------------------------------------------------------------------
# main analysis
# --------------------------------------------------------------------

def analyze_sequence(raw_sequence: str, window: int = pa.DEFAULT_HYDROPATHY_WINDOW) -> dict[str, Any]:
    """Run the complete validated analysis of one sequence.

    Raises InvalidSequence for unusable input. Individual sections that cannot be
    computed report ``status`` = ``no_data`` or ``error`` rather than raising.
    """
    sequence = normalize_sequence(raw_sequence)
    length = len(sequence)
    window = int(window)
    if window < 3 or window > 101:
        raise InvalidSequence("Hydropathy window must be between 3 and 101 residues.")

    counts = pa.amino_acid_counts(sequence)
    unknown = pa.unrecognised_count(sequence, counts)
    if unknown > 0:
        strange = sorted(set(sequence) - set(pa.KNOWN_CODES))
        raise InvalidSequence(
            f"Contains {unknown} unexpected character(s): {', '.join(strange)}. "
            "Use single-letter amino-acid codes only (20 standard plus U, O, X, B, Z, J).")

    result: dict[str, Any] = {
        "schema": "atlas.analysis/3",
        "engine_version": ENGINE_VERSION,
        "status": "success",
    }

    # ---- sequence -----------------------------------------------------
    extras = [
        {"code": c, "name": pa.AMINO_ACID_NAMES[c], "three": THREE_LETTER[c], "count": counts[c]}
        for c in pa.NON_STANDARD_CODES if counts[c] > 0]
    result["sequence"] = {
        "status": "success",
        "residues": sequence,
        "length": length,
        "validation": {
            "valid": True,
            "message": (
                "Valid sequence including non-standard code(s): "
                + ", ".join(f"{e['code']} x{e['count']}" for e in extras)
                if extras else "Valid protein sequence (20 standard amino acids)."),
            "non_standard": extras,
        },
    }

    # ---- composition ---------------------------------------------------
    rows = []
    for code in pa.AMINO_ACIDS:
        rows.append({
            "code": code, "three": THREE_LETTER[code], "name": pa.AMINO_ACID_NAMES[code],
            "family": pa.AMINO_ACID_CLASSES[code],
            "count": counts[code], "percent": _pct(counts[code], length),
            "residue_mass_da": round(pa.AMINO_ACID_MASSES[code] - pa.WATER_MASS, 3),
        })
    for e in extras:
        rows.append({
            "code": e["code"], "three": e["three"], "name": e["name"],
            "family": "Non-standard", "count": e["count"],
            "percent": _pct(e["count"], length),
            "residue_mass_da": round(pa.AMINO_ACID_MASSES[e["code"]] - pa.WATER_MASS, 3),
        })

    standard_rows = [r for r in rows if r["family"] != "Non-standard"]
    max_count = max((r["count"] for r in standard_rows), default=0)
    present = [r for r in standard_rows if r["count"] > 0]
    min_count = min((r["count"] for r in present), default=0)
    most = ({
        "count": max_count, "percent": _pct(max_count, length),
        "residues": [_residue_ref(r["code"]) for r in standard_rows if r["count"] == max_count],
    } if max_count > 0 else None)
    least = ({
        "count": min_count, "percent": _pct(min_count, length),
        "residues": [_residue_ref(r["code"]) for r in present if r["count"] == min_count],
        "only_one_type_present": len(present) == 1,
    } if present else None)

    families = []
    for cls in pa.CLASS_ORDER:
        members = [c for c in pa.AMINO_ACIDS if pa.AMINO_ACID_CLASSES[c] == cls]
        total = sum(counts[c] for c in members)
        families.append({"family": cls, "residues": "".join(members),
                         "count": total, "percent": _pct(total, length)})
    ns_total = sum(e["count"] for e in extras)
    if ns_total:
        families.append({"family": "Non-standard",
                         "residues": "".join(e["code"] for e in extras),
                         "count": ns_total, "percent": _pct(ns_total, length)})

    result["composition"] = {
        "status": "success",
        "rows": rows,
        "most_abundant": most,
        "least_abundant": least,
        "absent": [r["code"] for r in standard_rows if r["count"] == 0],
        "families": families,
        "families_total": sum(f["count"] for f in families),
        "classification_note": CLASSIFICATION_NOTE,
        "distinct_standard_residues": len(present),
    }

    # ---- physicochemical -----------------------------------------------
    weight = pa.molecular_weight_details(sequence)
    phys: dict[str, Any] = {"status": "success"}
    if weight is None:
        phys["molecular_weight"] = {"status": "no_data",
                                    "message": "No residues with a known mass."}
        phys["mass_breakdown"] = {"status": "no_data", "rows": []}
    else:
        brk = []
        for code in pa.KNOWN_CODES:
            n = counts[code]
            if n == 0:
                continue
            residue_mass = pa.AMINO_ACID_MASSES[code] - pa.WATER_MASS
            brk.append({
                "code": code, "three": THREE_LETTER[code], "name": pa.AMINO_ACID_NAMES[code],
                "count": n, "residue_mass_da": residue_mass,
                "mass_da": n * residue_mass,
                "percent": _pct(n * residue_mass, weight["weight"]),
                "approximate": code in pa.APPROXIMATE_CODES,
            })
        brk.sort(key=lambda r: r["mass_da"], reverse=True)
        summed = sum(r["mass_da"] for r in brk)
        total = summed + pa.WATER_MASS
        phys["molecular_weight"] = {
            "status": "success",
            "value_da": weight["weight"],
            "value_kda": weight["weight"] / 1000.0,
            "residue_sum_da": weight["residue_subtotal"],
            "water_da": pa.WATER_MASS,
            "residues_counted": weight["residues_used"],
            "approximated_residues": weight["approximated"],
            "skipped_residues": weight["skipped"],
            "convention": MASS_CONVENTION,
        }
        phys["mass_breakdown"] = {
            "status": "success",
            "rows": brk,
            "residue_sum_da": summed,
            "water_da": pa.WATER_MASS,
            "total_da": total,
            "reconciles": abs(total - weight["weight"]) < 0.005,
            "difference_da": total - weight["weight"],
            "convention": MASS_CONVENTION,
        }

    pi = pa.isoelectric_point(sequence)
    phys["isoelectric_point"] = {"status": "success", "value": pi,
                                 "method": "Bisection on the Henderson-Hasselbalch charge equation",
                                 "pka_set": "EMBOSS"}
    phys["net_charge_ph7"] = pa.net_charge_at_ph(sequence, 7.0)
    acidic = sum(counts[c] for c in pa.ACIDIC_RESIDUES)
    basic = sum(counts[c] for c in pa.BASIC_RESIDUES)
    uncharged = sum(counts[c] for c in pa.UNCHARGED_RESIDUES)
    phys["charge_composition"] = {
        "acidic_count": acidic, "basic_count": basic,
        "charged_count": acidic + basic, "uncharged_count": uncharged,
        "other_count": length - acidic - basic - uncharged,
        "acidic_percent": _pct(acidic, length), "basic_percent": _pct(basic, length),
        "charged_percent": _pct(acidic + basic, length),
        "uncharged_percent": _pct(uncharged, length),
        "acidic_residues": "".join(pa.ACIDIC_RESIDUES),
        "basic_residues": "".join(pa.BASIC_RESIDUES),
    }
    phys["charge_curve"] = [
        {"ph": round(float(ph), 2), "charge": pa.net_charge_at_ph(sequence, float(ph))}
        for ph in np.arange(0.0, 14.01, 0.25)]
    ext = pa.extinction_coefficient(sequence)
    phys["extinction_coefficient"] = {
        "reduced": ext["reduced"], "oxidised": ext["oxidised"],
        "cystines": ext["cystines"], "unit": "M-1 cm-1 at 280 nm"}
    phys["aromaticity_percent"] = pa.aromaticity(sequence)
    ala, val, ile, leu = (counts[c] / length for c in "AVIL")
    phys["aliphatic_index"] = (ala + 2.9 * val + 3.9 * (ile + leu)) * 100.0
    phys["gravy"] = pa.gravy_score(sequence)
    result["physicochemical"] = phys

    # ---- hydropathy ----------------------------------------------------
    per_residue = [pa.HYDROPATHY_SCALE.get(c, 0.0) for c in sequence]
    centres, averages = pa.sliding_hydropathy(sequence, window)
    hydro: dict[str, Any] = {
        "scale": "Kyte-Doolittle", "window": window, "gravy": phys["gravy"],
        "membrane_window": pa.TRANSMEMBRANE_WINDOW,
        "membrane_threshold": pa.TRANSMEMBRANE_THRESHOLD,
        "per_residue": per_residue,
    }
    if averages is None:
        hydro["status"] = "no_data"
        hydro["message"] = (f"The sequence ({length} aa) is shorter than the window "
                            f"({window} aa), so no sliding profile can be drawn. "
                            "The per-residue values are shown instead.")
        hydro["profile"] = None
    else:
        hydro["status"] = "success"
        hydro["profile"] = {"positions": centres.tolist(), "values": averages.tolist()}
    segments = []
    for start, end in pa.find_hydrophobic_segments(
            sequence, pa.TRANSMEMBRANE_WINDOW, pa.TRANSMEMBRANE_THRESHOLD):
        mean = float(np.mean(per_residue[start - 1:end]))
        segments.append({"start": start, "end": end, "length": end - start + 1, "mean": mean})
    hydro["membrane_segments"] = segments
    hydro["membrane_status"] = (
        "success" if length >= pa.TRANSMEMBRANE_WINDOW else "no_data")
    if length < pa.TRANSMEMBRANE_WINDOW:
        hydro["membrane_message"] = (
            f"Membrane-segment detection needs at least {pa.TRANSMEMBRANE_WINDOW} residues.")
    result["hydropathy"] = hydro

    # ---- features (repeats + residue positions) ----------------------------
    result["features"] = {
        "status": "success",
        "repeats": repeat_analysis(sequence, 3),
        "residue_positions": {
            code: [i + 1 for i, ch in enumerate(sequence) if ch == code]
            for code in pa.KNOWN_CODES if counts[code] > 0},
    }
    return clean_json(result)


def repeat_analysis(sequence: str, size: int = 3) -> dict[str, Any]:
    if len(sequence) < size:
        return {"status": "no_data", "size": size, "top_kmers": [], "longest_run": None,
                "message": f"The sequence is shorter than {size} residues."}
    counts = Counter(sequence[i:i + size] for i in range(len(sequence) - size + 1))
    top = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:10]
    best_len, best_code, best_start = 1, sequence[0], 0
    run_len, run_start = 1, 0
    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i - 1]:
            run_len += 1
        else:
            run_len, run_start = 1, i
        if run_len > best_len:
            best_len, best_code, best_start = run_len, sequence[i], run_start
    return {
        "status": "success", "size": size,
        "total_windows": len(sequence) - size + 1,
        "distinct_windows": len(counts),
        "top_kmers": [{"unit": u, "occurrences": n} for u, n in top],
        "longest_run": {"residue": best_code, "length": best_len,
                        "start": best_start + 1, "end": best_start + best_len},
    }


# --------------------------------------------------------------------
# custom pattern search (literal fragment or PROSITE-style pattern)
# --------------------------------------------------------------------

def compile_prosite(pattern: str) -> Optional[re.Pattern]:
    """Compile a PROSITE pattern to a look-ahead regex (finds overlapping hits).

    Supports x, [ABC], {ABC}, repeats x(2), x(2,4), N-terminal '<' and C-terminal '>',
    including '>' inside a class such as [ILV>] (residue OR the C-terminus).
    """
    text = pattern.strip().replace(" ", "")
    if text.endswith("."):
        text = text[:-1]
    if not text:
        return None
    parts: list[str] = []
    n_term = text.startswith("<")
    if n_term:
        text = text[1:]
    # a trailing '>' outside brackets is the C-terminal anchor ('[ILV>]' ends with ']')
    c_term = text.endswith(">")
    if c_term:
        text = text[:-1]
    for token in text.split("-"):
        if not token:
            return None
        rep = ""
        m = re.fullmatch(r"(.+?)\((\d+)(?:,(\d+))?\)", token)
        if m:
            token = m.group(1)
            rep = "{%s}" % m.group(2) if m.group(3) is None else "{%s,%s}" % (m.group(2), m.group(3))
        token = token.upper()
        if token == "X":
            parts.append("." + rep)
        elif re.fullmatch(r"\[[A-Z<>]+\]", token):
            body = token[1:-1]
            letters = body.replace("<", "").replace(">", "")
            alts = []
            if letters:
                alts.append("[" + letters + "]")
            if ">" in body:
                alts.append("$")
            if "<" in body:
                alts.append("^")
            if not alts:
                return None
            piece = alts[0] if len(alts) == 1 else "(?:" + "|".join(alts) + ")"
            parts.append(piece + rep)
        elif re.fullmatch(r"\{[A-Z]+\}", token):
            parts.append("[^" + token[1:-1] + "]" + rep)
        elif re.fullmatch(r"[A-Z]", token):
            parts.append(token + rep)
        else:
            return None
    regex = "".join(parts)
    if n_term:
        regex = "^" + regex
    if c_term:
        regex += "$"
    try:
        return re.compile("(?=(" + regex + "))")
    except re.error:
        return None


def _find_hits(sequence: str, regex: re.Pattern) -> list[dict[str, Any]]:
    hits = []
    for m in regex.finditer(sequence):
        text = m.group(1)
        if not text:
            continue
        hits.append({"start": m.start() + 1, "end": m.start() + len(text), "sequence": text})
    return hits


def search_pattern(raw_sequence: str, pattern: str) -> dict[str, Any]:
    """Search one sequence for a literal fragment or a PROSITE-style pattern."""
    sequence = normalize_sequence(raw_sequence)
    pattern_text = (pattern or "").strip().upper()
    if not pattern_text:
        raise InvalidSequence("Enter a pattern to search for.")
    is_prosite = any(t in pattern_text for t in ("-", "[", "{", "(", "<", ">"))
    if is_prosite:
        regex = compile_prosite(pattern_text)
        if regex is None:
            raise InvalidSequence(
                "That is not a valid PROSITE pattern. Example: N-{P}-[ST]-{P}. "
                "Use x for any residue, [ST] for a choice, {P} for 'not P', x(2,4) for repeats.")
        matches = _find_hits(sequence, regex)
        kind = "prosite"
    else:
        if not re.fullmatch(r"[A-Z]+", pattern_text):
            raise InvalidSequence("A literal search may only contain letters.")
        matches = []
        pos = sequence.find(pattern_text)
        while pos != -1:
            matches.append({"start": pos + 1, "end": pos + len(pattern_text),
                            "sequence": pattern_text})
            pos = sequence.find(pattern_text, pos + 1)
        kind = "literal"
    for i, m in enumerate(matches, 1):
        m["match"] = i
    return {"status": "success", "pattern": pattern_text, "kind": kind,
            "sequence_length": len(sequence), "count": len(matches), "matches": matches}


# --------------------------------------------------------------------
# PROSITE database + scan
# --------------------------------------------------------------------

class PrositeDatabase:
    """Parsed PROSITE pattern database with pre-compiled regular expressions."""

    def __init__(self, source: str = "") -> None:
        self.source = source
        self.entries: list[dict[str, Any]] = []
        self.total_records = 0
        self.pattern_records = 0
        self.unsupported = 0
        self.release = ""

    @classmethod
    def from_text(cls, text: str, source: str = "") -> "PrositeDatabase":
        db = cls(source)
        record: list[str] = []
        for line in text.splitlines():
            if line.startswith("//"):
                db._add_record(record)
                record = []
            else:
                record.append(line)
        db._add_record(record)
        m = re.search(r"Release\s+(\d{4}_\d+|[0-9.]+)", text[:4000])
        db.release = m.group(1) if m else ""
        return db

    def _add_record(self, lines: list[str]) -> None:
        if not lines:
            return
        entry: dict[str, Any] = {"id": None, "accession": None, "type": None,
                                 "description": "", "pattern": "", "doc": None,
                                 "frequent": False}
        desc: list[str] = []
        pat: list[str] = []
        for line in lines:
            tag = line[:2]
            body = line[5:].strip()
            if tag == "ID":
                fields = body.split(";")
                entry["id"] = fields[0].strip()
                if len(fields) > 1:
                    entry["type"] = fields[1].strip().rstrip(".")
            elif tag == "AC":
                entry["accession"] = body.split(";")[0].strip()
            elif tag == "DE":
                desc.append(body)
            elif tag == "PA":
                pat.append(body)
            elif tag == "DO":
                entry["doc"] = body.rstrip(";").strip()
            elif tag == "CC" and "/SKIP-FLAG=TRUE" in body.upper():
                entry["frequent"] = True
        if not entry["id"]:
            return
        self.total_records += 1
        entry["description"] = " ".join(desc).strip().rstrip(".")
        entry["pattern"] = "".join(pat).strip()
        if not entry["pattern"] or not entry["accession"]:
            return          # profile (MATRIX) records have no pattern; not scanned here
        self.pattern_records += 1
        regex = compile_prosite(entry["pattern"])
        if regex is None:
            self.unsupported += 1
            return
        entry["regex"] = regex
        self.entries.append(entry)

    @property
    def patterns_available(self) -> int:
        return len(self.entries)

    def info(self) -> dict[str, Any]:
        return {"source": self.source, "release": self.release,
                "records": self.total_records, "pattern_records": self.pattern_records,
                "patterns_usable": len(self.entries), "patterns_unsupported": self.unsupported}

    def scan(self, raw_sequence: str, exclude_frequent: bool = False) -> dict[str, Any]:
        sequence = normalize_sequence(raw_sequence)
        tested = 0
        found = []
        total_hits = 0
        for entry in self.entries:
            if exclude_frequent and entry["frequent"]:
                continue
            tested += 1
            hits = _find_hits(sequence, entry["regex"])
            if hits:
                total_hits += len(hits)
                found.append({
                    "id": entry["id"], "accession": entry["accession"],
                    "description": entry["description"], "pattern": entry["pattern"],
                    "frequent_pattern": entry["frequent"],
                    "url": f"https://prosite.expasy.org/{entry['accession']}",
                    "hits": hits})
        found.sort(key=lambda m: (-len(m["hits"]), m["accession"]))
        return {
            "status": "completed",
            "database": self.info(),
            "sequence_length": len(sequence),
            "exclude_frequent": exclude_frequent,
            "patterns_tested": tested,
            "patterns_matched": len(found),
            "total_hits": total_hits,
            "matches": found,
            "message": (
                f"Scan completed: {tested:,} PROSITE patterns tested, "
                f"{len(found):,} matched ({total_hits:,} occurrence"
                f"{'' if total_hits == 1 else 's'})."),
        }


def motif_status(state: str, message: str, **extra: Any) -> dict[str, Any]:
    """Explicit non-success scan states (unavailable / loading / failed)."""
    return {"status": state, "message": message, "patterns_tested": 0,
            "patterns_matched": 0, "total_hits": 0, "matches": [], **extra}


# --------------------------------------------------------------------
# comparison
# --------------------------------------------------------------------

def _kmers(sequence: str, k: int) -> set[str]:
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1)}


def compare_sequences(first_raw: str, second_raw: str, k: int = pa.KMER_SIZE) -> dict[str, Any]:
    a = normalize_sequence(first_raw)
    b = normalize_sequence(second_raw)
    try:
        k = int(k)
    except (TypeError, ValueError):
        raise InvalidSequence("k-mer size must be a whole number.")
    if k < 1:
        raise InvalidSequence("k-mer size must be at least 1.")
    if k > KMER_LIMIT:
        raise InvalidSequence(
            f"k-mer size is limited to {KMER_LIMIT} residues (you asked for {k}). "
            "Longer words are almost never shared between different proteins.")
    shortest = min(len(a), len(b))
    if k > shortest:
        raise InvalidSequence(
            f"k = {k} is longer than the shorter sequence ({shortest} aa), so it has no "
            f"{k}-mers to compare. Choose k of {shortest} or less.")

    set_a, set_b = _kmers(a, k), _kmers(b, k)
    shared = set_a & set_b
    union = set_a | set_b
    jaccard = len(shared) / len(union) * 100.0 if union else 0.0

    # how many identical windows would two random proteins with these compositions share?
    pa_comp = Counter(a)
    pb_comp = Counter(b)
    match_prob = sum((pa_comp[c] / len(a)) * (pb_comp[c] / len(b)) for c in pa_comp)
    windows_a, windows_b = len(a) - k + 1, len(b) - k + 1
    expected_pairs = windows_a * windows_b * (match_prob ** k)

    profile = []
    for kk in range(1, min(shortest, 10) + 1):
        sa, sb = _kmers(a, kk), _kmers(b, kk)
        un = sa | sb
        profile.append({"k": kk, "shared": len(sa & sb),
                        "similarity_percent": len(sa & sb) / len(un) * 100.0 if un else 0.0})

    if shared:
        note = (f"{len(shared):,} distinct {k}-mer{'s' if len(shared) != 1 else ''} occur in both "
                f"sequences ({jaccard:.2f}% of all distinct {k}-mers).")
    else:
        note = (f"No {k}-mer is shared between the two sequences. This is a genuine result, "
                f"not a failed calculation: unrelated proteins rarely share identical stretches "
                f"of {k} residues (about {expected_pairs:.2f} identical windows would be expected "
                f"by chance for these compositions). Try a smaller k (2 or 3) to look for weaker similarity.")

    return clean_json({
        "status": "success", "k": k, "kmer_limit": KMER_LIMIT,
        "lengths": [len(a), len(b)],
        "windows": [windows_a, windows_b],
        "distinct_kmers": [len(set_a), len(set_b)],
        "shared_kmers": len(shared), "union_kmers": len(union),
        "word_similarity_percent": jaccard,
        "expected_identical_windows_by_chance": expected_pairs,
        "shared_examples": sorted(shared)[:20],
        "composition_distance": pa.composition_distance(a, b),
        "profile": profile,
        "interpretation": note,
        "method": "Jaccard index of the distinct k-mer sets (shared / union). "
                  "A fast alignment-free screen, not a substitute for BLAST.",
    })


def dotplot_points(first_raw: str, second_raw: str, k: int = 3, limit: int = 6000) -> dict[str, Any]:
    """Shared k-mer positions for a dot plot (x = second sequence, y = first).

    If the plot would be too dense the word length is increased (never above
    the k-mer limit) until it is drawable; the word length actually used is returned.
    """
    a = normalize_sequence(first_raw)
    b = normalize_sequence(second_raw)
    k_max = min(KMER_LIMIT, len(a), len(b))
    if k < 1 or k > k_max:
        raise InvalidSequence(f"Dot-plot word length must be between 1 and {k_max}.")
    used = k
    while True:
        index: dict[str, list[int]] = {}
        for j in range(len(b) - used + 1):
            index.setdefault(b[j:j + used], []).append(j + 1)
        occurrences: dict[str, list[int]] = {}
        for i in range(len(a) - used + 1):
            word = a[i:i + used]
            if word in index:
                occurrences.setdefault(word, []).append(i + 1)
        pairs = sum(len(pos) * len(index[word]) for word, pos in occurrences.items())
        if pairs <= limit or used >= k_max:
            break
        used += 1
    points: list[dict[str, int]] = []
    for word, positions in occurrences.items():
        for i in positions:
            for j in index[word]:
                if len(points) < limit:
                    points.append({"x": j, "y": i})
    points.sort(key=lambda p: (p["y"], p["x"]))
    return {"k": used, "requested_k": k, "points": points, "pairs_total": pairs,
            "truncated": pairs > limit, "length_first": len(a), "length_second": len(b)}


def _pick(section: dict[str, Any], *path: str) -> Any:
    node: Any = section
    for key in path:
        if not isinstance(node, dict) or key not in node:
            return None
        node = node[key]
    return node


COMPARISON_PROPERTIES = (
    ("length", "Length", "aa", ("sequence", "length")),
    ("molecular_weight", "Molecular weight", "Da", ("physicochemical", "molecular_weight", "value_da")),
    ("isoelectric_point", "Isoelectric point", "pH", ("physicochemical", "isoelectric_point", "value")),
    ("net_charge_ph7", "Net charge at pH 7", "e", ("physicochemical", "net_charge_ph7")),
    ("gravy", "GRAVY", "", ("physicochemical", "gravy")),
    ("aromaticity", "Aromaticity", "%", ("physicochemical", "aromaticity_percent")),
    ("aliphatic_index", "Aliphatic index", "", ("physicochemical", "aliphatic_index")),
)


def compare_proteins(first_raw: str, second_raw: str, k: int = pa.KMER_SIZE) -> dict[str, Any]:
    """Everything the comparison screen shows, from the same engine as every other screen."""
    a = normalize_sequence(first_raw)
    b = normalize_sequence(second_raw)
    similarity = compare_sequences(a, b, k)       # validates k with a precise message
    analysis_a, analysis_b = analyze_sequence(a), analyze_sequence(b)
    properties = []
    for key, label, unit, path in COMPARISON_PROPERTIES:
        properties.append({"key": key, "label": label, "unit": unit,
                           "first": _pick(analysis_a, *path), "second": _pick(analysis_b, *path)})
    rows_a = {r["code"]: r for r in analysis_a["composition"]["rows"]}
    rows_b = {r["code"]: r for r in analysis_b["composition"]["rows"]}
    composition = [{"code": c, "name": pa.AMINO_ACID_NAMES[c],
                    "first_percent": rows_a[c]["percent"], "second_percent": rows_b[c]["percent"],
                    "first_count": rows_a[c]["count"], "second_count": rows_b[c]["count"]}
                   for c in pa.AMINO_ACIDS]
    return clean_json({
        "status": "success", "properties": properties, "composition": composition,
        "similarity": similarity, "dotplot": dotplot_points(a, b, min(k, 5) if k > 5 else k),
    })


# --------------------------------------------------------------------
# dataset-wide statistics (streamed, never held as sequences)
# --------------------------------------------------------------------

_AA_ORDER = pa.AMINO_ACIDS
_MASS_VEC = np.array([pa.AMINO_ACID_MASSES[c] - pa.WATER_MASS for c in _AA_ORDER])
_HYDRO_VEC = np.array([pa.HYDROPATHY_SCALE[c] for c in _AA_ORDER])
_IDX = {c: i for i, c in enumerate(_AA_ORDER)}


def iter_fasta(lines: Iterable[str]) -> Iterator[tuple[str, str]]:
    header, chunks = None, []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:], []
        elif header is not None:
            chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


class DatasetAggregator:
    """Accumulates dataset-wide statistics from a stream of sequences."""

    def __init__(self) -> None:
        self.n = 0
        self._lengths: list[int] = []
        self._counts: list[list[int]] = []
        self.nonstandard_totals: Counter = Counter()
        self.proteins_with_nonstandard = 0
        self.proteins_invalid = 0

    def add(self, sequence: str) -> None:
        sequence = sequence.strip().upper()
        if not sequence:
            self.proteins_invalid += 1
            return
        row = [sequence.count(c) for c in _AA_ORDER]
        known = sum(row)
        ns = 0
        for c in pa.NON_STANDARD_CODES:
            k = sequence.count(c)
            if k:
                self.nonstandard_totals[c] += k
                ns += k
        if ns:
            self.proteins_with_nonstandard += 1
        if known + ns != len(sequence):
            self.proteins_invalid += 1
        self.n += 1
        self._lengths.append(len(sequence))
        self._counts.append(row)

    @staticmethod
    def _describe(values: np.ndarray) -> dict[str, Any]:
        if values.size == 0:
            return {"count": 0}
        return {
            "count": int(values.size), "min": float(values.min()), "max": float(values.max()),
            "mean": float(values.mean()), "median": float(np.median(values)),
            "std": float(values.std()),
            "p05": float(np.percentile(values, 5)), "p25": float(np.percentile(values, 25)),
            "p75": float(np.percentile(values, 75)), "p95": float(np.percentile(values, 95)),
        }

    @staticmethod
    def _histogram(values: np.ndarray, edges: list[float], overflow_label: str,
                   underflow_label: Optional[str] = None) -> dict[str, Any]:
        counts, _ = np.histogram(np.clip(values, edges[0], edges[-1] - 1e-9), bins=edges)
        labels = []
        for i in range(len(edges) - 1):
            lo, hi = edges[i], edges[i + 1]
            if i == len(edges) - 2:
                labels.append(overflow_label)
            elif i == 0 and underflow_label:
                labels.append(underflow_label)
            else:
                labels.append(f"{lo:g} to {hi:g}")
        return {"edges": edges, "labels": labels, "counts": counts.tolist(),
                "centres": [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]}

    def finish(self) -> dict[str, Any]:
        if self.n == 0:
            return {"status": "no_data", "message": "No sequences were processed.", "total_proteins": 0}
        lengths = np.array(self._lengths, dtype=np.float64)
        counts = np.array(self._counts, dtype=np.float64)
        std_total = counts.sum(axis=1)
        total_residues = float(lengths.sum())

        mass = counts @ _MASS_VEC + pa.WATER_MASS
        gravy = (counts @ _HYDRO_VEC) / np.maximum(lengths, 1)
        arom_idx = [_IDX[c] for c in pa.AROMATIC_RESIDUES]
        aromaticity = counts[:, arom_idx].sum(axis=1) / np.maximum(lengths, 1) * 100.0

        # isoelectric point, vectorised bisection using the same pKa set as the engine
        def col(code: str) -> np.ndarray:
            return counts[:, _IDX[code]]
        low = np.zeros(self.n)
        high = np.full(self.n, 14.0)
        for _ in range(50):
            mid = (low + high) / 2.0
            charge = 1.0 / (1.0 + 10 ** (mid - pa.PKA_N_TERMINUS))
            for code, pka in pa.PKA_POSITIVE.items():
                charge += col(code) / (1.0 + 10 ** (mid - pka))
            charge -= 1.0 / (1.0 + 10 ** (pa.PKA_C_TERMINUS - mid))
            for code, pka in pa.PKA_NEGATIVE.items():
                charge -= col(code) / (1.0 + 10 ** (pka - mid))
            positive = charge > 0
            low = np.where(positive, mid, low)
            high = np.where(positive, high, mid)
        pi = (low + high) / 2.0

        aa_totals = counts.sum(axis=0)
        aa_rows = []
        for code in _AA_ORDER:
            n = float(aa_totals[_IDX[code]])
            aa_rows.append({"code": code, "name": pa.AMINO_ACID_NAMES[code],
                            "three": THREE_LETTER[code], "family": pa.AMINO_ACID_CLASSES[code],
                            "count": int(n), "percent": n / total_residues * 100.0})
        family_rows = []
        for cls in pa.CLASS_ORDER:
            n = sum(r["count"] for r in aa_rows if r["family"] == cls)
            family_rows.append({"family": cls, "count": n, "percent": n / total_residues * 100.0})

        length_edges = [float(x) for x in list(range(0, 2001, 50)) + [1e9]]
        length_edges[-1] = 2050.0
        mass_edges = [float(x) for x in range(0, 200001, 10000)]
        mass_edges.append(210000.0)

        out = {
            "status": "success",
            "total_proteins": int(self.n),
            "total_residues": int(total_residues),
            "length": {**self._describe(lengths),
                       "histogram": self._histogram(lengths, length_edges, ">2000")},
            "molecular_weight": {**self._describe(mass),
                                 "unit": "Da",
                                 "note": "Standard residues only; non-standard codes are not counted.",
                                 "histogram": self._histogram(mass, mass_edges, ">200000")},
            "isoelectric_point": {**self._describe(pi),
                                  "histogram": self._histogram(pi, [float(x) for x in range(0, 15)], "13-14")},
            "gravy": {**self._describe(gravy),
                      "histogram": self._histogram(gravy, [round(-2.0 + 0.2 * i, 1) for i in range(0, 21)], ">1.8", "<-1.8")},
            "aromaticity": {**self._describe(aromaticity),
                            "histogram": self._histogram(aromaticity, [float(x) for x in range(0, 22, 2)], ">20")},
            "amino_acids": aa_rows,
            "families": family_rows,
            "validation": {
                "proteins_checked": int(self.n),
                "proteins_with_nonstandard_codes": int(self.proteins_with_nonstandard),
                "proteins_with_unexpected_characters": int(self.proteins_invalid),
                "nonstandard_residue_totals": {k: int(v) for k, v in self.nonstandard_totals.items()},
                "standard_residue_share_percent": float(std_total.sum() / total_residues * 100.0),
            },
        }
        return clean_json(out)
