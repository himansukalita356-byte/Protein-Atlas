"""Generate 3D geometries for the 20 canonical L-amino acids.

This is a build-time helper (needs numpy, scipy, networkx). Its output,
data/amino_acid_structures.json, ships with the application, so users never
need to run it.

Method
------
1. Parse an isomeric SMILES string (Kekule form) with a small built-in parser.
2. Add hydrogens from valence rules.
3. Embed in 3D by minimising a restraint force field (bond lengths, angles,
   ring/planar-group templates, torsion preferences, van der Waals repulsion
   and chirality volumes) starting from a graph layout; several seeds are
   tried and the best geometry with correct stereochemistry is kept.
4. The result is verified independently (bond lengths, angles, ring planarity,
   clashes, and the CORN rule for L-alpha-carbon chirality).

The geometry is a clean, chemically sensible model of the neutral free amino
acid -- not a quantum-chemical or crystallographic structure.
"""
from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.optimize import least_squares

# L-amino acids, Kekule SMILES, stereo as in PubChem isomeric SMILES
SMILES = {
    "G": "NCC(=O)O",
    "A": "N[C@@H](C)C(=O)O",
    "V": "N[C@@H](C(C)C)C(=O)O",
    "L": "N[C@@H](CC(C)C)C(=O)O",
    "I": "CC[C@H](C)[C@@H](C(=O)O)N",
    "M": "N[C@@H](CCSC)C(=O)O",
    "P": "OC(=O)[C@@H]1CCCN1",
    "F": "N[C@@H](CC1=CC=CC=C1)C(=O)O",
    "W": "N[C@@H](CC1=CNC2=CC=CC=C12)C(=O)O",
    "Y": "N[C@@H](CC1=CC=C(O)C=C1)C(=O)O",
    "S": "N[C@@H](CO)C(=O)O",
    "T": "C[C@H]([C@@H](C(=O)O)N)O",
    "C": "N[C@@H](CS)C(=O)O",
    "N": "N[C@@H](CC(N)=O)C(=O)O",
    "Q": "N[C@@H](CCC(N)=O)C(=O)O",
    "D": "N[C@@H](CC(O)=O)C(=O)O",
    "E": "N[C@@H](CCC(O)=O)C(=O)O",
    "K": "N[C@@H](CCCCN)C(=O)O",
    "R": "N[C@@H](CCCN=C(N)N)C(=O)O",
    "H": "N[C@@H](CC1=CNC=N1)C(=O)O",
}

VALENCE = {"C": 4, "N": 3, "O": 2, "S": 2, "H": 1}
VDW_MIN = {"HH": 2.0, "HX": 2.4, "XX": 3.0}


# ----------------------------------------------------------------------------
# SMILES parsing (only what the 20 amino acids need)
# ----------------------------------------------------------------------------
def parse_smiles(smiles: str):
    atoms = []          # dict(el, chir, nh)
    bonds = {}          # (i,j) -> order
    nbr_order = {}      # atom index -> ordered neighbour list ("H" marks implicit H)
    stack, prev, pending_order = [], None, 1
    rings = {}          # digit -> (atom index, order)
    i = 0

    def add_neighbor(a, b):
        nbr_order.setdefault(a, []).append(b)

    while i < len(smiles):
        ch = smiles[i]
        if ch == "(":
            stack.append(prev)
            i += 1
        elif ch == ")":
            prev = stack.pop()
            i += 1
        elif ch == "=":
            pending_order = 2
            i += 1
        elif ch.isdigit():
            digit = int(ch)
            if digit in rings:
                other, order = rings.pop(digit)
                order = max(order, pending_order)
                pending_order = 1
                bonds[tuple(sorted((prev, other)))] = order
                # closing atom: partner appears at the position of the digit
                add_neighbor(prev, other)
                # opening atom: reserved slot filled now
                idx = nbr_order[other].index(("ring", digit))
                nbr_order[other][idx] = prev
            else:
                rings[digit] = (prev, pending_order)
                pending_order = 1
                add_neighbor(prev, ("ring", digit))
            i += 1
        else:
            if ch == "[":
                j = smiles.index("]", i)
                body = smiles[i + 1:j]
                el = body[0]
                chir = None
                rest = body[1:]
                if rest.startswith("@@"):
                    chir, rest = "@@", rest[2:]
                elif rest.startswith("@"):
                    chir, rest = "@", rest[1:]
                nh = 0
                if rest.startswith("H"):
                    nh = int(rest[1:]) if len(rest) > 1 else 1
                i = j + 1
                atoms.append({"el": el, "chir": chir, "nh": nh, "bracket": True})
            else:
                atoms.append({"el": ch, "chir": None, "nh": None, "bracket": False})
                i += 1
            idx = len(atoms) - 1
            if prev is not None:
                bonds[tuple(sorted((prev, idx)))] = pending_order
                add_neighbor(idx, prev)
                add_neighbor(prev, idx)
            pending_order = 1
            if atoms[idx]["bracket"] and atoms[idx]["nh"]:
                add_neighbor(idx, "H")
            prev = idx
    # bracket atom H must come right after the preceding atom: fix ordering
    for idx, a in enumerate(atoms):
        if a["chir"]:
            order = nbr_order[idx]
            if "H" in order:
                order.remove("H")
                pos = 1 if (idx > 0 and order and order[0] != "H" and _has_preceding(idx, order)) else 0
                order.insert(pos, "H")
    return atoms, bonds, nbr_order


def _has_preceding(idx, order):
    # the preceding atom (the one that introduced idx) is always first in the list
    # unless idx is the first atom of the string
    return idx != 0


def add_hydrogens(atoms, bonds, nbr_order):
    n_heavy = len(atoms)
    valence_used = [0] * n_heavy
    for (a, b), order in bonds.items():
        valence_used[a] += order
        valence_used[b] += order
    elements = [a["el"] for a in atoms]
    coords_index = list(range(n_heavy))
    h_of = {}
    for i, a in enumerate(atoms):
        nh = a["nh"] if a["nh"] is not None else VALENCE[a["el"]] - valence_used[i]
        h_of[i] = []
        for _ in range(nh):
            elements.append("H")
            h_idx = len(elements) - 1
            bonds[tuple(sorted((i, h_idx)))] = 1
            h_of[i].append(h_idx)
    # replace the "H" placeholder in neighbour order with the first hydrogen
    for i, order in nbr_order.items():
        if "H" in order:
            order[order.index("H")] = h_of[i][0]
    return elements, bonds, h_of


# ----------------------------------------------------------------------------
# Topology helpers
# ----------------------------------------------------------------------------
def analyse(elements, bonds):
    n = len(elements)
    g = nx.Graph()
    g.add_nodes_from(range(n))
    for (a, b), o in bonds.items():
        g.add_edge(a, b, order=o)
    dist = dict(nx.all_pairs_shortest_path_length(g))
    rings = [c for c in nx.cycle_basis(g)]
    ring_of = {}
    for r in rings:
        for a in r:
            ring_of.setdefault(a, []).append(tuple(r))
    double = {a for (a, b), o in bonds.items() if o == 2} | {b for (a, b), o in bonds.items() if o == 2}
    sp2 = set()
    for a in range(n):
        if elements[a] == "H":
            continue
        if a in double:
            sp2.add(a)
    # aromatic ring members are sp2 (pyrrole-type N-H included)
    for r in rings:
        if any(x in double for x in r) and len(r) in (5, 6):
            # ring counts as unsaturated only if it contains a double bond
            unsat = sum(1 for x in r if x in double)
            if unsat >= 2:
                sp2.update(r)
    # conjugated heteroatoms attached to an sp2 carbon (amide N, guanidine N, carboxyl OH)
    for a in range(n):
        if elements[a] in ("N", "O") and a not in sp2:
            for nb in g.neighbors(a):
                if elements[nb] == "C" and nb in sp2 and nb in double:
                    # carbon with an exocyclic/any double bond -> conjugation
                    sp2.add(a)
    return g, dist, rings, ring_of, sp2, double


def regular_polygon_template(ring, bond_len):
    """2D coordinates of a regular polygon for ring (ordered as a cycle)."""
    n = len(ring)
    radius = bond_len / (2 * math.sin(math.pi / n))
    return {a: np.array([radius * math.cos(2 * math.pi * k / n),
                         radius * math.sin(2 * math.pi * k / n)])
            for k, a in enumerate(ring)}


def ring_cycle_order(g, ring):
    sub = g.subgraph(ring)
    start = ring[0]
    order = [start]
    prev, cur = None, start
    while True:
        nxt = [x for x in sub.neighbors(cur) if x != prev and x not in order[1:]]
        if not nxt:
            break
        prev, cur = cur, nxt[0]
        if cur == start:
            break
        order.append(cur)
        if len(order) == len(ring):
            break
    return order


def planar_systems(g, rings, sp2, elements, bl):
    """Rigid planar templates: unsaturated rings (fused rings merged) plus exocyclic atoms."""
    aromatic_rings = [ring_cycle_order(g, r) for r in rings
                      if len(r) in (5, 6) and all(a in sp2 for a in r)]
    systems, used = [], set()
    for i, r in enumerate(aromatic_rings):
        if i in used:
            continue
        group, queue = {i}, [i]
        while queue:
            k = queue.pop()
            for j, r2 in enumerate(aromatic_rings):
                if j not in group and len(set(aromatic_rings[k]) & set(r2)) >= 2:
                    group.add(j)
                    queue.append(j)
        used |= group
        systems.append([aromatic_rings[k] for k in sorted(group)])
    templates = []
    for system in systems:
        coords = {}
        first = system[0]
        coords.update(regular_polygon_template(first, 1.40))
        placed = {0}
        progress = True
        while len(placed) < len(system) and progress:
            progress = False
            for k, ring in enumerate(system):
                if k in placed:
                    continue
                shared = [a for a in ring if a in coords]
                if len(shared) >= 2:
                    p, q = coords[shared[0]], coords[shared[1]]
                    existing_centre = np.mean([coords[a] for a in coords], axis=0)
                    n = len(ring)
                    side = np.linalg.norm(q - p)
                    radius = side / (2 * math.sin(math.pi / n))
                    mid = (p + q) / 2
                    perp = np.array([-(q - p)[1], (q - p)[0]])
                    perp /= np.linalg.norm(perp)
                    if np.dot(perp, mid - existing_centre) < 0:
                        perp = -perp
                    centre = mid + perp * math.sqrt(max(radius ** 2 - (side / 2) ** 2, 0.0))
                    for a in ring:
                        if a in coords:
                            continue
                        # place by angular position around the new centre
                        pass
                    ang_p = math.atan2(*(p - centre)[::-1])
                    ang_q = math.atan2(*(q - centre)[::-1])
                    order = ring_cycle_order(g, ring)
                    ip, iq = order.index(shared[0]), order.index(shared[1])
                    step = (ang_q - ang_p)
                    step = (step + math.pi) % (2 * math.pi) - math.pi
                    direction = 1 if step > 0 else -1
                    if abs(ip - iq) not in (1, n - 1):
                        direction = -direction
                    # walk around the ring starting from shared[0]
                    seq = order[ip:] + order[:ip]
                    if (order[(ip + 1) % n] == shared[1]):
                        walk = seq
                        walk_dir = direction
                    else:
                        walk = [seq[0]] + seq[:0:-1]
                        walk_dir = direction
                    for m, a in enumerate(walk):
                        ang = ang_p + walk_dir * 2 * math.pi * m / n
                        pos = centre + radius * np.array([math.cos(ang), math.sin(ang)])
                        if a not in coords:
                            coords[a] = pos
                    placed.add(k)
                    progress = True
        # exocyclic substituents (one ring atom each)
        ring_atoms = set(coords)
        centre = np.mean([coords[a] for a in ring_atoms], axis=0)
        for a in list(ring_atoms):
            ext = [x for x in g.neighbors(a) if x not in ring_atoms]
            if not ext:
                continue
            ring_nbrs = [x for x in g.neighbors(a) if x in ring_atoms]
            direction = coords[a] - np.mean([coords[x] for x in ring_nbrs], axis=0)
            direction /= np.linalg.norm(direction)
            for x in ext:
                length = bl(a, x)
                coords[x] = coords[a] + direction * length
        templates.append(coords)
    return templates


# ----------------------------------------------------------------------------
# Force field
# ----------------------------------------------------------------------------
def build_model(elements, bonds, atoms_smiles, nbr_order):
    n = len(elements)
    g, dist, rings, ring_of, sp2, double = analyse(elements, bonds)
    order_of = {tuple(sorted(k)): v for k, v in bonds.items()}

    def bl(a, b):
        ea, eb = sorted((elements[a], elements[b]))
        o = order_of[tuple(sorted((a, b)))]
        pair = ea + eb
        in_ring = any(a in r and b in r for r in rings)
        both_sp2 = a in sp2 and b in sp2
        table = {"CH": 1.09, "HN": 1.01, "HO": 0.96, "HS": 1.34}
        if pair in table:
            return table[pair]
        if pair == "CC":
            if o == 2:
                return 1.40 if in_ring else 1.34
            if both_sp2 and in_ring:
                return 1.40
            if a in sp2 or b in sp2:
                return 1.51
            return 1.53
        if pair == "CN":
            if o == 2:
                return 1.34 if in_ring else 1.29
            if in_ring and both_sp2:
                return 1.37
            if both_sp2:
                return 1.34
            return 1.47
        if pair == "CO":
            if o == 2:
                return 1.23
            return 1.34 if (a in sp2 and b in sp2) else 1.43
        if pair == "CS":
            return 1.81
        raise ValueError(pair)

    bond_list = sorted(order_of)
    b_idx = np.array(bond_list)
    b_len = np.array([bl(a, b) for a, b in bond_list])

    # angle restraints -> 1-3 distances
    angle_terms = []
    templates = planar_systems(g, rings, sp2, elements, bl)
    template_atoms = [set(t) for t in templates]
    for centre in range(n):
        nbrs = sorted(g.neighbors(centre))
        if len(nbrs) < 2:
            continue
        pairs = list(itertools.combinations(nbrs, 2))
        targets = {}
        if centre in sp2 and len(nbrs) in (2, 3):
            ring_pairs = {}
            for (a, b) in pairs:
                for r in ring_of.get(centre, []):
                    if a in r and b in r and len(r) in (5, 6) and all(x in sp2 for x in r):
                        ring_pairs[(a, b)] = (len(r) - 2) * 180.0 / len(r)
            free = [p for p in pairs if p not in ring_pairs]
            assigned = sum(ring_pairs.values())
            if len(nbrs) == 3:
                remaining = 360.0 - assigned
                for p in free:
                    targets[p] = remaining / max(len(free), 1) if ring_pairs else 120.0
            else:
                for p in free:
                    targets[p] = 120.0
            targets.update(ring_pairs)
        else:
            for p in pairs:
                in_r = [len(r) for r in ring_of.get(centre, []) if p[0] in r and p[1] in r]
                targets[p] = 104.0 if (in_r and in_r[0] == 5) else 109.5
        for (a, b), ang in targets.items():
            la, lb = bl(centre, a), bl(centre, b)
            d = math.sqrt(la * la + lb * lb - 2 * la * lb * math.cos(math.radians(ang)))
            angle_terms.append((a, b, d))
    a_idx = np.array([(a, b) for a, b, _ in angle_terms])
    a_len = np.array([d for _, _, d in angle_terms])

    # rigid planar templates: all pairwise distances inside each template
    tmpl_pairs, tmpl_len = [], []
    for coords in templates:
        keys = sorted(coords)
        for x, y in itertools.combinations(keys, 2):
            tmpl_pairs.append((x, y))
            tmpl_len.append(float(np.linalg.norm(coords[x] - coords[y])))
    # hydrogens on template atoms: keep them in plane as exocyclic substituents
    t_idx = np.array(tmpl_pairs) if tmpl_pairs else np.zeros((0, 2), int)
    t_len = np.array(tmpl_len)

    # torsion terms
    planar_t, stagger_t, perp_t = [], [], []
    for (b, c) in bond_list:
        if elements[b] == "H" or elements[c] == "H":
            continue
        in_ring_bc = any(b in r and c in r for r in rings)
        nb_b = [x for x in g.neighbors(b) if x != c]
        nb_c = [x for x in g.neighbors(c) if x != b]
        if not nb_b or not nb_c:
            continue
        both_sp2 = b in sp2 and c in sp2
        both_sp3 = b not in sp2 and c not in sp2
        for a in nb_b:
            for d in nb_c:
                if in_ring_bc and both_sp2:
                    continue  # covered by rigid template
                if both_sp2 and not in_ring_bc:
                    planar_t.append((a, b, c, d))
                elif both_sp3 and not in_ring_bc:
                    if elements[a] != "H" and elements[d] != "H":
                        stagger_t.append((a, b, c, d))
                elif (b in sp2) != (c in sp2) and not in_ring_bc:
                    if elements[a] != "H" and elements[d] != "H":
                        perp_t.append((a, b, c, d))
    planar_t = np.array(planar_t) if planar_t else np.zeros((0, 4), int)
    stagger_t = np.array(stagger_t) if stagger_t else np.zeros((0, 4), int)
    perp_t = np.array(perp_t) if perp_t else np.zeros((0, 4), int)

    # nonbonded pairs (graph distance >= 3; 1-4 pairs weighted less)
    nb_pairs, nb_min = [], []
    for i in range(n):
        for j in range(i + 1, n):
            d = dist[i][j]
            if d < 3:
                continue
            if any(i in ts and j in ts for ts in template_atoms):
                continue
            hi, hj = elements[i] == "H", elements[j] == "H"
            key = "HH" if hi and hj else ("HX" if hi or hj else "XX")
            scale = 0.82 if d == 3 else 1.0
            nb_pairs.append((i, j))
            nb_min.append(VDW_MIN[key] * scale)
    nb_idx = np.array(nb_pairs)
    nb_min = np.array(nb_min)

    # chirality: neighbour order from SMILES
    chir = []
    for idx, a in enumerate(atoms_smiles):
        if not a["chir"]:
            continue
        order = nbr_order[idx]
        n1, n2, n3, n4 = order
        sign = -1.0 if a["chir"] == "@" else 1.0
        chir.append((idx, n2, n3, n4, sign))
    return dict(n=n, g=g, sp2=sp2, rings=rings, b_idx=b_idx, b_len=b_len, a_idx=a_idx,
                a_len=a_len, t_idx=t_idx, t_len=t_len, planar_t=planar_t, stagger_t=stagger_t, perp_t=perp_t,
                nb_idx=nb_idx, nb_min=nb_min, chir=chir, bl=bl, templates=templates,
                elements=elements, order_of=order_of)


def dihedral(p, idx):
    if len(idx) == 0:
        return np.zeros(0)
    a, b, c, d = (p[idx[:, k]] for k in range(4))
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / np.linalg.norm(b1, axis=1)[:, None]
    v = b0 - (b0 * b1n).sum(1)[:, None] * b1n
    w = b2 - (b2 * b1n).sum(1)[:, None] * b1n
    x = (v * w).sum(1)
    y = (np.cross(b1n, v) * w).sum(1)
    return np.arctan2(y, x)


def chirality_volumes(p, chir):
    vols = []
    for c, n2, n3, n4, _ in chir:
        vols.append(np.linalg.det(np.array([p[n2] - p[c], p[n3] - p[c], p[n4] - p[c]])))
    return np.array(vols)


def residuals_factory(m, stage):
    b_idx, a_idx, t_idx = m["b_idx"], m["a_idx"], m["t_idx"]
    n = m["n"]

    def res(x):
        p = x.reshape(n, 3)
        out = []
        out.append((np.linalg.norm(p[b_idx[:, 0]] - p[b_idx[:, 1]], axis=1) - m["b_len"]) * 12.0)
        out.append((np.linalg.norm(p[a_idx[:, 0]] - p[a_idx[:, 1]], axis=1) - m["a_len"]) * 6.0)
        if len(t_idx):
            out.append((np.linalg.norm(p[t_idx[:, 0]] - p[t_idx[:, 1]], axis=1) - m["t_len"]) * 8.0)
        if len(m["planar_t"]):
            out.append(np.sin(dihedral(p, m["planar_t"])) * 2.0)
        if stage >= 2:
            if len(m["stagger_t"]):
                # heavy-atom chains prefer the extended (anti) arrangement
                out.append(0.30 * (1 + np.cos(dihedral(p, m["stagger_t"]))))
                out.append(0.35 * (1 + np.cos(3 * dihedral(p, m["stagger_t"]))))
            if len(m["perp_t"]):
                out.append(0.30 * np.cos(dihedral(p, m["perp_t"])))
            d = np.linalg.norm(p[m["nb_idx"][:, 0]] - p[m["nb_idx"][:, 1]], axis=1)
            out.append(np.maximum(0.0, m["nb_min"] - d) * 3.0)
        if m["chir"]:
            vols = chirality_volumes(p, m["chir"])
            targets = np.array([c[4] * 2.3 for c in m["chir"]])
            out.append((vols - targets) * 3.0)
        return np.concatenate(out)
    return res


def embed(m, seed):
    n = m["n"]
    rng = np.random.default_rng(seed)
    pos = nx.spring_layout(m["g"], dim=3, seed=seed, iterations=200)
    x0 = np.array([pos[i] for i in range(n)]) * 3.5 + rng.normal(scale=0.3, size=(n, 3))
    x = x0.ravel()
    for stage in (1, 2):
        sol = least_squares(residuals_factory(m, stage), x, method="trf", max_nfev=120,
                            xtol=1e-6, ftol=1e-8)
        x = sol.x
    p = x.reshape(n, 3)
    return p, float(np.sum(residuals_factory(m, 2)(x) ** 2))


# ----------------------------------------------------------------------------
# Independent verification
# ----------------------------------------------------------------------------
def verify(m, p, smiles_atoms, nbr_order, code):
    problems = []
    b_idx = m["b_idx"]
    lengths = np.linalg.norm(p[b_idx[:, 0]] - p[b_idx[:, 1]], axis=1)
    worst = np.max(np.abs(lengths - m["b_len"]))
    if worst > 0.06:
        problems.append(f"bond length deviation {worst:.3f}")
    # angle deviation
    a_idx = m["a_idx"]
    d13 = np.linalg.norm(p[a_idx[:, 0]] - p[a_idx[:, 1]], axis=1)
    if np.max(np.abs(d13 - m["a_len"])) > 0.12:
        problems.append(f"1-3 deviation {np.max(np.abs(d13 - m['a_len'])):.3f}")
    # planarity of rigid templates
    for coords in m["templates"]:
        pts = p[sorted(coords)]
        centred = pts - pts.mean(0)
        _, s, _ = np.linalg.svd(centred)
        if s[-1] / math.sqrt(len(pts)) > 0.05:
            problems.append(f"template not planar ({s[-1] / math.sqrt(len(pts)):.3f})")
    # clashes
    nb_idx, nb_min = m["nb_idx"], m["nb_min"]
    d = np.linalg.norm(p[nb_idx[:, 0]] - p[nb_idx[:, 1]], axis=1)
    if np.any(d < 0.72 * nb_min):
        problems.append(f"clash: min ratio {np.min(d / nb_min):.2f}")
    # chirality signs
    for (c, n2, n3, n4, sign) in m["chir"]:
        vol = np.linalg.det(np.array([p[n2] - p[c], p[n3] - p[c], p[n4] - p[c]]))
        if vol * sign <= 0:
            problems.append("chirality violated")
    return problems


def corn_is_L(p, elements, g, code):
    """Independent L-configuration check via the CORN rule at the alpha carbon.

    Alpha carbon = carbon bonded to N (amine) and to the carboxyl carbon.
    """
    if code == "G":
        return True
    for c in range(len(elements)):
        if elements[c] != "C":
            continue
        nb = list(g.neighbors(c))
        n_atoms = [x for x in nb if elements[x] == "N"]
        carboxyl = [x for x in nb if elements[x] == "C" and any(
            elements[y] == "O" and m_order(g, x, y) == 2 for y in g.neighbors(x))
            and sum(1 for y in g.neighbors(x) if elements[y] == "O") == 2]
        if len(n_atoms) == 1 and len(carboxyl) == 1:
            h = [x for x in nb if elements[x] == "H"][0]
            r = [x for x in nb if x not in (n_atoms[0], carboxyl[0], h)][0]
            v = [p[carboxyl[0]] - p[c], p[r] - p[c], p[n_atoms[0]] - p[c]]
            # viewed from H (H towards viewer): CO -> R -> N clockwise
            # clockwise seen from the H side <=> positive triple product here
            hv = p[h] - p[c]
            triple = np.linalg.det(np.array(v))
            # orientation of (CO, R, N) relative to the H direction
            return np.sign(triple) * np.sign(np.dot(np.cross(v[0], v[1]), hv)) != 0 and \
                (np.dot(np.cross(v[0], v[1]), v[2]) * 1.0 > 0)
    return None


def m_order(g, a, b):
    return g[a][b]["order"]


def sdf_block(code, elements, p, bonds_dict):
    lines = [f"{code}", "  atlas-embed", ""]
    order = sorted(bonds_dict)
    lines.append(f"{len(elements):3d}{len(order):3d}  0  0  1  0  0  0  0  0999 V2000")
    for el, (x, y, z) in zip(elements, p):
        lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {el:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
    for (a, b) in order:
        lines.append(f"{a + 1:3d}{b + 1:3d}{bonds_dict[(a, b)]:3d}  0")
    lines.append("M  END")
    return "\n".join(lines) + "\n$$$$\n"


def build(code, attempts=14):
    atoms, bonds, nbr_order = parse_smiles(SMILES[code])
    elements, bonds, h_of = add_hydrogens(atoms, dict(bonds), nbr_order)
    m = build_model(elements, bonds, atoms, nbr_order)
    best = None
    for seed in range(attempts):
        p, cost = embed(m, seed)
        problems = verify(m, p, atoms, nbr_order, code)
        score = (len(problems), cost)
        if best is None or score < best[0]:
            best = (score, p, problems)
        if not problems and cost < 5.0:
            break
    (score, p, problems) = best
    return m, elements, bonds, p, problems


def main(out_path: Path):
    result, report = {}, []
    for code in SMILES:
        m, elements, bonds, p, problems = build(code)
        p = p - p.mean(0)
        result[code] = {
            "smiles": SMILES[code],
            "atoms": [[el, round(float(x), 4), round(float(y), 4), round(float(z), 4)]
                      for el, (x, y, z) in zip(elements, p)],
            "bonds": [[a, b, int(o)] for (a, b), o in sorted(bonds.items())],
        }
        report.append((code, len(elements), problems))
        print(code, len(elements), "OK" if not problems else problems)
    out_path.write_text(json.dumps(result, separators=(",", ":")), encoding="utf-8")
    print("wrote", out_path)
    return report


if __name__ == "__main__":
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / "data" / "amino_acid_structures.json"
    main(target)
