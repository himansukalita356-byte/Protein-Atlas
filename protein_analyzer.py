"""
============================================================
            PROTEIN SEQUENCE ANALYZER
============================================================

A menu-driven command-line tool for exploring UniProt protein records
and PROSITE patterns fetched from public services.

Features
--------
- Basic composition, side-chain class, and validation analysis
- Physicochemical properties: molecular weight, isoelectric point,
  net charge vs pH, extinction coefficient, aromaticity
- PROSITE motif detection (single protein or whole-dataset scan),
  custom pattern search, and repeated-sequence (k-mer) analysis
- Dataset-wide amino acid composition and statistics
- Kyte-Doolittle hydropathy profiling and membrane-segment detection
- Pairwise sequence comparison, including a dot plot
- CSV / FASTA / text-report export

Requirements
------------
Python 3.10+, NumPy, pandas, Matplotlib, FastAPI, and Pydantic.
See requirements.txt.

Usage
-----
    python protein_analyzer.py
    python protein_analyzer.py --protein P04637
    python protein_analyzer.py --source local --fasta path/to/file.fasta --prosite path/to/prosite.dat

The default mode fetches one protein and PROSITE definitions from the
internet, so no dataset files are required. Use --source local only
when you intentionally want to use local files.

License: MIT (see LICENSE)
============================================================
"""

from __future__ import annotations

import argparse
import json
import os
import re
from typing import Any, Optional
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request as UrlRequest, urlopen

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.rule import Rule
from rich.table import Table
from rich.text import Text

__version__ = "2.0.0"

console = Console()

# Type aliases documenting the two core data shapes used throughout
# the program. Both are plain dicts at runtime (kept lightweight and
# dependency-free); these aliases exist purely so function
# signatures communicate intent.
#
# ProteinRecord keys: id, header, name, gene, organism, sequence,
#                      search_key
# MotifRecord keys:    id, accession, description, pattern, type,
#                      search_key
ProteinRecord = dict[str, Any]
MotifRecord = dict[str, Any]


# ============================================================
# FILE LOCATIONS
# ============================================================

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

# Use paths relative to the project directory first so the app works
# on any machine without a hard-coded user-specific path.
DEFAULT_FASTA_CANDIDATES = [
    os.path.join(PROJECT_ROOT, "uniprotkb_HUMAN_2026_09_14.fasta"),
    os.path.join(PROJECT_ROOT, "data", "uniprotkb_HUMAN_2026_09_14.fasta"),
    os.path.join(PROJECT_ROOT, "datasets", "uniprotkb_HUMAN_2026_09_14.fasta"),
]

DEFAULT_PROSITE_CANDIDATES = [
    os.path.join(PROJECT_ROOT, "prosite.DAT"),
    os.path.join(PROJECT_ROOT, "data", "prosite.DAT"),
    os.path.join(PROJECT_ROOT, "datasets", "prosite.DAT"),
]

FASTA_FILE = DEFAULT_FASTA_CANDIDATES[0]
PROSITE_FILE = DEFAULT_PROSITE_CANDIDATES[0]


# ============================================================
# AMINO ACID CONSTANTS
# ============================================================

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Non-standard codes that really do appear in UniProt files.
NON_STANDARD_CODES = "UOXBZJ"

KNOWN_CODES = AMINO_ACIDS + NON_STANDARD_CODES

# Codes whose mass is an approximation rather than an exact value.
APPROXIMATE_CODES = "XBZJ"

AMINO_ACID_NAMES = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine",
    "D": "Aspartic acid", "C": "Cysteine", "E": "Glutamic acid",
    "Q": "Glutamine", "G": "Glycine", "H": "Histidine",
    "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
    "M": "Methionine", "F": "Phenylalanine", "P": "Proline",
    "S": "Serine", "T": "Threonine", "W": "Tryptophan",
    "Y": "Tyrosine", "V": "Valine",
    "U": "Selenocysteine", "O": "Pyrrolysine",
    "X": "Unknown residue", "B": "Asn or Asp",
    "Z": "Gln or Glu", "J": "Leu or Ile"}

# Free amino acid masses in Daltons (average, not monoisotopic).
AMINO_ACID_MASSES = {
    "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10,
    "C": 121.16, "E": 147.13, "Q": 146.15, "G": 75.07,
    "H": 155.16, "I": 131.17, "L": 131.17, "K": 146.19,
    "M": 149.21, "F": 165.19, "P": 115.13, "S": 105.09,
    "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
    "U": 168.05, "O": 255.31,
    "B": 132.61, "Z": 146.64, "J": 131.17, "X": 128.16}

WATER_MASS = 18.015

# Kyte & Doolittle hydropathy scale.
HYDROPATHY_SCALE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5,
    "C": 2.5, "E": -3.5, "Q": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
    "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
    "U": 2.5, "O": -3.9}

# Side-chain classes used for grouped composition tables.
AMINO_ACID_CLASSES = {
    "A": "Nonpolar", "V": "Nonpolar", "L": "Nonpolar",
    "I": "Nonpolar", "M": "Nonpolar", "G": "Nonpolar",
    "P": "Nonpolar",
    "F": "Aromatic", "W": "Aromatic", "Y": "Aromatic",
    "S": "Polar", "T": "Polar", "C": "Polar",
    "N": "Polar", "Q": "Polar",
    "D": "Acidic", "E": "Acidic",
    "K": "Basic", "R": "Basic", "H": "Basic"}

CLASS_ORDER = ["Nonpolar", "Aromatic", "Polar", "Acidic", "Basic"]

CLASS_COLORS = {
    "Nonpolar": "#D97706",
    "Aromatic": "#7C3AED",
    "Polar": "#0F766E",
    "Acidic": "#C2410C",
    "Basic": "#1D4ED8"}

PLOT_COLORS = {
    "ink": "#263238",
    "muted": "#607D8B",
    "accent": "#D97706",
    "blue": "#1D4ED8",
    "teal": "#0F766E",
    "red": "#C2410C",
    "green": "#15803D",
    "surface": "#F7F9FA",
    "track": "#CFD8DC"}

ACIDIC_RESIDUES = ["D", "E"]
BASIC_RESIDUES = ["K", "R", "H"]
CHARGED_RESIDUES = ["D", "E", "K", "R", "H"]

UNCHARGED_RESIDUES = [
    "A", "C", "F", "G", "I", "L", "M", "N",
    "P", "Q", "S", "T", "V", "W", "Y"]

AROMATIC_RESIDUES = ["F", "W", "Y"]

# pKa values used for net charge and isoelectric point (EMBOSS set).
PKA_N_TERMINUS = 8.6
PKA_C_TERMINUS = 3.6

PKA_POSITIVE = {"K": 10.8, "R": 12.5, "H": 6.5}
PKA_NEGATIVE = {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}

# Molar extinction coefficients at 280 nm (M-1 cm-1).
EXTINCTION_TRYPTOPHAN = 5500
EXTINCTION_TYROSINE = 1490
EXTINCTION_CYSTINE = 125


# ============================================================
# PROGRAM SETTINGS
# ============================================================

MAX_DISPLAY_RESULTS = 20

# Page size specifically for browsing a whole-dataset PROSITE motif
# scan, which can return results for a very large share of the
# proteome. Kept separate from MAX_DISPLAY_RESULTS so that search
# screens elsewhere in the program are unaffected.
SCAN_RESULTS_PAGE_SIZE = 100

DEFAULT_HYDROPATHY_WINDOW = 9

TRANSMEMBRANE_WINDOW = 19
TRANSMEMBRANE_THRESHOLD = 1.6

DOT_PLOT_MAX_LENGTH = 1200

KMER_SIZE = 5

EXPORT_FOLDER = "protein_analyzer_exports"

# Dataset-wide results are slow to compute, so they are stored here
# after the first time and reused afterwards.
DATASET_CACHE = {}

# Set once in main(), from whichever FASTA path was actually loaded
# (the configured one, or a substitute the user typed in if that
# path was missing). Output titles use this instead of repeating
# the word "Dataset" with no indication of which dataset it is.
DATASET_SOURCE_PATH = None
REMOTE_QUERY = "*"
REMOTE_PAGE_SIZE = 100
REMOTE_OFFSET = 0
REMOTE_TOTAL = 0


def dataset_label() -> str:
    """A short, human-readable name for the loaded FASTA dataset."""

    if not DATASET_SOURCE_PATH:
        return f"UniProt query: {REMOTE_QUERY}"

    name = os.path.splitext(os.path.basename(DATASET_SOURCE_PATH))[0]

    return f"UniProt Human ({name})"


# ============================================================
# DISPLAY HELPERS
# ============================================================

def print_header(title: str) -> None:
    console.print()
    console.print(
        Panel(
            Text(title.upper(), justify="center", style="bold white"),
            border_style="bright_cyan",
            box=box.ROUNDED,
            padding=(0, 2),
        )
    )


def print_section(title: str) -> None:
    console.print(Rule(title, style="bright_cyan", characters="─"))


def format_residue(code: str) -> str:
    name = AMINO_ACID_NAMES.get(code)

    if name is None:
        return code

    return f"{code} ({name})"


def format_residue_list(codes: list[str]) -> str:
    if not codes:
        return "None"

    return ", ".join(format_residue(code) for code in codes)


def shorten(text: str, width: int) -> str:
    if len(text) <= width:
        return text

    return text[:width - 3] + "..."


def ask_yes_no(question: str) -> bool:
    answer = input(f"{question} (y/n): ").strip().lower()

    return answer in ("y", "yes")


def ask_number(question: str, minimum: int, maximum: int, default: int) -> int:
    answer = input(
        f"{question} [{minimum}-{maximum}, "
        f"Enter for {default}]: ").strip()

    if not answer:
        return default

    if not answer.isdigit():
        print("Not a number. Using the default value.")
        return default

    number = int(answer)

    if number < minimum or number > maximum:
        print("Out of range. Using the default value.")
        return default

    return number


def show_menu(title: str, options: list[str]) -> int:
    """Print a numbered menu and return the chosen number."""

    while True:
        print_header(title)

        menu = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
        menu.add_column("", style="bold bright_cyan", width=4, justify="right")
        menu.add_column("Option", style="white")
        for number, label in enumerate(options, start=1):
            menu.add_row(str(number), label)
        console.print(menu)
        console.print()

        choice = input("Enter your choice: ").strip()

        if choice.isdigit():
            number = int(choice)

            if 1 <= number <= len(options):
                return number

        print(
            f"\nInvalid choice. "
            f"Please enter a number from 1 to {len(options)}.")


def exit_program() -> None:
    print("\nExiting Protein Sequence Analyzer.")
    raise SystemExit


def prepare_plot() -> tuple[plt.Figure, plt.Axes]:
    """Create a figure with one set of axes and a white background."""

    figure, axes = plt.subplots(figsize=(11, 6))
    figure.patch.set_facecolor("white")
    axes.set_facecolor(PLOT_COLORS["surface"])

    return figure, axes


def finish_plot(axes: plt.Axes, title: str, x_label: str, y_label: str) -> None:
    figure = axes.figure

    axes.set_title(title, fontsize=14, fontweight="bold", pad=14)
    axes.set_xlabel(x_label, fontsize=11)
    axes.set_ylabel(y_label, fontsize=11)
    axes.grid(True, color=PLOT_COLORS["track"], linestyle=":", alpha=0.8)
    axes.set_axisbelow(True)

    plt.tight_layout()
    figure.text(
        0.99, 0.012,
        f"Protein Sequence Analyzer v{__version__}",
        ha="right", va="bottom", fontsize=8,
        color=PLOT_COLORS["muted"])
    show_and_close(figure)


def show_and_close(figure: plt.Figure) -> None:
    """Display a figure and guarantee the program keeps running.

    On some systems, closing the chart window raises an exception
    from inside Matplotlib's event loop. Left uncaught, that
    exception was propagating all the way up and terminating the
    whole program instead of just closing the chart. Catching it
    here, and always closing the figure afterwards, keeps chart
    display self-contained so the calling menu resumes normally.
    """

    try:
        plt.show()

    except Exception as error:
        print(f"\n(Chart window closed: {error})")

    finally:
        plt.close(figure)


# ============================================================
# FASTA LOADING
# ============================================================

def extract_protein_id(header: str) -> str:
    parts = header.split("|")

    if len(parts) >= 2 and parts[1].strip():
        return parts[1].strip()

    header_parts = header.split()

    if header_parts:
        return header_parts[0]

    return "Unknown_ID"


def extract_header_field(header: str, key: str) -> str:
    """Return the value of a UniProt header field such as OS= or GN=."""

    match = re.search(
        rf"{key}=(.*?)(?=\s+[A-Z]{{2}}=|$)",
        header)

    if match:
        return match.group(1).strip()

    return ""


def extract_protein_name(header: str) -> str:
    text = header

    for key in ("OS=", "OX=", "GN=", "PE=", "SV="):
        if " " + key in text:
            text = text.split(" " + key)[0]

    parts = text.split("|")

    if len(parts) >= 3 and parts[2].strip():
        text = parts[2].strip()

    # "P53_HUMAN Cellular tumour antigen p53" -> drop the entry name.
    words = text.split(" ", 1)

    if len(words) == 2 and "_" in words[0] and words[0].isupper():
        text = words[1]

    return text.strip()


def build_protein_record(header: str, sequence: str) -> ProteinRecord:
    return {
        "id": extract_protein_id(header),
        "header": header,
        "name": extract_protein_name(header),
        "gene": extract_header_field(header, "GN"),
        "organism": extract_header_field(header, "OS"),
        "sequence": sequence,
        # Pre-computed once so searching stays fast.
        "search_key": (
            extract_protein_id(header) + " " + header).lower()}


def load_fasta(filename: str) -> list[ProteinRecord]:
    proteins = []
    current_header = None
    current_sequence = []

    try:
        with open(filename, "r", encoding="utf-8") as file:
            for line in file:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    if current_header is not None:
                        sequence = "".join(current_sequence).upper()

                        if sequence:
                            proteins.append(
                                build_protein_record(
                                    current_header, sequence))

                    current_header = line[1:].strip()
                    current_sequence = []

                else:
                    current_sequence.append(line)

        if current_header is not None:
            sequence = "".join(current_sequence).upper()

            if sequence:
                proteins.append(
                    build_protein_record(current_header, sequence))

        print(
            f"\nDataset loaded successfully: "
            f"{len(proteins):,} protein records")

        return proteins

    except FileNotFoundError:
        print("\nERROR: FASTA file not found:")
        print(filename)
        return []

    except PermissionError:
        print("\nERROR: Permission denied while reading the FASTA file.")
        return []

    except UnicodeDecodeError:
        print("\nERROR: Unable to decode the FASTA file as UTF-8.")
        return []

    except OSError as error:
        print(f"\nERROR while loading FASTA file: {error}")
        return []


# ============================================================
# PROSITE LOADING
# ============================================================

def parse_prosite_record(record: str) -> MotifRecord:
    motif = {
        "id": None,
        "accession": None,
        "description": None,
        "pattern": None,
        "type": None}

    description_parts = []
    pattern_parts = []

    for line in record.splitlines():
        if line.startswith("ID"):
            body = line[5:].strip()
            fields = body.split(";")

            motif["id"] = fields[0].strip()

            if len(fields) >= 2:
                motif["type"] = fields[1].strip().rstrip(".")

        elif line.startswith("AC"):
            motif["accession"] = line[5:].split(";")[0].strip()

        elif line.startswith("DE"):
            # Descriptions can wrap across several lines.
            description_parts.append(line[5:].strip())

        elif line.startswith("PA"):
            pattern_parts.append(line[5:].strip())

    motif["description"] = " ".join(description_parts).strip()
    motif["pattern"] = "".join(pattern_parts).strip()

    return motif


def load_prosite(filename: str) -> list[MotifRecord]:
    motifs = []
    current_record = []

    try:
        with open(filename, "r", encoding="utf-8",
                  errors="replace") as file:
            for line in file:
                line = line.rstrip("\n")

                if line.strip() == "//":
                    if current_record:
                        motif = parse_prosite_record(
                            "\n".join(current_record))

                        if (
                            motif["id"]
                            and motif["accession"]
                            and motif["pattern"]
                        ):
                            motif["search_key"] = (
                                f"{motif['id']} "
                                f"{motif['accession']} "
                                f"{motif['description']}").lower()

                            motifs.append(motif)

                        current_record = []

                else:
                    current_record.append(line)

        print(
            f"\nPROSITE database loaded successfully: "
            f"{len(motifs):,} pattern records")

        return motifs

    except FileNotFoundError:
        print("\nERROR: PROSITE database file not found:")
        print(filename)
        return []

    except PermissionError:
        print("\nERROR: Permission denied while reading the PROSITE file.")
        return []

    except OSError as error:
        print(f"\nERROR while loading PROSITE database: {error}")
        return []


def load_prosite_text(content: str) -> list[MotifRecord]:
    """Parse PROSITE records supplied from memory instead of a local file."""

    motifs = []

    for record in content.split("//"):
        if not record.strip():
            continue

        motif = parse_prosite_record(record)

        if motif["id"] and motif["accession"] and motif["pattern"]:
            motif["search_key"] = (
                f"{motif['id']} {motif['accession']} "
                f"{motif['description']}").lower()
            motifs.append(motif)

    return motifs


def fetch_remote_protein(protein_id: str) -> ProteinRecord:
    """Fetch one UniProt record without requiring a local FASTA file."""

    base_url = os.getenv(
        "UNIPROT_API_BASE_URL", "https://rest.uniprot.org").rstrip("/")
    url = f"{base_url}/uniprotkb/{quote(protein_id.strip())}.json"
    request = UrlRequest(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "protein-sequence-analyzer/2.0",
        })

    try:
        with urlopen(request, timeout=20) as response:
            entry = json.loads(response.read().decode("utf-8"))
    except HTTPError as error:
        raise RuntimeError(
            f"UniProt could not find '{protein_id}' (HTTP {error.code}).") from error
    except (URLError, TimeoutError, json.JSONDecodeError) as error:
        raise RuntimeError(f"UniProt is unavailable: {error}") from error

    accession = entry.get("primaryAccession", protein_id.strip())
    description = entry.get("proteinDescription", {})
    recommended = description.get("recommendedName", {})
    full_name = recommended.get("fullName", {})
    genes = entry.get("genes", [])
    gene = genes[0].get("geneName", {}).get("value", "") if genes else ""
    organism = entry.get("organism", {}).get("scientificName", "")
    sequence = entry.get("sequence", {}).get("value", "")
    name = full_name.get("value", "")
    header = (
        f"sp|{accession}|{entry.get('uniProtkbId', '')} {name} "
        f"OS={organism} GN={gene}")

    return build_protein_record(header, sequence)


def fetch_remote_proteins(
        query: str = "*", offset: int = 0, size: int = 100
        ) -> tuple[list[ProteinRecord], int]:
    """Fetch one page from the complete UniProt protein universe.

    UniProt search results are paginated so the analyzer can work with
    unrestricted queries without downloading millions of records at startup.
    The returned total is the current number of records matching the query.
    """

    if size < 1 or size > 500:
        raise ValueError("Remote page size must be between 1 and 500.")

    base_url = os.getenv(
        "UNIPROT_API_BASE_URL", "https://rest.uniprot.org").rstrip("/")
    url = f"{base_url}/uniprotkb/search?{urlencode({
        'query': query or '*',
        'format': 'json',
        'size': size,
        'offset': offset})}"
    request = UrlRequest(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "protein-sequence-analyzer/2.0",
        })

    try:
        with urlopen(request, timeout=20) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except HTTPError as error:
        raise RuntimeError(
            f"UniProt search failed (HTTP {error.code}).") from error
    except (URLError, TimeoutError, json.JSONDecodeError) as error:
        raise RuntimeError(f"UniProt search is unavailable: {error}") from error

    records = []
    for entry in payload.get("results", []):
        accession = entry.get("primaryAccession", "")
        description = entry.get("proteinDescription", {})
        recommended = description.get("recommendedName", {})
        full_name = recommended.get("fullName", {})
        genes = entry.get("genes", [])
        gene = genes[0].get("geneName", {}).get("value", "") if genes else ""
        organism = entry.get("organism", {}).get("scientificName", "")
        sequence = entry.get("sequence", {}).get("value", "")
        name = full_name.get("value", "")
        header = (
            f"sp|{accession}|{entry.get('uniProtkbId', '')} {name} "
            f"OS={organism} GN={gene}")

        if accession and sequence:
            records.append(build_protein_record(header, sequence))

    return records, int(payload.get("totalResults", len(records)))


def load_more_remote_proteins(proteins: list[ProteinRecord]) -> None:
    """Append the next remote UniProt page to the current session."""

    global REMOTE_OFFSET, REMOTE_TOTAL

    if REMOTE_OFFSET >= REMOTE_TOTAL:
        print("\nAll matching remote UniProt records are already loaded.")
        return

    try:
        print(
            f"\nFetching remote records "
            f"{REMOTE_OFFSET + 1:,}-"
            f"{min(REMOTE_OFFSET + REMOTE_PAGE_SIZE, REMOTE_TOTAL):,}...")
        page, total = fetch_remote_proteins(
            REMOTE_QUERY, REMOTE_OFFSET, REMOTE_PAGE_SIZE)
    except RuntimeError as error:
        print(f"\nUnable to load more remote proteins: {error}")
        return

    known_ids = {protein["id"] for protein in proteins}
    new_records = [
        protein for protein in page if protein["id"] not in known_ids]
    proteins.extend(new_records)
    REMOTE_OFFSET += len(page)
    REMOTE_TOTAL = total
    DATASET_CACHE.clear()

    print(
        f"Loaded {len(new_records):,} new protein(s). "
        f"Session contains {len(proteins):,} of {REMOTE_TOTAL:,} matches.")


def load_remote_prosite() -> list[MotifRecord]:
    """Download PROSITE definitions into memory without writing a file."""

    url = os.getenv(
        "PROSITE_URL", "https://ftp.expasy.org/databases/prosite/prosite.dat")
    request = UrlRequest(
        url,
        headers={
            "Accept": "text/plain",
            "User-Agent": "protein-sequence-analyzer/2.0",
        })

    try:
        with urlopen(request, timeout=20) as response:
            content = response.read().decode("utf-8", errors="replace")
    except (HTTPError, URLError, TimeoutError) as error:
        raise RuntimeError(f"Remote PROSITE source is unavailable: {error}") from error

    return load_prosite_text(content)


def resolve_file(path: str, description: str) -> Optional[str]:
    """Use the configured path, or find a nearby copy in the project.

    This is intentionally tolerant of machine-specific install locations.
    It checks the explicit path first, then common files next to the script,
    and only then prompts the user for a replacement.
    """

    candidate_paths = []

    if path:
        candidate_paths.append(path)
        candidate_paths.append(os.path.expanduser(path))

    if description.lower().startswith("fasta"):
        candidate_paths.extend(DEFAULT_FASTA_CANDIDATES)
    elif description.lower().startswith("prosite"):
        candidate_paths.extend(DEFAULT_PROSITE_CANDIDATES)

    seen = set()

    for candidate in candidate_paths:
        normalized = os.path.normpath(candidate)

        if normalized not in seen and os.path.exists(normalized):
            seen.add(normalized)
            return normalized

    print(f"\nThe configured {description} was not found.")
    for candidate in candidate_paths:
        if candidate:
            print(f"- {candidate}")

    typed = input(
        f"Enter the full path to the {description} "
        f"(or press Enter to cancel): ").strip().strip('"')

    if typed and os.path.exists(typed):
        return os.path.normpath(typed)

    if typed:
        print("That path does not exist either.")

    return None


# ============================================================
# SEQUENCE BASICS
# ============================================================

def amino_acid_counts(sequence: str) -> dict[str, int]:
    """Count every known residue code in one pass per code."""

    counts = {}

    for code in KNOWN_CODES:
        counts[code] = sequence.count(code)

    return counts


def standard_counts(counts: dict[str, int]) -> dict[str, int]:
    return {code: counts[code] for code in AMINO_ACIDS}


def non_standard_counts(counts: dict[str, int]) -> dict[str, int]:
    return {
        code: counts[code]
        for code in NON_STANDARD_CODES
        if counts[code] > 0}


def unrecognised_count(sequence: str, counts: dict[str, int]) -> int:
    known_total = sum(counts.values())

    return len(sequence) - known_total


def amino_acid_composition(sequence: str) -> dict[str, float]:
    counts = amino_acid_counts(sequence)
    length = len(sequence)

    if length == 0:
        return {code: 0.0 for code in AMINO_ACIDS}

    return {
        code: (counts[code] / length) * 100
        for code in AMINO_ACIDS}


def most_abundant_amino_acids(counts: dict[str, int]) -> list[str]:
    standard = standard_counts(counts)
    maximum = max(standard.values())

    if maximum == 0:
        return []

    return [code for code in AMINO_ACIDS if standard[code] == maximum]


def least_abundant_amino_acids(counts: dict[str, int]) -> list[str]:
    standard = standard_counts(counts)

    present = [count for count in standard.values() if count > 0]

    if not present:
        return []

    minimum = min(present)

    return [code for code in AMINO_ACIDS if standard[code] == minimum]


def validate_protein(sequence: str) -> tuple[bool, str]:
    """Return (is_usable, message).

    Non-standard UniProt codes such as U or X are accepted, because
    they really occur in the database. Only unexpected characters
    make a sequence unusable.
    """

    if not sequence:
        return False, "Sequence is empty."

    counts = amino_acid_counts(sequence)
    unknown = unrecognised_count(sequence, counts)

    if unknown > 0:
        strange = sorted(set(sequence) - set(KNOWN_CODES))

        return False, (
            f"Contains {unknown} unexpected character(s): "
            f"{', '.join(strange)}")

    extras = non_standard_counts(counts)

    if extras:
        listed = ", ".join(
            f"{format_residue(code)} x{count}"
            for code, count in extras.items())

        return True, (
            f"Valid sequence, including non-standard code(s): {listed}")

    return True, "Valid protein sequence (20 standard amino acids)."


# ============================================================
# SEARCH AND SELECTION
# ============================================================

def search_proteins(proteins: list[ProteinRecord], query: str) -> list[ProteinRecord]:
    query = query.strip().lower()

    if not query:
        return []

    return [
        protein for protein in proteins
        if query in protein["search_key"]]


def display_search_results(matches: list[dict[str, Any]], start: int = 0) -> None:
    total = len(matches)
    end = min(start + MAX_DISPLAY_RESULTS, total)

    print(f"\nSearch Results (showing {start + 1}-{end} of {total:,}):")

    print(
        f"{'No.':>4}  {'UniProt ID':<12}  "
        f"{'Gene':<10}  {'Protein Name':<45}  {'Length':>7}")

    print(
        f"{'-' * 4}  {'-' * 12}  "
        f"{'-' * 10}  {'-' * 45}  {'-' * 7}")

    for offset in range(start, end):
        protein = matches[offset]

        print(
            f"{offset + 1:>4}  "
            f"{shorten(protein['id'], 12):<12}  "
            f"{shorten(protein['gene'], 10):<10}  "
            f"{shorten(protein['name'], 45):<45}  "
            f"{len(protein['sequence']):>7,}")

    if total > end:
        print(
            f"\n{total - end:,} more result(s) not shown. "
            f"Type 'n' for the next page.")

    if total > 1000:
        print(
            "Tip: searching by the exact UniProt ID or gene name "
            "gives a far shorter list.")


def select_from_search_results(matches: list[dict[str, Any]]) -> Optional[ProteinRecord]:
    if not matches:
        return None

    if len(matches) == 1:
        protein = matches[0]

        print(f"\nAutomatically selected: {protein['id']}")

        return protein

    start = 0

    while True:
        display_search_results(matches, start)

        choice = input(
            "\nEnter result number, UniProt ID, 'n' for next page, "
            "'p' for previous page, or 0 to cancel: ").strip()

        if choice == "0":
            return None

        if choice.lower() == "n":
            if start + MAX_DISPLAY_RESULTS < len(matches):
                start += MAX_DISPLAY_RESULTS
            else:
                print("\nYou are already on the last page.")

            continue

        if choice.lower() == "p":
            if start >= MAX_DISPLAY_RESULTS:
                start -= MAX_DISPLAY_RESULTS
            else:
                print("\nYou are already on the first page.")

            continue

        if choice.isdigit():
            number = int(choice)

            if 1 <= number <= len(matches):
                return matches[number - 1]

            print(
                f"Invalid selection. "
                f"Enter a number from 1 to {len(matches):,}.")

            continue

        for protein in matches:
            if protein["id"].lower() == choice.lower():
                return protein

        print("UniProt ID not found in the search results.")


def search_and_select_protein(proteins: list[ProteinRecord], prompt: Optional[str] = None) -> Optional[ProteinRecord]:
    if prompt is None:
        prompt = "\nEnter UniProt ID, gene name or keyword: "

    query = input(prompt).strip()

    if not query:
        print("Search query cannot be empty.")
        return None

    matches = search_proteins(proteins, query)

    if not matches:
        if DATASET_SOURCE_PATH is None:
            try:
                print(f"\nSearching all UniProt proteins for '{query}'...")
                remote_matches, total = fetch_remote_proteins(
                    query, 0, REMOTE_PAGE_SIZE)

                known_ids = {protein["id"] for protein in proteins}
                proteins.extend(
                    protein for protein in remote_matches
                    if protein["id"] not in known_ids)
                DATASET_CACHE.clear()
                print(
                    f"Loaded {len(remote_matches):,} result(s) from "
                    f"{total:,} matching UniProt record(s).")

                if remote_matches:
                    return select_from_search_results(remote_matches)

                if " " not in query:
                    protein = fetch_remote_protein(query)
                    proteins.append(protein)
                    DATASET_CACHE.clear()
                    return protein
            except RuntimeError as error:
                print(f"\nUnable to fetch '{query}': {error}")

        print("\nNo matching proteins found.")
        return None

    return select_from_search_results(matches)


def search_prosite_motifs(motifs: list[MotifRecord], query: str) -> list[MotifRecord]:
    query = query.strip().lower()

    if not query:
        return []

    return [motif for motif in motifs if query in motif["search_key"]]


def display_motif_search_results(matches: list[dict[str, Any]], start: int = 0) -> None:
    total = len(matches)
    end = min(start + MAX_DISPLAY_RESULTS, total)

    print(f"\nPROSITE Motifs (showing {start + 1}-{end} of {total:,}):")

    print(
        f"{'No.':>4}  {'PROSITE ID':<26}  "
        f"{'Accession':<10}  {'Description':<50}")

    print(
        f"{'-' * 4}  {'-' * 26}  {'-' * 10}  {'-' * 50}")

    for offset in range(start, end):
        motif = matches[offset]

        print(
            f"{offset + 1:>4}  "
            f"{shorten(motif['id'], 26):<26}  "
            f"{shorten(motif['accession'], 10):<10}  "
            f"{shorten(motif['description'], 50):<50}")

    if total > end:
        print(
            f"\n{total - end:,} more motif(s) not shown. "
            f"Type 'n' for the next page.")


def select_prosite_motif(matches: list[dict[str, Any]]) -> Optional[MotifRecord]:
    if not matches:
        return None

    if len(matches) == 1:
        motif = matches[0]

        print(f"\nAutomatically selected PROSITE motif: {motif['id']}")

        return motif

    start = 0

    while True:
        display_motif_search_results(matches, start)

        choice = input(
            "\nEnter motif number, PROSITE ID, accession, "
            "'n' for next page, 'p' for previous, or 0 to cancel: "
        ).strip()

        if choice == "0":
            return None

        if choice.lower() == "n":
            if start + MAX_DISPLAY_RESULTS < len(matches):
                start += MAX_DISPLAY_RESULTS
            else:
                print("\nYou are already on the last page.")

            continue

        if choice.lower() == "p":
            if start >= MAX_DISPLAY_RESULTS:
                start -= MAX_DISPLAY_RESULTS
            else:
                print("\nYou are already on the first page.")

            continue

        if choice.isdigit():
            number = int(choice)

            if 1 <= number <= len(matches):
                return matches[number - 1]

            print(
                f"Invalid selection. "
                f"Enter a number from 1 to {len(matches):,}.")

            continue

        for motif in matches:
            if (
                motif["id"].lower() == choice.lower()
                or motif["accession"].lower() == choice.lower()
            ):
                return motif

        print("PROSITE ID or accession not found.")


def search_and_select_motif(motifs: list[MotifRecord]) -> Optional[MotifRecord]:
    query = input(
        "\nEnter PROSITE ID, accession or keyword: ").strip()

    if not query:
        print("Motif search query cannot be empty.")
        return None

    matches = search_prosite_motifs(motifs, query)

    if not matches:
        print("\nNo matching PROSITE motifs found.")
        return None

    return select_prosite_motif(matches)


# ============================================================
# BASIC PROTEIN ANALYSIS
# ============================================================

def composition_table(sequence: str) -> pd.DataFrame:
    counts = amino_acid_counts(sequence)
    length = len(sequence)

    if length == 0:
        return pd.DataFrame(
            columns=["Code", "Amino Acid", "Class", "Count", "Percentage"])

    rows = []

    for code in KNOWN_CODES:
        if code not in AMINO_ACIDS and counts[code] == 0:
            continue

        rows.append({
            "Code": code,
            "Amino Acid": AMINO_ACID_NAMES[code],
            "Class": AMINO_ACID_CLASSES.get(code, "Non-standard"),
            "Count": counts[code],
            "Percentage": round(counts[code] / length * 100, 2)})

    table = pd.DataFrame(rows)

    return table.sort_values(
        by="Count", ascending=False).reset_index(drop=True)


def class_summary_table(sequence: str) -> pd.DataFrame:
    counts = amino_acid_counts(sequence)
    length = len(sequence)

    if length == 0:
        return pd.DataFrame(
            columns=["Class", "Residues", "Count", "Percentage"])

    rows = []

    for class_name in CLASS_ORDER:
        members = [
            code for code in AMINO_ACIDS
            if AMINO_ACID_CLASSES[code] == class_name]

        total = sum(counts[code] for code in members)

        rows.append({
            "Class": class_name,
            "Residues": "".join(members),
            "Count": total,
            "Percentage": round(total / length * 100, 2)})

    return pd.DataFrame(rows)


def display_composition_chart(sequence: str, title: str) -> None:
    counts = standard_counts(amino_acid_counts(sequence))

    series = pd.Series(counts).sort_values(ascending=False)
    series = series[series > 0]

    if series.empty:
        print("\nNo composition data available for plotting.")
        return

    length = len(sequence)
    percentages = series.to_numpy(dtype=float) / length * 100

    codes = list(series.index)

    colors = [
        CLASS_COLORS[AMINO_ACID_CLASSES[code]] for code in codes]

    legend_labels = [
        f"{code} - {AMINO_ACID_NAMES[code]} ({percentages[i]:.1f}%)"
        for i, code in enumerate(codes)]

    figure, axes = plt.subplots(figsize=(11, 8))
    figure.patch.set_facecolor("white")

    wedges, _, autotexts = axes.pie(
        percentages,
        labels=codes,
        colors=colors,
        autopct=lambda pct: f"{pct:.1f}%" if pct >= 3 else "",
        pctdistance=0.78,
        startangle=90,
        wedgeprops={"edgecolor": "white", "linewidth": 1.0},
        textprops={"fontsize": 10, "fontweight": "bold"})

    for autotext in autotexts:
        autotext.set_color("white")
        autotext.set_fontsize(8)

    axes.set_title(title, fontsize=14, fontweight="bold", pad=16)
    axes.axis("equal")

    class_handles = [
        plt.Rectangle((0, 0), 1, 1, color=CLASS_COLORS[name])
        for name in CLASS_ORDER]

    class_legend = axes.legend(
        class_handles, CLASS_ORDER,
        title="Side-chain class", fontsize=9, title_fontsize=9,
        loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)

    axes.add_artist(class_legend)

    axes.legend(
        wedges, legend_labels,
        title="Amino Acid (Full Name) - Share",
        loc="lower left", bbox_to_anchor=(1.02, 0.0),
        fontsize=8, title_fontsize=9, frameon=False)

    print("\nLegend:", "; ".join(legend_labels))

    plt.tight_layout()
    show_and_close(figure)


def display_class_chart(sequence: str, title: str) -> None:
    table = class_summary_table(sequence)

    values = table["Percentage"].to_numpy(dtype=float)
    colors = [CLASS_COLORS[name] for name in table["Class"]]

    figure, axes = plt.subplots(figsize=(8, 8))
    figure.patch.set_facecolor("white")

    wedges, _, autotexts = axes.pie(
        values,
        labels=table["Class"],
        colors=colors,
        explode=[0.03] * len(table),
        autopct="%1.1f%%",
        pctdistance=0.74,
        startangle=90,
        wedgeprops={"edgecolor": "white", "linewidth": 1.5},
        textprops={"fontsize": 11, "fontweight": "bold"})

    for autotext in autotexts:
        autotext.set_color("white")
        autotext.set_fontsize(11)

    print(
        "\nLegend:",
        "; ".join(
            f"{row['Class']} ({row['Residues']}) - "
            f"{row['Count']:,} residues, {row['Percentage']:.2f}%"
            for _, row in table.iterrows()))

    axes.set_title(title, fontsize=14, fontweight="bold", pad=16)
    axes.axis("equal")

    plt.tight_layout()
    show_and_close(figure)


def display_basic_analysis(protein: ProteinRecord) -> None:
    sequence = protein["sequence"]

    usable, message = validate_protein(sequence)

    print_header("BASIC PROTEIN ANALYSIS")

    print(f"UniProt ID      : {protein['id']}")
    print(f"Protein Name    : {protein['name']}")
    print(f"Gene            : {protein['gene'] or 'Not listed'}")
    print(f"Organism        : {protein['organism'] or 'Not listed'}")
    print(f"Sequence Length : {len(sequence):,} amino acids")
    print(f"Validation      : {message}")

    if not usable:
        return

    counts = amino_acid_counts(sequence)

    print(
        "\nMost Abundant Amino Acid(s) :",
        format_residue_list(most_abundant_amino_acids(counts)))

    print(
        "Least Abundant Amino Acid(s):",
        format_residue_list(least_abundant_amino_acids(counts)))

    print_section("Amino Acid Composition (sorted by abundance)")
    print(composition_table(sequence).to_string(index=False))

    print_section("Composition by Side-Chain Class")
    print(class_summary_table(sequence).to_string(index=False))

    if ask_yes_no("\nDisplay composition pie chart?"):
        display_composition_chart(
            sequence,
            f"Amino Acid Composition - {protein['id']}")

    if ask_yes_no("Display side-chain class chart?"):
        display_class_chart(
            sequence,
            f"Side-Chain Classes - {protein['id']}")


# ============================================================
# MOLECULAR WEIGHT
# ============================================================

def molecular_weight_details(sequence: str) -> Optional[dict[str, Any]]:
    """Return the molecular weight together with quality information."""

    if not sequence:
        return None

    counts = amino_acid_counts(sequence)

    total = 0.0
    residues_used = 0
    approximated = 0

    for code in KNOWN_CODES:
        count = counts[code]

        if count == 0:
            continue

        residue_mass = AMINO_ACID_MASSES[code] - WATER_MASS

        total += count * residue_mass
        residues_used += count

        if code in APPROXIMATE_CODES:
            approximated += count

    skipped = unrecognised_count(sequence, counts)

    if residues_used == 0:
        return None

    return {
        "weight": total + WATER_MASS,
        "residue_subtotal": total,
        "residues_used": residues_used,
        "approximated": approximated,
        "skipped": skipped}


def protein_molecular_weight(sequence: str) -> Optional[float]:
    details = molecular_weight_details(sequence)

    if details is None:
        return None

    return details["weight"]


def average_residue_mass(sequence: str) -> Optional[float]:
    details = molecular_weight_details(sequence)

    if details is None:
        return None

    return details["weight"] / details["residues_used"]


def molecular_weight_breakdown(sequence: str) -> Optional[dict[str, Any]]:
    details = molecular_weight_details(sequence)

    if details is None:
        return None

    counts = amino_acid_counts(sequence)
    total_weight = details["weight"]

    rows = []

    for code in KNOWN_CODES:
        count = counts[code]

        if count == 0:
            continue

        residue_mass = AMINO_ACID_MASSES[code] - WATER_MASS
        contribution = count * residue_mass

        rows.append({
            "Code": code,
            "Amino Acid": AMINO_ACID_NAMES[code],
            "Count": count,
            "Residue Mass (Da)": round(residue_mass, 2),
            "Weight Contributed (Da)": round(contribution, 2),
            "Percentage": round(contribution / total_weight * 100, 2)})

    table = pd.DataFrame(rows).sort_values(
        by="Weight Contributed (Da)",
        ascending=False).reset_index(drop=True)

    return {"details": details, "table": table}


def display_molecular_weight(protein: ProteinRecord) -> None:
    sequence = protein["sequence"]
    details = molecular_weight_details(sequence)

    if details is None:
        print("\nUnable to calculate molecular weight.")
        return

    print_header("MOLECULAR WEIGHT")

    print(f"UniProt ID      : {protein['id']}")
    print(f"Protein Name    : {protein['name']}")
    print(f"Sequence Length : {len(sequence):,} amino acids")

    print(
        f"\nMolecular Weight     : {details['weight']:,.2f} Da")

    print(
        f"Molecular Weight     : "
        f"{details['weight'] / 1000:,.2f} kDa")

    print(
        f"Average Residue Mass : "
        f"{details['weight'] / details['residues_used']:,.2f} Da")

    if details["approximated"]:
        print(
            f"\nNote: {details['approximated']} residue(s) use an "
            f"average mass (codes X, B, Z or J), so the result is "
            f"an estimate.")

    if details["skipped"]:
        print(
            f"Note: {details['skipped']} unexpected character(s) "
            f"were ignored.")


def display_molecular_weight_breakdown(protein: ProteinRecord) -> None:
    breakdown = molecular_weight_breakdown(protein["sequence"])

    if breakdown is None:
        print("\nUnable to calculate a molecular weight breakdown.")
        return

    details = breakdown["details"]

    print_header("MOLECULAR WEIGHT BREAKDOWN")

    print(f"UniProt ID      : {protein['id']}")
    print(f"Protein Name    : {protein['name']}")
    print(f"Sequence Length : {len(protein['sequence']):,} amino acids")

    print_section("Weight Contributed by Each Amino Acid Type")
    print(breakdown["table"].to_string(index=False))

    print(
        f"\nSum of residue contributions        : "
        f"{details['residue_subtotal']:,.2f} Da")

    print(
        f"+ One water molecule (chain ends)    : "
        f"{WATER_MASS:.3f} Da")

    print(
        f"= Total molecular weight             : "
        f"{details['weight']:,.2f} Da "
        f"({details['weight'] / 1000:,.2f} kDa)")


def dataset_weights(proteins: list[ProteinRecord]) -> np.ndarray:
    """Molecular weight of every protein, cached after the first run."""

    if "weights" in DATASET_CACHE:
        return DATASET_CACHE["weights"]

    print("\nCalculating molecular weights for the whole dataset...")

    weights = []

    for index, protein in enumerate(proteins, start=1):
        weight = protein_molecular_weight(protein["sequence"])

        if weight is not None:
            weights.append(weight)

        if index % 20000 == 0:
            print(f"  processed {index:,} proteins...")

    array = np.array(weights, dtype=np.float64)

    DATASET_CACHE["weights"] = array

    return array


def describe_array(array: np.ndarray, unit: str) -> dict[str, str]:
    return {
        "Count": f"{len(array):,}",
        "Mean": f"{np.mean(array):,.2f} {unit}",
        "Median": f"{np.median(array):,.2f} {unit}",
        "Minimum": f"{np.min(array):,.2f} {unit}",
        "Maximum": f"{np.max(array):,.2f} {unit}",
        "Std. deviation (sample)": (
            f"{np.std(array, ddof=1):,.2f} {unit}"),
        "25th percentile": f"{np.percentile(array, 25):,.2f} {unit}",
        "75th percentile": f"{np.percentile(array, 75):,.2f} {unit}"}


def print_statistics(title: str, array: np.ndarray, unit: str) -> None:
    if array.size == 0:
        print("\nNo data available.")
        return

    print_section(title)

    for label, value in describe_array(array, unit).items():
        print(f"{label:<24}: {value}")


def display_histogram(array: np.ndarray, title: str, x_label: str, color: str, log_x: bool = False) -> None:
    data = np.asarray(array, dtype=np.float64)

    if log_x:
        data = data[data > 0]

    if data.size == 0:
        print("\nNo positive data available for this chart.")
        return

    figure, axes = prepare_plot()

    if log_x:
        minimum = float(np.min(data))
        maximum = float(np.max(data))
        if minimum == maximum:
            minimum *= 0.9
            maximum *= 1.1
        bins = np.logspace(np.log10(minimum), np.log10(maximum), 60)
    else:
        bins = 60

    axes.hist(
        data, bins=bins, color=color,
        edgecolor="white", linewidth=0.5)

    if log_x:
        axes.set_xscale("log")

    mean_value = float(np.mean(data))
    median_value = float(np.median(data))

    axes.axvline(
        mean_value, color=PLOT_COLORS["red"], linestyle="--", linewidth=1.6,
        label=f"Mean = {mean_value:,.1f}")

    axes.axvline(
        median_value, color=PLOT_COLORS["ink"], linestyle=":", linewidth=1.6,
        label=f"Median = {median_value:,.1f}")

    axes.legend(frameon=False, fontsize=9)

    finish_plot(axes, title, x_label, "Number of proteins")


def molecular_weight_menu(proteins: list[ProteinRecord]) -> None:
    while True:
        choice = show_menu(
            "MOLECULAR WEIGHT ANALYSIS",
            ["Single Protein Molecular Weight",
             "Molecular Weight Breakdown by Amino Acid",
             "Dataset Molecular Weight Statistics",
             "Molecular Weight Distribution (chart)",
             "Go Back",
             "Exit"])

        if choice == 1:
            protein = search_and_select_protein(proteins)

            if protein is not None:
                display_molecular_weight(protein)

        elif choice == 2:
            protein = search_and_select_protein(proteins)

            if protein is not None:
                display_molecular_weight_breakdown(protein)

        elif choice == 3:
            weights = dataset_weights(proteins)
            print_statistics(
                "Molecular Weight Statistics", weights, "Da")

        elif choice == 4:
            weights = dataset_weights(proteins)

            if weights.size:
                display_histogram(
                    weights,
                    "Molecular Weight Distribution",
                    "Molecular weight (Da, log scale)",
                    "#8e44ad", log_x=True)

        elif choice == 5:
            return

        elif choice == 6:
            exit_program()


# ============================================================
# CHARGE, NET CHARGE AND ISOELECTRIC POINT
# ============================================================

def charge_composition(sequence: str) -> Optional[dict[str, Any]]:
    if not sequence:
        return None

    counts = amino_acid_counts(sequence)
    length = len(sequence)

    total_acidic = sum(counts[code] for code in ACIDIC_RESIDUES)
    total_basic = sum(counts[code] for code in BASIC_RESIDUES)
    total_uncharged = sum(counts[code] for code in UNCHARGED_RESIDUES)
    other = length - total_acidic - total_basic - total_uncharged

    return {
        "length": length,
        "counts": counts,
        "total_acidic": total_acidic,
        "total_basic": total_basic,
        "total_charged": total_acidic + total_basic,
        "total_uncharged": total_uncharged,
        "other": other,
        "acidic_percentage": total_acidic / length * 100,
        "basic_percentage": total_basic / length * 100,
        "charged_percentage": (
            (total_acidic + total_basic) / length * 100),
        "uncharged_percentage": total_uncharged / length * 100}


def net_charge_at_ph(sequence: str, ph: float) -> float:
    counts = amino_acid_counts(sequence)

    charge = 1.0 / (1.0 + 10 ** (ph - PKA_N_TERMINUS))

    for code, pka in PKA_POSITIVE.items():
        charge += counts[code] / (1.0 + 10 ** (ph - pka))

    charge -= 1.0 / (1.0 + 10 ** (PKA_C_TERMINUS - ph))

    for code, pka in PKA_NEGATIVE.items():
        charge -= counts[code] / (1.0 + 10 ** (pka - ph))

    return charge


def isoelectric_point(sequence: str) -> Optional[float]:
    """Find the pH where the net charge is zero, by bisection."""

    if not sequence:
        return None

    low = 0.0
    high = 14.0

    for _ in range(100):
        middle = (low + high) / 2.0
        charge = net_charge_at_ph(sequence, middle)

        if charge > 0:
            low = middle
        else:
            high = middle

    return (low + high) / 2.0


def display_charge_curve(protein: ProteinRecord) -> None:
    sequence = protein["sequence"]

    ph_values = np.arange(0.0, 14.01, 0.1)

    charges = np.array(
        [net_charge_at_ph(sequence, float(ph)) for ph in ph_values])

    pi_value = isoelectric_point(sequence)

    figure, axes = prepare_plot()

    axes.plot(ph_values, charges, color=PLOT_COLORS["blue"], linewidth=2.2)

    axes.fill_between(
        ph_values, charges, 0,
        where=(charges > 0), color=PLOT_COLORS["blue"], alpha=0.25)

    axes.fill_between(
        ph_values, charges, 0,
        where=(charges < 0), color=PLOT_COLORS["red"], alpha=0.25)

    axes.axhline(0, color=PLOT_COLORS["ink"], linewidth=1.2)

    axes.axvline(
        pi_value, color=PLOT_COLORS["green"], linestyle="--", linewidth=1.8,
        label=f"pI = {pi_value:.2f}")

    axes.axvline(
        7.0, color=PLOT_COLORS["muted"], linestyle=":", linewidth=1.4,
        label="pH 7.0")

    axes.legend(frameon=False, fontsize=10)

    finish_plot(
        axes,
        f"Net Charge vs pH - {protein['id']}",
        "pH", "Net charge (elementary charges)")


def display_charge_analysis(protein: ProteinRecord) -> None:
    sequence = protein["sequence"]
    analysis = charge_composition(sequence)

    if analysis is None:
        print("\nUnable to perform charge analysis.")
        return

    counts = analysis["counts"]

    print_header("CHARGE AND IONISATION ANALYSIS")

    print(f"UniProt ID      : {protein['id']}")
    print(f"Protein Name    : {protein['name']}")
    print(f"Sequence Length : {analysis['length']:,} amino acids")

    print_section("Acidic Residues")

    for code in ACIDIC_RESIDUES:
        print(f"{format_residue(code):<22}: {counts[code]:,}")

    print(
        f"{'Total acidic':<22}: {analysis['total_acidic']:,} "
        f"({analysis['acidic_percentage']:.2f}%)")

    print_section("Basic / Ionisable Residues")

    for code in BASIC_RESIDUES:
        print(f"{format_residue(code):<22}: {counts[code]:,}")

    print(
        f"{'Total basic':<22}: {analysis['total_basic']:,} "
        f"({analysis['basic_percentage']:.2f}%)")

    print_section("Overall Composition")

    print(
        f"{'Charged / ionisable':<22}: "
        f"{analysis['total_charged']:,} "
        f"({analysis['charged_percentage']:.2f}%)")

    print(
        f"{'Uncharged':<22}: {analysis['total_uncharged']:,} "
        f"({analysis['uncharged_percentage']:.2f}%)")

    if analysis["other"]:
        print(
            f"{'Non-standard / other':<22}: {analysis['other']:,}")

    charge_difference = (
        analysis["total_basic"] - analysis["total_acidic"])

    print(
        f"{'Basic minus acidic':<22}: {charge_difference:+,}")

    print_section("Charge Calculated from pKa Values")

    pi_value = isoelectric_point(sequence)

    print(f"{'Isoelectric point (pI)':<26}: {pi_value:.2f}")

    print(
        f"{'Net charge at pH 7.0':<26}: "
        f"{net_charge_at_ph(sequence, 7.0):+.2f}")

    print(
        f"{'Net charge at pH 5.0':<26}: "
        f"{net_charge_at_ph(sequence, 5.0):+.2f}")

    print(
        f"{'Net charge at pH 9.0':<26}: "
        f"{net_charge_at_ph(sequence, 9.0):+.2f}")

    if pi_value < 7.0:
        print(
            "\nInterpretation: the protein is acidic overall and "
            "carries a negative net charge at neutral pH.")
    else:
        print(
            "\nInterpretation: the protein is basic overall and "
            "carries a positive net charge at neutral pH.")

    if ask_yes_no("\nDisplay the net charge vs pH curve?"):
        display_charge_curve(protein)

    if ask_yes_no("Display charge composition chart?"):
        display_class_chart(
            sequence,
            f"Side-Chain Classes - {protein['id']}")


# ============================================================
# OTHER PHYSICOCHEMICAL PROPERTIES
# ============================================================

def extinction_coefficient(sequence: str) -> dict[str, int]:
    counts = amino_acid_counts(sequence)

    reduced = (
        counts["W"] * EXTINCTION_TRYPTOPHAN
        + counts["Y"] * EXTINCTION_TYROSINE)

    cystines = counts["C"] // 2

    oxidised = reduced + cystines * EXTINCTION_CYSTINE

    return {
        "reduced": reduced,
        "oxidised": oxidised,
        "cystines": cystines}


def aromaticity(sequence: str) -> float:
    if not sequence:
        return 0.0

    counts = amino_acid_counts(sequence)

    total = sum(counts[code] for code in AROMATIC_RESIDUES)

    return total / len(sequence) * 100


def display_general_properties(protein: ProteinRecord) -> None:
    sequence = protein["sequence"]

    weight = protein_molecular_weight(sequence)
    coefficients = extinction_coefficient(sequence)

    print_header("GENERAL PHYSICOCHEMICAL SUMMARY")

    print(f"UniProt ID      : {protein['id']}")
    print(f"Protein Name    : {protein['name']}")

    rows = [
        ("Length", f"{len(sequence):,} amino acids"),
        ("Molecular weight", f"{weight:,.2f} Da"),
        ("Molecular weight", f"{weight / 1000:,.2f} kDa"),
        ("Isoelectric point (pI)", f"{isoelectric_point(sequence):.2f}"),
        ("Net charge at pH 7.0",
         f"{net_charge_at_ph(sequence, 7.0):+.2f}"),
        ("GRAVY (Kyte-Doolittle)", f"{gravy_score(sequence):.3f}"),
        ("Aromaticity", f"{aromaticity(sequence):.2f}%"),
        ("Extinction coeff. (reduced)",
         f"{coefficients['reduced']:,} M-1 cm-1"),
        ("Extinction coeff. (cystines)",
         f"{coefficients['oxidised']:,} M-1 cm-1")]

    print_section("Calculated Properties")

    for label, value in rows:
        print(f"{label:<30}: {value}")

    if coefficients["reduced"] == 0:
        print(
            "\nNote: no tryptophan or tyrosine is present, so this "
            "protein absorbs very little light at 280 nm.")


def physicochemical_menu(proteins: list[ProteinRecord]) -> None:
    while True:
        choice = show_menu(
            "PHYSICOCHEMICAL PROPERTIES",
            ["General Property Summary",
             "Molecular Weight Analysis",
             "Charge and Isoelectric Point",
             "Go Back",
             "Exit"])

        if choice == 1:
            protein = search_and_select_protein(proteins)

            if protein is not None:
                display_general_properties(protein)

        elif choice == 2:
            molecular_weight_menu(proteins)

        elif choice == 3:
            protein = search_and_select_protein(proteins)

            if protein is not None:
                display_charge_analysis(protein)

        elif choice == 4:
            return

        elif choice == 5:
            exit_program()


# ============================================================
# PROSITE PATTERN CONVERSION
# ============================================================

def prosite_token_to_regex(token: str) -> Optional[str]:
    """Convert one hyphen-separated PROSITE element into regex text."""

    repeat = ""

    repeat_match = re.fullmatch(
        r"(.+?)\((\d+)(?:,(\d+))?\)", token)

    if repeat_match:
        element_text = repeat_match.group(1)
        minimum = repeat_match.group(2)
        maximum = repeat_match.group(3)

        if maximum is None:
            repeat = "{" + minimum + "}"
        else:
            repeat = "{" + minimum + "," + maximum + "}"

    else:
        element_text = token

    element_text = element_text.upper()

    if element_text == "X":
        return "." + repeat

    if re.fullmatch(r"\[[A-Z<>]+\]", element_text):
        # Terminal markers inside a class are ignored.
        members = element_text[1:-1].replace("<", "").replace(">", "")

        if not members:
            return None

        return "[" + members + "]" + repeat

    if re.fullmatch(r"\{[A-Z]+\}", element_text):
        return "[^" + element_text[1:-1] + "]" + repeat

    if re.fullmatch(r"[A-Z]", element_text):
        return element_text + repeat

    return None


def prosite_pattern_to_regex(pattern: str) -> Optional[str]:
    """Convert a full PROSITE pattern into a regular expression."""

    pattern = pattern.strip().replace(" ", "")

    if not pattern:
        return None

    if pattern.endswith("."):
        pattern = pattern[:-1]

    n_terminal = pattern.startswith("<")
    c_terminal = pattern.endswith(">")

    if n_terminal:
        pattern = pattern[1:]

    if c_terminal:
        pattern = pattern[:-1]

    parts = []

    for token in pattern.split("-"):
        if not token:
            continue

        converted = prosite_token_to_regex(token)

        if converted is None:
            return None

        parts.append(converted)

    if not parts:
        return None

    regex = "".join(parts)

    if n_terminal:
        regex = "^" + regex

    if c_terminal:
        regex = regex + "$"

    return regex


def compile_prosite_pattern(pattern: str) -> Optional[re.Pattern]:
    """Compile a PROSITE pattern so that overlapping hits are found.

    The pattern is wrapped in a look-ahead group. A plain search
    consumes each match, which would hide motifs that overlap, and
    overlapping motifs are common (N-glycosylation sites, for
    example).
    """

    regex_text = prosite_pattern_to_regex(pattern)

    if regex_text is None:
        return None

    try:
        return re.compile("(?=(" + regex_text + "))")

    except re.error:
        return None


def find_pattern_occurrences(sequence: str, regex: re.Pattern) -> list[dict[str, Any]]:
    matches = []

    for match in regex.finditer(sequence):
        text = match.group(1)

        if not text:
            continue

        matches.append({
            "start": match.start() + 1,
            "end": match.start() + len(text),
            "sequence": text})

    return matches


def find_motif_matches(sequence: str, motif: MotifRecord) -> Optional[list[dict[str, Any]]]:
    regex = compile_prosite_pattern(motif["pattern"])

    if regex is None:
        return None

    return find_pattern_occurrences(sequence, regex)


# ============================================================
# MOTIF DETECTION
# ============================================================

def display_motif_information(motif: MotifRecord) -> None:
    print_header("PROSITE MOTIF INFORMATION")

    print(f"PROSITE ID  : {motif['id']}")
    print(f"Accession   : {motif['accession']}")
    print(f"Type        : {motif['type'] or 'PATTERN'}")
    print(f"Description : {motif['description']}")
    print(f"Pattern     : {motif['pattern']}")

    regex_text = prosite_pattern_to_regex(motif["pattern"])

    if regex_text is None:
        print(
            "\nThis pattern uses syntax the converter does not "
            "support, so it cannot be searched for.")
    else:
        print(f"\nRegular expression: {regex_text}")


def display_motif_matches(protein: ProteinRecord, motif: MotifRecord, matches: list[dict[str, Any]]) -> None:
    print_header("MOTIF DETECTION RESULT")

    print(f"PROSITE ID  : {motif['id']}")
    print(f"Accession   : {motif['accession']}")
    print(f"Description : {motif['description']}")
    print(f"Pattern     : {motif['pattern']}")

    print_section("Protein")

    print(f"UniProt ID   : {protein['id']}")
    print(f"Protein Name : {protein['name']}")
    print(
        f"Length       : "
        f"{len(protein['sequence']):,} amino acids")

    if not matches:
        print("\nMatches found: 0")
        print("This motif does not occur in this protein.")
        return

    print(f"\nMatches found: {len(matches):,}")

    table = pd.DataFrame([
        {"No.": index,
         "Start": match["start"],
         "End": match["end"],
         "Length": len(match["sequence"]),
         "Matched Sequence": shorten(match["sequence"], 45)}
        for index, match in enumerate(matches, start=1)])

    print_section("Match Positions")

    if len(table) > MAX_DISPLAY_RESULTS:
        print(table.head(MAX_DISPLAY_RESULTS).to_string(index=False))

        print(
            f"\n{len(table) - MAX_DISPLAY_RESULTS:,} further "
            f"match(es) not shown.")
    else:
        print(table.to_string(index=False))

    if len(matches) > 1 and ask_yes_no(
            "\nDisplay a map of match positions?"):
        display_motif_position_map(protein, motif, matches)


def display_motif_position_map(protein: ProteinRecord, motif: MotifRecord, matches: list[dict[str, Any]]) -> None:
    length = len(protein["sequence"])

    figure, axes = prepare_plot()

    axes.hlines(
        1, 0, length, color=PLOT_COLORS["track"], linewidth=10)

    starts = [match["start"] for match in matches]
    widths = [
        match["end"] - match["start"] + 1 for match in matches]

    axes.barh(
        [1] * len(matches), widths, left=starts,
        height=0.42, color=PLOT_COLORS["accent"],
        edgecolor=PLOT_COLORS["red"])

    axes.set_ylim(0.5, 1.5)
    axes.set_yticks([])
    axes.set_xlim(0, length)

    finish_plot(
        axes,
        f"{motif['id']} positions in {protein['id']}",
        "Residue position", "")


def detect_motif_in_protein(proteins: list[ProteinRecord], motif: MotifRecord) -> None:
    protein = search_and_select_protein(proteins)

    if protein is None:
        return

    matches = find_motif_matches(protein["sequence"], motif)

    if matches is None:
        print(
            "\nThis PROSITE pattern uses syntax the converter does "
            "not support.")

        print(f"Pattern: {motif['pattern']}")
        return

    display_motif_matches(protein, motif, matches)


def column_width(header: str, values: list[str]) -> int:
    """Widest string needed for one printed table column."""

    widest_value = max((len(value) for value in values), default=0)

    return max(len(header), widest_value)


def print_scan_results_page(results: list[dict[str, Any]], start: int, end: int) -> None:
    """Print one page of dataset motif-scan results.

    UniProt ID, Gene and Occurrences use fixed-width columns sized
    to what's on the current page, so the ID is always shown in
    full. Protein Name is printed last and is never truncated -
    since nothing follows it, an unusually long name cannot break
    the alignment of the columns before it.
    """

    page = results[start:end]

    ids = [item["protein"]["id"] for item in page]
    genes = [item["protein"]["gene"] or "-" for item in page]
    occurrence_text = [f"{len(item['matches']):,}" for item in page]

    number_width = len(f"{end:,}")
    id_width = column_width("UniProt ID", ids)
    gene_width = column_width("Gene", genes)
    occurrence_width = column_width("Occurrences", occurrence_text)

    print(
        f"{'No.':>{number_width}}  "
        f"{'UniProt ID':<{id_width}}  "
        f"{'Gene':<{gene_width}}  "
        f"{'Occurrences':>{occurrence_width}}  "
        f"Protein Name")

    print(
        f"{'-' * number_width}  "
        f"{'-' * id_width}  "
        f"{'-' * gene_width}  "
        f"{'-' * occurrence_width}  "
        f"{'-' * 12}")

    for offset, item in enumerate(page):
        position = start + offset + 1
        protein = item["protein"]

        print(
            f"{position:>{number_width},}  "
            f"{protein['id']:<{id_width}}  "
            f"{(protein['gene'] or '-'):<{gene_width}}  "
            f"{len(item['matches']):>{occurrence_width},}  "
            f"{protein['name']}")


def browse_scan_results(results: list[dict[str, Any]], motif: MotifRecord) -> None:
    """Let the user page through the complete motif-scan result list.

    The full list stays in memory and only one page (at most
    SCAN_RESULTS_PAGE_SIZE rows) is ever printed at a time, so a
    match count in the millions never floods the terminal. Every
    result remains reachable by paging forward, and a result can be
    opened for full match-position detail from any page.
    """

    total = len(results)
    start = 0

    while True:
        end = min(start + SCAN_RESULTS_PAGE_SIZE, total)

        print_section(f"Showing results {start + 1:,}-{end:,} of {total:,}")

        print_scan_results_page(results, start, end)

        options = [f"a result number ({start + 1:,}-{end:,}) for details"]

        if end < total:
            options.append("'n' for the next page")

        if start > 0:
            options.append("'p' for the previous page")

        options.append("'q' to return to the previous menu")

        choice = input("\nEnter " + ", ".join(options) + ": ").strip().lower()

        if choice in ("q", ""):
            return

        if choice == "n":
            if end < total:
                start = end
            else:
                print("\nYou are already on the last page.")

            continue

        if choice == "p":
            if start > 0:
                start = max(0, start - SCAN_RESULTS_PAGE_SIZE)
            else:
                print("\nYou are already on the first page.")

            continue

        if choice.isdigit():
            number = int(choice)

            if start + 1 <= number <= end:
                item = results[number - 1]

                display_motif_matches(
                    item["protein"], motif, item["matches"])

                continue

            print(
                f"\nThat number is not on the current page. "
                f"Enter a number from {start + 1:,} to {end:,}, "
                f"or use n/p/q.")

            continue

        print("\nInvalid input. Enter a result number, 'n', 'p', or 'q'.")


def scan_dataset_for_motif(proteins: list[ProteinRecord], motif: MotifRecord) -> None:
    regex = compile_prosite_pattern(motif["pattern"])

    if regex is None:
        print(
            "\nThis PROSITE pattern cannot be converted into a "
            "searchable expression.")

        print(f"Pattern: {motif['pattern']}")
        return

    print_header("DATASET MOTIF SCAN")

    print(f"PROSITE ID  : {motif['id']}")
    print(f"Description : {motif['description']}")
    print(f"Pattern     : {motif['pattern']}")

    print(f"\nScanning {len(proteins):,} protein records...")

    results = []
    total_matches = 0

    for index, protein in enumerate(proteins, start=1):
        matches = find_pattern_occurrences(protein["sequence"], regex)

        if matches:
            # Each matching record keeps its own distinct UniProt ID,
            # gene name, protein name and match list. Records are
            # never merged, even when several share the same gene or
            # protein name (as happens with isoforms and paralogs) -
            # the UniProt ID is what identifies a record.
            results.append({"protein": protein, "matches": matches})
            total_matches += len(matches)

        if index % 20000 == 0:
            print(f"  scanned {index:,} proteins...")

    print("\nScan complete.")

    if not results:
        print("\nNo proteins in the dataset contain this motif.")
        return

    # Most heavily-featured proteins first.
    results.sort(key=lambda item: len(item["matches"]), reverse=True)

    print(f"\nProteins containing the motif : {len(results):,}")
    print(f"Total motif occurrences       : {total_matches:,}")

    print(
        f"Share of dataset              : "
        f"{len(results) / len(proteins) * 100:.2f}%")

    print(
        f"\nResults are sorted by motif occurrence count (highest "
        f"first) and shown {SCAN_RESULTS_PAGE_SIZE} at a time below.")

    browse_scan_results(results, motif)


def motif_detection_menu(proteins: list[ProteinRecord], motifs: list[MotifRecord]) -> None:
    while True:
        choice = show_menu(
            "MOTIF DETECTION",
            ["View Motif Information",
             "Detect Motif in a Selected Protein",
             "Scan the Whole Dataset for a Motif",
             "Go Back",
             "Exit"])

        if choice == 1:
            motif = search_and_select_motif(motifs)

            if motif is not None:
                display_motif_information(motif)

        elif choice == 2:
            motif = search_and_select_motif(motifs)

            if motif is not None:
                detect_motif_in_protein(proteins, motif)

        elif choice == 3:
            motif = search_and_select_motif(motifs)

            if motif is not None:
                scan_dataset_for_motif(proteins, motif)

        elif choice == 4:
            return

        elif choice == 5:
            exit_program()


# ============================================================
# SEQUENCE PATTERN ANALYSIS
# ============================================================

def find_subsequence(sequence: str, query: str) -> list[int]:
    matches = []
    position = sequence.find(query)

    while position != -1:
        matches.append(position + 1)
        position = sequence.find(query, position + 1)

    return matches


def custom_pattern_search(proteins: list[ProteinRecord]) -> None:
    print_header("SEQUENCE PATTERN SEARCH")

    print(
        "Enter a short sequence such as GGSGG, or a PROSITE-style "
        "pattern such as N-{P}-[ST]-{P}.")

    query = input("\nPattern: ").strip().upper()

    if not query:
        print("Pattern cannot be empty.")
        return

    protein = search_and_select_protein(proteins)

    if protein is None:
        return

    sequence = protein["sequence"]

    if "-" in query or "[" in query or "{" in query or "(" in query:
        regex = compile_prosite_pattern(query)

        if regex is None:
            print("\nThat pattern could not be understood.")
            return

        matches = find_pattern_occurrences(sequence, regex)
    else:
        matches = [
            {"start": position,
             "end": position + len(query) - 1,
             "sequence": query}
            for position in find_subsequence(sequence, query)]

    print_header("PATTERN SEARCH RESULT")

    print(f"Pattern      : {query}")
    print(f"UniProt ID   : {protein['id']}")
    print(f"Protein Name : {protein['name']}")
    print(f"Matches      : {len(matches):,}")

    if not matches:
        print("\nThis pattern does not occur in this protein.")
        return

    table = pd.DataFrame([
        {"No.": index,
         "Start": match["start"],
         "End": match["end"],
         "Matched Sequence": match["sequence"]}
        for index, match in enumerate(matches, start=1)])

    print_section("Match Positions")
    print(table.to_string(index=False))


def repeat_analysis(proteins: list[ProteinRecord]) -> None:
    protein = search_and_select_protein(proteins)

    if protein is None:
        return

    sequence = protein["sequence"]
    length = len(sequence)

    size = ask_number(
        "\nRepeat unit length (the size of the sliding sequence "
        "window, in residues)", 2, 10, 3)

    if length < size:
        print("\nSequence is too short for this repeat length.")
        return

    # Overlapping windows: a window starts at every position, so for
    # ABCDEFG with size 4 the windows are ABCD, BCDE, CDEF, DEFG -
    # each window overlaps the previous one by (size - 1) residues.
    counts = {}

    for start in range(length - size + 1):
        piece = sequence[start:start + size]
        counts[piece] = counts.get(piece, 0) + 1

    series = pd.Series(counts).sort_values(ascending=False)
    top = series.head(10)

    print_header("REPEATED SEQUENCE (K-MER) ANALYSIS")

    print(f"UniProt ID      : {protein['id']}")
    print(f"Protein Name    : {protein['name']}")
    print(f"Sequence Length : {length:,} amino acids")

    print(
        f"\nA {size}-residue unit means a window of {size} "
        f"consecutive amino acids read from the sequence. Windows "
        f"overlap and slide forward one residue at a time - for "
        f"example, with a unit length of 4 the sequence ABCDEFG "
        f"produces the overlapping windows ABCD, BCDE, CDEF, DEFG.")

    print(
        f"\nThis sequence contains {len(counts):,} overlapping "
        f"{size}-residue windows in total.")

    print(
        f"Of those, {len(series):,} are distinct - that is, "
        f"{len(series):,} unique {size}-residue sequences occur "
        f"one or more times.")

    print_section(
        f"Top 10 Most Frequent {size}-Residue Units "
        f"(out of {len(series):,} distinct units)")

    table = pd.DataFrame({
        "Rank": range(1, len(top) + 1),
        "Unit": top.index,
        "Occurrences": top.to_numpy(),
        "Share of Sequence (%)": (
            top.to_numpy() / length * 100).round(3)})

    print(table.to_string(index=False))

    if len(top) < 10:
        print(
            f"\nOnly {len(top):,} distinct {size}-residue unit(s) "
            f"exist for this protein, so fewer than 10 are listed.")
    else:
        print(
            f"\nThese are only the top 10 by frequency, not the "
            f"complete set of {len(series):,} distinct units.")

    print(
        "\nNote: this counts how often a short sequence pattern "
        "recurs in this one protein. A high count is a "
        "sequence-pattern statistic only - it does not by itself "
        "mean the unit is a functional motif, a disease-associated "
        "repeat, or evidence of a disordered region unless that is "
        "confirmed independently.")

    # Longest single-residue run: a separate analysis from the
    # k-mer counts above. This looks for the longest unbroken
    # stretch of one amino acid repeated in a row (e.g. LLLL is a
    # run of 4 consecutive leucines), not a repeated multi-residue
    # unit.
    longest_run = 1
    longest_code = sequence[0]
    longest_start = 0
    current_run = 1
    current_start = 0

    for index in range(1, length):
        if sequence[index] == sequence[index - 1]:
            current_run += 1
        else:
            current_start = index
            current_run = 1

        if current_run > longest_run:
            longest_run = current_run
            longest_code = sequence[index]
            longest_start = current_start

    longest_end = longest_start + longest_run - 1

    print_section("Longest Single-Residue Run")

    print(
        "This is a separate statistic from the k-mer analysis "
        "above: the longest unbroken stretch of one amino acid "
        "repeated consecutively, anywhere in the sequence.")

    if longest_run == 1:
        print(
            "\nNo amino acid repeats immediately after itself "
            "anywhere in this sequence, so the longest possible "
            "run is a single residue.")

    print(f"\nResidue        : {format_residue(longest_code)}")
    print(f"Run length     : {longest_run}")
    print(f"Sequence       : {longest_code * longest_run}")
    print(f"Start position : {longest_start + 1:,}")
    print(f"End position   : {longest_end + 1:,}")

    print(
        "\nNote: a long single-residue run is a composition "
        "observation about this sequence. It does not by itself "
        "indicate a functional or pathological feature unless "
        "confirmed independently.")


def sequence_features_menu(proteins: list[ProteinRecord], motifs: list[MotifRecord]) -> None:
    while True:
        choice = show_menu(
            "SEQUENCE FEATURES",
            ["Motif Detection (PROSITE)",
             "Custom Pattern Search",
             "Repeated Sequence Analysis",
             "Go Back",
             "Exit"])

        if choice == 1:
            motif_detection_menu(proteins, motifs)

        elif choice == 2:
            custom_pattern_search(proteins)

        elif choice == 3:
            repeat_analysis(proteins)

        elif choice == 4:
            return

        elif choice == 5:
            exit_program()


# ============================================================
# AMINO ACID ANALYSIS
# ============================================================

def display_amino_acid_reference() -> None:
    rows = []

    for code in AMINO_ACIDS:
        rows.append({
            "Code": code,
            "Amino Acid": AMINO_ACID_NAMES[code],
            "Class": AMINO_ACID_CLASSES[code],
            "Free Mass (Da)": AMINO_ACID_MASSES[code],
            "Residue Mass (Da)": round(
                AMINO_ACID_MASSES[code] - WATER_MASS, 2),
            "Hydropathy": HYDROPATHY_SCALE[code]})

    table = pd.DataFrame(rows)

    print_header("AMINO ACID REFERENCE TABLE")
    print(table.to_string(index=False))

    print(
        "\nHydropathy values are from the Kyte-Doolittle scale: "
        "positive means water-repelling, negative means water-liking.")


def dataset_composition(proteins: list[ProteinRecord]) -> dict[str, Any]:
    """Average amino acid composition across the whole dataset."""

    if "composition" in DATASET_CACHE:
        return DATASET_CACHE["composition"]

    print("\nCalculating dataset composition...")

    totals = {code: 0 for code in AMINO_ACIDS}
    total_residues = 0

    for index, protein in enumerate(proteins, start=1):
        sequence = protein["sequence"]

        for code in AMINO_ACIDS:
            totals[code] += sequence.count(code)

        total_residues += len(sequence)

        if index % 20000 == 0:
            print(f"  processed {index:,} proteins...")

    percentages = ({
        code: totals[code] / total_residues * 100
        for code in AMINO_ACIDS} if total_residues else {
        code: 0.0 for code in AMINO_ACIDS})

    result = {
        "totals": totals,
        "percentages": percentages,
        "total_residues": total_residues}

    DATASET_CACHE["composition"] = result

    return result


def display_dataset_composition(proteins: list[ProteinRecord]) -> None:
    result = dataset_composition(proteins)

    table = pd.DataFrame([
        {"Code": code,
         "Amino Acid": AMINO_ACID_NAMES[code],
         "Class": AMINO_ACID_CLASSES[code],
         "Total Count": result["totals"][code],
         "Percentage": round(result["percentages"][code], 3)}
        for code in AMINO_ACIDS])

    table = table.sort_values(
        by="Percentage", ascending=False).reset_index(drop=True)

    print_header(f"{dataset_label()} - AMINO ACID COMPOSITION".upper())

    print(f"Proteins        : {len(proteins):,}")
    print(f"Total residues  : {result['total_residues']:,}")

    print_section("Composition Across the Whole Dataset")
    print(table.to_string(index=False))

    if not ask_yes_no("\nDisplay dataset composition chart?"):
        return

    display_dataset_composition_chart(result)


def display_dataset_composition_chart(result: dict[str, Any]) -> None:
    """Pie chart of amino-acid composition across the whole dataset.

    Every amino acid's share here is a portion of one whole - 100%
    of all residues in the dataset - so a pie chart is the correct
    chart type, for the same reason a single protein's composition
    chart is a pie chart. Uses the exact same style (side-chain
    class colors, a class-color legend, and a full-name/percentage
    legend) so the two composition views look and read the same way.
    """

    percentages = np.array(
        [result["percentages"][code] for code in AMINO_ACIDS])

    order = np.argsort(percentages)[::-1]
    codes = [AMINO_ACIDS[index] for index in order]
    percentages = percentages[order]

    colors = [CLASS_COLORS[AMINO_ACID_CLASSES[code]] for code in codes]

    legend_labels = [
        f"{code} - {AMINO_ACID_NAMES[code]} ({percentages[i]:.2f}%)"
        for i, code in enumerate(codes)]

    figure, axes = plt.subplots(figsize=(11, 8))
    figure.patch.set_facecolor("white")

    wedges, _, autotexts = axes.pie(
        percentages,
        labels=codes,
        colors=colors,
        autopct=lambda pct: f"{pct:.1f}%" if pct >= 3 else "",
        pctdistance=0.78,
        startangle=90,
        wedgeprops={"edgecolor": "white", "linewidth": 1.0},
        textprops={"fontsize": 10, "fontweight": "bold"})

    for autotext in autotexts:
        autotext.set_color("white")
        autotext.set_fontsize(8)

    axes.set_title(
        f"Amino Acid Composition - {dataset_label()}",
        fontsize=14, fontweight="bold", pad=16)

    axes.axis("equal")

    class_handles = [
        plt.Rectangle((0, 0), 1, 1, color=CLASS_COLORS[name])
        for name in CLASS_ORDER]

    class_legend = axes.legend(
        class_handles, CLASS_ORDER,
        title="Side-chain class", fontsize=9, title_fontsize=9,
        loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)

    axes.add_artist(class_legend)

    axes.legend(
        wedges, legend_labels,
        title="Amino Acid (Full Name) - Share",
        loc="lower left", bbox_to_anchor=(1.02, 0.0),
        fontsize=8, title_fontsize=9, frameon=False)

    print("\nLegend:", "; ".join(legend_labels))

    plt.tight_layout()
    show_and_close(figure)


def compare_protein_to_dataset(proteins: list[ProteinRecord]) -> None:
    protein = search_and_select_protein(proteins)

    if protein is None:
        return

    sequence = protein["sequence"]
    dataset = dataset_composition(proteins)

    protein_percentages = amino_acid_composition(sequence)

    rows = []

    for code in AMINO_ACIDS:
        protein_value = protein_percentages[code]
        dataset_value = dataset["percentages"][code]

        rows.append({
            "Code": code,
            "Amino Acid": AMINO_ACID_NAMES[code],
            "Protein (%)": round(protein_value, 2),
            "Dataset (%)": round(dataset_value, 2),
            "Difference": round(protein_value - dataset_value, 2)})

    table = pd.DataFrame(rows)

    table = table.sort_values(
        by="Difference", ascending=False).reset_index(drop=True)

    print_header("COMPOSITION COMPARED WITH THE DATASET")

    print(f"UniProt ID   : {protein['id']}")
    print(f"Protein Name : {protein['name']}")

    print_section("Enrichment and Depletion (percentage points)")
    print(table.to_string(index=False))

    enriched = table.head(3)["Code"].tolist()
    depleted = table.tail(3)["Code"].tolist()

    print(
        "\nMost enriched :",
        format_residue_list(enriched))

    print(
        "Most depleted :",
        format_residue_list(list(reversed(depleted))))

    if not ask_yes_no("\nDisplay comparison chart?"):
        return

    ordered = table.set_index("Code").loc[list(AMINO_ACIDS)]

    figure, axes = prepare_plot()

    positions = np.arange(len(AMINO_ACIDS))
    width = 0.4

    axes.bar(
        positions - width / 2, ordered["Protein (%)"],
        width, label=protein["id"], color=PLOT_COLORS["accent"],
        edgecolor="white")

    axes.bar(
        positions + width / 2, ordered["Dataset (%)"],
        width, label="Dataset average", color=PLOT_COLORS["ink"],
        edgecolor="white")

    axes.set_xticks(positions)
    axes.set_xticklabels(list(AMINO_ACIDS), fontweight="bold")
    axes.legend(frameon=False, fontsize=10)

    finish_plot(
        axes,
        f"Composition of {protein['id']} vs Dataset Average",
        "Amino acid", "Share of sequence (%)")


def residue_position_analysis(proteins: list[ProteinRecord]) -> None:
    protein = search_and_select_protein(proteins)

    if protein is None:
        return

    code = input(
        "\nEnter a one-letter amino acid code: ").strip().upper()

    if code not in AMINO_ACID_NAMES:
        print("That is not a recognised amino acid code.")
        return

    sequence = protein["sequence"]

    positions = [
        index + 1
        for index, letter in enumerate(sequence)
        if letter == code]

    print_header("RESIDUE POSITION ANALYSIS")

    print(f"UniProt ID   : {protein['id']}")
    print(f"Protein Name : {protein['name']}")
    print(f"Residue      : {format_residue(code)}")
    print(f"Occurrences  : {len(positions):,}")

    if not positions:
        print("\nThis residue does not occur in this protein.")
        return

    print(
        f"Share        : "
        f"{len(positions) / len(sequence) * 100:.2f}%")

    print(
        "\nFirst positions:",
        ", ".join(str(value) for value in positions[:25]))

    if len(positions) > 1:
        gaps = np.diff(np.array(positions))

        print(
            f"\nAverage gap between occurrences: "
            f"{np.mean(gaps):.1f} residues")

        print(f"Smallest gap: {np.min(gaps)}")
        print(f"Largest gap : {np.max(gaps)}")

    if not ask_yes_no("\nDisplay a position map?"):
        return

    figure, axes = prepare_plot()

    axes.hlines(1, 0, len(sequence), color=PLOT_COLORS["track"], linewidth=12)

    axes.vlines(
        positions, 0.72, 1.28,
        color=CLASS_COLORS.get(
            AMINO_ACID_CLASSES.get(code, "Polar"), "#2d3436"),
        linewidth=1.2)

    axes.set_ylim(0.5, 1.5)
    axes.set_yticks([])
    axes.set_xlim(0, len(sequence))

    finish_plot(
        axes,
        f"Positions of {format_residue(code)} in {protein['id']}",
        "Residue position", "")


def amino_acid_menu(proteins: list[ProteinRecord]) -> None:
    while True:
        choice = show_menu(
            "AMINO ACID ANALYSIS",
            ["Amino Acid Reference Table",
             "Dataset Amino Acid Composition",
             "Compare a Protein with the Dataset",
             "Positions of a Single Residue",
             "Go Back",
             "Exit"])

        if choice == 1:
            display_amino_acid_reference()

        elif choice == 2:
            display_dataset_composition(proteins)

        elif choice == 3:
            compare_protein_to_dataset(proteins)

        elif choice == 4:
            residue_position_analysis(proteins)

        elif choice == 5:
            return

        elif choice == 6:
            exit_program()


# ============================================================
# HYDROPATHY ANALYSIS
# ============================================================

def hydropathy_values(sequence: str) -> np.ndarray:
    """Kyte-Doolittle value for every residue, as a numpy array."""

    return np.array(
        [HYDROPATHY_SCALE.get(letter, 0.0) for letter in sequence],
        dtype=np.float64)


def gravy_score(sequence: str) -> float:
    """Average hydropathy across the whole sequence.

    Uses the same rule as hydropathy_values(): a residue with no
    published Kyte-Doolittle value (the non-standard codes X, B, Z,
    J) contributes 0.0 rather than being dropped from the average.
    Keeping both functions on this one rule means the GRAVY summary
    number and the sliding-window profile always agree, instead of
    quietly disagreeing whenever a sequence contains one of those
    codes.
    """

    if not sequence:
        return 0.0

    total = sum(HYDROPATHY_SCALE.get(code, 0.0) for code in sequence)

    return total / len(sequence)


def sliding_hydropathy(sequence: str, window: int) -> tuple:
    values = hydropathy_values(sequence)

    if values.size < window:
        return None, None

    weights = np.ones(window) / window
    averages = np.convolve(values, weights, mode="valid")

    # Position of the middle residue of each window.
    centres = np.arange(len(averages)) + window // 2 + 1

    return centres, averages


def find_hydrophobic_segments(sequence: str, window: int, threshold: float) -> list[tuple[int, int]]:
    centres, averages = sliding_hydropathy(sequence, window)

    if averages is None:
        return []

    above = averages > threshold

    segments = []
    start = None

    for index, flag in enumerate(above):
        if flag and start is None:
            start = index

        elif not flag and start is not None:
            segments.append((int(centres[start]), int(centres[index - 1])))
            start = None

    if start is not None:
        segments.append((int(centres[start]), int(centres[-1])))

    return segments


def display_hydropathy_analysis(protein: ProteinRecord) -> None:
    sequence = protein["sequence"]

    window = ask_number(
        "\nSliding window size", 3, 51, DEFAULT_HYDROPATHY_WINDOW)

    centres, averages = sliding_hydropathy(sequence, window)

    print_header("HYDROPATHY ANALYSIS")

    print(f"UniProt ID   : {protein['id']}")
    print(f"Protein Name : {protein['name']}")
    print(f"Length       : {len(sequence):,} amino acids")
    print(f"Window size  : {window}")

    gravy = gravy_score(sequence)

    print(f"\nGRAVY score  : {gravy:.3f}")

    if gravy > 0:
        print(
            "A positive GRAVY score suggests an overall "
            "hydrophobic, water-repelling protein.")
    else:
        print(
            "A negative GRAVY score suggests an overall "
            "hydrophilic, water-soluble protein.")

    if averages is None:
        print(
            "\nThe sequence is shorter than the window, so no "
            "sliding profile can be drawn.")

        return

    print(f"Most hydrophobic window average : {np.max(averages):.2f}")
    print(f"Most hydrophilic window average : {np.min(averages):.2f}")

    segments = find_hydrophobic_segments(
        sequence, TRANSMEMBRANE_WINDOW, TRANSMEMBRANE_THRESHOLD)

    print_section("Possible Membrane-Spanning Segments")

    print(
        f"Windows of {TRANSMEMBRANE_WINDOW} residues with an average "
        f"above {TRANSMEMBRANE_THRESHOLD}:")

    if not segments:
        print("None found. This protein is probably not membrane-bound.")
    else:
        table = pd.DataFrame([
            {"No.": index,
             "Start": start,
             "End": end,
             "Length": end - start + 1}
            for index, (start, end) in enumerate(segments, start=1)])

        print(table.to_string(index=False))

        print(
            "\nThis is a simple indication only, not a prediction "
            "from a trained model.")

    if not ask_yes_no("\nDisplay the hydropathy plot?"):
        return

    figure, axes = prepare_plot()

    axes.plot(centres, averages, color=PLOT_COLORS["ink"], linewidth=1.4)

    axes.fill_between(
        centres, averages, 0,
        where=(averages > 0), color=PLOT_COLORS["accent"], alpha=0.45,
        label="Hydrophobic")

    axes.fill_between(
        centres, averages, 0,
        where=(averages <= 0), color=PLOT_COLORS["blue"], alpha=0.45,
        label="Hydrophilic")

    axes.axhline(0, color=PLOT_COLORS["ink"], linewidth=1.0)

    axes.axhline(
        TRANSMEMBRANE_THRESHOLD, color=PLOT_COLORS["red"],
        linestyle="--", linewidth=1.2,
        label=f"Threshold {TRANSMEMBRANE_THRESHOLD}")

    axes.legend(frameon=False, fontsize=9, loc="upper right")

    finish_plot(
        axes,
        f"Kyte-Doolittle Hydropathy - {protein['id']} "
        f"(window {window})",
        "Residue position", "Average hydropathy")


def dataset_gravy(proteins: list[ProteinRecord]) -> np.ndarray:
    if "gravy" in DATASET_CACHE:
        return DATASET_CACHE["gravy"]

    print("\nCalculating GRAVY scores for the whole dataset...")

    scores = []

    for index, protein in enumerate(proteins, start=1):
        scores.append(gravy_score(protein["sequence"]))

        if index % 20000 == 0:
            print(f"  processed {index:,} proteins...")

    array = np.array(scores, dtype=np.float64)

    DATASET_CACHE["gravy"] = array

    return array


def hydropathy_menu(proteins: list[ProteinRecord]) -> None:
    while True:
        choice = show_menu(
            "HYDROPATHY ANALYSIS",
            ["Hydropathy Profile of a Protein",
             "GRAVY Scores Across the Dataset",
             "Go Back",
             "Exit"])

        if choice == 1:
            protein = search_and_select_protein(proteins)

            if protein is not None:
                display_hydropathy_analysis(protein)

        elif choice == 2:
            scores = dataset_gravy(proteins)

            print_statistics(
                "GRAVY Statistics", scores, "")

            hydrophobic = int(np.sum(scores > 0))

            print(
                f"\nProteins with a positive GRAVY score: "
                f"{hydrophobic:,} "
                f"({hydrophobic / len(scores) * 100:.1f}%)")

            if ask_yes_no("\nDisplay the GRAVY distribution?"):
                display_histogram(
                    scores,
                    "Distribution of GRAVY Scores",
                    "GRAVY score", "#16a085")

        elif choice == 3:
            return

        elif choice == 4:
            exit_program()


# ============================================================
# PROTEIN FUNCTION ANALYSIS
# ============================================================

def compile_all_motifs(motifs: list[MotifRecord]) -> list[dict[str, Any]]:
    """Compile every PROSITE pattern once and keep the result."""

    cache_key = (
        "compiled_motifs",
        tuple((motif["accession"], motif["pattern"]) for motif in motifs))

    if cache_key in DATASET_CACHE:
        return DATASET_CACHE[cache_key]

    print("\nPreparing PROSITE patterns...")

    compiled = []
    unsupported = 0

    for motif in motifs:
        regex = compile_prosite_pattern(motif["pattern"])

        if regex is None:
            unsupported += 1
            continue

        compiled.append({"motif": motif, "regex": regex})

    print(
        f"Patterns ready: {len(compiled):,} "
        f"(unsupported syntax: {unsupported:,})")

    DATASET_CACHE[cache_key] = compiled

    return compiled


def motif_profile(protein: ProteinRecord, motifs: list[MotifRecord]) -> list[dict[str, Any]]:
    compiled = compile_all_motifs(motifs)
    sequence = protein["sequence"]

    found = []

    for item in compiled:
        matches = find_pattern_occurrences(sequence, item["regex"])

        if matches:
            found.append({
                "motif": item["motif"],
                "count": len(matches),
                "first": matches[0]["start"]})

    return found


def display_function_analysis(protein: ProteinRecord, motifs: list[MotifRecord]) -> None:
    print_header("PROTEIN FUNCTION ANALYSIS")

    print(f"UniProt ID   : {protein['id']}")
    print(f"Protein Name : {protein['name']}")
    print(f"Gene         : {protein['gene'] or 'Not listed'}")
    print(f"Organism     : {protein['organism'] or 'Not listed'}")
    print(f"Length       : {len(protein['sequence']):,} amino acids")

    evidence = extract_header_field(protein["header"], "PE")

    if evidence:
        print(f"Evidence     : {evidence}")

    print("\nSearching every PROSITE pattern against this protein...")

    found = motif_profile(protein, motifs)

    if not found:
        print("\nNo PROSITE patterns were found in this protein.")
        return

    found.sort(key=lambda item: item["count"], reverse=True)

    print(f"\nPatterns found: {len(found):,}")

    table = pd.DataFrame([
        {"No.": index,
         "PROSITE ID": shorten(item["motif"]["id"], 24),
         "Accession": item["motif"]["accession"],
         "Description": shorten(item["motif"]["description"], 44),
         "Hits": item["count"],
         "First": item["first"]}
        for index, item in enumerate(
            found[:MAX_DISPLAY_RESULTS], start=1)])

    print_section("Functional Patterns Detected")
    print(table.to_string(index=False))

    if len(found) > MAX_DISPLAY_RESULTS:
        print(
            f"\n{len(found) - MAX_DISPLAY_RESULTS:,} further "
            f"pattern(s) not shown.")

    print(
        "\nNote: short PROSITE patterns occur very often by chance. "
        "Longer, more specific patterns carry far more meaning.")


def function_menu(proteins: list[ProteinRecord], motifs: list[MotifRecord]) -> None:
    while True:
        choice = show_menu(
            "PROTEIN FUNCTION ANALYSIS",
            ["Motif Profile of a Protein",
             "Go Back",
             "Exit"])

        if choice == 1:
            protein = search_and_select_protein(proteins)

            if protein is not None:
                display_function_analysis(protein, motifs)

        elif choice == 2:
            return

        elif choice == 3:
            exit_program()


# ============================================================
# SEQUENCE COMPARISON
# ============================================================

def shared_kmers(first: str, second: str, size: int) -> tuple[int, float]:
    first_set = {
        first[index:index + size]
        for index in range(len(first) - size + 1)}

    second_set = {
        second[index:index + size]
        for index in range(len(second) - size + 1)}

    shared = first_set & second_set
    union = first_set | second_set

    if not union:
        return 0, 0.0

    return len(shared), len(shared) / len(union) * 100


def composition_distance(first: str, second: str) -> float:
    """How different two compositions are, in percentage points."""

    first_values = np.array(
        [amino_acid_composition(first)[code] for code in AMINO_ACIDS])

    second_values = np.array(
        [amino_acid_composition(second)[code] for code in AMINO_ACIDS])

    return float(np.sum(np.abs(first_values - second_values)) / 2)


def display_dot_plot(first_protein: ProteinRecord, second_protein: ProteinRecord, window: int) -> None:
    first = first_protein["sequence"][:DOT_PLOT_MAX_LENGTH]
    second = second_protein["sequence"][:DOT_PLOT_MAX_LENGTH]

    first_array = np.array(list(first))
    second_array = np.array(list(second))

    # True wherever two residues are identical.
    matrix = first_array[:, None] == second_array[None, :]

    if window > 1:
        # Keep only diagonal runs of the requested length.
        smoothed = np.ones_like(matrix, dtype=bool)

        for shift in range(window):
            smoothed[
                : matrix.shape[0] - window + 1,
                : matrix.shape[1] - window + 1] &= matrix[
                    shift: matrix.shape[0] - window + 1 + shift,
                    shift: matrix.shape[1] - window + 1 + shift]

        smoothed[matrix.shape[0] - window + 1:, :] = False
        smoothed[:, matrix.shape[1] - window + 1:] = False

        matrix = smoothed

    positions = np.argwhere(matrix)

    figure, axes = prepare_plot()

    axes.scatter(
        positions[:, 1] + 1, positions[:, 0] + 1,
        s=1.2, color=PLOT_COLORS["ink"], alpha=0.7)

    axes.set_xlim(0, len(second))
    axes.set_ylim(0, len(first))
    axes.invert_yaxis()

    finish_plot(
        axes,
        f"Dot Plot (window {window})",
        f"{second_protein['id']} position",
        f"{first_protein['id']} position")


def display_sequence_comparison(first_protein: ProteinRecord, second_protein: ProteinRecord, motifs: list[MotifRecord]) -> None:
    first = first_protein["sequence"]
    second = second_protein["sequence"]

    print_header("SEQUENCE COMPARISON")

    rows = [
        ("UniProt ID", first_protein["id"], second_protein["id"]),
        ("Gene",
         first_protein["gene"] or "-", second_protein["gene"] or "-"),
        ("Length", f"{len(first):,}", f"{len(second):,}"),
        ("Molecular weight (Da)",
         f"{protein_molecular_weight(first):,.0f}",
         f"{protein_molecular_weight(second):,.0f}"),
        ("Isoelectric point",
         f"{isoelectric_point(first):.2f}",
         f"{isoelectric_point(second):.2f}"),
        ("Net charge at pH 7",
         f"{net_charge_at_ph(first, 7.0):+.2f}",
         f"{net_charge_at_ph(second, 7.0):+.2f}"),
        ("GRAVY",
         f"{gravy_score(first):.3f}",
         f"{gravy_score(second):.3f}"),
        ("Aromaticity (%)",
         f"{aromaticity(first):.2f}",
         f"{aromaticity(second):.2f}")]

    table = pd.DataFrame(
        rows, columns=["Property", "Protein A", "Protein B"])

    print_section("Side-by-Side Properties")
    print(table.to_string(index=False))

    shared_count, similarity = shared_kmers(first, second, KMER_SIZE)

    print_section("Sequence Similarity")

    print(
        f"Shared {KMER_SIZE}-residue words : {shared_count:,}")

    print(
        f"Word-based similarity      : {similarity:.2f}%")

    print(
        f"Composition difference     : "
        f"{composition_distance(first, second):.2f} percentage points")

    if len(first) == len(second):
        identical = sum(
            1 for index in range(len(first))
            if first[index] == second[index])

        print(
            f"Position-by-position match : "
            f"{identical / len(first) * 100:.2f}%")

    if ask_yes_no("\nCompare amino acid composition on a chart?"):
        display_comparison_chart(first_protein, second_protein)

    if ask_yes_no("Display a dot plot?"):
        if (
            len(first) > DOT_PLOT_MAX_LENGTH
            or len(second) > DOT_PLOT_MAX_LENGTH
        ):
            print(
                f"\nBoth sequences will be trimmed to the first "
                f"{DOT_PLOT_MAX_LENGTH:,} residues to keep the plot "
                f"readable.")

        window = ask_number("Dot plot window", 1, 15, 4)

        display_dot_plot(first_protein, second_protein, window)

    if ask_yes_no("Compare PROSITE motif profiles? (slow)"):
        display_motif_comparison(first_protein, second_protein, motifs)


def display_comparison_chart(first_protein: ProteinRecord, second_protein: ProteinRecord) -> None:
    first = amino_acid_composition(first_protein["sequence"])
    second = amino_acid_composition(second_protein["sequence"])

    figure, axes = prepare_plot()

    positions = np.arange(len(AMINO_ACIDS))
    width = 0.4

    axes.bar(
        positions - width / 2,
        [first[code] for code in AMINO_ACIDS],
        width, label=first_protein["id"],
        color=PLOT_COLORS["accent"], edgecolor="white")

    axes.bar(
        positions + width / 2,
        [second[code] for code in AMINO_ACIDS],
        width, label=second_protein["id"],
        color=PLOT_COLORS["blue"], edgecolor="white")

    axes.set_xticks(positions)
    axes.set_xticklabels(list(AMINO_ACIDS), fontweight="bold")
    axes.legend(frameon=False, fontsize=10)

    finish_plot(
        axes, "Amino Acid Composition Comparison",
        "Amino acid", "Share of sequence (%)")


def display_motif_comparison(first_protein: ProteinRecord, second_protein: ProteinRecord, motifs: list[MotifRecord]) -> None:
    print("\nBuilding motif profiles for both proteins...")

    first_profile = motif_profile(first_protein, motifs)
    second_profile = motif_profile(second_protein, motifs)

    first_ids = {item["motif"]["id"] for item in first_profile}
    second_ids = {item["motif"]["id"] for item in second_profile}

    shared = sorted(first_ids & second_ids)
    only_first = sorted(first_ids - second_ids)
    only_second = sorted(second_ids - first_ids)

    print_section("Motif Profile Comparison")

    print(f"Patterns in {first_protein['id']:<12}: {len(first_ids):,}")
    print(f"Patterns in {second_protein['id']:<12}: {len(second_ids):,}")
    print(f"{'Shared patterns':<24}: {len(shared):,}")

    if first_ids or second_ids:
        overlap = len(shared) / len(first_ids | second_ids) * 100
        print(f"{'Profile overlap':<24}: {overlap:.1f}%")

    print("\nShared (first 10):")
    print(", ".join(shared[:10]) or "None")

    print(f"\nOnly in {first_protein['id']} (first 10):")
    print(", ".join(only_first[:10]) or "None")

    print(f"\nOnly in {second_protein['id']} (first 10):")
    print(", ".join(only_second[:10]) or "None")


def comparison_menu(proteins: list[ProteinRecord], motifs: list[MotifRecord]) -> None:
    while True:
        choice = show_menu(
            "SEQUENCE COMPARISON",
            ["Compare Two Proteins",
             "Go Back",
             "Exit"])

        if choice == 1:
            print("\nSelect the FIRST protein.")
            first = search_and_select_protein(proteins)

            if first is None:
                continue

            print("\nSelect the SECOND protein.")
            second = search_and_select_protein(proteins)

            if second is None:
                continue

            display_sequence_comparison(first, second, motifs)

        elif choice == 2:
            return

        elif choice == 3:
            exit_program()


# ============================================================
# DATASET STATISTICS
# ============================================================

def dataset_lengths(proteins: list[ProteinRecord]) -> np.ndarray:
    if "lengths" in DATASET_CACHE:
        return DATASET_CACHE["lengths"]

    lengths = np.array(
        [len(protein["sequence"]) for protein in proteins],
        dtype=np.float64)

    DATASET_CACHE["lengths"] = lengths

    return lengths


def display_dataset_overview(proteins: list[ProteinRecord]) -> None:
    lengths = dataset_lengths(proteins)

    print_header("DATASET OVERVIEW")

    print(f"Protein records : {len(proteins):,}")
    print(f"Total residues  : {int(np.sum(lengths)):,}")

    genes = sum(1 for protein in proteins if protein["gene"])

    print(f"Records with a gene name: {genes:,}")

    print_statistics(
        "Sequence Length Statistics", lengths, "residues")

    order = np.argsort(lengths)[::-1]

    longest = pd.DataFrame([
        {"UniProt ID": proteins[index]["id"],
         "Gene": shorten(proteins[index]["gene"], 10),
         "Protein Name": shorten(proteins[index]["name"], 42),
         "Length": len(proteins[index]["sequence"])}
        for index in order[:10]])

    shortest = pd.DataFrame([
        {"UniProt ID": proteins[index]["id"],
         "Gene": shorten(proteins[index]["gene"], 10),
         "Protein Name": shorten(proteins[index]["name"], 42),
         "Length": len(proteins[index]["sequence"])}
        for index in order[::-1][:10]])

    print_section("Ten Longest Proteins")
    print(longest.to_string(index=False))

    print_section("Ten Shortest Proteins")
    print(shortest.to_string(index=False))


def display_length_distribution(proteins: list[ProteinRecord]) -> None:
    lengths = dataset_lengths(proteins)

    display_histogram(
        lengths,
        "Distribution of Protein Lengths",
        "Sequence length (residues, log scale)",
        "#2980b9", log_x=True)


def display_length_vs_weight(proteins: list[ProteinRecord]) -> None:
    lengths = dataset_lengths(proteins)
    weights = dataset_weights(proteins)

    size = min(len(lengths), len(weights))

    figure, axes = prepare_plot()

    axes.scatter(
        lengths[:size], weights[:size],
        s=4, alpha=0.25, color=PLOT_COLORS["teal"], edgecolors="none")

    correlation = float(
        np.corrcoef(lengths[:size], weights[:size])[0, 1])

    print(
        f"\nCorrelation between length and molecular weight: "
        f"{correlation:.4f}")

    print(
        "A value close to 1 means the two measures rise together "
        "almost perfectly, which is expected.")

    finish_plot(
        axes,
        f"Length vs Molecular Weight (r = {correlation:.4f})",
        "Sequence length (residues)", "Molecular weight (Da)")


def display_class_distribution(proteins: list[ProteinRecord]) -> None:
    result = dataset_composition(proteins)

    totals = {name: 0 for name in CLASS_ORDER}

    for code in AMINO_ACIDS:
        totals[AMINO_ACID_CLASSES[code]] += result["totals"][code]

    values = np.array(
        [totals[name] for name in CLASS_ORDER], dtype=np.float64)

    percentages = values / np.sum(values) * 100

    table = pd.DataFrame({
        "Class": CLASS_ORDER,
        "Total Residues": values.astype(np.int64),
        "Percentage": percentages.round(2)})

    print_header("SIDE-CHAIN CLASSES ACROSS THE DATASET")
    print(table.to_string(index=False))

    if not ask_yes_no("\nDisplay the class chart?"):
        return

    figure, axes = plt.subplots(figsize=(8, 8))
    figure.patch.set_facecolor("white")

    wedges, _, autotexts = axes.pie(
        percentages,
        labels=CLASS_ORDER,
        colors=[CLASS_COLORS[name] for name in CLASS_ORDER],
        explode=[0.03] * len(CLASS_ORDER),
        autopct="%1.1f%%",
        pctdistance=0.74,
        startangle=90,
        wedgeprops={"edgecolor": "white", "linewidth": 1.5},
        textprops={"fontsize": 11, "fontweight": "bold"})

    for autotext in autotexts:
        autotext.set_color("white")
        autotext.set_fontsize(11)

    axes.set_title(
        "Side-Chain Classes Across the Dataset",
        fontsize=14, fontweight="bold", pad=16)

    axes.axis("equal")

    plt.tight_layout()
    show_and_close(figure)


def dataset_statistics_menu(proteins: list[ProteinRecord]) -> None:
    while True:
        choice = show_menu(
            "DATASET STATISTICS",
            ["Dataset Overview and Length Statistics",
             "Length Distribution (chart)",
             "Molecular Weight Statistics",
             "Length vs Molecular Weight (chart)",
             "Side-Chain Classes Across the Dataset",
             "Go Back",
             "Exit"])

        if choice == 1:
            display_dataset_overview(proteins)

        elif choice == 2:
            display_length_distribution(proteins)

        elif choice == 3:
            print_statistics(
                "Molecular Weight Statistics",
                dataset_weights(proteins), "Da")

        elif choice == 4:
            display_length_vs_weight(proteins)

        elif choice == 5:
            display_class_distribution(proteins)

        elif choice == 6:
            return

        elif choice == 7:
            exit_program()


# ============================================================
# EXPORT
# ============================================================

def ensure_export_folder() -> str:
    if not os.path.isdir(EXPORT_FOLDER):
        os.makedirs(EXPORT_FOLDER)

    return EXPORT_FOLDER


def save_table(table: pd.DataFrame, filename: str) -> None:
    folder = ensure_export_folder()
    path = os.path.join(folder, filename)

    try:
        table.to_csv(path, index=False)

        print(f"\nSaved: {os.path.abspath(path)}")

    except OSError as error:
        print(f"\nCould not save the file: {error}")


def export_protein_report(protein: ProteinRecord) -> None:
    sequence = protein["sequence"]

    folder = ensure_export_folder()
    path = os.path.join(folder, f"{protein['id']}_report.txt")

    weight = protein_molecular_weight(sequence)
    coefficients = extinction_coefficient(sequence)

    lines = [
        "PROTEIN SEQUENCE ANALYZER - REPORT",
        "=" * 60,
        f"UniProt ID      : {protein['id']}",
        f"Protein Name    : {protein['name']}",
        f"Gene            : {protein['gene'] or 'Not listed'}",
        f"Organism        : {protein['organism'] or 'Not listed'}",
        f"Sequence Length : {len(sequence):,} amino acids",
        "",
        "CALCULATED PROPERTIES",
        "-" * 60,
        f"Molecular weight       : {weight:,.2f} Da "
        f"({weight / 1000:,.2f} kDa)",
        f"Isoelectric point (pI) : {isoelectric_point(sequence):.2f}",
        f"Net charge at pH 7.0   : "
        f"{net_charge_at_ph(sequence, 7.0):+.2f}",
        f"GRAVY score            : {gravy_score(sequence):.3f}",
        f"Aromaticity            : {aromaticity(sequence):.2f}%",
        f"Extinction (reduced)   : {coefficients['reduced']:,}",
        "",
        "AMINO ACID COMPOSITION",
        "-" * 60,
        composition_table(sequence).to_string(index=False),
        "",
        "SIDE-CHAIN CLASSES",
        "-" * 60,
        class_summary_table(sequence).to_string(index=False),
        "",
        "SEQUENCE",
        "-" * 60]

    for start in range(0, len(sequence), 60):
        lines.append(
            f"{start + 1:>7}  {sequence[start:start + 60]}")

    try:
        with open(path, "w", encoding="utf-8") as file:
            file.write("\n".join(lines))

        print(f"\nSaved: {os.path.abspath(path)}")

    except OSError as error:
        print(f"\nCould not save the report: {error}")


def export_protein_fasta(protein: ProteinRecord) -> None:
    folder = ensure_export_folder()
    path = os.path.join(folder, f"{protein['id']}.fasta")

    sequence = protein["sequence"]

    lines = [">" + protein["header"]]

    for start in range(0, len(sequence), 60):
        lines.append(sequence[start:start + 60])

    try:
        with open(path, "w", encoding="utf-8") as file:
            file.write("\n".join(lines) + "\n")

        print(f"\nSaved: {os.path.abspath(path)}")

    except OSError as error:
        print(f"\nCould not save the FASTA file: {error}")


def export_dataset_summary(proteins: list[ProteinRecord]) -> None:
    print("\nBuilding the dataset summary table...")

    rows = []

    for index, protein in enumerate(proteins, start=1):
        sequence = protein["sequence"]
        weight = protein_molecular_weight(sequence)

        rows.append({
            "UniProt ID": protein["id"],
            "Gene": protein["gene"],
            "Protein Name": protein["name"],
            "Length": len(sequence),
            "Molecular Weight (Da)": (
                round(weight, 2) if weight is not None else ""),
            "GRAVY": round(gravy_score(sequence), 4)})

        if index % 20000 == 0:
            print(f"  processed {index:,} proteins...")

    save_table(pd.DataFrame(rows), "dataset_summary.csv")


def export_menu(proteins: list[ProteinRecord]) -> None:
    while True:
        choice = show_menu(
            "EXPORT ANALYSIS",
            ["Export a Full Protein Report (text)",
             "Export a Protein's Composition (CSV)",
             "Export a Protein as FASTA",
             "Export a Dataset Summary (CSV, slow)",
             "Go Back",
             "Exit"])

        if choice == 1:
            protein = search_and_select_protein(proteins)

            if protein is not None:
                export_protein_report(protein)

        elif choice == 2:
            protein = search_and_select_protein(proteins)

            if protein is not None:
                save_table(
                    composition_table(protein["sequence"]),
                    f"{protein['id']}_composition.csv")

        elif choice == 3:
            protein = search_and_select_protein(proteins)

            if protein is not None:
                export_protein_fasta(protein)

        elif choice == 4:
            export_dataset_summary(proteins)

        elif choice == 5:
            return

        elif choice == 6:
            exit_program()


# ============================================================
# VIEW A PROTEIN SEQUENCE
# ============================================================

def display_protein_sequence(proteins: list[ProteinRecord]) -> None:
    protein = search_and_select_protein(proteins)

    if protein is None:
        return

    sequence = protein["sequence"]

    print_header("PROTEIN RECORD")

    print(f"UniProt ID   : {protein['id']}")
    print(f"Protein Name : {protein['name']}")
    print(f"Gene         : {protein['gene'] or 'Not listed'}")
    print(f"Organism     : {protein['organism'] or 'Not listed'}")
    print(f"Length       : {len(sequence):,} amino acids")

    print_section("Sequence (60 residues per line)")

    for start in range(0, len(sequence), 60):
        print(f"{start + 1:>7}  {sequence[start:start + 60]}")


# ============================================================
# MAIN PROGRAM
# ============================================================

def print_banner(source: str, protein: Optional[str] = None) -> None:
    source_label = "Remote UniProt" if source == "uniprot" else "Local files"
    detail = f" | Protein: {protein}" if protein else ""
    subtitle = Text("UniProt + PROSITE protein analysis", style="bright_cyan")
    metadata = Text(
        f"v{__version__}  |  Source: {source_label}{detail}",
        style="dim white",
    )
    body = Text.assemble(
        ("PROTEIN SEQUENCE ANALYZER\n", "bold bright_white"),
        subtitle,
        "\n",
        metadata,
    )
    console.print()
    console.print(
        Panel(
            body,
            title="[bold bright_green]● READY[/bold bright_green]",
            title_align="left",
            border_style="bright_blue",
            box=box.DOUBLE,
            padding=(1, 3),
        )
    )


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments.

    Both arguments are optional. When omitted, the hardcoded
    defaults under FILE LOCATIONS are used exactly as before, so
    running the program with no arguments behaves identically to
    earlier versions that had no CLI at all.
    """

    parser = argparse.ArgumentParser(
        prog="protein_analyzer.py",
        description=(
            "Protein Sequence Analyzer - explore a UniProt FASTA "
            "dataset together with the PROSITE pattern database."))

    parser.add_argument(
        "--source", choices=("uniprot", "local"), default="uniprot",
        help="Data source for the CLI (default: uniprot)")

    parser.add_argument(
        "--protein", metavar="ACCESSION",
        help="Optional UniProt accession to load first in remote mode")

    parser.add_argument(
        "--query", default="*",
        help=(
            "UniProt query for remote mode; '*' means all proteins "
            "(default: *)"))

    parser.add_argument(
        "--page-size", type=int, default=REMOTE_PAGE_SIZE,
        help="Number of remote UniProt records to load per page (1-500)")

    parser.add_argument(
        "--fasta", metavar="PATH", default=FASTA_FILE,
        help="Path to the local UniProt FASTA file (only with --source local)")

    parser.add_argument(
        "--prosite", metavar="PATH", default=PROSITE_FILE,
        help="Path to the local PROSITE file (only with --source local)")

    parser.add_argument(
        "--version", action="version",
        version=f"Protein Sequence Analyzer {__version__}")

    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()

    if not 1 <= arguments.page_size <= 500:
        raise SystemExit("--page-size must be between 1 and 500.")

    print_banner(
        arguments.source,
        arguments.protein if arguments.source == "uniprot" else None)

    DATASET_CACHE.clear()

    global DATASET_SOURCE_PATH, REMOTE_QUERY, REMOTE_PAGE_SIZE
    global REMOTE_OFFSET, REMOTE_TOTAL

    if arguments.source == "uniprot":
        REMOTE_QUERY = arguments.query
        REMOTE_PAGE_SIZE = arguments.page_size
        if arguments.protein:
            print(f"\nFetching UniProt record {arguments.protein}...")
            proteins = [fetch_remote_protein(arguments.protein)]
            REMOTE_OFFSET = 1
            REMOTE_TOTAL = 1
        else:
            print(
                f"\nFetching the first {REMOTE_PAGE_SIZE:,} records "
                f"matching UniProt query '{REMOTE_QUERY}'...")
            proteins, total = fetch_remote_proteins(
                REMOTE_QUERY, 0, REMOTE_PAGE_SIZE)
            REMOTE_OFFSET = len(proteins)
            REMOTE_TOTAL = total
            print(f"Remote scope contains {total:,} matching records.")
        print("Downloading PROSITE definitions...")
        motifs = load_remote_prosite()
        DATASET_SOURCE_PATH = None
    else:
        fasta_path = resolve_file(arguments.fasta, "FASTA file")

        if fasta_path is None:
            print("\nNo FASTA file available. Use the default remote mode instead.")
            return

        DATASET_SOURCE_PATH = fasta_path
        print("\nLoading protein dataset...")
        proteins = load_fasta(fasta_path)

        if not proteins:
            print("\nNo protein records available. Program terminated.")
            return

        prosite_path = resolve_file(arguments.prosite, "PROSITE database file")

        if prosite_path is None:
            print("\nNo PROSITE file available. Use the default remote mode instead.")
            return

        print("\nLoading PROSITE database...")
        motifs = load_prosite(prosite_path)

    if not motifs:
        print("\nNo PROSITE motif records available. Program terminated.")
        return

    while True:
        menu_options = [
            "Basic Protein Analysis",
            "Physicochemical Properties",
            "Sequence Features",
            "Amino Acid Analysis",
            "Hydropathy Analysis",
            "Protein Function Analysis",
            "Sequence Comparison",
            "Dataset Statistics",
            "Export Analysis",
            "View a Protein Sequence"]

        if arguments.source == "uniprot":
            menu_options.append("Load More Remote Proteins")

        menu_options.append("Exit")
        choice = show_menu("MAIN MENU", menu_options)

        if choice == 1:
            protein = search_and_select_protein(proteins)

            if protein is not None:
                display_basic_analysis(protein)

        elif choice == 2:
            physicochemical_menu(proteins)

        elif choice == 3:
            sequence_features_menu(proteins, motifs)

        elif choice == 4:
            amino_acid_menu(proteins)

        elif choice == 5:
            hydropathy_menu(proteins)

        elif choice == 6:
            function_menu(proteins, motifs)

        elif choice == 7:
            comparison_menu(proteins, motifs)

        elif choice == 8:
            dataset_statistics_menu(proteins)

        elif choice == 9:
            export_menu(proteins)

        elif choice == 10:
            display_protein_sequence(proteins)

        elif arguments.source == "uniprot" and choice == 11:
            load_more_remote_proteins(proteins)

        elif choice == len(menu_options):
            print("\nProtein Sequence Analyzer closed.")
            break


if __name__ == "__main__":
    try:
        main()

    except KeyboardInterrupt:
        print("\n\nInterrupted. Protein Sequence Analyzer closed.")

    except EOFError:
        print("\n\nNo more input. Protein Sequence Analyzer closed.")

    except SystemExit:
        pass

    except Exception as error:  # noqa: BLE001 - last-resort safety net
        print(
            "\n\nAn unexpected error stopped the program: "
            f"{type(error).__name__}: {error}")

        print("This is not expected to happen. If it does, please "
            "report the exact steps that led to it.")
