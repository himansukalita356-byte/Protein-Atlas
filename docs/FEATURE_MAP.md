# Backend calculation → Frontend analysis → Visualization

Every number on screen comes from **one** engine (`atlas_engine.py`). The same `analyze_sequence()`
result feeds the protein pages, Sequence Lab, Compare, Export and (through `DatasetAggregator`) the
dataset statistics; nothing is recalculated in the browser.

| # | Original module | Engine / service function | API endpoint | Front-end destination | Visualization |
|---|---|---|---|---|---|
| 1 | Basic Protein Analysis | `analyze_sequence` → `composition` (counts, %, most/least abundant, families) | `GET /proteins/{id}/analysis`, `POST /analyze` | Protein → **Overview**; Sequence Lab | Stat tiles, family table, interactive doughnut, interactive composition bar chart |
| 2 | Physicochemical Properties | `physicochemical` (MW, mass breakdown + reconciliation, pI, charge curve, extinction, aromaticity, aliphatic index, GRAVY) | same | Protein → **Physicochemical** | Mass-breakdown table + contribution chart, net-charge-vs-pH line chart with pI marker |
| 3 | Sequence Features | `search_pattern`, `PrositeDatabase.scan`, `repeat_analysis` | `POST /pattern`, `GET /proteins/{id}/prosite`, `POST /prosite/scan` | Protein → **Sequence Features**; **Pattern Lab** | Match table (Match/Start/End/Sequence), motif map, highlighted sequence, k-mer table |
| 4 | Amino Acid Analysis | `composition.rows`, `features.residue_positions`, dataset amino-acid frequencies | analysis endpoints, `GET /dataset/stats` | Protein → **Amino Acids** | Coloured bar chart (count / percent), 20-row table, residue position finder + track, enrichment vs dataset |
| 5 | Hydropathy Analysis | `hydropathy` (sliding window, GRAVY, membrane-like segments) | analysis endpoints (`hydropathy_window`) | Protein → **Hydropathy** | Zoom/pan line chart with shaded membrane bands and threshold; per-residue bars for short sequences |
| 6 | Protein Function Analysis | `PrositeDatabase.scan` (all scannable patterns) + UniProt annotation | prosite endpoints | Protein → **Function Scan** | Scan-status card, ranked hits table, motif map |
| 7 | Sequence Comparison | `compare_sequences`, `compare_proteins`, `dotplot_points` | `POST /compare` | **Compare** | Property table, relative bar chart, composition bars, similarity-vs-k line chart, dot plot |
| 8 | Dataset Statistics | `DatasetAggregator` via `DatasetStatsManager` (streams the full dataset once, caches) | `GET /dataset/stats`, `/dataset/universe` | **Dataset** | Length / MW / pI / GRAVY / aromaticity histograms, amino-acid bars, family doughnut, validation tiles, CSV |
| 9 | Export Analysis | exporters over the same analysis object | (client side) | Protein → **Export**; **Export** page | FASTA, composition CSV, mass CSV, hydropathy CSV, PROSITE CSV, JSON, text report |
| – | Molecular structure | `StructureService`: PDB → AlphaFold → SWISS-MODEL → ESMFold | `GET /proteins/{id}/structures`, `GET /structures/{source}/{id}/file` | Protein → 3D workspace (`?view=3d`) | 3Dmol.js: cartoon / stick / ball & stick / spacefill / surface, colour by N→C, chain, structure, confidence; residue inspector; deep-zoom detail lens |
| – | Amino-acid reference | `AMINO_ACID_REFERENCE` + `atlas_aminoacids` (offline 3D models) | `GET /reference/amino-acids[/{code}]` | **Reference** | Class-coloured grid, full property cards, interactive atom/bond 3D model |

## States (never a bare "undefined")

| State | What the user sees |
|---|---|
| Success | The value (`Count: 7`) |
| Genuinely unavailable | `Data unavailable` |
| External database unavailable | `Structure database unavailable` / `PROSITE database unavailable — no scan was performed` |
| Invalid input | `Check your input` + the specific reason (e.g. which k is allowed) |
| Calculation failure | `Analysis unavailable — calculation error` |
| Loading | Spinner with what is being fetched |
