# =============================================================================
# ENGINE: Bio Security Auditing Platform (BSAP)
# VERSION: 0.9.2
# AUTHOR:  BSAP contributors
# LICENSE: MIT
# STATUS:  OPEN SOURCE / RESEARCH TOOL
# =============================================================================

**Release 0.9.2** — forensic genomic **triage** engine for DNA, RNA, and protein
sequences. Scans user-supplied databases for lab-associated markers (cloning
scars, vector backbones, selection cassettes, reporters, promoters, and related
signals), extracts sequence flanks, cross-links shared flanks across Clinical /
Industrial / Reference sets, and produces colour-coded heatmaps, a multi-file
waterfall view, evidence logs, optional audio sonification, and an Excel report.

**Not a medical diagnostic.** It does not prove disease etiology, legal origin, or
that any experimental physical model is correct. Hypothesis-labelled fields in
logs (for example charge-density RF heuristics) are configuration aids, not
established science. Validate critical findings with independent laboratory methods.

All core modules ship as editable Python source under the MIT License so researchers
can inspect, fork, and adapt the pipeline (including `Targets.json` panels).

---

## Quick start (local)

```bash
chmod +x Installer.sh
./Installer.sh
streamlit run APP.py
```

Open **http://localhost:8501**, pick a **target panel** if desired, then **Run Forensic Audit**.

| Entry | Command |
|-------|---------|
| **Dashboard (primary)** | `streamlit run APP.py` |
| **Headless audit** | `python Alpha4.py` |
| **Standalone website** | After an audit: `python -m http.server 8765` → http://localhost:8765/index.html |

The website (`index.html`) is optional. It reads `Forensic_Results/Global_Evidence/portal_snapshot.json`
(and `Audit_Library/` when present). You can also open it from the Streamlit sidebar button.

---

## Installation

### 1. Clone

```bash
git clone https://github.com/YourUsername/Genomic_Audit_Engine.git
cd Genomic_Audit_Engine
```

### 2. Installer

`Installer.sh` checks Python 3.10+, creates a venv, installs `Requirements.txt`,
and runs `Setup_env.py` (directory tree; optional NCBI sync when `TEST_MODE=0`).

```bash
chmod +x Installer.sh
./Installer.sh
```

### 3. Docker (optional)

```bash
bash Run_Engine.sh
```

Dashboard: **http://localhost:8501**

---

## Running an audit

### Web UI

1. `streamlit run APP.py`
2. Sidebar: pick **Target panel** → **Save panel setting**
3. **Run Forensic Audit**
4. Review waterfall, per-file artifacts, cross-similarity, triage, Excel download

The engine module is always **`Alpha4.py`** (open source).

### CLI

```bash
python Alpha4.py
```

On first run, if input folders or FASTA are missing, Alpha4 can call `Setup_env`
automatically (`env_needs_bootstrap`).

---

## What the pipeline does

### Phase 1 — Scan & evidence

- Load the active panel from `Configuration.txt` (`TARGETS_PANEL`) via `Constants.resolve_targets_path`
- Master marker list remains in **`Targets.json`**; panel JSON files are optional subsets
- Dual-strand Aho-Corasick (FWD + reverse-complement patterns)
- Circular bridge window for plasmids
- Large files: mmap / chunked scan
- **Forensic** markers: evidence blocks, flanks, optional biophysics / circuit / codon bias
- **Composition** markers (short density motifs): counts in
  `Forensic_Results/Global_Evidence/composition/*.TARGET_COUNTS.json` (and optional hit indexes)
- Per-file heatmaps under `*_Heatmaps/FWD|REV/`

### Phase 2 — Flank cross-link

- Harvested LEFT/RIGHT flanks searched across the corpus (exact match, both strands)
- `CROSS_MATCH_LOG.txt`, similarity tracker, pair support (`PAIR_WINDOW_BP`), same-accession flags
- Flank-cross heatmaps when Phase 2 hits exist

### Optional Phase 2.5 — Flank-of-flank (`ENABLE_FLANK_OF_FLANK=1`)

- For each Phase-2 FOUND site, write outer sequence context to
  `Global_Evidence/FLANK_OF_FLANK_LOG.txt` (`FLANK_OF_FLANK_BP` bases each side)
- Lists nearby marker *names* in that window when the automaton is still loaded
- **Off by default.** Research navigation only — not a third proof of origin

### Phase 3 — Report & triage

- `Compiler.py` → `FINAL_GENOMIC_REPORT.xlsx` (forensic log, biophysics summary, cross-links, similarity)
- `WATERFALL_SUMMARY.png`, `TRIAGE_SUMMARY.json`, optional BED under `Global_Evidence/coordinates/`
- `portal_snapshot.json` for the HTML portal

---

## Target panels

| Panel | File | Use |
|-------|------|-----|
| full | `Targets.json` | All markers |
| lab_backbone | `targets_lab_backbone.json` | Chassis, tags, scars, promoters |
| viral_vector | `targets_viral_vector.json` | LTR, RRE, WPRE, CMV-style |
| composition | `targets_composition.json` | Short motifs / density only |
| custom | path in `TARGETS_CUSTOM_PATH` | Your own JSON |

Markers use `role`: `"forensic"` or `"composition"`. The UI can hide composition
columns in the waterfall until you enable **Show composition density**.

### Per-panel config snippets

Copy into `Configuration.txt` (or change only the panel line), then **Save panel**
in the UI or restart Streamlit, and re-run the audit.

**Full manifest (default research):**

```text
TARGETS_PANEL=full
TARGETS_CUSTOM_PATH=
COMPOSITION_FAST_ONLY=1
```

**Lab backbone only (AbR, ori, T7, tags, scars):**

```text
TARGETS_PANEL=lab_backbone
COMPOSITION_FAST_ONLY=1
ENABLE_HEATMAPS=1
```

**Viral / vector focus:**

```text
TARGETS_PANEL=viral_vector
COMPOSITION_FAST_ONLY=1
ENABLE_HEATMAPS=1
```

**Composition density study (short motifs; counts + optional positions):**

```text
TARGETS_PANEL=composition
COMPOSITION_FAST_ONLY=1
STORE_COMPOSITION_POSITIONS=1
ENABLE_HEATMAPS=0
```

**Custom JSON:**

```text
TARGETS_PANEL=custom
TARGETS_CUSTOM_PATH=path/to/my_markers.json
```

Edit the live manifest under **Manage Markers**, or edit JSON on disk and re-run.

---

## Configuration (`Configuration.txt`)

| Key | Role |
|-----|------|
| `TEST_MODE` | `1` = calibration FASTA via Setup_env; `0` = production |
| `TARGETS_PANEL` | `full` \| `lab_backbone` \| `viral_vector` \| `composition` \| `custom` |
| `FLANK_AUTO` / `MIN_FLANK_LEN` | Flank length for harvest and Phase 2 |
| `PAIR_WINDOW_BP` | LEFT+RIGHT pair proximity on destination |
| `MAX_CORES` | Worker processes |
| `ENABLE_HEATMAPS` | Per-file maps |
| `ENABLE_SONIFICATION` | WAV (default off for bulk speed) |
| `ENABLE_BIOPHYSICS` | Per-hit biophysical text |
| `ENABLE_CIRCUIT` | Local neighbor window |
| `ENABLE_CODON_BIAS` | RSCU-style signature text |
| `ENABLE_RF_HEURISTICS` | Explicitly experimental; prefer `0` for public defaults |
| `ENABLE_FLANK_OF_FLANK` | Optional outer context after Phase 2 (`0` default) |
| `FLANK_OF_FLANK_BP` | Outer window size when flank-of-flank is on |
| `COMPOSITION_FAST_ONLY` | Short motifs via fast count, not automaton |
| `STORE_COMPOSITION_POSITIONS` | Optional full position JSON for short motifs |
| `ENABLE_COORDINATE_EXPORT` | BED/TSV for long hits |
| `ARCHIVE_PREVIOUS_RUNS` | Archive last outputs under `_archive/` |

Suggested **ship defaults for bulk triage:**

```text
COMPOSITION_FAST_ONLY=1
STORE_COMPOSITION_POSITIONS=0
ENABLE_RF_HEURISTICS=0
ENABLE_SONIFICATION=0
ENABLE_FLANK_OF_FLANK=0
```

Tune `ENABLE_*` upward for deep single-construct reviews.

---

## Directory layout (summary)

```
Clinical_Samples/          # INPUT
Industrial_Sources/        # INPUT
Reference_Archive/         # INPUT
Forensic_Results/          # OUTPUT
  {Clinical|Industrial|Reference}_Evidence/{FWD,REV,FLANKS,WAV}/
  {Clinical|Industrial|Reference}_Heatmaps/{FWD,REV,FLANKS}/
  {Clinical|Industrial|Reference}_Spectrograms/{FWD,REV,FLANKS}/  # feature spectrograms
  Global_Heatmaps/         # waterfall, legend
  Global_Evidence/         # logs, composition/, hits/, coordinates/, triage, portal_*
  FINAL_GENOMIC_REPORT.xlsx
Alpha4.py                  # Engine
APP.py                     # Streamlit dashboard
Compiler.py
Constants.py / Configuration.txt
Targets.json + targets_*.json
Setup_env.py / DL.py
Installer.sh / Dockerfile / docker-compose.yml
index.html                 # Portal (snapshot-driven)
```

---

## Requirements

See `Requirements.txt`: streamlit, numpy, matplotlib, pandas, xlsxwriter,
pyahocorasick, tqdm, psutil, requests.

---


### Interpreting large genomes (T2T / multi‑GB RNA)

With `WATERFALL_INCLUDE_COMPOSITION=1`, portal `total_hits` and tag columns
**include short density motifs** (CG, NGG, 6‑cutters, …). A human T2T row can
show hundreds of millions of composition hits — that is density counting, not
“millions of lab constructs.”

For **decision triage**, trust:
- `TRIAGE_SUMMARY.json` lab‑marker rankings (composition excluded where designed)
- Phase‑2 flank edges / incorporation warnings
- Long forensic markers (ori, AbR, T7, loxP, reporters)

Set `WATERFALL_INCLUDE_COMPOSITION=0` if the portal matrix should show forensic
tags only.

Stress‑test defaults (large genomes, bounded runtime/disk):

```text
COMPOSITION_AC_MIN_LEN=6
COMPOSITION_FAST_ONLY=1
STORE_COMPOSITION_POSITIONS=0
ENABLE_BIOPHYSICS=0
ENABLE_CIRCUIT=0
ENABLE_CODON_BIAS=0
ENABLE_RF_HEURISTICS=0
ENABLE_HEATMAPS=1
```

CLI without Streamlit: `python Alpha4.py`. Run bundles land in
`Audit_Library/<stamp>/` for static hosting with `index.html` (sibling of Forensic_Results).

## Known limitations (2.6.0)

- **Standalone website** (`index.html`) is optional; open via Streamlit sidebar or `python -m http.server`. Streamlit does not embed the portal.
- **Composition vs forensic:** short motifs under `COMPOSITION_FAST_ONLY=1` appear as
  counts, not full evidence TXT / dense heatmaps — by design.
- **RF / experimental physics lines** are configuration heuristics, not validated biophysics.
- **Flank-of-flank** skips very large single FASTAs (>~80 MB) to avoid loading whole chromosomes into RAM for outer-context rows.
- NCBI accessions in `DL.py` are commented out by default; drop your own FASTAs into the input folders.

---

## Design notes for contributors

- **Open by design:** source is meant to be readable and changeable. Specialize by
  editing panels and config, not by locking modules.
- **Triage, not verdict of truth:** ranked files/pairs and shared flanks point at
  regions to inspect in IGV or the lab; they do not establish chain of custody alone.
- **Honest labels:** keep hypothesis / HEURISTIC wording on experimental physics
  fields; do not present them as validated biophysics.
- **Chain of custody helpers:** SHA-256 on sources, archive of previous runs, BED export.

---

## License

MIT — see `LICENSE`.

---

## Disclaimer

This software is a research and forensic **analysis aid**. It does not diagnose
disease, prove laboratory origin in a legal sense by itself, or certify that any
heuristic model is scientifically established. Use at your own risk; critical
claims require independent experimental validation.
