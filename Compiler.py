import json
import os
import pandas as pd
import re
import io
from functools import lru_cache
import Constants  # <--- Crucial Sync Point
import BiophysicalProfile

# --- SYNCED FROM SINGLE SOURCE OF TRUTH ---
BASE_OUT = Constants.BASE_OUT

# Dynamically gather all evidence and map directories from Constants
_CAT = {"CLN": "Clinical", "IND": "Industrial", "REF": "Reference"}
EVIDENCE_DIRS = []
EVIDENCE_META = {}   # evidence dir -> (category, strand label)
for _key, _path in Constants.DIRS.items():
    _m = re.fullmatch(r"OUT_(CLN|IND|REF)_EVI_(FWD|REV)", _key)
    if _m:
        EVIDENCE_DIRS.append(_path)
        EVIDENCE_META[_path] = (_CAT[_m.group(1)], "FORWARD" if _m.group(2) == "FWD" else "REVERSE")
MAP_DIRS = [v for k, v in Constants.DIRS.items() if k.startswith("OUT_") and "_MAP" in k]


@lru_cache(maxsize=None)
def _read_text(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()


def _open_text(path):
    """File-like view of a cached read, so each evidence file is read from disk once per compile."""
    return io.StringIO(_read_text(path))


def _cap(text, limit=32000):
    """Excel cells hold at most 32,767 characters."""
    return text if len(text) <= limit else text[:limit - 20] + " ...[truncated]"

def extract_source_basename(evidence_filename):
    """
    Pulls the original source fasta filename out of an evidence filename of
    the exact form Alpha4.py writes: '{Category}_Evidence_{MOL}_{STRAND}_{base_name}.txt'.
    This is a straight split on the known delimiter Alpha4.py itself uses,
    which is far more reliable than trying to algorithmically strip
    prefixes/suffixes for a naming scheme clean_name() was never designed for.
    """
    if "_Evidence_" not in evidence_filename:
        return evidence_filename
    after = evidence_filename.split("_Evidence_", 1)[1]      # "{MOL}_{STRAND}_{base_name}.txt"
    parts = after.split("_", 2)
    if len(parts) < 3:
        return after
    base_with_ext = parts[2]
    if base_with_ext.endswith(".txt"):
        base_with_ext = base_with_ext[:-4]
    return base_with_ext  # e.g. "Clinical_DNA_Test.fasta"

_MAP_DIR_STRAND = {}
for _key, _path in Constants.DIRS.items():
    _m = re.fullmatch(r"OUT_(?:CLN|IND|REF)_MAP_(FWD|REV)", _key)
    if _m:
        _MAP_DIR_STRAND[_path] = "FORWARD" if _m.group(1) == "FWD" else "REVERSE"

_MAP_INDEX = None


def check_map_exists(source_basename, strand_label=None):
    """Path of '{base_name}_HEATMAP.png' (exact name Alpha4.py writes) for this strand, or NOT_FOUND."""
    global _MAP_INDEX
    if _MAP_INDEX is None:
        _MAP_INDEX = {}
        for root_dir in MAP_DIRS:
            strand = _MAP_DIR_STRAND.get(root_dir)
            for root, _, files in os.walk(root_dir):
                for f in files:
                    if f.endswith("_HEATMAP.png"):
                        _MAP_INDEX.setdefault((f, strand), os.path.join(root, f))
                        _MAP_INDEX.setdefault((f, None), os.path.join(root, f))
    name = f"{source_basename}_HEATMAP.png"
    return _MAP_INDEX.get((name, strand_label)) or _MAP_INDEX.get((name, None), "NOT_FOUND")

def get_latest_evidence_time():
    """Finds the timestamp of the most recently modified evidence file."""
    latest_time = 0
    for d in EVIDENCE_DIRS:
        if not os.path.exists(d): continue
        for root, _, files in os.walk(d):
            for f in files:
                if f.endswith(".txt"):
                    t = os.path.getmtime(os.path.join(root, f))
                    if t > latest_time: latest_time = t
    return latest_time

def load_targets_map():
    """Helper to map marker names to their specific Threat Tags."""
    try:
        with open(Constants.TARGETS_PATH, "r", encoding="utf-8") as f:
            data = json.load(f)
            # Flatten the nested dictionary to just name -> tag
            return {k: v.get('tag', 'UNKNOWN') for k, v in data.get('FORENSIC_MARKERS', {}).items()}
    except Exception as e:
        print(f"[!] Error loading Targets.json for tagging: {e}")
        return {}

def parse_evidence_logs():
    all_records = []
    tag_map = load_targets_map()
    
    for d in EVIDENCE_DIRS:
        if not os.path.exists(d): continue
        
        category, strand_label = EVIDENCE_META[d]

        for root, _, files in os.walk(d):
            for f in files:
                if "_Evidence_" in f and f.endswith(".txt"):
                    file_path = os.path.join(root, f)
                    subject_name = extract_source_basename(f)
                    map_link = check_map_exists(subject_name, strand_label)
                    
                    current_mol = "Unknown"
                    current_hash = "Unknown"
                    file_results = {}

                    # Molecule type is embedded in the filename itself (Alpha4.py's naming pattern)
                    mol_from_name = re.search(r"_Evidence_(DNA|RNA|PROT)_", f)
                    if mol_from_name: current_mol = mol_from_name.group(1)

                    with _open_text(file_path) as log:
                        content = log.read()
                        # Extract Metadata
                        hash_match = re.search(r"FILE HASH\s+:\s+(\w+)", content)
                        if hash_match: current_hash = hash_match.group(1)

                        # Alpha4.py separates hits with a dashed line, not an equals line
                        blocks = content.split("-" * 60)
                        
                        for block in blocks:
                            # Marker NAME lives on the "GENOMIC HIT DETECTED" line.
                            # (TARGET: holds the raw aligned sequence, not the marker name.)
                            target_match = re.search(r"GENOMIC HIT DETECTED:\s+([A-Za-z0-9_-]+)", block)
                            pos_match = re.search(r"LOCATION:\s+(\d+)", block)
                            
                            if target_match and pos_match:
                                marker_name = target_match.group(1).strip()
                                pos = pos_match.group(1).strip()
                                
                                group_key = (subject_name, strand_label)
                                if group_key not in file_results:
                                    file_results[group_key] = {
                                        "Subject File": subject_name,
                                        "SHA-256 Hash": current_hash,
                                        "Category": category,
                                        "Molecule": current_mol,
                                        "Markers": [],
                                        "Positions": [],
                                        "Threat Tag": "UNTAGGED",
                                        "Strand": strand_label,
                                        "Visual Map": map_link
                                    }
                                
                                file_results[group_key]["Markers"].append(marker_name)
                                file_results[group_key]["Positions"].append(pos)
                                
                                # Update tag if not UNTAGGED
                                tag = tag_map.get(marker_name, "UNTAGGED")
                                if tag != "UNTAGGED":
                                    file_results[group_key]["Threat Tag"] = tag

                    # Flatten results for this file
                    for data in file_results.values():
                        all_records.append({
                            "Subject File": data["Subject File"],
                            "SHA-256 Hash": data["SHA-256 Hash"],
                            "Category": data["Category"],
                            "Molecule": data["Molecule"],
                            "Marker Detected": ", ".join(sorted(set(data["Markers"]))),
                            "Threat Tag": data["Threat Tag"],
                            "Genomic Position": _cap(", ".join(data["Positions"])),
                            "Strand": data["Strand"],
                            "Visual Map": data["Visual Map"]
                        })
    
    return pd.DataFrame(all_records)

def parse_evidence_hits_detailed():
    """One row per individual hit (not aggregated by file/strand like the
    main log), used to compute real per-hit biophysical stats."""
    rows = []
    for d in EVIDENCE_DIRS:
        if not os.path.exists(d): continue
        category, strand_label = EVIDENCE_META[d]

        for root, _, files in os.walk(d):
            for f in files:
                if "_Evidence_" not in f or not f.endswith(".txt"): continue
                mol_from_name = re.search(r"_Evidence_(DNA|RNA|PROT)_", f)
                mol_type = mol_from_name.group(1) if mol_from_name else "Unknown"
                subject_name = extract_source_basename(f)
                with _open_text(os.path.join(root, f)) as log:
                    content = log.read()
                for block in content.split("-" * 60):
                    marker_m = re.search(r"GENOMIC HIT DETECTED:\s+([A-Za-z0-9_-]+)", block)
                    pos_m = re.search(r"LOCATION:\s+(\d+)", block)
                    match_m = re.search(r"^MATCH:\s+([A-Za-z*]+)", block, re.MULTILINE)
                    if marker_m and pos_m and match_m:
                        rows.append({
                            "Subject File": subject_name, "Category": category, "Molecule": mol_type,
                            "Marker": marker_m.group(1).strip(), "Position": pos_m.group(1).strip(),
                            "Strand": strand_label, "Match Sequence": match_m.group(1).strip()
                        })
    return pd.DataFrame(rows)

def build_biophysical_summary(detailed_df):
    """Computes BiophysicalProfile.build_profile() for every real hit and
    aggregates (mean/min/max/count) grouped by Category + Molecule."""
    if detailed_df.empty:
        return pd.DataFrame()

    records = []
    for _, row in detailed_df.iterrows():
        if row["Molecule"] not in ("DNA", "RNA", "PROT"):
            continue
        try:
            profile = BiophysicalProfile.build_profile(row["Match Sequence"], mol_type=row["Molecule"])
            records.append({
                "Category": row["Category"], "Molecule": row["Molecule"], "Marker": row["Marker"],
                "Mass_Da": profile["molecular_mass_da"],
                "Net_Charge": profile["net_charge_at_ph"] if profile["net_charge_at_ph"] is not None else float("nan"),
                "Coulomb_Force_Vacuum_N": profile["coulomb_force_vacuum_newtons"],
                "Coulomb_Force_Water_Screened_N": profile["coulomb_force_water_screened_newtons"],
                "Vibrational_THz": (profile.get("sample_vibrational_hz") or float("nan")) / 1e12,
                "Photon_Energy_eV": profile.get("sample_photon_energy_ev", float("nan")),
            })
        except Exception as e:
            print(f"[!] Biophysical summary skipped a hit ({row['Marker']}): {e}")

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)
    numeric_cols = ["Mass_Da", "Net_Charge", "Coulomb_Force_Vacuum_N",
                     "Coulomb_Force_Water_Screened_N", "Vibrational_THz", "Photon_Energy_eV"]
    summary = df.groupby(["Category", "Molecule"])[numeric_cols].agg(["mean", "min", "max", "count"])
    summary.columns = ["_".join(col).strip() for col in summary.columns.values]
    return summary.reset_index()

def parse_cross_similarity_tracker():
    """Load Phase 2 CROSS_SIMILARITY_TRACKER.json into a flat DataFrame for Excel."""
    path = os.path.join(Constants.DIRS["OUT_LOGS"], "CROSS_SIMILARITY_TRACKER.json")
    if not os.path.exists(path):
        return pd.DataFrame()
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        return pd.DataFrame()
    rows = []
    for e in data.get("edges") or []:
        rows.append({
            "File A": e.get("file_a", ""),
            "File B": e.get("file_b", ""),
            "Shared Flank Hits": e.get("shared_flank_hits", 0),
            "Pair Support Count": e.get("pair_support_count", 0),
            "Similarity Score": e.get("similarity_score", 0),
            "Same Accession": e.get("same_accession", False),
            "Accession": e.get("accession", ""),
            "Markers": "; ".join(e.get("markers") or []),
            "Top Label": e.get("top_label", ""),
            "Note": e.get("note", ""),
        })
    return pd.DataFrame(rows)

def parse_cross_match_log():
    """Parses Phase 2's Global_Evidence/CROSS_MATCH_LOG.txt into a DataFrame."""
    log_path = os.path.join(Constants.DIRS["OUT_LOGS"], "CROSS_MATCH_LOG.txt")
    if not os.path.exists(log_path):
        return pd.DataFrame()

    content = _read_text(log_path)
    rows = []
    for block in content.split("=" * 40):
        if "CROSS-LINK" not in block: continue
        def grab(pattern, default=""):
            m = re.search(pattern, block)
            return m.group(1).strip() if m else default

        rows.append({
            "Flank Sequence": grab(r"CROSS-LINK:\s+(\S+)"),
            "Marker": grab(r"MARKER:\s+(\S+)"),
            "Side": grab(r"SIDE:\s+(\S+)"),
            "Strand": grab(r"STRAND:\s+(\S+)"),
            "Identity %": grab(r"IDENTITY:\s+([\d.]+)"),
            "Found In": grab(r"FOUND:\s+(\S+)"),
            "Origin File": grab(r"ORIGIN:\s+(\S+)"),
            "Verdict": grab(r"VERDICT:\s+(.+?)\s+\(Risk"),
            "Risk Score": grab(r"Risk Score:\s+(\d+)"),
            "Nearby Context": grab(r"NEARBY CONTEXT:\s+(.+)"),
        })
    return pd.DataFrame(rows)


def build_start_here_sheet():
    """One-page orientation: config stamp + triage top files + multi-domain + warnings."""
    import json
    rows = []
    logs = Constants.DIRS.get("OUT_LOGS", "")
    stamp_path = os.path.join(logs, "RUN_CONFIG_STAMP.json")
    triage_path = os.path.join(logs, "TRIAGE_SUMMARY.json")
    warn_path = os.path.join(logs, "FLANK_INCORPORATION_WARNINGS.json")

    stamp = {}
    if os.path.isfile(stamp_path):
        try:
            with open(stamp_path, encoding="utf-8") as f:
                stamp = json.load(f)
        except Exception:
            pass
    elif os.path.isfile(triage_path):
        try:
            with open(triage_path, encoding="utf-8") as f:
                stamp = (json.load(f) or {}).get("config_stamp") or {}
        except Exception:
            pass

    rows.append({"section": "CONFIG_STAMP", "item": "config_sha256_16", "value": stamp.get("config_sha256_16", "")})
    rows.append({"section": "CONFIG_STAMP", "item": "created_utc", "value": stamp.get("created_utc", "")})
    rows.append({"section": "CONFIG_STAMP", "item": "targets_path", "value": stamp.get("targets_path", "")})
    for k, v in (stamp.get("gates") or {}).items():
        rows.append({"section": "CONFIG_GATE", "item": k, "value": str(v)})

    triage = {}
    if os.path.isfile(triage_path):
        try:
            with open(triage_path, encoding="utf-8") as f:
                triage = json.load(f)
        except Exception:
            triage = {}

    for i, r in enumerate((triage.get("top_files_by_lab_markers") or [])[:10], 1):
        rows.append({
            "section": "TOP_FILE",
            "item": f"rank_{i}",
            "value": r.get("file"),
            "triage_score": r.get("triage_score", r.get("lab_marker_hits")),
            "lab_marker_hits": r.get("lab_marker_hits"),
            "unique_long_tags": r.get("unique_long_tags"),
            "phase2_flank_support": r.get("phase2_flank_support"),
            "note": r.get("note", ""),
        })
    for i, r in enumerate((triage.get("multi_domain_files") or [])[:10], 1):
        rows.append({
            "section": "MULTI_DOMAIN",
            "item": f"multi_{i}",
            "value": r.get("file"),
            "triage_score": r.get("triage_score"),
            "note": f"domains={r.get('domains')}",
        })
    for i, r in enumerate((triage.get("top_file_pairs_by_flank_score") or [])[:10], 1):
        rows.append({
            "section": "TOP_PAIR",
            "item": f"pair_{i}",
            "value": f"{r.get('file_a')} <-> {r.get('file_b')}",
            "triage_score": r.get("similarity_score"),
            "note": f"shared_flank_hits={r.get('shared_flank_hits')}",
        })

    if os.path.isfile(warn_path):
        try:
            with open(warn_path, encoding="utf-8") as f:
                w = json.load(f)
            warnings = w if isinstance(w, list) else (w.get("warnings") or w.get("rows") or [])
            for i, r in enumerate(warnings[:10], 1):
                if not isinstance(r, dict):
                    continue
                rows.append({
                    "section": "INCORPORATION_WARNING",
                    "item": f"warn_{i}",
                    "value": r.get("message") or r.get("origin_file"),
                    "note": r.get("verdict", ""),
                    "triage_score": r.get("risk_score"),
                })
        except Exception:
            pass

    if not rows:
        rows = [{"section": "NOTE", "item": "empty", "value": "No triage/stamp yet — run a full audit."}]
    return pd.DataFrame(rows)


def run_compiler():
    """Main entry point: Generates the FINAL_GENOMIC_REPORT.xlsx."""
    global _MAP_INDEX
    _MAP_INDEX = None
    _read_text.cache_clear()
    output_path = os.path.join(BASE_OUT, "FINAL_GENOMIC_REPORT.xlsx")
    
    # --- IDEMPOTENCY CHECK (Don't waste CPU/RAM if nothing changed) ---
    if os.path.exists(output_path):
        report_time = os.path.getmtime(output_path)
        evidence_time = get_latest_evidence_time()
        _cross = os.path.join(Constants.DIRS["OUT_LOGS"], "CROSS_MATCH_LOG.txt")
        if os.path.exists(_cross):
            evidence_time = max(evidence_time, os.path.getmtime(_cross))
        _sim = os.path.join(Constants.DIRS["OUT_LOGS"], "CROSS_SIMILARITY_TRACKER.json")
        if os.path.exists(_sim):
            evidence_time = max(evidence_time, os.path.getmtime(_sim))
        _triage = os.path.join(Constants.DIRS["OUT_LOGS"], "TRIAGE_SUMMARY.json")
        if os.path.exists(_triage):
            evidence_time = max(evidence_time, os.path.getmtime(_triage))    
        if report_time > evidence_time:
            print(f"[-] Skip: {output_path} is already up to date.")
            return

    df = parse_evidence_logs()
    if df.empty:
        print("[!] No forensic data found in evidence folders.")
        return

    # Use 'xlsxwriter' to allow for large-scale reporting on low-ram systems
    with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
        try:
            start_df = build_start_here_sheet()
            start_df.to_excel(writer, sheet_name='Start_Here', index=False)
        except Exception as _sh_err:
            pd.DataFrame([{"Note": f"Start_Here failed: {_sh_err}"}]).to_excel(
                writer, sheet_name='Start_Here', index=False
            )
        df.to_excel(writer, sheet_name='Forensic_Log', index=False)
        workbook = writer.book
        worksheet = writer.sheets['Forensic_Log']
        
        # Forensic Formatting (Highlights missing evidence maps in Red)
        red_format = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})
        worksheet.conditional_format(1, 8, len(df), 8, {
            'type': 'cell',
            'criteria': 'equal to',
            'value': '"NOT_FOUND"',
            'format': red_format
        })

        # --- Sheet 2: Biophysical Summary (real, aggregated mass/charge/vibrational stats) ---
        detailed_hits = parse_evidence_hits_detailed()
        bio_summary = build_biophysical_summary(detailed_hits)
        if not bio_summary.empty:
            bio_summary.to_excel(writer, sheet_name='Biophysical_Summary', index=False)
        else:
            pd.DataFrame([{"Note": "No biophysical data available (ENABLE_BIOPHYSICS may be off, or no hits found)"}]).to_excel(
                writer, sheet_name='Biophysical_Summary', index=False)

        # --- Sheet 3: Cross-Links (Phase 2 flank cross-matching, folded into the same report) ---
        cross_df = parse_cross_match_log()
        if not cross_df.empty:
            cross_df.to_excel(writer, sheet_name='Cross_Links', index=False)
        else:
            pd.DataFrame([{"Note": "No cross-links found (Phase 2 found no shared flanks)"}]).to_excel(
                writer, sheet_name='Cross_Links', index=False)

        # --- Sheet 4: Cross-Similarity (file-pair flank reuse graph) ---
        sim_df = parse_cross_similarity_tracker()
        if not sim_df.empty:
            sim_df.to_excel(writer, sheet_name="Cross_Similarity", index=False)
        else:
            pd.DataFrame([{
                "Note": "No CROSS_SIMILARITY_TRACKER.json yet (run Phase 2 / full audit)"
            }]).to_excel(writer, sheet_name="Cross_Similarity", index=False)

    print(f"[+] Compiled {len(df)} forensic hits into Audit Matrix.")
    print(f"[+] Report: {output_path}")