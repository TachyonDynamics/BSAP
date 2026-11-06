# =============================================================================
# REPLACE entire Constants.py with this file content (or merge carefully).
# Anchors TARGETS_PATH resolution to active panel without deleting master list.
# =============================================================================
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
BASE_OUT = os.path.join(ROOT, "Forensic_Results")
# Durable snapshot library lives NEXT TO Forensic_Results (not inside it),
# so wiping Forensic_Results does not erase past audit packages.
AUDIT_LIBRARY = os.path.join(ROOT, "Audit_Library")
CONFIG_PATH = os.path.join(ROOT, "Configuration.txt")
# Master marker library (never delete markers from this file)
TARGETS_MASTER_PATH = os.path.join(ROOT, "Targets.json")
TARGETS_PATH = TARGETS_MASTER_PATH  # default; may be overridden by panel


def load_config(path=CONFIG_PATH):
    config = {}
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                k, v = line.split("=", 1)
                config[k.strip()] = v.strip()
    return config


CONFIG = load_config()


def cfg_int(key, default):
    try:
        return int(CONFIG.get(key, default))
    except (TypeError, ValueError):
        return default


def cfg_flag(key, default=False):
    return str(CONFIG.get(key, "1" if default else "0")).strip() == "1"


def resolve_targets_path(config=None):
    """
    TARGETS_PANEL values:
      full | Targets.json | (empty)  -> master Targets.json (all markers)
      lab_backbone                   -> targets_lab_backbone.json
      viral_vector                   -> targets_viral_vector.json
      composition                    -> targets_composition.json
      custom                         -> TARGETS_CUSTOM_PATH file
    Master Targets.json is never reduced; panels are optional search subsets.
    """
    cfg = config if config is not None else CONFIG
    panel = str(cfg.get("TARGETS_PANEL", "full")).strip().lower()
    custom = str(cfg.get("TARGETS_CUSTOM_PATH", "")).strip()
    aliases = {
        "full": TARGETS_MASTER_PATH,
        "all": TARGETS_MASTER_PATH,
        "targets.json": TARGETS_MASTER_PATH,
        "": TARGETS_MASTER_PATH,
        "lab_backbone": os.path.join(ROOT, "targets_lab_backbone.json"),
        "lab": os.path.join(ROOT, "targets_lab_backbone.json"),
        "viral_vector": os.path.join(ROOT, "targets_viral_vector.json"),
        "viral": os.path.join(ROOT, "targets_viral_vector.json"),
        "composition": os.path.join(ROOT, "targets_composition.json"),
        "comp": os.path.join(ROOT, "targets_composition.json"),
    }
    if panel == "custom" and custom:
        path = custom if os.path.isabs(custom) else os.path.join(ROOT, custom)
    else:
        path = aliases.get(panel, TARGETS_MASTER_PATH)
        if panel not in aliases and panel.endswith(".json"):
            path = panel if os.path.isabs(panel) else os.path.join(ROOT, panel)
    if not os.path.exists(path):
        path = TARGETS_MASTER_PATH
    return path


# Active panel path used by Alpha4.initialize_globals
TARGETS_PATH = resolve_targets_path(CONFIG)

_CATS = {"CLN": "Clinical", "IND": "Industrial", "REF": "Reference"}

DIRS = {
    "CLN": os.path.join(ROOT, "Clinical_Samples"),
    "IND": os.path.join(ROOT, "Industrial_Sources"),
    "REF": os.path.join(ROOT, "Reference_Archive"),
}
for _k, _label in _CATS.items():
    for _strand in ("FWD", "REV"):
        DIRS[f"OUT_{_k}_EVI_{_strand}"] = os.path.join(BASE_OUT, f"{_label}_Evidence", _strand)
        DIRS[f"OUT_{_k}_MAP_{_strand}"] = os.path.join(BASE_OUT, f"{_label}_Heatmaps", _strand)
    DIRS[f"OUT_{_k}_AUDIO"] = os.path.join(BASE_OUT, f"{_label}_Evidence", "WAV")
    DIRS[f"OUT_{_k}_FLANKS"] = os.path.join(BASE_OUT, f"{_label}_Evidence", "FLANKS")
    DIRS[f"OUT_{_k}_MAP_FLANKS"] = os.path.join(BASE_OUT, f"{_label}_Heatmaps", "FLANKS")
    # Sequence-feature spectrograms (composition display; NOT measured RF)
    # Layout mirrors heatmaps: FWD (sense), REV (RC composition), FLANKS (harvest strip)
    _spec_root = os.path.join(BASE_OUT, f"{_label}_Spectrograms")
    DIRS[f"OUT_{_k}_SPECTROGRAM"] = _spec_root  # root (legacy / gallery)
    DIRS[f"OUT_{_k}_SPECTROGRAM_FWD"] = os.path.join(_spec_root, "FWD")
    DIRS[f"OUT_{_k}_SPECTROGRAM_REV"] = os.path.join(_spec_root, "REV")
    DIRS[f"OUT_{_k}_SPECTROGRAM_FLANKS"] = os.path.join(_spec_root, "FLANKS")

DIRS.update({
    "OUT_MAPS": os.path.join(BASE_OUT, "Global_Heatmaps"),
    "OUT_LOGS": os.path.join(BASE_OUT, "Global_Evidence"),
    "CACHE": os.path.join(BASE_OUT, ".sequence_cache"),
    "AUDIT_LIBRARY": AUDIT_LIBRARY,
    "REPORT": os.path.join(BASE_OUT, "FINAL_GENOMIC_REPORT.xlsx"),
})
