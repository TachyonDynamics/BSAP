# BSAP local UI — primary entry point:
#   streamlit run APP.py
# Headless audit (no UI):
#   python Alpha4.py
# Standalone website (optional, after an audit):
#   python -m http.server 8765   then open http://localhost:8765/index.html
#   or use the sidebar button "Open website (index.html)"
import os
import sys

_ROOT = os.path.dirname(os.path.abspath(__file__))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
os.chdir(_ROOT)

# --- Windows: suppress Proactor ConnectionResetError when the browser closes ---
# Closing the tab only drops the WebSocket; Streamlit itself keeps running until
# Ctrl+C or the Stop button. Without this patch, WinError 10054 can hang the
# server on "Stopping..." and leave the console stuck.
if sys.platform == "win32":
    import asyncio
    try:
        from asyncio.proactor_events import _ProactorBasePipeTransport

        _orig_call_connection_lost = _ProactorBasePipeTransport._call_connection_lost

        def _silent_call_connection_lost(self, exc):
            try:
                _orig_call_connection_lost(self, exc)
            except (ConnectionResetError, ConnectionAbortedError, BrokenPipeError, OSError):
                pass

        _ProactorBasePipeTransport._call_connection_lost = _silent_call_connection_lost
    except Exception:
        pass

    def _ignore_conn_reset(loop, context):
        exc = context.get("exception")
        if isinstance(exc, (ConnectionResetError, ConnectionAbortedError, BrokenPipeError)):
            return
        msg = context.get("message", "")
        if "ConnectionResetError" in str(msg) or "forcibly closed" in str(msg):
            return
        loop.default_exception_handler(context)

    try:
        asyncio.get_event_loop_policy().get_event_loop().set_exception_handler(_ignore_conn_reset)
    except Exception:
        pass

    # Quiet asyncio logger noise from closed sockets
    import logging as _logging
    _logging.getLogger("asyncio").setLevel(_logging.CRITICAL)

import streamlit as st
import pandas as pd
import os
import json
import glob
import Constants
import logging

class _DropStreamlitNoise(logging.Filter):
    def filter(self, record):
        msg = record.getMessage()
        if "ScriptRunContext" in msg:
            return False
        if "Session state does not function" in msg:
            return False
        if "missing ScriptRunContext" in msg:
            return False
        if "ConnectionResetError" in msg or "forcibly closed" in msg:
            return False
        if "_call_connection_lost" in msg:
            return False
        return True

for _name in (
    "streamlit",
    "streamlit.runtime.scriptrunner_utils.script_run_context",
    "streamlit.runtime.state.session_state_proxy",
    "streamlit.runtime.scriptrunner_utils",
):
    _log = logging.getLogger(_name)
    _log.setLevel(logging.ERROR)
    _log.addFilter(_DropStreamlitNoise())

logging.getLogger().addFilter(_DropStreamlitNoise())  
    
# --- ENGINE AGNOSTIC IMPORT ---
import Alpha4 as Engine
import Compiler  

# --- UNIFIED FORENSIC HIERARCHY ---
DIRS = Constants.DIRS
BASE_OUT = Constants.BASE_OUT
AUDIT_LIBRARY = getattr(Constants, "AUDIT_LIBRARY", None) or os.path.join(_ROOT, "Audit_Library")


def _list_audit_runs():
    """Return [(id, abs_path), ...] newest first from Audit_Library."""
    runs = []
    root = AUDIT_LIBRARY
    if not os.path.isdir(root):
        return runs
    for name in os.listdir(root):
        if name in ("latest", "index.json"):
            continue
        path = os.path.join(root, name)
        if os.path.isdir(path) and (
            os.path.isfile(os.path.join(path, "TRIAGE_SUMMARY.json"))
            or os.path.isfile(os.path.join(path, "portal_snapshot.json"))
        ):
            runs.append((name, path))
    runs.sort(key=lambda x: x[0], reverse=True)
    return runs


def evidence_logs_dir():
    """Live Global_Evidence or selected Audit_Library run folder."""
    review = st.session_state.get("review_run_path")
    if review and os.path.isdir(review):
        return review
    return DIRS.get("OUT_LOGS") or os.path.join(BASE_OUT, "Global_Evidence")


def evidence_maps_dir():
    review = st.session_state.get("review_run_path")
    if review and os.path.isdir(review):
        # packages store waterfall/grid next to triage
        return review
    return DIRS.get("OUT_MAPS") or os.path.join(BASE_OUT, "Global_Heatmaps")


st.set_page_config(page_title="BSAP — Bio Security Auditing Platform", layout="wide")

# Auto-bootstrap once per session: folders, config, optional baseline FASTAs (DL non-interactive skip)
if "env_bootstrapped" not in st.session_state:
    try:
        Engine.ensure_environment()
        st.session_state["env_bootstrapped"] = True
        st.session_state["env_boot_error"] = None
    except Exception as _env_err:
        st.session_state["env_bootstrapped"] = False
        st.session_state["env_boot_error"] = str(_env_err)

# Optional standalone website server (index.html) — not embedded in Streamlit
def _portal_state():
    return st.session_state.setdefault("_portal", {"httpd": None, "url": None, "port": None})

def _start_portal_server():
    import socket
    import threading
    import webbrowser
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

    ps = _portal_state()
    if ps.get("httpd") is not None and ps.get("url"):
        try:
            webbrowser.open(ps["url"])
        except Exception:
            pass
        return ps["url"]

    def _free_port(start=8765):
        for port in range(start, start + 40):
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                try:
                    s.bind(("127.0.0.1", port))
                    return port
                except OSError:
                    continue
        return start

    class _Handler(SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=_ROOT, **kwargs)
        def log_message(self, fmt, *args):
            pass  # quiet in Streamlit console

    port = _free_port()
    httpd = ThreadingHTTPServer(("127.0.0.1", port), _Handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{port}/index.html"
    ps["httpd"] = httpd
    ps["url"] = url
    ps["port"] = port
    try:
        webbrowser.open(url)
    except Exception:
        pass
    return url

# --- HELPER: RECURSIVE FILE FINDER ---
def get_all_files(directory, extension):
    search_path = os.path.join(directory, "**", f"*{extension}")
    return glob.glob(search_path, recursive=True)

# --- HELPER: VALIDATION STATUS ---
def get_validation_status():
    """Returns (environment_initialized, audit_passed).
    TEST_MODE=1: the 9 calibration files exist (written as '<Label>_<name>') and the cross-link log links them.
    Production: every input folder holds FASTA files, and a report has been built."""
    test_mode = str(Constants.CONFIG.get("TEST_MODE", "0")) == "1"
    if test_mode:
        labels = {"CLN": "Clinical", "IND": "Industrial", "REF": "Reference"}
        names = ("DNA_Test.fasta", "RNA_Test.fasta", "Prot_Test.fasta")
        env_ok = all(os.path.exists(os.path.join(DIRS[c], f"{lbl}_{n}"))
                     for c, lbl in labels.items() for n in names)
        passed = False
        cross_log = os.path.join(DIRS["OUT_LOGS"], "CROSS_MATCH_LOG.txt")
        if os.path.exists(cross_log):
            with open(cross_log, "r", errors="replace") as f:
                content = f.read()
            passed = "DNA_Test.fasta" in content and "RNA_Test.fasta" in content
        return env_ok, passed

    def has_fasta(folder):
        if not os.path.isdir(folder):
            return False
        try:
            return any(
                f.lower().endswith((".fa", ".fasta", ".fna", ".fas"))
                for f in os.listdir(folder)
            )
        except OSError:
            return False

    # Real runs only need FASTA in at least one sector (same idea as env_needs_bootstrap).
    env_ok = any(has_fasta(DIRS[c]) for c in ("CLN", "IND", "REF"))
    return env_ok, os.path.exists(DIRS["REPORT"])

# --- SIDEBAR ---
st.sidebar.header("BSAP Control")

# Visual Badge Placement
env_initialized, passed_audit = get_validation_status()

if not env_initialized:
    status_text, status_color = "Environment Not Initialized", "red"
elif not passed_audit:
    status_text, status_color = "Initialized / Audit Not Yet Run", "orange"
else:
    status_text, status_color = "Validated", "green"

st.sidebar.caption("CLI (no Streamlit): `python Alpha4.py`")
st.sidebar.subheader("Standalone website")
st.sidebar.caption(
    "Optional viewer only — not required for audits. "
    "`index.html` loads `Forensic_Results/Global_Evidence/portal_snapshot.json` "
    "(and `Audit_Library/` history when present). "
    "Before the first audit the page may show demo data."
)
if st.sidebar.button("Open website (index.html)", key="btn_open_portal"):
    try:
        url = _start_portal_server()
        st.sidebar.success(f"Portal: {url}")
        snap = os.path.join(evidence_logs_dir(), "portal_snapshot.json")
        if not os.path.exists(snap):
            st.sidebar.warning("No portal_snapshot.json yet — run an audit first.")
    except Exception as e:
        st.sidebar.error(f"Portal failed: {e}")
st.sidebar.caption(
    "Or from a terminal in the BSAP folder: "
    "`python -m http.server 8765` → http://localhost:8765/index.html"
)
st.sidebar.divider()

st.sidebar.markdown(f"**System Status:** :{status_color}[{status_text}]")
st.sidebar.divider()
if "audit_msg" in st.session_state:
    st.sidebar.success(st.session_state.pop("audit_msg"))

st.sidebar.subheader("Target panel")
st.sidebar.caption(
    "Search subset only. Master Targets.json keeps the full list. "
    "Save panel → re-run audit."
)
_panel_opts = {
    "full (all markers)": "full",
    "lab_backbone": "lab_backbone",
    "viral_vector": "viral_vector",
    "composition (density only)": "composition",
    "custom path": "custom",
}
_cur = str(Constants.CONFIG.get("TARGETS_PANEL", "full")).strip().lower()
_labels = list(_panel_opts.keys())
_default_i = 0
for i, (_lab, _val) in enumerate(_panel_opts.items()):
    if _val == _cur:
        _default_i = i
        break
_sel_label = st.sidebar.selectbox("Active panel", _labels, index=_default_i, key="panel_sel")
_sel_val = _panel_opts[_sel_label]
_custom = st.sidebar.text_input(
    "Custom JSON path (if custom)",
    value=str(Constants.CONFIG.get("TARGETS_CUSTOM_PATH", "")),
    key="panel_custom",
)
if st.sidebar.button("💾 Save panel setting", key="btn_save_panel"):
    cfg_path = Constants.CONFIG_PATH
    try:
        lines = []
        if os.path.exists(cfg_path):
            with open(cfg_path, "r", encoding="utf-8") as cf:
                lines = cf.readlines()
        keys_done = set()
        out_lines = []
        for line in lines:
            if line.strip().startswith("TARGETS_PANEL="):
                out_lines.append(f"TARGETS_PANEL={_sel_val}\n")
                keys_done.add("TARGETS_PANEL")
            elif line.strip().startswith("TARGETS_CUSTOM_PATH="):
                out_lines.append(f"TARGETS_CUSTOM_PATH={_custom.strip()}\n")
                keys_done.add("TARGETS_CUSTOM_PATH")
            else:
                out_lines.append(line)
        if "TARGETS_PANEL" not in keys_done:
            out_lines.append(f"TARGETS_PANEL={_sel_val}\n")
        if "TARGETS_CUSTOM_PATH" not in keys_done:
            out_lines.append(f"TARGETS_CUSTOM_PATH={_custom.strip()}\n")
        with open(cfg_path, "w", encoding="utf-8") as cf:
            cf.writelines(out_lines)
        Constants.CONFIG = Constants.load_config()
        if hasattr(Constants, "resolve_targets_path"):
            Constants.TARGETS_PATH = Constants.resolve_targets_path(Constants.CONFIG)
        st.sidebar.success(f"Panel set to {_sel_val}. Re-run audit to apply.")
    except Exception as e:
        st.sidebar.error(f"Could not write Configuration.txt: {e}")
st.sidebar.caption(
    f"Current: `{os.path.basename(getattr(Constants, 'TARGETS_PATH', 'Targets.json'))}`"
)


st.sidebar.markdown("---")
st.sidebar.subheader("Review past run")
_runs = _list_audit_runs()
_run_labels = ["(live Forensic_Results)"] + [r[0] for r in _runs]
_run_paths = {r[0]: r[1] for r in _runs}
_default_idx = 0
_cur_review = st.session_state.get("review_run_id")
if _cur_review and _cur_review in _run_paths:
    try:
        _default_idx = _run_labels.index(_cur_review)
    except ValueError:
        _default_idx = 0
_pick = st.sidebar.selectbox(
    "Audit_Library package",
    _run_labels,
    index=min(_default_idx, max(0, len(_run_labels) - 1)),
    key="audit_run_picker",
    help="Load TRIAGE / portal / report from a stamped package without re-scanning.",
)
if _pick == "(live Forensic_Results)":
    st.session_state["review_run_path"] = None
    st.session_state["review_run_id"] = None
else:
    st.session_state["review_run_path"] = _run_paths.get(_pick)
    st.session_state["review_run_id"] = _pick
if st.session_state.get("review_run_path"):
    st.sidebar.caption(f"Reviewing: `{st.session_state['review_run_id']}`")
    _stamp_p = os.path.join(st.session_state["review_run_path"], "RUN_CONFIG_STAMP.json")
    if os.path.isfile(_stamp_p):
        try:
            import json as _json
            with open(_stamp_p, encoding="utf-8") as _sf:
                _st = _json.load(_sf)
            st.sidebar.caption(
                f"config `{_st.get('config_sha256_16', '')}` · {_st.get('created_utc', '')}"
            )
        except Exception:
            pass
    _xlsx = os.path.join(st.session_state["review_run_path"], "FINAL_GENOMIC_REPORT.xlsx")
    if os.path.isfile(_xlsx):
        with open(_xlsx, "rb") as _xf:
            st.sidebar.download_button(
                "Download package report",
                data=_xf.read(),
                file_name=f"BSAP_{st.session_state['review_run_id']}_report.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                key="dl_review_xlsx",
            )
else:
    st.sidebar.caption("Showing live Forensic_Results outputs")

if st.sidebar.button("🚀 Run Forensic Audit", key="btn_run_audit"):
    # Run as a separate process so Windows multiprocessing does NOT re-import
    # this Streamlit script (that spam: "use streamlit run [FILE_NAME]").
    with st.status("Scanning Sequences...", expanded=True) as status:
        st.write("Starting Alpha4 engine in a background process...")
        st.caption("Console progress is in the terminal. UI stays responsive after finish.")
        try:
            import subprocess
            env = os.environ.copy()
            env["PYTHONUNBUFFERED"] = "1"
            # Avoid Streamlit trying to attach to worker contexts
            env.pop("STREAMLIT_SERVER_PORT", None)
            proc = subprocess.run(
                [sys.executable, os.path.join(_ROOT, "Alpha4.py")],
                cwd=_ROOT,
                env=env,
                capture_output=False,  # live logs in the same terminal
            )
            if proc.returncode != 0:
                status.update(label="Audit failed", state="error")
                st.error(f"Alpha4 exited with code {proc.returncode}. See terminal for details.")
                st.stop()
        except Exception as exc:
            status.update(label="Audit failed", state="error")
            st.error(f"Audit failed: {exc}")
            st.stop()

        status.update(label="Audit Complete!", state="complete", expanded=False)

    st.session_state["audit_msg"] = "Analysis finished via Alpha4 (subprocess)."
    st.rerun()

# Download Report Button
if os.path.exists(DIRS["REPORT"]):
    with open(DIRS["REPORT"], "rb") as f:
        st.sidebar.download_button(
            label="📊 Download Excel Report",
            data=f,
            file_name="Forensic_Audit_Report.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

# Stop server cleanly (preferred over closing the browser alone)
st.sidebar.markdown("---")
st.sidebar.caption(
    "Closing the browser tab does **not** stop Streamlit. "
    "Use **Stop server** below or Ctrl+C in this terminal."
)
if st.sidebar.button("⏹ Stop server", key="btn_stop_server", type="primary"):
    st.sidebar.warning("Shutting down… return to the terminal.")
    # Hard exit so Windows does not hang on asyncio proactor / "Stopping..."
    os._exit(0)

# --- MAIN INTERFACE ---
st.title("BSAP — Bio Security Auditing Platform")

tab1, tab2, tab3 = st.tabs(["📊 Forensic Audit", "🧭 Snapshot triage", "⚙️ Manage Markers"])

with tab1:
    st.subheader("Visual Genetic Signatures")

    # --- WATERFALL (multi-file colour bands) ---
    st.markdown("### Multi-File Waterfall")
    st.caption(
        "Triage view: one row per file, columns = marker-tag categories. "
        "Not a medical diagnostic. Colour intensity = hit count."
    )
    show_comp = st.checkbox(
        "Show composition density (short / epigenetic motifs when stored as composition=true)",
        value=False,
        key="wf_show_comp",
    )
    waterfall_path = os.path.join(evidence_maps_dir(), "WATERFALL_SUMMARY.png")
    if os.path.exists(waterfall_path):
        st.image(waterfall_path, caption="Forensic Waterfall Summary", width="stretch")
    else:
        st.info("No waterfall yet — run an audit to generate colour bands across all scanned files.")

    try:
        names, cats, matrix, colors = Engine.build_waterfall_matrix(
            max_files=2000, include_comp=bool(show_comp)
        )
        if names and cats:
            if not cats:
                st.warning("All columns filtered out — enable composition density to see them.")
            else:
                page_size = st.slider("Files per page", 10, 100, 25, key="wf_page_size")
                n_pages = max(1, (len(names) + page_size - 1) // page_size)
                page = st.number_input("Page", 1, n_pages, 1, key="wf_page") - 1
                start, end = page * page_size, min((page + 1) * page_size, len(names))
                rows = []
                for i in range(start, end):
                    row = {"File": names[i]}
                    for j, cat in enumerate(cats):
                        row[cat] = matrix[i][j]
                    rows.append(row)
                st.dataframe(pd.DataFrame(rows), width="stretch", height=min(400, 40 + 28 * len(rows)))
                st.caption(
                    f"Showing files {start + 1}–{end} of {len(names)} · "
                    f"composition tags {'ON' if show_comp else 'HIDDEN'}"
                )
    except Exception as e:
        st.caption(f"Waterfall matrix unavailable: {e}")

    st.divider()
    st.subheader("Per-file outputs (Evidence · Flanks · Audio · Composition · Heatmaps)")
    st.caption(
        "Reads Forensic_Results on every page load. "
        "A new audit is only needed to create/update files — not to view existing ones."
    )

    subject_map = {}

    def _register(subject, kind, path):
        if not subject:
            return
        subject_map.setdefault(subject, {})[kind] = path

    for key, kind in [
        ("OUT_CLN_FLANKS", "flanks"),
        ("OUT_IND_FLANKS", "flanks"),
        ("OUT_REF_FLANKS", "flanks"),
        ("OUT_CLN_AUDIO", "wav"),
        ("OUT_IND_AUDIO", "wav"),
        ("OUT_REF_AUDIO", "wav"),
        ("OUT_CLN_EVI_FWD", "evi_fwd"),
        ("OUT_CLN_EVI_REV", "evi_rev"),
        ("OUT_IND_EVI_FWD", "evi_fwd"),
        ("OUT_IND_EVI_REV", "evi_rev"),
        ("OUT_REF_EVI_FWD", "evi_fwd"),
        ("OUT_REF_EVI_REV", "evi_rev"),
        ("OUT_CLN_MAP_FWD", "map_fwd"),
        ("OUT_CLN_MAP_REV", "map_rev"),
        ("OUT_CLN_MAP_FLANKS", "map_flanks"),
        ("OUT_IND_MAP_FWD", "map_fwd"),
        ("OUT_IND_MAP_REV", "map_rev"),
        ("OUT_IND_MAP_FLANKS", "map_flanks"),
        ("OUT_REF_MAP_FWD", "map_fwd"),
        ("OUT_REF_MAP_REV", "map_rev"),
        ("OUT_REF_MAP_FLANKS", "map_flanks"),
        ("OUT_CLN_SPECTROGRAM_FWD", "spectrogram_fwd"),
        ("OUT_IND_SPECTROGRAM_FWD", "spectrogram_fwd"),
        ("OUT_REF_SPECTROGRAM_FWD", "spectrogram_fwd"),
        ("OUT_CLN_SPECTROGRAM_REV", "spectrogram_rev"),
        ("OUT_IND_SPECTROGRAM_REV", "spectrogram_rev"),
        ("OUT_REF_SPECTROGRAM_REV", "spectrogram_rev"),
        ("OUT_CLN_SPECTROGRAM_FLANKS", "spectrogram_flanks"),
        ("OUT_IND_SPECTROGRAM_FLANKS", "spectrogram_flanks"),
        ("OUT_REF_SPECTROGRAM_FLANKS", "spectrogram_flanks"),
    ]:
        
        d = DIRS.get(key)
        if not d or not os.path.isdir(d):
            continue
        for root, _, files in os.walk(d):
            for fname in files:
                path = os.path.join(root, fname)
                if kind == "flanks" and fname.endswith("_flanks.txt"):
                    _register(fname.replace("_flanks.txt", ""), "flanks", path)
                elif kind == "wav" and fname.lower().endswith(".wav"):
                    subj = fname.replace("_AUDIO.wav", "").replace(".wav", "")
                    _register(subj, "wav", path)
                elif kind.startswith("evi") and "_Evidence_" in fname and fname.endswith(".txt"):
                    after = fname.split("_Evidence_", 1)[1]
                    parts = after.split("_", 2)
                    if len(parts) >= 3:
                        _register(parts[2].replace(".txt", ""), kind, path)
                elif kind == "map_flanks" and fname.endswith("_FLANKS_HEATMAP.png"):
                    _register(fname.replace("_FLANKS_HEATMAP.png", ""), "map_flanks", path)
                elif kind.startswith("map") and fname.endswith("_HEATMAP.png"):
                    _register(fname.replace("_HEATMAP.png", ""), kind, path)
                elif kind == "spectrogram_fwd" and fname.endswith("_FEATURE_SPECTROGRAM.png") and "_RC" not in fname and "_FLANKS" not in fname:
                    _register(fname.replace("_FEATURE_SPECTROGRAM.png", ""), "spectrogram_fwd", path)
                elif kind == "spectrogram_rev" and "_FEATURE_SPECTROGRAM_RC.png" in fname:
                    _register(fname.replace("_FEATURE_SPECTROGRAM_RC.png", ""), "spectrogram_rev", path)
                elif kind == "spectrogram_flanks" and "_FEATURE_SPECTROGRAM_FLANKS.png" in fname:
                    _register(fname.replace("_FEATURE_SPECTROGRAM_FLANKS.png", ""), "spectrogram_flanks", path)

    for comp_root in [
        os.path.join(DIRS.get("OUT_LOGS", BASE_OUT), "composition"),
        os.path.join(BASE_OUT, "Global_Evidence", "composition"),
    ]:
        if not os.path.isdir(comp_root):
            continue
        for fname in os.listdir(comp_root):
            path = os.path.join(comp_root, fname)
            if fname.endswith(".TARGET_COUNTS.json") or fname.endswith(".composition.json"):
                subj = fname.replace(".TARGET_COUNTS.json", "").replace(".composition.json", "")
                _register(subj, "composition_json", path)
            elif fname.endswith(".TARGET_COUNTS.txt"):
                _register(fname.replace(".TARGET_COUNTS.txt", ""), "composition_txt", path)

    subjects = sorted(subject_map.keys())
    st.caption(f"{len(subjects)} files with at least one artifact under Forensic_Results")

    if not subjects:
        st.info("No outputs found. Run an audit once, or confirm Forensic_Results is next to APP.py.")
    else:
        fpage_size = st.slider("Files per page", 10, 100, 25, key="file_page_size")
        fn_pages = max(1, (len(subjects) + fpage_size - 1) // fpage_size)
        fpage = st.number_input("File list page", 1, fn_pages, 1, key="file_page") - 1
        fstart, fend = fpage * fpage_size, min((fpage + 1) * fpage_size, len(subjects))
        page_subjects = subjects[fstart:fend]

        rows = []
        for s in page_subjects:
            info = subject_map[s]
            rows.append({
                "File": s[:70],
                "FWD": "Y" if info.get("evi_fwd") else "",
                "REV": "Y" if info.get("evi_rev") else "",
                "Flanks": "Y" if info.get("flanks") else "",
                "WAV": "Y" if info.get("wav") else "",
                "Comp": "Y" if info.get("composition_json") or info.get("composition_txt") else "",
                "Map": "Y" if info.get("map_fwd") or info.get("map_rev") else "",
                "FlankMap": "Y" if info.get("map_flanks") else "",
                "Spec": "Y" if info.get("spectrogram_fwd") or info.get("spectrogram_rev") or info.get("spectrogram_flanks") else "",
            })
        st.dataframe(pd.DataFrame(rows), width="stretch", height=min(420, 40 + 28 * len(rows)))
        st.caption(f"List {fstart + 1}–{fend} of {len(subjects)}")

        selected = st.selectbox("Open file detail", page_subjects, key="file_detail_select")
        if selected:
            info = subject_map[selected]
            c1, c2 = st.columns(2)
            with c1:
                st.markdown(f"**{selected}**")
                if info.get("composition_txt"):
                    st.markdown("**Target counts**")
                    with open(info["composition_txt"], "r", errors="replace") as f:
                        st.text(f.read())
                elif info.get("composition_json"):
                    with open(info["composition_json"], "r", errors="replace") as f:
                        st.json(json.load(f))
                if info.get("flanks"):
                    st.markdown("**Flanks** (marker · pos · strand · left · match · right)")
                    with open(info["flanks"], "r", errors="replace") as f:
                        st.text(f.read())
            with c2:
                if info.get("wav") and os.path.exists(info["wav"]):
                    st.markdown("**Sonification (WAV)**")
                    st.audio(info["wav"])
                for label, k in [
                    ("Heatmap FWD", "map_fwd"),
                    ("Heatmap REV", "map_rev"),
                    ("Flanks heatmap (FWD+REV)", "map_flanks"),
                    ("Spectrogram FWD (NOT RF)", "spectrogram_fwd"),
                    ("Spectrogram RC (NOT RF)", "spectrogram_rev"),
                    ("Spectrogram FLANKS strip (NOT RF)", "spectrogram_flanks"),
                ]:
                    if info.get(k) and os.path.exists(info[k]):
                        st.image(info[k], caption=label, width="stretch")
                for label, k in [("Evidence FWD", "evi_fwd"), ("Evidence REV", "evi_rev")]:
                    if info.get(k):
                        with st.expander(label):
                            with open(info[k], "r", errors="replace") as f:
                                st.text(f.read())

    st.divider()
    st.subheader("Feature spectrograms (composition display — not measured RF)")
    st.caption(
        "Per-file maps under Clinical|Industrial|Reference_Spectrograms/{FWD,REV,FLANKS}. "
        "ENABLE_RF_HEURISTICS=1; SPECTROGRAM_RC / SPECTROGRAM_FLANKS toggle RC and flanks strip."
    )
    spec_maps = []
    for key, label in (
        ("OUT_CLN_SPECTROGRAM_FWD", "FWD"),
        ("OUT_IND_SPECTROGRAM_FWD", "FWD"),
        ("OUT_REF_SPECTROGRAM_FWD", "FWD"),
        ("OUT_CLN_SPECTROGRAM_REV", "RC"),
        ("OUT_IND_SPECTROGRAM_REV", "RC"),
        ("OUT_REF_SPECTROGRAM_REV", "RC"),
        ("OUT_CLN_SPECTROGRAM_FLANKS", "FLANKS"),
        ("OUT_IND_SPECTROGRAM_FLANKS", "FLANKS"),
        ("OUT_REF_SPECTROGRAM_FLANKS", "FLANKS"),
    ):
        d = DIRS.get(key)
        if d and os.path.isdir(d):
            for root, _, files in os.walk(d):
                for fname in files:
                    if fname.endswith(".png") and "SPECTROGRAM" in fname.upper():
                        spec_maps.append((label, os.path.join(root, fname)))
    spec_maps = sorted(spec_maps, key=lambda x: x[1])
    if not spec_maps:
        st.info(
            "No feature spectrograms yet — set ENABLE_RF_HEURISTICS=1 "
            "(and optionally SPECTROGRAM_RC=1 / SPECTROGRAM_FLANKS=1) then re-run audit."
        )
    else:
        sm_size = st.slider("Spectrograms per page", 5, 30, 10, key="spec_page_size")
        sm_pages = max(1, (len(spec_maps) + sm_size - 1) // sm_size)
        sm_page = st.number_input("Spectrogram page", 1, sm_pages, 1, key="spec_page") - 1
        sa, sb = sm_page * sm_size, min((sm_page + 1) * sm_size, len(spec_maps))
        for label, path in spec_maps[sa:sb]:
            base = os.path.basename(path)
            with st.expander(f"📡 [{label}] {base}", expanded=False):
                st.image(path, caption=f"{label} — NOT RF measurement", width="stretch")
        st.caption(f"Spectrograms {sa + 1}–{sb} of {len(spec_maps)}")

    st.divider()
    st.subheader("Flanks heatmaps (scroll — open for full FWD/REV)")
    st.caption("Combined FWD+REV marker map per file under */Heatmaps/FLANKS/")

    flanks_maps = []
    for key in ("OUT_CLN_MAP_FLANKS", "OUT_IND_MAP_FLANKS", "OUT_REF_MAP_FLANKS"):
        d = DIRS.get(key)
        if d and os.path.isdir(d):
            for root, _, files in os.walk(d):
                for fname in files:
                    if fname.endswith("_FLANKS_HEATMAP.png"):
                        flanks_maps.append(os.path.join(root, fname))
    flanks_maps = sorted(flanks_maps)

    if not flanks_maps:
        st.info("No FLANKS heatmaps yet — run an audit after Alpha4 MAP_FLANKS + Constants paths are in place.")
    else:
        fm_size = st.slider("Flanks maps per page", 5, 30, 10, key="fm_page_size")
        fm_pages = max(1, (len(flanks_maps) + fm_size - 1) // fm_size)
        fm_page = st.number_input("Flanks map page", 1, fm_pages, 1, key="fm_page") - 1
        a, b = fm_page * fm_size, min((fm_page + 1) * fm_size, len(flanks_maps))
        for path in flanks_maps[a:b]:
            base = os.path.basename(path).replace("_FLANKS_HEATMAP.png", "")
            with st.expander(f"🧬 {base}", expanded=False):
                st.image(path, caption="Flanks map (FWD+REV)", width="stretch")
                info = subject_map.get(base, {})
                cols = st.columns(2)
                if info.get("map_fwd"):
                    cols[0].image(info["map_fwd"], caption="FWD", width="stretch")
                if info.get("map_rev"):
                    cols[1].image(info["map_rev"], caption="REV", width="stretch")
        st.caption(f"Flanks maps {a + 1}–{b} of {len(flanks_maps)}")

    st.divider()
    st.subheader("Cross-file similarity tracker")
    st.caption(
        "Exact shared LEFT/RIGHT flanks between files (Phase 2). "
        "Not whole-genome identity — ranks genomes that reuse the same marker context."
    )
    sim_json = os.path.join(evidence_logs_dir(), "CROSS_SIMILARITY_TRACKER.json")
    if os.path.exists(sim_json):
        try:
            with open(sim_json, "r", encoding="utf-8") as f:
                sim = json.load(f)
            st.caption(sim.get("description", ""))
            edges = sim.get("edges") or []
            if edges:
                table_rows = []
                for e in edges:
                    table_rows.append({
                        "file_a": e.get("file_a", ""),
                        "file_b": e.get("file_b", ""),
                        "shared_flank_hits": e.get("shared_flank_hits", 0),
                        "pair_support_count": e.get("pair_support_count", 0),
                        "similarity_score": e.get("similarity_score", 0),
                        "same_accession": e.get("same_accession", False),
                        "markers": "; ".join(e.get("markers") or []),
                        "note": e.get("note", ""),
                    })
                st.dataframe(
                    pd.DataFrame(table_rows),
                    width="stretch",
                    height=min(360, 40 + 28 * min(len(table_rows), 12)),
                )
                st.caption(
                    f"{len(table_rows)} file-pair edges · "
                    f"{sim.get('n_raw_hits', '?')} raw Phase 2 hits"
                )
            else:
                st.info("Tracker exists but has no edges yet.")
        except Exception as e:
            st.warning(f"Could not read similarity tracker: {e}")
    else:
        st.info("No CROSS_SIMILARITY_TRACKER.json yet — run a full audit (Phase 2).")

    st.divider()
    st.subheader("Triage — where to look next")
    st.caption(
        "Ranked files and pairs from lab markers + flank graph. "
        "Research navigation only — not a diagnosis or origin proof."
    )
    triage_path = os.path.join(evidence_logs_dir(), "TRIAGE_SUMMARY.json")
    if os.path.exists(triage_path):
        try:
            with open(triage_path, "r", encoding="utf-8") as f:
                triage = json.load(f)
            st.caption(triage.get("description", ""))
            c1, c2 = st.columns(2)
            with c1:
                st.markdown("**Top files by lab-marker count** (composition excluded)")
                top_f = triage.get("top_files_by_lab_markers") or []
                _tmin = triage.get("triage_min_len")
                _stamp = triage.get("config_stamp") or {}
                if _stamp.get("config_sha256_16"):
                    st.caption(
                        f"Run config `{_stamp.get('config_sha256_16')}` · {_stamp.get('created_utc', '')}"
                    )
                if _tmin is not None:
                    st.caption(
                        f"triage_score = long_hits(≥{_tmin} bp) + 5×unique_tags + Phase-2 flank support. "
                        f"Short motifs excluded from the score."
                    )
                if top_f:
                    _df = pd.DataFrame(top_f)
                    _pref = [
                        c for c in (
                            "file",
                            "triage_score",
                            "lab_marker_hits",
                            "unique_long_tags",
                            "phase2_flank_support",
                            "short_forensic_hits_excluded",
                            "triage_min_len",
                            "note",
                        )
                        if c in _df.columns
                    ]
                    _rest = [c for c in _df.columns if c not in _pref]
                    st.dataframe(
                        _df[_pref + _rest] if _pref else _df,
                        width="stretch",
                        height=min(280, 40 + 28 * min(len(top_f), 8)),
                    )
                else:
                    st.info("No lab-marker tallies yet.")
            with c2:
                st.markdown("**Multi-domain files** (viral + lab/chassis + tag)")
                multi = triage.get("multi_domain_files") or []
                if multi:
                    mrows = [{
                        "file": m.get("file"),
                        "domains": m.get("domains"),
                        "viral": m.get("viral"),
                        "bacterial_or_lab": m.get("bacterial_or_lab"),
                        "tag_or_reporter": m.get("tag_or_reporter"),
                    } for m in multi]
                    st.dataframe(
                        pd.DataFrame(mrows), width="stretch",
                        height=min(280, 40 + 28 * min(len(mrows), 8)),
                    )
                else:
                    st.info("No multi-domain files flagged.")
            st.markdown("**Top file pairs by flank score**")
            pairs = triage.get("top_file_pairs_by_flank_score") or []
            if pairs:
                prows = [{
                    "file_a": e.get("file_a"),
                    "file_b": e.get("file_b"),
                    "score": e.get("similarity_score"),
                    "shared_flanks": e.get("shared_flank_hits"),
                    "pair_support": e.get("pair_support_count"),
                    "same_accession": e.get("same_accession"),
                    "markers": "; ".join(e.get("markers") or []),
                } for e in pairs]
                st.dataframe(
                    pd.DataFrame(prows), width="stretch",
                    height=min(320, 40 + 28 * min(len(prows), 10)),
                )
            else:
                st.info("No pair rankings yet — need Phase 2.")
            bed_dir = os.path.join(evidence_logs_dir(), "coordinates")
            if os.path.isdir(bed_dir) and any(x.endswith(".bed") for x in os.listdir(bed_dir)):
                st.caption(f"BED tracks (long hits): `{bed_dir}` — open in IGV / genome browsers.")
        except Exception as e:
            st.warning(f"Could not read triage summary: {e}")
    else:
        st.info("No TRIAGE_SUMMARY.json yet — run a full audit.")

    st.divider()
    st.subheader("Raw log search")
    all_logs = get_all_files(BASE_OUT, ".txt")
    if all_logs:
        log_map = {os.path.basename(p): p for p in sorted(all_logs)}
        selected_log = st.selectbox("Select log file", list(log_map.keys()), key="raw_log_select")
        if selected_log:
            with open(log_map[selected_log], "r", errors="replace") as f:
                st.text(f.read())
    else:
        st.warning("No .txt logs under Forensic_Results yet.")    


with tab2:
    st.subheader("Snapshot triage")
    st.caption(
        "Reads `Forensic_Results/Global_Evidence/portal_snapshot.json` from the last audit. "
        "Same JSON powers the standalone `index.html` website (sidebar → Open website, "
        "or `python -m http.server` from the project root)."
    )
    snap_path = os.path.join(evidence_logs_dir(), "portal_snapshot.json")
    snap = None
    if os.path.exists(snap_path):
        try:
            with open(snap_path, "r", encoding="utf-8") as sf:
                snap = json.load(sf)
        except Exception as e:
            st.error(f"Could not read portal_snapshot.json: {e}")
    if not snap:
        st.info(
            "No portal_snapshot.json yet. Run a full forensic audit first. "
            "The website (`index.html`) uses the same file when served from the project root."
        )
    else:
        files = snap.get("files") or []
        edges = snap.get("edges") or []
        cats = snap.get("categories") or []
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Files", f"{len(files):,}")
        m2.metric("Tags", f"{len(cats):,}")
        m3.metric("Cross edges", f"{len(edges):,}")
        m4.metric("Incorp. links", f"{len(snap.get('incorporation') or []):,}")
        st.caption(f"Live from `{snap_path}`")

        view = st.radio(
            "View",
            ["Waterfall matrix", "Per-file inventory", "Cross-similarity", "Incorporation map", "Inc summary"],
            horizontal=True,
            key="snap_view_radio",
        )

        if view == "Waterfall matrix":
            if not files:
                st.warning("Snapshot has no files.")
            else:
                rows = []
                for f in files:
                    row = {"file": f.get("name", ""), "total": int(f.get("total_hits") or 0)}
                    by = f.get("by_tag") or {}
                    for c in cats:
                        row[c] = int(by.get(c, 0))
                    rows.append(row)
                st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

        elif view == "Per-file inventory":
            inv = []
            for f in files:
                inv.append({
                    "file": f.get("name", ""),
                    "category": f.get("category", ""),
                    "total_hits": int(f.get("total_hits") or 0),
                    "tags": ", ".join(sorted((f.get("by_tag") or {}).keys())),
                })
            if inv:
                st.dataframe(pd.DataFrame(inv), width="stretch", hide_index=True)
            else:
                st.info("No files in snapshot.")

        elif view == "Cross-similarity":
            if not edges:
                st.info("No cross-similarity edges in this snapshot.")
            else:
                erows = []
                for e in edges:
                    if not isinstance(e, dict):
                        continue
                    erows.append({
                        "file_a": e.get("a") or e.get("file_a") or e.get("source", ""),
                        "file_b": e.get("b") or e.get("file_b") or e.get("target", ""),
                        "score": e.get("similarity_score") or e.get("score") or e.get("shared_flank_hits") or "",
                        "label": e.get("label") or e.get("top_label") or "",
                        "markers": e.get("markers") or e.get("top_markers") or "",
                    })
                st.dataframe(pd.DataFrame(erows), width="stretch", hide_index=True)

        elif view == "Incorporation map":
            st.caption(
                "100% flank identity from marker A next to a different forensic marker B. "
                "Not whole-genome identity."
            )
            grid_png = os.path.join(evidence_maps_dir(), "INCORPORATION_GRID.png")
            if os.path.exists(grid_png):
                st.image(grid_png, caption="Incorporation grid", width="stretch")
            inc = snap.get("incorporation") or []
            hide_dup = st.checkbox("Hide same-accession pairs", value=True, key="inc_hide_dup")
            min_c = st.number_input("Min context count", min_value=0, value=1, key="inc_min_c")
            rows = []
            for r in inc:
                if not isinstance(r, dict):
                    continue
                if hide_dup and r.get("same_accession"):
                    continue
                if int(r.get("context_count") or 0) < int(min_c):
                    continue
                rows.append({
                    "origin": r.get("origin_file") or r.get("origin", ""),
                    "target": r.get("target_file") or r.get("target", ""),
                    "context_count": r.get("context_count", 0),
                    "max_risk": r.get("max_risk", 0),
                    "same_accession": r.get("same_accession", False),
                    "foreign_markers": "; ".join(r.get("foreign_markers") or [])
                    if isinstance(r.get("foreign_markers"), list)
                    else (r.get("foreign_markers") or ""),
                })
            if not rows:
                st.info("No incorporation links (or all filtered out).")
            else:
                st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

        elif view == "Inc summary":
            summary = snap.get("inc_summary") or {}
            fam = summary.get("family_file_counts") or {}
            if fam:
                st.write("**Family presence (file counts)**")
                st.dataframe(
                    pd.DataFrame(
                        [{"family": k, "files": v} for k, v in sorted(fam.items(), key=lambda x: -x[1])]
                    ),
                    width="stretch",
                    hide_index=True,
                )
            rows = summary.get("files") or []
            if not rows:
                st.info("No Inc summary yet — needs Phase-2 edges from a full audit.")
            else:
                table = []
                for r in rows:
                    table.append({
                        "file": r.get("file", ""),
                        "families": "; ".join(r.get("families") or []),
                        "inc_markers": "; ".join(r.get("inc_markers") or []),
                        "chassis": r.get("chassis_origin_hits", 0),
                        "backbone": r.get("industrial_backbone_hits", 0),
                        "seamless": r.get("seamless_assembly_hits", 0),
                        "AmpR-class": r.get("antibiotic_res_hits", 0),
                        "reporter": r.get("reporter_hits", 0),
                    })
                st.dataframe(pd.DataFrame(table), width="stretch", hide_index=True)
                if summary.get("description"):
                    st.caption(summary["description"])


with tab3:
    st.header("Live Manifest Editor")
    try:
        with open(Constants.TARGETS_PATH, "r", encoding="utf-8") as f:
            manifest_data = json.load(f)
        
        new_json_str = st.text_area(
            "Targets.json Raw Configuration", 
            value=json.dumps(manifest_data, indent=4), 
            height=500
        )
        
        if st.button("💾 Save"):
            try:
                updated_data = json.loads(new_json_str)
                with open(Constants.TARGETS_PATH, "w", encoding="utf-8") as f:
                    json.dump(updated_data, f, indent=4)
                st.toast("Manifest Saved!", icon="✅")
            except Exception as e:
                st.error(f"Format Error: {e}")
    except FileNotFoundError:
        st.error("Targets.json not found.")