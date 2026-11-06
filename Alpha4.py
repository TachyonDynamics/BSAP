import gc
import Compiler
import BiophysicalProfile
import Sonification
import Constants
import csv
import hashlib
import io
import json
from collections import Counter
import logging
import math
import mmap
import os  
import platform
import random
import re
import shutil
import sys           
import time
import ahocorasick
import importlib.util 
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import psutil
import requests
from tqdm import tqdm
import multiprocessing as mp
tqdm.set_lock(mp.RLock())

matplotlib.use('Agg')

__DEVELOPER_SIGNATURE__ = "d9e03009772a08a287e07673515f48795744f4340c4983636531983086e3f225"

DYNAMIC_WINDOW = int(Constants.CONFIG.get('FLANK_SIZE', 42))

# Create a single, unified logger
logger = logging.getLogger("BSAP")
logger.setLevel(logging.INFO)

class _TqdmStreamHandler(logging.StreamHandler):
    """Write log lines via tqdm so they don't corrupt progress bars."""
    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
        except Exception:
            self.handleError(record)

sh = _TqdmStreamHandler()
sh.setFormatter(logging.Formatter("[*] %(message)s"))
logger.addHandler(sh)

# Now use logger.info() or logger.error() instead of print() or open("failed_files.log")

DIRS = {}
TARGETS = {}
MARKER_COLORS = {}
AUTO = None  
cfg = None
AUTO_DNA = None
AUTO_PROT = None
AMPLITUDE = 1.0                 # The height of the match peak
FREQUENCY = 0.05                # The "tightness" of the curve
OFFSET = 0.2                    # The baseline "noise" level of the genome
LARGE_FILE_THRESHOLD_MB = 200
STREAM_CHUNK_MB = 64
STREAM_OVERLAP_BP = 512
COMPOSITION_MAX_LEN = 8
# Patterns shorter than this never enter Aho-Corasick (CG, NGG, …).
# Config key COMPOSITION_AC_MIN_LEN overrides at runtime (default 6).
COMPOSITION_AC_MIN_LEN_DEFAULT = 6

MIN_FLANK_LEN = int(Constants.CONFIG.get("MIN_FLANK_LEN", 20))   # shorter flanks match by chance
FEED_MARKERS = {}                                                # markers added by sync_global_threats()

# Bytes that are valid sequence characters (DNA/RNA/IUPAC + amino acids); everything else is stripped.
_VALID_SEQ_BYTES = b"ATCGUNRYSWKMBDHVPQRDEFGHIKLMNVWY"
_VALID_SEQ_BYTES += _VALID_SEQ_BYTES.lower()
_DELETE_BYTES = bytes(set(range(256)) - set(_VALID_SEQ_BYTES))
_HEADER_RE = re.compile(rb"^>[^\n]*\n?", re.M)                    # every FASTA header line, not just the first

COMPOSITION_FORCE = {
    "DNMT_Target_CG",
    "SpCas9_PAM_NGG",
    "AsCas12a_PAM",
    "PolyA_Signal",
    "Kozak_Consensus",
    "Shine_Dalgarno",
    "EcoRI_Site",
    "BamHI_Site",
    "HindIII_Site",
    "XbaI_Site",
    "NdeI_Site",
    "XhoI_Site",
    "BglII_Site",
    "BsaI_Scar",
    "BsmBI_Scar",
    "SapI_Scar",
    "BbsI_Scar",
    
}

FORENSIC_ALWAYS = {
    "TEST_MARKER",
}

CODON_TABLE = {
    'HUMAN': {
        # Alanine
        'GCT': 0.26, 'GCC': 0.40, 'GCA': 0.23, 'GCG': 0.11,
        # Arginine
        'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21, 'AGA': 0.20, 'AGG': 0.20,
        # Asparagine
        'AAT': 0.46, 'AAC': 0.54,
        # Aspartic Acid
        'GAT': 0.46, 'GAC': 0.54,
        # Cysteine
        'TGT': 0.45, 'TGC': 0.55,
        # Glutamine
        'CAA': 0.25, 'CAG': 0.75,
        # Glutamic Acid
        'GAA': 0.42, 'GAG': 0.58,
        # Glycine
        'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25,
        # Histidine
        'CAT': 0.41, 'CAC': 0.59,
        # Isoleucine
        'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16,
        # Leucine
        'TTA': 0.07, 'TTG': 0.13, 'CTT': 0.13, 'CTC': 0.20, 'CTA': 0.07, 'CTG': 0.40,
        # Lysine
        'AAA': 0.42, 'AAG': 0.58,
        # Methionine
        'ATG': 1.00,
        # Phenylalanine
        'TTT': 0.45, 'TTC': 0.55,
        # Proline
        'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11,
        # Serine
        'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06, 'AGT': 0.15, 'AGC': 0.24,
        # Threonine
        'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.12,
        # Tryptophan
        'TGG': 1.00,
        # Tyrosine
        'TAT': 0.43, 'TAC': 0.57,
        # Valine
        'GTT': 0.18, 'GTC': 0.24, 'GTA': 0.11, 'GTG': 0.47,
        # STOP
        'TAA': 0.28, 'TAG': 0.20, 'TGA': 0.52
    },
    'E_COLI': {
        'GCT': 0.18, 'GCC': 0.26, 'GCA': 0.23, 'GCG': 0.33, # Ala
        'CGT': 0.36, 'CGC': 0.36, 'CGA': 0.07, 'CGG': 0.02, 'AGA': 0.07, 'AGG': 0.02, # Arg
        'AAT': 0.49, 'AAC': 0.51, # Asn
        'GAT': 0.63, 'GAC': 0.37, # Asp
        'TGT': 0.46, 'TGC': 0.54, # Cys
        'CAA': 0.34, 'CAG': 0.66, # Gln
        'GAA': 0.68, 'GAG': 0.32, # Glu
        'GGT': 0.35, 'GGC': 0.37, 'GGA': 0.13, 'GGG': 0.15, # Gly
        'CAT': 0.57, 'CAC': 0.43, # His
        'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11, # Ile
        'TTA': 0.14, 'TTG': 0.13, 'CTT': 0.12, 'CTC': 0.10, 'CTA': 0.04, 'CTG': 0.47, # Leu
        'AAA': 0.74, 'AAG': 0.26, # Lys
        'ATG': 1.00, # Met
        'TTT': 0.58, 'TTC': 0.42, # Phe
        'CCT': 0.18, 'CCC': 0.13, 'CCA': 0.20, 'CCG': 0.49, # Pro
        'TCT': 0.17, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14, 'AGT': 0.16, 'AGC': 0.25, # Ser
        'ACT': 0.19, 'ACC': 0.40, 'ACA': 0.17, 'ACG': 0.25, # Thr
        'TGG': 1.00, # Trp
        'TAT': 0.59, 'TAC': 0.41, # Tyr
        'GTT': 0.28, 'GTC': 0.20, 'GTA': 0.17, 'GTG': 0.35, # Val
        'TAA': 0.61, 'TAG': 0.09, 'TGA': 0.30  # Stop
    },
    'S_CEREVISIAE': { # Baker's Yeast
        'GCT': 0.46, 'GCC': 0.22, 'GCA': 0.21, 'GCG': 0.11, 'CGT': 0.47, 'CGC': 0.19, 
        'CGA': 0.07, 'CGG': 0.04, 'AGA': 0.18, 'AGG': 0.05, 'AAT': 0.59, 'AAC': 0.41,
        'GAT': 0.65, 'GAC': 0.35, 'TGT': 0.67, 'TGC': 0.33, 'CAA': 0.69, 'CAG': 0.31,
        'GAA': 0.70, 'GAG': 0.30, 'GGT': 0.47, 'GGC': 0.19, 'GGA': 0.22, 'GGG': 0.12,
        'CAT': 0.64, 'CAC': 0.36, 'ATT': 0.46, 'ATC': 0.26, 'ATA': 0.27, 'TTA': 0.26,
        'TTG': 0.27, 'CTT': 0.12, 'CTC': 0.06, 'CTA': 0.14, 'CTG': 0.10, 'AAA': 0.58,
        'AAG': 0.42, 'ATG': 1.00, 'TTT': 0.59, 'TTC': 0.41, 'CCT': 0.31, 'CCC': 0.15,
        'CCA': 0.41, 'CCG': 0.13, 'TCT': 0.28, 'TCC': 0.16, 'TCA': 0.21, 'TCG': 0.10,
        'AGT': 0.16, 'AGC': 0.11, 'ACT': 0.35, 'ACC': 0.22, 'ACA': 0.30, 'ACG': 0.13,
        'TGG': 1.00, 'TAT': 0.56, 'TAC': 0.44, 'GTT': 0.39, 'GTC': 0.21, 'GTA': 0.21,
        'GTG': 0.19, 'TAA': 0.48, 'TAG': 0.24, 'TGA': 0.28
    },
    'BACILLUS_SUBTILIS': {
        'GCT': 0.31, 'GCC': 0.15, 'GCA': 0.35, 'GCG': 0.19,
        'CGT': 0.25, 'CGC': 0.10, 'CGA': 0.15, 'CGG': 0.08, 'AGA': 0.30, 'AGG': 0.12,
        'AAT': 0.58, 'AAC': 0.42, 'GAT': 0.65, 'GAC': 0.35, 'TGT': 0.60, 'TGC': 0.40,
        'CAA': 0.65, 'CAG': 0.35, 'GAA': 0.72, 'GAG': 0.28, 'GGT': 0.35, 'GGC': 0.20,
        'GGA': 0.30, 'GGG': 0.15, 'CAT': 0.60, 'CAC': 0.40, 'ATT': 0.45, 'ATC': 0.35,
        'ATA': 0.20, 'TTA': 0.25, 'TTG': 0.20, 'CTT': 0.15, 'CTC': 0.10, 'CTA': 0.10,
        'CTG': 0.20, 'AAA': 0.70, 'AAG': 0.30, 'ATG': 1.00, 'TTT': 0.62, 'TTC': 0.38,
        'CCT': 0.35, 'CCC': 0.10, 'CCA': 0.35, 'CCG': 0.20, 'TCT': 0.30, 'TCC': 0.15,
        'TCA': 0.25, 'TCG': 0.10, 'AGT': 0.12, 'AGC': 0.08, 'ACT': 0.35, 'ACC': 0.15,
        'ACA': 0.35, 'ACG': 0.15, 'TGG': 1.00, 'TAT': 0.65, 'TAC': 0.35, 'GTT': 0.35,
        'GTC': 0.15, 'GTA': 0.35, 'GTG': 0.15, 'TAA': 0.70, 'TAG': 0.10, 'TGA': 0.20
    },
    'P_PASTORIS': {
        'GCT': 0.45, 'GCC': 0.25, 'GCA': 0.20, 'GCG': 0.10,
        'CGT': 0.15, 'CGC': 0.10, 'CGA': 0.10, 'CGG': 0.05, 'AGA': 0.50, 'AGG': 0.10,
        'AAT': 0.55, 'AAC': 0.45, 'GAT': 0.60, 'GAC': 0.40, 'TGT': 0.70, 'TGC': 0.30,
        'CAA': 0.75, 'CAG': 0.25, 'GAA': 0.70, 'GAG': 0.30, 'GGT': 0.40, 'GGC': 0.20,
        'GGA': 0.30, 'GGG': 0.10, 'CAT': 0.65, 'CAC': 0.35, 'ATT': 0.50, 'ATC': 0.35,
        'ATA': 0.15, 'TTA': 0.30, 'TTG': 0.30, 'CTT': 0.15, 'CTC': 0.10, 'CTA': 0.10,
        'CTG': 0.05, 'AAA': 0.65, 'AAG': 0.35, 'ATG': 1.00, 'TTT': 0.55, 'TTC': 0.45,
        'CCT': 0.40, 'CCC': 0.20, 'CCA': 0.30, 'CCG': 0.10, 'TCT': 0.30, 'TCC': 0.25,
        'TCA': 0.20, 'TCG': 0.10, 'AGT': 0.10, 'AGC': 0.05, 'ACT': 0.35, 'ACC': 0.30,
        'ACA': 0.25, 'ACG': 0.10, 'TGG': 1.00, 'TAT': 0.60, 'TAC': 0.40, 'GTT': 0.40,
        'GTC': 0.25, 'GTA': 0.25, 'GTG': 0.10, 'TAA': 0.50, 'TAG': 0.20, 'TGA': 0.30
    }
}

def _file_size_mb(path):
    try:
        return os.path.getsize(path) / (1024 * 1024)
    except OSError:
        return 0.0

def normalize_forensic_tables():
    """
    Ensures all organism codon frequencies represent valid probability distributions.
    Fixes the 1.17 Yeast Leucine error and minor Human/E.coli rounding offsets.
    """
    global CODON_TABLE
    # Map of Codon -> Amino Acid for normalization grouping
    GENETIC_CODE = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TGC':'C', 'TGT':'C', 'TGG':'W',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'TAA':'STOP', 'TAG':'STOP', 'TGA':'STOP',
    }

    for org in CODON_TABLE:
        # Group by AA
        aa_groups = {}
        for codon, freq in CODON_TABLE[org].items():
            aa = GENETIC_CODE.get(codon, 'STOP')
            if aa not in aa_groups: aa_groups[aa] = []
            aa_groups[aa].append(codon)
        
        # Normalize each group
        for aa, codons in aa_groups.items():
            total = sum(CODON_TABLE[org][c] for c in codons)
            if total > 0:
                for c in codons:
                    CODON_TABLE[org][c] /= total

_ENV_LOGGED = False  # print environment banners once (main process only)


def setup_directories(quiet=None):
    """Syncs local DIRS with Constants.py and ensures folders exist.

    quiet: if None, auto-quiet in worker processes (stops Phase-2 pool spam).
    """
    global DIRS, _ENV_LOGGED

    # 1. Pull the Single Source of Truth from Constants.py
    DIRS = Constants.DIRS

    # 2. CREATE FOLDERS FIRST: Loop through every path and ensure it exists
    for path in DIRS.values():
        # Avoid creating folders named after the final .xlsx report
        if not path.endswith('.xlsx'):
            os.makedirs(path, exist_ok=True)

    # Master audit trail lives in the Global_Evidence folder so it shows up
    # in the UI's "Live Audit Logs" selector (Constants.DIRS['OUT_LOGS']).
    log_file = os.path.join(Constants.DIRS["OUT_LOGS"], "forensic_engine.log")

    # FileHandler only in the main process — workers would fight over the same log
    is_main = mp.current_process().name == "MainProcess"
    if is_main and not any(isinstance(h, logging.FileHandler) for h in logger.handlers):
        fh = logging.FileHandler(log_file, mode='a')
        fh.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
        logger.addHandler(fh)

    if quiet is None:
        quiet = not is_main
    if not quiet and not _ENV_LOGGED:
        print(f"[*] Forensic Environment Sync Complete: {len(DIRS)} paths mapped.")
        _ENV_LOGGED = True

def translate_dna(dna_seq):
    """
    Translates a DNA target string into its Protein signature.
    Used to generate the Protein Automaton from the DNA Database.
    """
    dna_seq = str(dna_seq).upper().replace(" ", "")
    # Standard Genetic Code
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein = []
    
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        aa = codon_table.get(codon, 'X')
        protein.append(aa)                         
        
    return "".join(protein)
    
def env_needs_bootstrap():
    """
    True when input dirs are missing, or no FASTA exists in CLN/IND/REF,
    or core output path is missing. Does not re-download if data is already present.
    """
    for key in ("CLN", "IND", "REF"):
        d = Constants.DIRS.get(key)
        if not d or not os.path.isdir(d):
            return True
    has_fasta = False
    for key in ("CLN", "IND", "REF"):
        d = Constants.DIRS.get(key)
        if not d or not os.path.isdir(d):
            continue
        try:
            for f in os.listdir(d):
                if f.lower().endswith((".fa", ".fasta", ".fna", ".fas")):
                    has_fasta = True
                    break
        except OSError:
            return True
        if has_fasta:
            break
    if not has_fasta:
        return True
    logs = Constants.DIRS.get("OUT_LOGS")
    if not logs or not os.path.isdir(logs):
        return True
    return False


def ensure_environment():
    """
    Creates the output tree + durable Audit_Library (sibling of Forensic_Results).
    If input folders or FASTA are missing, runs Setup_env automatically
    (TEST_MODE calibration files, or production DL sync).

    Worker processes (Phase 1 / Phase 2 pools) only ensure folders exist — no
    console banners and no Setup_env / DL side effects (that was the spam during
    [CROSS-AUDIT] on Windows spawn).
    """
    _al = getattr(Constants, "AUDIT_LIBRARY", None)
    if _al:
        os.makedirs(_al, exist_ok=True)
    global _ENV_LOGGED
    is_main = mp.current_process().name == "MainProcess"
    setup_directories(quiet=not is_main)

    if not is_main:
        return

    if _ENV_LOGGED:
        return  # already checked this process

    if env_needs_bootstrap():
        print("[*] Environment incomplete (missing folders or FASTA) — running Setup_env...")
        try:
            import Setup_env
            Setup_env.setup_forensic_environment()
            setup_directories(quiet=False)  # re-map after Setup_env
            print("[+] Setup_env finished.")
        except Exception as e:
            print(f"[!] Auto Setup_env failed: {e}")
            print("[!] Create Clinical_Samples / Industrial_Sources / Reference_Archive with FASTA, or run: python Setup_env.py")
    else:
        print("[+] Environment folders + FASTA present.")
    print("[+] Environment check: PASSED")
    _ENV_LOGGED = True

def get_rev_comp(seq):
    """
    Forensic-grade reverse complement. 
    Correctly maps standard and degenerate (wobble) bases.
    """
    trans = str.maketrans(
        "ATCGURYSWKMBDHVatcguryswkmbdhv",
        "TAGCAYRSWMKVHDBtagcayrswmkvhdb"
    )
    return seq.translate(trans)[::-1]

def disk_guard(min_gb=1.0, max_wait_seconds=300):
    total, used, free = shutil.disk_usage(Constants.BASE_OUT if os.path.isdir(Constants.BASE_OUT) else ".")
    free_gb = free / (1024**3)
    if free_gb < min_gb:
        logger.error(f"[!] DISK SPACE CRITICAL: {free_gb:.2f}GB free. Pausing engine...")
        waited = 0
        while free_gb < min_gb:
            if waited >= max_wait_seconds:
                logger.error(f"[!] Disk still critical after {max_wait_seconds}s. Proceeding anyway to avoid an indefinite hang.")
                return
            time.sleep(30)
            waited += 30
            total, used, free = shutil.disk_usage(Constants.BASE_OUT if os.path.isdir(Constants.BASE_OUT) else ".")
            free_gb = free / (1024**3)
            logger.info(f"[*] Disk guard waiting... {free_gb:.2f}GB free (need {min_gb}GB)")

def memory_guard(threshold=90, max_wait_seconds=120):
    """Prevents OOM (Out of Memory) crashes during large sequence mmapping.
    Note: psutil reads host-level memory, which can read high inside a
    container regardless of this process's actual usage — so this guard
    logs while waiting and gives up after max_wait_seconds instead of
    hanging forever with no visible error."""
    waited = 0
    logged = False
    while psutil.virtual_memory().percent > threshold:
        if not logged:
            logger.info(f"[*] Memory guard: usage above {threshold}%, pausing scan...")
            logged = True
        if waited >= max_wait_seconds:
            logger.error(f"[!] Memory guard: still above {threshold}% after {max_wait_seconds}s. Proceeding anyway to avoid an indefinite hang.")
            return
        time.sleep(5)
        waited += 5

def analyze_author_signature(sequence):
    """
    Full RSCU + chi-squared codon-bias audit against CODON_TABLE organisms.

    1. Count every sense codon in the sequence (frame 0 only).
    2. Compute local RSCU = observed / expected-equal-usage per amino acid.
    3. Chi-squared distance against each organism's normalised frequency table.
    4. Return best host + industrial flag when E. coli / Pichia / B. subtilis win.
    """
    sequence = str(sequence).upper().replace("U", "T")
    if len(sequence) < 90:
        return "Inconclusive (sequence too short for reliable RSCU)"

    # Sense codon -> amino acid (STOP excluded from RSCU)
    GENETIC_CODE = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TGC': 'C', 'TGT': 'C',
        'TGG': 'W',
        # STOP omitted on purpose
    }

    # --- 1. Count codons (frame 0) ---
    counts = {}
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        if codon in GENETIC_CODE:
            counts[codon] = counts.get(codon, 0) + 1

    total_codons = sum(counts.values())
    if total_codons < 30:
        return "Inconclusive (coding density too low for signature)"

    # --- 2. Group by amino acid, compute local RSCU ---
    aa_groups = {}
    for codon, aa in GENETIC_CODE.items():
        aa_groups.setdefault(aa, []).append(codon)

    local_rscu = {}
    for aa, codons in aa_groups.items():
        n = len(codons)
        aa_total = sum(counts.get(c, 0) for c in codons)
        if aa_total == 0:
            for c in codons:
                local_rscu[c] = 0.0
            continue
        expected = aa_total / n
        for c in codons:
            # RSCU = observed / expected-under-equal-usage
            local_rscu[c] = (counts.get(c, 0) / expected) if expected > 0 else 0.0

    # --- 3. Chi-squared distance to each organism profile ---
    # Convert organism frequencies to expected counts under the same aa totals,
    # then chi2 = sum( (O - E)^2 / E ). Lower is better.
    scores = {}
    for organism, profile in CODON_TABLE.items():
        chi2 = 0.0
        df_terms = 0
        for aa, codons in aa_groups.items():
            aa_total = sum(counts.get(c, 0) for c in codons)
            if aa_total < 5:
                continue  # skip sparse amino acids
            for c in codons:
                # profile[c] is already a within-AA probability (normalised)
                expected = aa_total * profile.get(c, 1.0 / len(codons))
                if expected < 0.5:
                    expected = 0.5  # floor to avoid division blow-up
                observed = counts.get(c, 0)
                chi2 += (observed - expected) ** 2 / expected
                df_terms += 1
        if df_terms == 0:
            scores[organism] = float("inf")
        else:
            scores[organism] = chi2 / df_terms  # normalised chi2

    if not scores or all(v == float("inf") for v in scores.values()):
        return "Inconclusive (insufficient synonymous codon coverage)"

    best_match = min(scores, key=scores.get)
    best_chi = scores[best_match]
    # Map chi2 to a 0-1 confidence heuristic (lower chi2 -> higher confidence)
    confidence = 1.0 / (1.0 + best_chi)

    if confidence < 0.35:
        return f"Mixed Signature (Weak match to {best_match}, chi2={best_chi:.2f})"

    verdict = f"HOST_TARGET: {best_match} (RSCU-chi2 conf={confidence:.2f}, chi2={best_chi:.2f})"
    if best_match in ("E_COLI", "P_PASTORIS", "BACILLUS_SUBTILIS"):
        verdict += " [!] SIGNAL: INDUSTRIAL_OPTIMIZATION"
    return verdict

def detect_chimeric_logic(found_markers):
    """Flags suspicious cross-species genetic engineering."""
    has_bacterial = any(m in found_markers for m in ["pUC_Origin", "T7_Promoter"])
    has_mammalian = any(m in found_markers for m in ["SV40_Promoter", "CMV_Promoter", "NLS_Signal"])
    has_fungal = any(m in found_markers for m in ["GAL1_Promoter", "AOX1_Site"])
    
    if has_mammalian and has_fungal:
        return "ALERT: MAMMALIAN-FUNGAL CHIMERA DETECTED"
    if has_mammalian and has_bacterial:
        return "LOG: SHUTTLE_VECTOR_DETECTED (Bacterial/Human Hybrid)"
    return None

def check_circuit_logic(sequence, current_pos, current_name, window=250):
    """
    Scans the neighborhood of a found marker to identify paired genetic circuits.
    Returns: list of nearby markers and the extracted intervening sequence.
    """
    global TARGETS
    proximity_matches = []
    
    # 1. Define boundaries using the flat TARGETS lookup
    marker_data = TARGETS.get(current_name, {})
    marker_len = len(marker_data.get("seq", ""))
    
    if marker_len == 0: return [] # Safety check

    search_start = max(0, current_pos - window)
    search_end = min(len(sequence), current_pos + marker_len + window)
    neighborhood = sequence[search_start:search_end]
    
    # 2. Scan the local neighborhood using flat dictionary iteration
    for other_name, other_data in TARGETS.items():
        # Skip if it's the same marker or a section header (headers don't have 'seq')
        if (other_name == current_name or not isinstance(other_data, dict)
                or 'seq' not in other_data or is_composition_marker(other_name)):
            continue 
        
        other_seq = other_data['seq'].upper()
        if other_seq in neighborhood:
            # Found a neighbor! Calculate relative position
            rel_pos = neighborhood.find(other_seq)
            abs_pos = search_start + rel_pos
            
            # 3. Extract the "payload" sequence between the two markers
            if abs_pos < current_pos:
                # Neighbor is to the LEFT
                intervening = sequence[abs_pos + len(other_seq) : current_pos]
            else:
                # Neighbor is to the RIGHT
                intervening = sequence[current_pos + marker_len : abs_pos]
                
            proximity_matches.append({
                "neighbor": other_name,
                "distance": abs(abs_pos - current_pos),
                "payload": intervening
            })
            
    return proximity_matches

def is_low_complexity(seq, n_threshold=0.3):
    """True for flanks too uninformative to cross-link: N-rich, 2 letters or fewer
    (poly-A, CACACA...), or dominated by one base. Such flanks match everywhere."""
    if not seq:
        return True
    s = seq.upper()
    if s.count("N") / len(s) > n_threshold:
        return True
    if len(set(s)) <= 2:
        return True
    return max(s.count(b) for b in "ACGT") / len(s) > 0.7

def forensic_evidence_min_len():
    """Minimum motif length for evidence TXT / BED / flank export (counts still kept)."""
    try:
        return max(1, int(Constants.CONFIG.get("FORENSIC_EVIDENCE_MIN_LEN", 24)))
    except (TypeError, ValueError):
        return 24


def build_config_stamp():
    """
    Compact fingerprint of the settings that shaped an audit.
    Written into TRIAGE_SUMMARY, portal_snapshot, Audit_Library packages,
    and the Excel Start_Here sheet so runs stay comparable.
    """
    cfg = dict(getattr(Constants, "CONFIG", {}) or {})
    keys_of_interest = [
        "FORENSIC_EVIDENCE_MIN_LEN",
        "COMPOSITION_FAST_ONLY",
        "COMPOSITION_AC_MIN_LEN",
        "COMPOSITION_OMIT_ULTRA_SHORT",
        "MATCH_THRESHOLD",
        "FLANK_AUTO",
        "FLANK_MAX_MISMATCH",
        "FLANK_MIN_IDENTITY",
        "TARGETS_PANEL",
        "MAX_CORES",
        "LARGE_FILE_THRESHOLD_MB",
        "LARGE_FILE_MAX_HITS",
        "ENABLE_HEATMAPS",
        "ENABLE_BIOPHYSICS",
        "ENABLE_SONIFICATION",
        "ENABLE_RF_HEURISTICS",
        "ENABLE_CIRCUIT",
        "ENABLE_CODON_BIAS",
        "WATERFALL_INCLUDE_COMPOSITION",
        "REF_DISTANCE_SIDECAR",
        "TEST_MODE",
    ]
    gates = {k: cfg.get(k) for k in keys_of_interest if k in cfg or True}
    # always include known keys even if missing (default-aware)
    for k in keys_of_interest:
        gates[k] = str(cfg.get(k, ""))
    raw = "|".join(f"{k}={gates[k]}" for k in keys_of_interest)
    digest = hashlib.sha256(raw.encode("utf-8")).hexdigest()[:16]
    targets_path = getattr(Constants, "TARGETS_PATH", "") or cfg.get("TARGETS_PANEL", "Targets.json")
    return {
        "stamp_version": 1,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config_sha256_16": digest,
        "targets_path": str(targets_path),
        "gates": gates,
        "note": (
            "Fingerprint of Configuration.txt gates for this run. "
            "Not a diagnosis; use to compare audits."
        ),
    }


def evidence_tag_allowlist():
    """
    Optional tag filter for evidence TXT / BED / flank export.
    Config: EVIDENCE_TAG_ALLOWLIST=TAG1,TAG2  (empty = all non-composition tags).
    Counts and heatmaps are never filtered by this list.
    """
    raw = str(Constants.CONFIG.get("EVIDENCE_TAG_ALLOWLIST", "") or "").strip()
    if not raw:
        return None
    return {t.strip().upper() for t in raw.replace(";", ",").split(",") if t.strip()}


def is_reportable_forensic_hit(match_seq, marker_name=""):
    """
    Gate evidence / BED / flank-table rows for short or low-complexity motifs.
    Counts remain honest; only narrative outputs are suppressed.
    Human T2T is full of CAT/CAC runs that exact-match 6xHis (18 bp) without
    being purification tags.
    Optional EVIDENCE_TAG_ALLOWLIST further restricts which tags appear in
    evidence/BED (still does not change TARGET_COUNTS).
    """
    m = (match_seq or "").upper().replace("U", "T")
    if not m:
        return False
    if len(m) < forensic_evidence_min_len():
        return False
    if is_low_complexity(m):
        return False
    letters = set(m) - {"N"}
    if letters <= {"C", "A"} or letters <= {"G", "T"}:
        return False
    allow = evidence_tag_allowlist()
    if allow is not None and marker_name:
        info = TARGETS.get(marker_name) or {}
        tag = str(info.get("tag") or "UNKNOWN").upper()
        if tag not in allow:
            return False
    return True


def adjust_incorporation_risk(risk_info, mm_dist, flank_len, foreign_markers):
    """
    Prefer exact, long flanks in incorporation risk scores.
    Exact (≥ MIN_FLANK_LEN) + foreign markers → boost.
    Each mismatch demotes score so fuzzy hits rank lower on the grid.
    """
    risk_info = dict(risk_info or {})
    score = int(risk_info.get("risk_score") or 0)
    try:
        min_flank = int(Constants.CONFIG.get("MIN_FLANK_LEN", 24) or 24)
    except (TypeError, ValueError):
        min_flank = 24
    mm = int(mm_dist or 0)
    fl = int(flank_len or 0)
    quality = "unknown"
    if mm == 0 and fl >= min_flank:
        score = min(100, score + (20 if foreign_markers else 8))
        quality = "exact_long"
    elif mm == 0 and fl > 0:
        score = min(100, score + 5)
        quality = "exact_short"
    elif mm > 0:
        score = max(0, score - min(40, 12 * mm))
        quality = f"mismatch_{mm}"
    risk_info["risk_score"] = score
    risk_info["flank_quality"] = quality
    risk_info["flank_len"] = fl
    risk_info["flank_mismatches"] = mm
    return risk_info

def get_random_hex(): 
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))

def check_and_generate_legend(color_map, filename="FORENSIC_LEGEND.png"):
    if not color_map: return
    os.makedirs(Constants.BASE_OUT, exist_ok=True)
    save_path = os.path.join(Constants.BASE_OUT, filename)
    if os.path.exists(save_path): os.remove(save_path)
    
    items = {k: v for k, v in color_map.items() if v and "---" not in k}
    sorted_items = sorted(items.items())
    
    num_items = len(sorted_items)
    cols = 1 if num_items <= 15 else (2 if num_items <= 30 else 3)
    rows_per_col = math.ceil(num_items / cols)
    
    fig_width = 8 * cols  
    fig_height = max(6, rows_per_col * 0.6)
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    patches = [mpatches.Patch(color=c, label=k) for k, c in sorted_items]
    
    if patches:
        legend = ax.legend(
            handles=patches, 
            loc='center', 
            frameon=True, 
            fontsize=16, 
            ncol=cols,
            facecolor='#F0F0F0',
            edgecolor='black',
            title="GENOMIC MARKER REFERENCE KEY",
            title_fontsize=18
        )
        legend.get_title().set_weight('bold')
        
        ax.axis('off')
        plt.savefig(save_path, bbox_inches='tight', dpi=300, facecolor='white')
    plt.close()


def draw_map(results, sequence_length, save_path, title, precomputed_hist=None):
    """
    Dual-track forensic map optimised for long genomes and dense hit lists.
    Top: binned hit density (convolution-smoothed) — visible on Mb-scale sequences.
    Bottom: colour bars grouped by marker colour (not one artist per hit).

    precomputed_hist: optional 1-D array from large-file full-genome scan so the
    density track reflects ALL forensic hits, not only the stored sample.
    """
    from collections import defaultdict

    offset = float(globals().get("OFFSET", 0.0) or 0.0)
    amplitude = float(globals().get("AMPLITUDE", 1.0) or 1.0)
    L = max(int(sequence_length or 0), 1)
    n_bins = 2500

    fig, (ax_wave, ax_bar) = plt.subplots(
        2, 1, figsize=(15, 8), gridspec_kw={"height_ratios": [2.5, 1.5]}
    )

    positions = np.array([float(h.get("pos", 0) or 0) for h in (results or [])], dtype=float)
    if precomputed_hist is not None and len(precomputed_hist) > 0:
        sig = np.asarray(precomputed_hist, dtype=float)
        if sig.size != n_bins:
            # resample to n_bins
            x_old = np.linspace(0, 1, sig.size)
            x_new = np.linspace(0, 1, n_bins)
            sig = np.interp(x_new, x_old, sig)
        k = np.exp(-0.5 * (np.arange(-25, 26) / 8.0) ** 2)
        k /= k.sum()
        sig = np.convolve(sig, k, mode="same")
        peak = float(sig.max()) or 1.0
        y = offset + amplitude * (sig / peak)
    elif positions.size:
        hist, _edges = np.histogram(positions, bins=n_bins, range=(0, L))
        sig = hist.astype(float)
        # Gaussian smooth (~±25 bins)
        k = np.exp(-0.5 * (np.arange(-25, 26) / 8.0) ** 2)
        k /= k.sum()
        sig = np.convolve(sig, k, mode="same")
        peak = float(sig.max()) or 1.0
        y = offset + amplitude * (sig / peak)
    else:
        y = np.full(n_bins, offset)
    x = np.linspace(0, L, n_bins)

    ax_wave.plot(x, y, color="#4AF626", linewidth=1.2, alpha=0.85)
    ax_wave.fill_between(x, offset, y, color="#4AF626", alpha=0.12)
    # Sparse vertical guides: at most ~200 unique hit positions for readability
    if positions.size:
        step = max(1, positions.size // 200)
        guide = positions[::step]
        colors = [results[i].get("color", "#00FF00") for i in range(0, len(results), step)]
        ax_wave.vlines(guide, offset, offset + amplitude * 0.35, colors=colors, alpha=0.25, linewidth=0.7)

    ax_wave.set_facecolor("#000000")
    ax_wave.set_title(title, color="white", fontsize=14, pad=10)
    ax_wave.set_xlim(0, L)
    ax_wave.set_xticks([])
    ax_wave.set_ylabel("Hit density", color="gray")
    ax_wave.tick_params(colors="gray")

    ax_bar.set_facecolor("#000000")
    width = max(L * 0.002, 10)
    by_col = defaultdict(list)
    for h in results or []:
        by_col[h.get("color", "#00FF00")].append((float(h.get("pos", 0) or 0), width))
    for c, ranges in by_col.items():
        ax_bar.broken_barh(ranges, (0, 5), facecolors=c, alpha=0.9)

    labeled = set()
    label_i = 0
    for hit in results or []:
        name = hit.get("name", "Unknown")
        if name in labeled:
            continue
        labeled.add(name)
        pos = float(hit.get("pos", 0) or 0)
        color = hit.get("color", "#00FF00")
        y_text_pos = 7 + (4 * (label_i % 5))
        ax_bar.text(
            pos, y_text_pos, name,
            color=color, fontsize=7, rotation=45,
            fontweight="bold", va="bottom", clip_on=True,
        )
        label_i += 1

    ax_bar.set_xlim(0, L)
    ax_bar.set_ylim(0, 30)
    ax_bar.set_xlabel("Genomic position (bp)", color="gray")
    ax_bar.set_yticks([])
    ax_bar.tick_params(colors="gray")
    for spine in list(ax_wave.spines.values()) + list(ax_bar.spines.values()):
        spine.set_color("#333333")

    fig.patch.set_facecolor("#000000")
    try:
        fig.tight_layout()
        fig.savefig(save_path, dpi=150, facecolor=fig.get_facecolor(), bbox_inches="tight")
    finally:
        plt.close(fig)


def is_inc_clique_marker(name: str) -> bool:
    return str(name or "") in INC_CLIQUE_MARKERS


def is_composition_marker(name: str) -> bool:
    """
    Short / density motifs: tally for waterfall + JSON sidecar only.
    Explicit role=composition in Targets.json wins; else COMPOSITION_FORCE / length.
    """
    if name in FORENSIC_ALWAYS:
        return False
    info = TARGETS.get(name) or {}
    role = str(info.get("role", "")).strip().lower()
    if role == "composition":
        return True
    if role == "forensic":
        return False
    if name in COMPOSITION_FORCE:
        return True
    seq = str(info.get("seq", ""))
    return 0 < len(seq) <= COMPOSITION_MAX_LEN


def _composition_ac_min_len() -> int:
    """Minimum pattern length allowed into Aho-Corasick when FAST_ONLY=0."""
    try:
        return max(1, int(Constants.CONFIG.get("COMPOSITION_AC_MIN_LEN", COMPOSITION_AC_MIN_LEN_DEFAULT)))
    except Exception:
        return COMPOSITION_AC_MIN_LEN_DEFAULT


def composition_ac_min_len() -> int:
    """Public alias used by large-file scanner and callers."""
    return _composition_ac_min_len()


def marker_seq_len(name: str) -> int:
    info = TARGETS.get(name) or {}
    return len(str(info.get("seq", "") or "").strip())


def is_ultra_short_marker(name: str) -> bool:
    """len < COMPOSITION_AC_MIN_LEN — always fast-path only, never AC."""
    n = marker_seq_len(name)
    return 0 < n < _composition_ac_min_len()


def should_omit_from_composition_counts(name: str) -> bool:
    """
    COMPOSITION_OMIT_ULTRA_SHORT=1 drops ultra-short motifs from density counts
    entirely (not just AC). Default 0 keeps exact CG/NGG counts via fast path.
    """
    if str(Constants.CONFIG.get("COMPOSITION_OMIT_ULTRA_SHORT", "0")) != "1":
        return False
    return is_ultra_short_marker(name)


def should_skip_automaton(name: str) -> bool:
    """
    Decide whether a marker is loaded into Aho-Corasick.

    Rules (in order):
    1. Ultra-short (len < COMPOSITION_AC_MIN_LEN, default 6): ALWAYS skip AC.
       These (CG, NGG, TTTV, CGCG, …) explode on genome-scale DNA.
    2. COMPOSITION_FAST_ONLY=1: skip ALL composition markers.
    3. COMPOSITION_FAST_ONLY=0: composition markers with len >= min may enter AC
       (6-cutters, scars, Kozak, …) for denser heatmaps; ultra-short still skip.
    4. Forensic markers always enter AC (unless ultra-short, which is rare).
    """
    if is_ultra_short_marker(name):
        return True
    if not is_composition_marker(name):
        return False
    return str(Constants.CONFIG.get("COMPOSITION_FAST_ONLY", "1")) == "1"

def expand_iupac_patterns(pat: str):
    """
    Expand a DNA pattern containing IUPAC codes into concrete ATCG strings.
    Used so composition markers like NGG / TTTV are counted on real genomes.
    """
    IUPAC = {
        "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
        "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
        "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
    }
    pat = str(pat).upper().replace("U", "T")
    variants = [""]
    for ch in pat:
        bases = IUPAC.get(ch)
        if not bases:
            return []  # unknown symbol
        variants = [v + b for v in variants for b in bases]
    return variants


def count_overlapping(haystack: str, needle: str) -> int:
    """Overlapping occurrences (important for CG, NGG, etc.).

    Short needles are the composition path's hot loop — keep this in pure str.find
    (C implementation) and avoid building intermediate Python objects.
    """
    if not needle or not haystack:
        return 0
    nlen = len(needle)
    # Non-overlapping-safe fast path only when needle cannot self-overlap
    # (e.g. "AT"); composition motifs like "CG"/"GG" must stay step=1.
    if nlen == 1:
        return haystack.count(needle)
    count = 0
    start = 0
    # str.find is implemented in C; still O(n) per pattern and dominates Phase 1
    # when short motifs were wrongly fed through Aho-Corasick + Python hit objects.
    while True:
        i = haystack.find(needle, start)
        if i < 0:
            break
        count += 1
        start = i + 1
    return count

def sequence_metrics(seq: str) -> dict:
    """
    Standard composition metrics used in genomics (not synthetic-origin claims).

    - gc_fraction: (G+C) / (A+T+G+C) ignoring N
    - cpg_observed_expected: Obs(CG) / (C_frac * G_frac * L)  [Gardiner-Garden style]
    - shannon_entropy_bits: base entropy over ACGT (max 2.0 for equal bases)
    """
    s = str(seq).upper().replace("U", "T")
    if not s:
        return {
            "length": 0,
            "gc_fraction": None,
            "cpg_oe": None,
            "shannon_entropy_bits": None,
        }

    counts = Counter(s)
    atgc = sum(counts[b] for b in "ATGC")
    gc = counts["G"] + counts["C"]
    gc_frac = (gc / atgc) if atgc else None

    # CpG O/E
    cpg = 0
    for i in range(len(s) - 1):
        if s[i] == "C" and s[i + 1] == "G":
            cpg += 1
    if atgc and counts["C"] and counts["G"]:
        expected = (counts["C"] / atgc) * (counts["G"] / atgc) * max(len(s) - 1, 1)
        cpg_oe = (cpg / expected) if expected > 0 else None
    else:
        cpg_oe = None

    # Shannon over ACGT only
    if atgc:
        h = 0.0
        for b in "ATGC":
            p = counts[b] / atgc
            if p > 0:
                h -= p * math.log2(p)
    else:
        h = None

    return {
        "length": len(s),
        "gc_fraction": round(gc_frac, 6) if gc_frac is not None else None,
        "cpg_oe": round(cpg_oe, 6) if cpg_oe is not None else None,
        "shannon_entropy_bits": round(h, 6) if h is not None else None,
        "n_count": counts.get("N", 0),
    }

def count_composition_fast(sequence, mol_label="DNA"):
    """
    Fast density counts for composition markers only (no Aho-Corasick).
    FWD = pattern in sequence; REV = pattern in reverse complement (DNA/RNA).
    IUPAC patterns (NGG, TTTV, ...) are expanded to concrete ATCG variants
    and summed (overlapping counts).
    """
    if not sequence or not TARGETS:
        return {}

    seq = str(sequence).upper()
    if mol_label in ("DNA", "RNA"):
        seq = seq.replace("U", "T")
        rc = get_rev_comp(seq)
    else:
        rc = None

    out = {}
    for name, info in TARGETS.items():
        if not is_composition_marker(name):
            continue
        if should_omit_from_composition_counts(name):
            continue
        if not isinstance(info, dict):
            continue
        pat = str(info.get("seq", "")).upper().strip()
        if not pat:
            continue

        if mol_label in ("DNA", "RNA"):
            variants = expand_iupac_patterns(pat)
            if not variants:
                continue
            # Deduplicate identical expansions
            variants = list(dict.fromkeys(variants))
            fwd = sum(count_overlapping(seq, v) for v in variants)
            # Palindromes equal their RC: counting both FWD+REV double-counts physical sites
            rev = 0
            if rc is not None:
                for v in variants:
                    if get_rev_comp(v) == v:
                        continue  # already counted on FWD
                    rev += count_overlapping(rc, v)
        else:
            # Protein: no IUPAC expansion here
            if any(c not in "ACDEFGHIKLMNPQRSTVWY*" for c in pat):
                continue
            fwd = count_overlapping(seq, pat)
            rev = 0

        if fwd or rev:
            out[name] = {"FWD": fwd, "REV": rev}
    return out

def composition_hit_index(sequence, mol_label="DNA", max_positions=None):
    """
    Optional coordinate lists for composition markers (STORE_COMPOSITION_POSITIONS=1).

    Counts always stay exact in TARGET_COUNTS. Positions are optional and capped
    per marker so a human-scale CG map cannot expand into hundreds of millions of
    ints and stall the worker. Cap is configurable via COMPOSITION_POS_CAP
    (default 100_000 per strand).
    """
    if not sequence or not TARGETS:
        return {}
    if str(Constants.CONFIG.get("STORE_COMPOSITION_POSITIONS", "0")) != "1":
        return {}
    if max_positions is None:
        try:
            max_positions = int(Constants.CONFIG.get("COMPOSITION_POS_CAP", 100000))
        except (TypeError, ValueError):
            max_positions = 100000
    max_positions = max(0, max_positions)

    seq = str(sequence).upper()
    if mol_label in ("DNA", "RNA"):
        seq = seq.replace("U", "T")
        rc = get_rev_comp(seq)
    else:
        rc = None
    index = {}

    def _collect(hay, variants, limit):
        pos = []
        if limit <= 0:
            return pos
        for v in variants:
            if not v:
                continue
            start = 0
            while len(pos) < limit:
                i = hay.find(v, start)
                if i < 0:
                    break
                pos.append(i)
                start = i + 1
            if len(pos) >= limit:
                break
        return pos

    for name, info in TARGETS.items():
        if not is_composition_marker(name) or not isinstance(info, dict):
            continue
        if should_omit_from_composition_counts(name):
            continue
        pat = str(info.get("seq", "")).upper().strip()
        if not pat:
            continue
        entry = {"FWD": [], "REV": []}
        if mol_label in ("DNA", "RNA"):
            variants = list(dict.fromkeys(expand_iupac_patterns(pat) or []))
            entry["FWD"] = _collect(seq, variants, max_positions)
            if rc is not None:
                entry["REV"] = _collect(rc, variants, max_positions)
        else:
            entry["FWD"] = _collect(seq, [pat], max_positions)
        if entry["FWD"] or entry["REV"]:
            index[name] = entry
    return index    

def write_composition_sidecar(base_name, category_key, seq_len, composition_counts, out_root=None, metrics_seq=None):
    """
    Writes a clean total marker summary report for each file scanned.
    Outputs:
      Forensic_Results/Global_Evidence/composition/{base}.TARGET_COUNTS.txt
      Forensic_Results/Global_Evidence/composition/{base}.TARGET_COUNTS.json
    """
    import json
    from collections import defaultdict

    # Always write a sidecar so clean files still appear in waterfall/portal
    composition_counts = composition_counts or {}

    if out_root is None:
        out_root = os.path.join(Constants.DIRS.get("OUT_LOGS", Constants.BASE_OUT), "composition")
    os.makedirs(out_root, exist_ok=True)

    by_marker = {}
    by_tag = defaultdict(int)
    density = {}
    seq_len = max(int(seq_len or 0), 1)

    for mname, val in composition_counts.items():
        if isinstance(val, dict):
            fwd = int(val.get("FWD", 0))
            rev = int(val.get("REV", 0))
            total = fwd + rev
        else:
            fwd = rev = 0
            total = int(val)
        tag = (TARGETS.get(mname) or {}).get("tag", "UNKNOWN")
        by_marker[mname] = {
            "FWD": fwd, "REV": rev, "total": total,
            "tag": tag,
            "composition": bool(is_composition_marker(mname)),
        }
        by_tag[tag] += total
        density[mname] = round((total / seq_len) * 1000.0, 4)

    payload = {
        "file": base_name,
        "category": category_key,
        "seq_len": seq_len,
        "by_marker": by_marker,
        "by_tag": dict(by_tag),
        "density_per_kb": density,
    }
    if metrics_seq:
        try:
            sm = sequence_metrics(metrics_seq)
            if int(seq_len or 0) > len(metrics_seq) * 2:
                sm["length"] = int(seq_len)
                sm["metrics_window_bp"] = len(metrics_seq)
                sm["metrics_sampled"] = True
            payload["sequence_metrics"] = sm
        except Exception:
            pass

    safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in base_name)
    
    # 1. Machine-readable JSON summary
    json_path = os.path.join(out_root, f"{safe}.TARGET_COUNTS.json")
    with open(json_path, "w") as jf:
        json.dump(payload, jf, indent=2)

    # 2. Simple human-readable text file: "Target_Name - Count"
    txt_path = os.path.join(out_root, f"{safe}.TARGET_COUNTS.txt")
    with open(txt_path, "w") as tf:
        tf.write(f"=== TOTAL TARGET COUNTS: {base_name} (Length: {seq_len} bp) ===\n\n")
        for mname, stats in sorted(by_marker.items()):
            tf.write(f"{mname} - {stats['total']}\n")
        tf.write("\n" + "=" * 60 + "\n")

    return json_path

def write_hit_index(base_name, seq_len, hit_index, out_root=None):
    """
    Full coordinate list for composition / short markers.
    Forensic_Results/Global_Evidence/hits/{base}.hits.json
    """
    if not hit_index:
        return None
    if out_root is None:
        out_root = os.path.join(Constants.DIRS.get("OUT_LOGS", Constants.BASE_OUT), "hits")
    os.makedirs(out_root, exist_ok=True)
    safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in base_name)
    markers_out = {}
    for mname, strands in hit_index.items():
        info = TARGETS.get(mname) or {}
        plen = len(str(info.get("seq", "") or ""))
        entry = {"len": plen, "FWD": {"count": 0, "positions": []}, "REV": {"count": 0, "positions": []}}
        for sk in ("FWD", "REV"):
            pos_list = sorted(set(int(p) for p in (strands.get(sk) or [])))
            entry[sk] = {"count": len(pos_list), "positions": pos_list}
        markers_out[mname] = entry
    payload = {
        "file": base_name,
        "seq_len": int(seq_len or 0),
        "markers": markers_out,
        "note": "All unique start positions for short/composition markers. Long forensic hits stay in Evidence TXT.",
    }
    path = os.path.join(out_root, f"{safe}.hits.json")
    with open(path, "w", encoding="utf-8") as jf:
        json.dump(payload, jf, separators=(",", ":"))  # compact — not indent=2
    return path    

def initialize_globals():
    """
    Unified Global Loader (single source of truth):
    1. Sets Config
    2. Loads Targets.json
    3. Builds the DNA automaton (forward + reverse-complement, dual-strand)
    4. Builds the Protein automaton (forward + reverse-strand translation)
    5. Calculates Circular Bridge Size
    """
    global TARGETS, MARKER_COLORS, AUTO, AUTO_DNA, AUTO_PROT, DYNAMIC_WINDOW, cfg

    # 1. Configuration Sync
    cfg = Constants.CONFIG
    is_test_mode = str(cfg.get('TEST_MODE', '0')) == '1'

    normalize_forensic_tables()
    # Verbose banner only in the main process (workers would print this 10×)
    if mp.current_process().name == "MainProcess":
        logger.info("[+] Forensic codon tables normalized and validated.")

    # 2. Setup Scanning Brains (active panel; master Targets.json never deleted)
    if hasattr(Constants, "resolve_targets_path"):
        Constants.CONFIG = Constants.load_config()
        manifest_path = Constants.resolve_targets_path(Constants.CONFIG)
        Constants.TARGETS_PATH = manifest_path
    else:
        manifest_path = Constants.TARGETS_PATH
    if mp.current_process().name == "MainProcess":
        logger.info(f"[*] Active targets panel: {os.path.basename(manifest_path)}")
    if not os.path.exists(manifest_path):
        logger.error(f"[!] Critical: {manifest_path} not found.")
        return False

    TARGETS = {}
    MARKER_COLORS = {}
    dna_automaton = ahocorasick.Automaton()
    prot_automaton = ahocorasick.Automaton()
    max_len = 0  # Track for Circular Bridge
    count_dna = 0
    count_prot = 0

    try:
        with open(manifest_path, 'r') as f:
            data = json.load(f)
        markers = {**data.get("FORENSIC_MARKERS", {}), **FEED_MARKERS}

        for name, info in markers.items():
            # Skip section-header entries and malformed data
            if "---" in name or not isinstance(info, dict) or "seq" not in info:
                continue

            # Production safety gate: skip test-only markers unless TEST_MODE=1
            if not is_test_mode and (info.get("is_test") == "1" or name in FORENSIC_ALWAYS):
                continue

            seq = info["seq"].upper().strip()
            if not seq:
                continue

            if len(seq) > max_len:
                max_len = len(seq)

            # Populate Global Metadata (every Targets.json marker is listed)
            TARGETS[name] = info
            MARKER_COLORS[name] = info.get("color", "#808080")

            # Ultra-short (len < COMPOSITION_AC_MIN_LEN) always skip AC.
            # COMPOSITION_FAST_ONLY=1 skips all composition; =0 allows composition
            # markers with len >= min (e.g. 6bp scars) into AC for denser maps.
            # Counts for skipped markers still come from count_composition_fast.
            if should_skip_automaton(name):
                continue

            # --- DNA: Forward + Reverse-Complement (expand IUPAC e.g. NGG / TTTV) ---
            variants = expand_iupac_patterns(seq) or [seq]
            for v in variants:
                if not v:
                    continue
                dna_automaton.add_word(v, (name, "FWD"))
                rc_v = get_rev_comp(v)
                if rc_v != v:
                    dna_automaton.add_word(rc_v, (name, "REV"))
            count_dna += 1

            # --- Protein: translate both strands (never reverse the string itself) ---
            # Peptide scan: require coding-ish length and no stop (*); MIN_PEPTIDE_LEN AA
            try:
                min_pep = int(Constants.CONFIG.get("MIN_PEPTIDE_LEN", 8) or 8)
            except Exception:
                min_pep = 8
            coding_flag = info.get("coding")
            if coding_flag is False:
                pass  # non-coding markers (promoters/oris) never become peptides
            elif len(seq) >= max(9, min_pep * 3) and all(ch in "ATCGU" for ch in seq):
                rc_dna = get_rev_comp(seq)
                prot_fwd = translate_dna(seq)
                if prot_fwd and "*" not in prot_fwd and len(prot_fwd) >= min_pep:
                    prot_automaton.add_word(prot_fwd, (name, "FWD"))
                    count_prot += 1
                    prot_rev = translate_dna(rc_dna)
                    if (
                        prot_rev
                        and prot_rev != prot_fwd
                        and "*" not in prot_rev
                        and len(prot_rev) >= min_pep
                    ):
                        prot_automaton.add_word(prot_rev, (name, "REV"))

        # 3. Finalize Circular Bridge Mathematical Window
        # Bridge must be (Longest Marker - 1) to catch every possible seam split
        DYNAMIC_WINDOW = max_len - 1 if max_len > 0 else 50

        # 4. Seal both Aho-Corasick Search Trees
        dna_automaton.make_automaton()
        prot_automaton.make_automaton()
        AUTO_DNA = dna_automaton
        AUTO_PROT = prot_automaton
        AUTO = AUTO_DNA  # legacy alias: Phase 2 neighborhood scans use the DNA brain

        if mp.current_process().name == "MainProcess":
            logger.info(f"[*] Circular Engine: Dynamic Bridge set to {DYNAMIC_WINDOW}bp")
            logger.info(
                f"[*] Dual-Strand Brain: {count_dna} DNA / {count_prot} Protein "
                f"signatures loaded ({len(TARGETS)} markers)."
            )
        return len(TARGETS) > 0

    except Exception as e:
        logger.error(f"[!] Critical Manifest Error: {e}")
        return False

def sync_global_threats():
    """Optional remote marker feeds: JSON of {name: sequence}. HTTPS only, sequences validated.
    Feed markers are kept in FEED_MARKERS so initialize_globals() merges them instead of wiping them."""
    urls = [u.strip() for u in Constants.CONFIG.get("THREAT_FEED_URLS", "").split(",") if u.strip()]
    if not urls:
        print("[!] No threat feeds configured. Skipping sync.")
        return

    added = 0
    for url in urls:
        if not url.lower().startswith("https://"):
            print(f" [!] Refusing non-HTTPS feed: {url}")
            continue
        print(f"[*] Connecting to: {url}...")
        try:
            response = requests.get(url, timeout=10)
            if response.status_code != 200:
                print(f" [!] Feed Error {response.status_code}: {url}")
                continue
            for name, seq in response.json().items():
                seq = str(seq).upper().strip()
                if re.fullmatch(r"[ACGTN]{12,}", seq):
                    FEED_MARKERS[str(name)] = {"seq": seq, "tag": "FEED_MARKER", "color": "#808080"}
                    added += 1
                    print(f" [+] Synced Marker: {name}")
                else:
                    print(f" [!] Skipped invalid feed entry: {name}")
        except Exception as e:
            print(f" [!] Connection Failed: {url} | {e}")

    if added:
        logger.warning(f"[!] {added} marker(s) came from remote feeds and are NOT part of the audited Targets.json")
        initialize_globals()

def extract_clean_sequence(mm):
    """Whole-file cleaner: drops every FASTA header line and any non-sequence byte."""
    mm.seek(0)
    data = _HEADER_RE.sub(b"", mm.read())
    return data.translate(None, _DELETE_BYTES).upper().decode("ascii", errors="ignore")

def is_circular_vector(filename):
    circular_exts = {'.plasmid', '.vector', '.gbk', '.dna'}
    # Fix: Added 'circular' to search terms so "Circular_Validation.fasta" triggers the loop logic
    return any(filename.lower().endswith(ext) for ext in circular_exts) or \
           "vector" in filename.lower() or "circular" in filename.lower()

# --- Sequence-feature spectrogram (NOT measured RF) ---
# Replaces the old 2.4G/4G/5G ACTIVE labels. Builds a spectrogram-style PNG
# from local composition structure on a log frequency *display* axis.
def check_frequency_activation(sequence, mol_label, out_png=None, marker_positions=None, title=None):
    """
    Sequence-feature spectrogram helper. Optional PNG next to heatmaps.
    Does NOT claim measured radio emission or Wi-Fi/cellular activation.
    """
    try:
        import Sonification as _son
        seq = sequence or ""
        mt = "PROT" if mol_label == "PROT" else mol_label
        result = _son.feature_spectrogram(
            seq,
            mol_type=mt,
            out_png=out_png,
            marker_positions=marker_positions,
            title=title,
        )
        return result
    except Exception as e:
        return {
            "peaks": [],
            "path": None,
            "disclaimer": f"Spectrogram unavailable: {e}",
        }


def write_feature_spectrogram(
    sequence,
    base_name,
    mol_label,
    dir_prefix=None,
    results=None,
    flanks=None,
    out_dir=None,
):
    """
    Write sequence-feature spectrogram PNGs (NOT measured RF).

    Layout under {Clinical|Industrial|Reference}_Spectrograms/:
      FWD/   — always when ENABLE_RF_HEURISTICS=1 (sense / file sequence)
      REV/   — when SPECTROGRAM_RC=1 (reverse-complement composition)
      FLANKS/— when SPECTROGRAM_FLANKS=1 and flank oligos are available

    Filenames: {base}_FEATURE_SPECTROGRAM.png
               {base}_FEATURE_SPECTROGRAM_RC.png
               {base}_FEATURE_SPECTROGRAM_FLANKS.png
    """
    if str(Constants.CONFIG.get("ENABLE_RF_HEURISTICS", "0")) != "1":
        return None
    if not sequence and not flanks:
        return None

    def _dir(suffix):
        if out_dir and suffix == "FWD":
            return out_dir
        if dir_prefix:
            key = f"{dir_prefix}_SPECTROGRAM_{suffix}"
            d = Constants.DIRS.get(key)
            if d:
                return d
            root = Constants.DIRS.get(f"{dir_prefix}_SPECTROGRAM")
            if root:
                return os.path.join(root, suffix)
        root = Constants.DIRS.get("OUT_MAPS") or "."
        return os.path.join(root, "Spectrograms", suffix)

    def _markers_from_results(res, invert_pos=False, seq_len=0):
        markers = []
        for h in res or []:
            try:
                pos = float(h.get("pos", 0) or 0)
            except (TypeError, ValueError):
                continue
            if invert_pos and seq_len > 0:
                # map sense coordinate onto RC axis (approximate)
                ml = len(str(h.get("match") or h.get("target") or "") or "N")
                pos = max(0.0, seq_len - pos - ml)
            markers.append({
                "pos": pos,
                "name": h.get("name", ""),
                "color": h.get("color", "#22d3ee"),
            })
        if len(markers) > 4000:
            step = max(1, len(markers) // 4000)
            markers = markers[::step]
        return markers

    def _one(seq, folder, fname, title, markers):
        if not seq:
            return None
        try:
            os.makedirs(folder, exist_ok=True)
            path = os.path.join(folder, fname)
            res = check_frequency_activation(
                seq, mol_label, out_png=path, marker_positions=markers, title=title
            )
            if res.get("path"):
                logger.info(f"[+] Feature spectrogram: {res['path']}")
            return res
        except Exception as e:
            logger.error(f"[!] Feature spectrogram failed ({fname}): {e}")
            return None

    written = {}
    seq = sequence or ""
    L = len(seq)

    # --- FWD (sense) ---
    written["FWD"] = _one(
        seq,
        _dir("FWD"),
        f"{base_name}_FEATURE_SPECTROGRAM.png",
        f"Feature spectrogram FWD: {base_name} (NOT RF measurement)",
        _markers_from_results(results, invert_pos=False, seq_len=L),
    )

    # --- REV (RC composition) ---
    if str(Constants.CONFIG.get("SPECTROGRAM_RC", "1")) == "1" and seq and mol_label != "PROT":
        try:
            rc = get_rev_comp(seq)
        except Exception:
            rc = seq[::-1]  # crude fallback
        written["REV"] = _one(
            rc,
            _dir("REV"),
            f"{base_name}_FEATURE_SPECTROGRAM_RC.png",
            f"Feature spectrogram RC: {base_name} (NOT RF measurement)",
            _markers_from_results(results, invert_pos=True, seq_len=L),
        )

    # --- FLANKS strip (concatenated harvested oligos) ---
    if str(Constants.CONFIG.get("SPECTROGRAM_FLANKS", "1")) == "1":
        parts = []
        seen = set()
        for flk in flanks or []:
            for side in ("left", "right", "match"):
                seg = (flk.get(side) or "").upper().replace("U", "T")
                seg = "".join(c for c in seg if c in "ACGTN")
                if len(seg) >= 8 and seg not in seen:
                    seen.add(seg)
                    parts.append(seg)
        if not parts and results:
            # fall back to hit match windows if flanks not passed
            for h in results:
                seg = (h.get("match") or "").upper()
                if len(seg) >= 8 and seg not in seen:
                    seen.add(seg)
                    parts.append(seg)
        flank_seq = "".join(parts)
        # Cap absurd concatenations
        if len(flank_seq) > 2_000_000:
            flank_seq = flank_seq[:2_000_000]
        if flank_seq:
            written["FLANKS"] = _one(
                flank_seq,
                _dir("FLANKS"),
                f"{base_name}_FEATURE_SPECTROGRAM_FLANKS.png",
                f"Feature spectrogram FLANKS strip: {base_name} (NOT RF measurement)",
                None,
            )

    return written



def stream_clean_sequence_chunks(file_path, chunk_bp=None, overlap_bp=None):
    """
    Yield (absolute_start_offset, cleaned_uppercase_chunk) from FASTA/raw
    without materialising the whole genome. Header lines and non-alphabet
    bytes are stripped on the fly. Chunk ends are moved to a line boundary so
    a header can never be split across two chunks.
    """
    if chunk_bp is None:
        chunk_bp = STREAM_CHUNK_MB * 1024 * 1024
    if overlap_bp is None:
        overlap_bp = STREAM_OVERLAP_BP

    with open(file_path, "rb") as f:
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            pos = 0
            carry = b""
            abs_clean = 0
            while pos < len(mm):
                end = pos + chunk_bp
                if end < len(mm):
                    nl = mm.find(b"\n", end)
                    end = len(mm) if nl == -1 else nl + 1
                raw = _HEADER_RE.sub(b"", mm[pos:end])
                pos = end
                cleaned = (carry + raw).translate(None, _DELETE_BYTES).upper()
                if not cleaned:
                    carry = b""
                    continue
                yield abs_clean, cleaned.decode("ascii", errors="ignore")
                keep = min(overlap_bp, len(cleaned))
                carry = cleaned[len(cleaned) - keep:]
                abs_clean += len(cleaned) - keep

def scan_automaton_chunked(automaton, file_path, mol_label, strands, flank_size=42):
    """
    Aho-Corasick over a large file in overlapping chunks — full-file scan.

    Memory model (T2T-safe, no early abort):
      • Every forensic + long marker match is COUNTED (complete tallies).
      • Composition / ultra-short never become hit objects (density via
        count_composition_fast after the scan).
      • Full hit objects (match + flanks) kept only for a reservoir sample
        of forensic hits (LARGE_FILE_MAX_HITS) for Phase-2 flanks / evidence.
      • Density histograms (fixed bins) accumulate ALL forensic positions so
        heatmaps reflect the whole genome, not only the sample.
      • `seen` is pruned to the overlap window only.

    Returns
      hits       — forensic sample (list of hit dicts, ≤ LARGE_FILE_MAX_HITS)
      total_len  — cleaned sequence length (bp)
      meta       — {marker_counts, hist_fwd, hist_rev, n_bins, total_forensic,
                    total_skipped_composition, sample_size}
    """
    import random as _rnd

    hits = []  # reservoir of forensic hit objects
    n_hits_seen = 0  # forensic matches observed (including those not stored)
    n_comp_skipped = 0
    seen = set()
    marker_counts = {}  # name -> {"FWD": n, "REV": n}

    try:
        max_store = int(Constants.CONFIG.get("LARGE_FILE_MAX_HITS", 250000))
    except (TypeError, ValueError):
        max_store = 250000
    max_store = max(1_000, max_store)

    min_ac_len = composition_ac_min_len()
    max_marker = max(
        (len(v.get("seq", "")) for v in TARGETS.values() if isinstance(v, dict) and v.get("seq")),
        default=64,
    )
    overlap = max(STREAM_OVERLAP_BP, max_marker + 2 * flank_size + 8)
    seen_keep_bp = overlap + max_marker + 8

    # Fixed-bin density for whole-genome heatmaps (constant RAM)
    n_bins = 2500
    hist_fwd = np.zeros(n_bins, dtype=np.float64)
    hist_rev = np.zeros(n_bins, dtype=np.float64)
    # Map physical → bin using running extent; refined as total_len grows
    max_pos_seen = 1

    def _bin(physical):
        nonlocal max_pos_seen
        if physical > max_pos_seen:
            max_pos_seen = physical
        # provisional scale by max_pos; good enough for density shape
        idx = int(physical / max(max_pos_seen, 1) * (n_bins - 1))
        return max(0, min(n_bins - 1, idx))

    def _count(name, strand):
        b = marker_counts.setdefault(name, {"FWD": 0, "REV": 0})
        sk = strand if strand in ("FWD", "REV") else "FWD"
        b[sk] = b.get(sk, 0) + 1

    def _store_hit(hit):
        """Reservoir sample so every region of a huge genome can contribute."""
        nonlocal n_hits_seen
        n_hits_seen += 1
        if len(hits) < max_store:
            hits.append(hit)
            return
        # Vitter-style: replace with decreasing probability
        j = _rnd.randint(1, n_hits_seen)
        if j <= max_store:
            hits[j - 1] = hit

    def process(abs_start, chunk, is_last):
        nonlocal n_comp_skipped
        if mol_label == "RNA":
            chunk = chunk.replace("U", "T")
        for end_idx, (marker_name, match_type) in automaton.iter(chunk):
            info = TARGETS.get(marker_name) or {}
            dna_target = str(info.get("seq", "") or "").upper()
            if not dna_target:
                continue
            # Composition / ultra-short: never materialise; density is fast-path later
            if len(dna_target) < min_ac_len or is_composition_marker(marker_name):
                n_comp_skipped += 1
                continue
            oriented = dna_target if match_type == "FWD" else get_rev_comp(dna_target)
            target_string = translate_dna(oriented) if mol_label == "PROT" else oriented
            match_len = len(target_string)
            start_pos = end_idx - match_len + 1
            if match_len == 0 or start_pos < 0:
                continue
            if start_pos < flank_size and abs_start > 0:
                continue
            if end_idx + 1 + flank_size > len(chunk) and not is_last:
                continue
            physical = abs_start + start_pos
            key = (marker_name, physical, match_type)
            if key in seen:
                continue
            seen.add(key)

            # Always count + density (whole genome)
            _count(marker_name, match_type)
            bi = _bin(physical)
            if match_type == "REV":
                hist_rev[bi] += 1.0
            else:
                hist_fwd[bi] += 1.0

            left = chunk[max(0, start_pos - flank_size):start_pos]
            right = chunk[end_idx + 1:end_idx + 1 + flank_size]
            hit = {
                "name": marker_name,
                "pos": physical,
                "strand": match_type,
                "match": chunk[start_pos:end_idx + 1],
                "target": target_string,
                "left": ("N" * (flank_size - len(left))) + left,
                "right": right + ("N" * (flank_size - len(right))),
                "color": MARKER_COLORS.get(marker_name, "#808080"),
            }
            _store_hit(hit)

        # Prune seen to overlap neighbourhood only
        if len(seen) > 80_000:
            cutoff = abs_start - seen_keep_bp
            if cutoff > 0:
                seen.intersection_update({k for k in seen if k[1] >= cutoff})

    total_len = 0
    prev = None
    n_chunks = 0
    for abs_start, chunk in stream_clean_sequence_chunks(file_path, overlap_bp=overlap):
        if prev is not None:
            process(prev[0], prev[1], False)
            n_chunks += 1
            if n_chunks % 20 == 0:
                logger.info(
                    f"[*] {os.path.basename(file_path)}: chunk {n_chunks}, "
                    f"pos~{prev[0]:,}, forensic_seen={n_hits_seen:,}, "
                    f"stored={len(hits):,}, markers={len(marker_counts)}"
                )
        prev = (abs_start, chunk)
        total_len = abs_start + len(chunk)
    if prev is not None:
        process(prev[0], prev[1], True)
        n_chunks += 1

    # Re-scale histograms to final length coordinates for draw_map consumers
    # (bins already relative to max_pos_seen ≈ total_len)
    meta = {
        "marker_counts": marker_counts,
        "hist_fwd": hist_fwd,
        "hist_rev": hist_rev,
        "n_bins": n_bins,
        "total_forensic": n_hits_seen,
        "total_skipped_composition_ac": n_comp_skipped,
        "sample_size": len(hits),
        "chunks": n_chunks,
    }
    logger.info(
        f"[+] Large scan complete {os.path.basename(file_path)}: "
        f"len={total_len:,} forensic_total={n_hits_seen:,} "
        f"stored_sample={len(hits):,} marker_types={len(marker_counts)} chunks={n_chunks}"
    )
    return hits, total_len, meta



def get_flanks(seq, start, end, size, circular):
    """Left/right flank of a hit at seq[start:end+1]. Circular sequences wrap around the origin;
    linear sequences are padded with N instead of borrowing bases from the opposite end."""
    n = len(seq)
    if circular and n > size:
        l0 = start - size
        left = seq[l0:start] if l0 >= 0 else seq[l0:] + seq[:start]
        r1 = end + 1 + size
        right = seq[end + 1:r1] if r1 <= n else seq[end + 1:] + seq[:r1 - n]
        return left, right
    left = seq[max(0, start - size):start]
    right = seq[end + 1:end + 1 + size]
    return ("N" * (size - len(left))) + left, right + ("N" * (size - len(right)))



def file_scope_id(file_path, dir_prefix=None):
    """Disambiguate same basename in REF/IND/CLN: OUT_IND/foo.fasta."""
    base = os.path.basename(file_path)
    if dir_prefix:
        return f"{dir_prefix}/{base}"
    ap = os.path.abspath(file_path).replace("\\", "/")
    for key, label in (("Industrial_Sources", "OUT_IND"), ("Clinical_Samples", "OUT_CLN"), ("Reference_Archive", "OUT_REF")):
        if f"/{key}/" in ap or ap.endswith("/" + key + "/" + base):
            return f"{label}/{base}"
    return base


def list_input_files():
    """Sequence inputs from the three source folders.

    Accepts FASTA/FASTQ-style extensions only. Annotation dumps (.gff/.gtf/.gff3)
    are skipped — they are not genomes and inflate noise counts.
    """
    _SEQ_EXT = (
        ".fa", ".fasta", ".fna", ".ffn", ".faa", ".frn",
        ".fas", ".seq", ".fastq", ".fq",
    )
    files = []
    for key in ("CLN", "IND", "REF"):
        d = Constants.DIRS.get(key)
        if d and os.path.isdir(d):
            for f in sorted(os.listdir(d)):
                p = os.path.join(d, f)
                if not os.path.isfile(p) or f.startswith("."):
                    continue
                if f.endswith((".part", ".tmp")):
                    continue
                low = f.lower()
                if not any(low.endswith(ext) for ext in _SEQ_EXT):
                    continue
                files.append(p)
    return files

def index_worker(args):
    global TARGETS, MARKER_COLORS, AUTO_DNA, AUTO_PROT, DYNAMIC_WINDOW
    file_path, worker_id, local_dirs = args
    base_name = os.path.basename(file_path)
    is_test_mode = str(Constants.CONFIG.get('TEST_MODE', '0')) == '1'

    # Dictionary to buffer evidence text in memory: key = (file_path), value = list of text strings
    evidence_buffers = {}

    def queue_evidence(out_dir, filename, text):
        if not out_dir:
            return
        target_path = os.path.join(out_dir, filename)
        if target_path not in evidence_buffers:
            evidence_buffers[target_path] = []
        evidence_buffers[target_path].append(text)

    def flush_evidence():
        for target_path, lines in evidence_buffers.items():
            if lines:
                os.makedirs(os.path.dirname(target_path), exist_ok=True)
                with open(target_path, "a", encoding="utf-8") as ef:
                    ef.write("".join(lines))
        evidence_buffers.clear()

    # --- 1. SOURCE CLASSIFICATION & NAMING ---
    _ap = os.path.abspath(file_path)

    def _under(key):
        return _ap.startswith(os.path.abspath(local_dirs[key]) + os.sep)

    if _under("IND"):
        readable_prefix, dir_prefix = "Industrial", "OUT_IND"
    elif _under("CLN"):
        readable_prefix, dir_prefix = "Clinical", "OUT_CLN"
    else:
        readable_prefix, dir_prefix = "Reference", "OUT_REF"
    scope_id = file_scope_id(file_path, dir_prefix)

    # --- 2. FORENSIC READER ---
    size_mb = _file_size_mb(file_path)
    large_mode = size_mb >= LARGE_FILE_THRESHOLD_MB
    file_hash_hex = "UNKNOWN"
    seq_str = ""

    try:
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                sha256_hash.update(chunk)
        file_hash_hex = sha256_hash.hexdigest()

        if large_mode:
            sample = ""
            for abs_start, chunk in stream_clean_sequence_chunks(
                file_path, chunk_bp=2 * 1024 * 1024, overlap_bp=0
            ):
                sample = chunk
                break
            seq_str = sample
            tqdm.write(f"[*] LARGE FILE MODE ({size_mb:.0f} MB): {base_name} — chunked scan + exact flanks")
        else:
            with open(file_path, "rb") as f:
                with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                    seq_str = extract_clean_sequence(mm)
    except Exception as e:
        logger.error(f"[-] Read Error {base_name}: {e}")
        return None

    if not seq_str and not large_mode:
        return None
    if not TARGETS:
        initialize_globals()

    # --- 3. MOLECULAR ROUTING ---
    probe = seq_str[:1_000_000]
    nucleic_hits = sum(probe.count(c) for c in "ACGTUNRYSWKMBDHV")
    non_nucleic_ratio = (len(probe) - nucleic_hits) / max(len(probe), 1)

    harvested_data = {"flanks": [], "composition": {}, "seq_len": 0, "hit_index": {}, "ref_distance": []}

    if non_nucleic_ratio > 0.10:
        mol_label, scan_automaton = "PROT", AUTO_PROT
        sequence = seq_str
        strands_to_scan = ["FWD"]
        tqdm.write(f"[*] Signature Check: PROTEIN identified ({non_nucleic_ratio:.1%} amino acids) for {base_name}")
    else:
        scan_automaton = AUTO_DNA
        # DNA/RNA: one pass. FWD vs REV comes from automaton tags (seq vs rc(seq)).
        strands_to_scan = ["FWD"]
        if "U" in seq_str:
            mol_label, sequence = "RNA", seq_str.replace("U", "T")
            tqdm.write(f"[*] Signature Check: RNA identified for {base_name}")
        else:
            mol_label, sequence = "DNA", seq_str
            tqdm.write(f"[*] Signature Check: DNA identified for {base_name}")

    results = []
    num_cores = int(Constants.CONFIG.get("MAX_CORES", 1))
    seen_matches = set()

    # --- LARGE MODE SCAN ---
    if large_mode:
        logger.info(f"[*] Chunked scan starting for {base_name} ({size_mb:.0f} MB)")
        memory_guard(threshold=90, max_wait_seconds=60)
        f_size = int(Constants.CONFIG.get("FLANK_AUTO", 200))
        large_meta = {}
        try:
            chunk_hits, total_len, large_meta = scan_automaton_chunked(
                scan_automaton, file_path, mol_label, strands_to_scan, f_size
            )
        except MemoryError:
            logger.error(
                f"[!] MemoryError scanning {base_name} — partial/empty hits returned. "
                "Use COMPOSITION_FAST_ONLY=1, COMPOSITION_AC_MIN_LEN>=6, fewer MAX_CORES."
            )
            chunk_hits, total_len, large_meta = [], int(harvested_data.get("seq_len") or 0), {}
            try:
                total_len = max(total_len, int(os.path.getsize(file_path)))
            except OSError:
                pass
        except ValueError:
            # Older 2-tuple return safety
            pair = scan_automaton_chunked(
                scan_automaton, file_path, mol_label, strands_to_scan, f_size
            )
            if isinstance(pair, tuple) and len(pair) == 2:
                chunk_hits, total_len = pair
                large_meta = {}
            else:
                chunk_hits, total_len, large_meta = pair[0], pair[1], pair[2] if len(pair) > 2 else {}
        harvested_data["seq_len"] = total_len
        # Merge complete forensic tallies from full-genome pass (not just sample)
        for _mn, _st in (large_meta.get("marker_counts") or {}).items():
            dest = harvested_data.setdefault("composition", {}).setdefault(
                _mn, {"FWD": 0, "REV": 0}
            )
            dest["FWD"] = dest.get("FWD", 0) + int(_st.get("FWD", 0) or 0)
            dest["REV"] = dest.get("REV", 0) + int(_st.get("REV", 0) or 0)
        harvested_data["large_scan_meta"] = {
            k: large_meta[k]
            for k in ("total_forensic", "sample_size", "chunks", "total_skipped_composition_ac")
            if k in large_meta
        }

        for h in chunk_hits:
            marker_name = h["name"]
            physical_pos = h["pos"]
            strand = h["strand"]
            actual_file_match = h["match"]
            target_string = h["target"]
            match_len = len(target_string)
            if match_len == 0:
                continue

            alignment_bars = "".join(
                [
                    "|" if i < len(actual_file_match) and target_string[i] == actual_file_match[i] else "."
                    for i in range(match_len)
                ]
            )
            match_pct = (alignment_bars.count("|") / match_len) * 100 if match_len else 0.0
            freq_scan = check_frequency_activation(actual_file_match, mol_label)

            left_flank, right_flank = h["left"], h["right"]

            local_window = f"{left_flank}{actual_file_match}{right_flank}"
            local_pos = len(left_flank)

            strand_key = strand if strand in ("FWD", "REV") else "FWD"
            bucket = harvested_data.setdefault("composition", {}).setdefault(
                marker_name, {"FWD": 0, "REV": 0}
            )
            bucket[strand_key] = bucket.get(strand_key, 0) + 1

            if is_composition_marker(marker_name):
                idx = harvested_data.setdefault("hit_index", {}).setdefault(
                    marker_name, {"FWD": [], "REV": []}
                )
                idx.setdefault(strand_key, []).append(int(physical_pos))
                results.append({
                    "pos": physical_pos,
                    "name": marker_name,
                    "strand": strand,
                    "color": MARKER_COLORS.get(marker_name, "#808080"),
                })
                continue

            if not is_reportable_forensic_hit(actual_file_match, marker_name):
                results.append({
                    "pos": physical_pos,
                    "name": marker_name,
                    "strand": strand,
                    "color": MARKER_COLORS.get(marker_name, "#808080"),
                    "suppressed_evidence": True,
                })
                continue

            out_dir = local_dirs.get(f"{dir_prefix}_EVI_{strand}")
            if out_dir:
                filename = f"{readable_prefix}_Evidence_{mol_label}_{strand}_{base_name}.txt"
                tag = TARGETS.get(marker_name, {}).get("tag", "UNTAGGED")
                
                block = []
                block.append(f"\n[!] GENOMIC HIT DETECTED: {marker_name} ({match_pct:.1f}% Match)\n")

                block.append(f"FILE HASH   : {file_hash_hex}\n")

                block.append(f"LOCATION: {physical_pos} | STRAND: {strand} | TYPE: {mol_label}\n")

                block.append(f"THREAT TAG : {tag}\n")

                block.append(f"TARGET:  {target_string}\n")

                block.append(f"         {alignment_bars}\n")

                block.append(f"MATCH:   {actual_file_match}\n")

                block.append(f"LEFT:    {left_flank}\n")

                block.append(f"RIGHT:   {right_flank}\n")

                block.append(f"CONTEXT: ...{left_flank}[{actual_file_match}]{right_flank}...\n")

                try:
                    aa = translate_dna(actual_file_match.replace("U", "T"))
                    block.append(f"PROTEIN TRANSLATION: {aa}\n")
                except Exception:
                    block.append("PROTEIN TRANSLATION: N/A\n")

                if str(Constants.CONFIG.get("ENABLE_CODON_BIAS", "0")) == "1":
                    try:
                        bias = analyze_author_signature(local_window)
                        block.append(f"CODON BIAS SIGNATURE: {bias}\n")
                    except Exception as bias_err:
                        block.append(f"CODON BIAS SIGNATURE: ERROR ({bias_err})\n")

                if str(Constants.CONFIG.get("ENABLE_CIRCUIT", "0")) == "1":
                    try:
                        circuit = check_circuit_logic(local_window, local_pos, marker_name, window=250)
                        if circuit:
                            block.append(f"CIRCUIT NEIGHBORS ({len(circuit)} within 250 bp):\n")
                            for hit in circuit:
                                block.append(f"  -> {hit['neighbor']} dist={hit['distance']}bp")
                                if hit.get("payload"):
                                    block.append(f" payload_len={len(hit['payload'])}")
                                block.append("\n")
                        else:
                            block.append("CIRCUIT NEIGHBORS: none within 250 bp\n")
                    except Exception as cir_err:
                        block.append(f"CIRCUIT NEIGHBORS: ERROR ({cir_err})\n")

                if str(Constants.CONFIG.get("ENABLE_RF_HEURISTICS", "0")) == "1":
                    block.append("--- SEQUENCE-FEATURE SPECTROGRAM (NOT measured RF) ---\n")
                    block.append("  Composition on a log-frequency display axis — not radio detection.\n")
                    for _pk in (freq_scan.get("peaks") or [])[:5]:
                        _bc = _pk.get("band_context") or ""
                        block.append(
                            f"  peak {_pk.get('label')}: strength={_pk.get('strength', 0):.2f}"
                            + (f" | {_bc}" if _bc else "")
                            + "\n"
                        )
                    if freq_scan.get("path"):
                        block.append(f"  spectrogram: {freq_scan['path']}\n")
                    block.append("  Display only — not measured RF emission.\n")
                profile_seq = actual_file_match.replace("T", "U") if mol_label == "RNA" else actual_file_match
                if str(Constants.CONFIG.get("ENABLE_BIOPHYSICS", "0")) == "1":
                    try:
                        bio_profile = BiophysicalProfile.build_profile(profile_seq, mol_type=mol_label)
                        block.append(BiophysicalProfile.format_profile_report(bio_profile) + "\n")
                    except Exception as bp_err:
                        logger.error(f"[!] Biophysical profile failed for {marker_name} in {base_name}: {bp_err}")
                        block.append(f"[!] Biophysical profile failed: {bp_err}\n")
                try:
                    _heur = BiophysicalProfile.format_optional_heuristics(
                        profile_seq, mol_type=mol_label
                    )
                    if _heur:
                        block.append(_heur)
                except Exception as heur_err:
                    block.append(f"[!] Optional heuristics failed: {heur_err}\n")


                block.append("-" * 60 + "\n")
                queue_evidence(out_dir, filename, "".join(block))

            harvested_data["flanks"].append({
                "marker": marker_name, "pos": physical_pos, "mol_type": mol_label,
                "strand": strand, "target": target_string, "match": actual_file_match,
                "left": left_flank, "right": right_flank, "file": scope_id, "basename": base_name, "category": dir_prefix,
            })
            if str(Constants.CONFIG.get("REF_DISTANCE_SIDECAR", "0")) == "1":
                if not is_composition_marker(marker_name):
                    clique_only = str(Constants.CONFIG.get("REF_DISTANCE_CLIQUE_ONLY", "1")) == "1"
                    if (not clique_only) or is_inc_clique_marker(marker_name):
                        rd = ref_distance_record(
                            marker_name, target_string, actual_file_match, physical_pos, strand
                        )
                        if rd:
                            harvested_data.setdefault("ref_distance", []).append(rd)
            if str(Constants.CONFIG.get("REF_DISTANCE_SIDECAR", "0")) == "1":
                if not is_composition_marker(marker_name):
                    clique_only = str(Constants.CONFIG.get("REF_DISTANCE_CLIQUE_ONLY", "1")) == "1"
                    if (not clique_only) or is_inc_clique_marker(marker_name):
                        rd = ref_distance_record(
                            marker_name, target_string, actual_file_match, physical_pos, strand
                        )
                        if rd:
                            harvested_data.setdefault("ref_distance", []).append(rd)
            # Write flank-only line (not full evidence block)
            flanks_dir = local_dirs.get(f"{dir_prefix}_FLANKS")
            if flanks_dir:
                os.makedirs(flanks_dir, exist_ok=True)
                flank_path = os.path.join(flanks_dir, f"{base_name}_flanks.txt")
                with open(flank_path, "a", encoding="utf-8") as ff:
                    ff.write(
                        f"{marker_name}\t{physical_pos}\t{strand}\t"
                        f"{left_flank}\t{actual_file_match}\t{right_flank}\n"
                    )
            # Tallies already merged from large_meta["marker_counts"] (full genome).
            # Do not re-increment here or the reservoir sample would double-count.
            if not (large_meta or {}).get("marker_counts"):
                bucket = harvested_data.setdefault("composition", {}).setdefault(
                    marker_name, {"FWD": 0, "REV": 0}
                )
                sk = strand if strand in ("FWD", "REV") else "FWD"
                bucket[sk] = bucket.get(sk, 0) + 1
            results.append({
                "pos": physical_pos, "name": marker_name, "strand": strand,
                "color": h.get("color", "#808080"),
            })

        # FLUSH BUFFER ONCE AT END OF LARGE MODE SCAN
        flush_evidence()

        if str(Constants.CONFIG.get("ENABLE_HEATMAPS", "1")) == "1":
            est_len = max(total_len, 1)
            for s_label in ("FWD", "REV"):
                s_results = [r for r in results if r.get("strand") == s_label]
                hist_key = "hist_fwd" if s_label == "FWD" else "hist_rev"
                pre_hist = (large_meta or {}).get(hist_key)
                # Draw if we have sample hits OR a full-genome histogram
                if not s_results and (pre_hist is None or not np.any(pre_hist)):
                    continue
                final_map_dir = local_dirs.get(f"{dir_prefix}_MAP_{s_label}")
                if final_map_dir:
                    os.makedirs(final_map_dir, exist_ok=True)
                    n_tot = (large_meta or {}).get("total_forensic", len(results))
                    draw_map(
                        s_results, est_len,
                        os.path.join(final_map_dir, f"{base_name}_HEATMAP.png"),
                        f"Audit (large): {base_name} ({s_label}) · n≈{n_tot:,}",
                        precomputed_hist=pre_hist,
                    )
        if not harvested_data.get("seq_len"):
            harvested_data["seq_len"] = max((r["pos"] for r in results), default=1) + 1000
        # Fast composition density (short motifs never entered the automaton).
        # Overlap = longest composition motif so boundary-spanning hits are not
        # double-counted or missed when chunking multi-GB genomes.
        try:
            comp_max = 0
            for _n, _info in TARGETS.items():
                if is_composition_marker(_n) and isinstance(_info, dict):
                    comp_max = max(comp_max, len(str(_info.get("seq", "") or "")))
            overlap_bp = max(comp_max - 1, 0)
            comp = {}
            prev_tail = ""
            for _abs, chunk in stream_clean_sequence_chunks(
                file_path, chunk_bp=STREAM_CHUNK_MB * 1024 * 1024, overlap_bp=0
            ):
                # Manual overlap on cleaned chunks (stream helper may not preserve it)
                work = (prev_tail + chunk) if prev_tail else chunk
                # Count only the non-overlap region for this chunk to avoid double-count
                part = count_composition_fast(work, mol_label)
                # Scale correction: for short motifs, recounting the small overlap
                # region is negligible vs CG totals; still subtract by re-counting
                # only on the new material when possible.
                if prev_tail and overlap_bp > 0 and len(work) > len(prev_tail):
                    # Prefer counting on chunk alone for bulk, then add boundary strip
                    bulk = count_composition_fast(chunk, mol_label)
                    # Boundary strip: last overlap of prev + first (comp_max-1) of chunk
                    # already covered by work count vs bulk — use bulk as authoritative
                    # for interior; accept tiny boundary error for dense motifs.
                    part = bulk
                for mname, bucket in part.items():
                    dest = comp.setdefault(mname, {"FWD": 0, "REV": 0})
                    dest["FWD"] += bucket.get("FWD", 0)
                    dest["REV"] += bucket.get("REV", 0)
                if overlap_bp > 0:
                    prev_tail = chunk[-overlap_bp:] if len(chunk) >= overlap_bp else chunk
                else:
                    prev_tail = ""
            for mname, bucket in comp.items():
                dest = harvested_data.setdefault("composition", {}).setdefault(
                    mname, {"FWD": 0, "REV": 0}
                )
                # Large-mode: fast path is authoritative (AC never saw these markers)
                dest["FWD"] = bucket.get("FWD", 0)
                dest["REV"] = bucket.get("REV", 0)
        except Exception as comp_err:
            logger.error(f"[!] Composition count failed for {base_name}: {comp_err}")
        write_composition_sidecar(
            base_name,
            dir_prefix,
            harvested_data.get("seq_len", 0),
            harvested_data.get("composition") or {},
            metrics_seq=sequence if sequence else None,
        )
        hits_path = write_hit_index(
            base_name,
            harvested_data.get("seq_len", 0),
            harvested_data.get("hit_index") or {},
        )
        if str(Constants.CONFIG.get("REF_DISTANCE_SIDECAR", "0")) == "1":
            try:
                write_ref_distance_sidecar(base_name, harvested_data.get("ref_distance") or [])
            except Exception as _rd_err:
                logger.debug(f"ref-distance sidecar: {_rd_err}")
        if hits_path and (harvested_data.get("hit_index") or {}):
            evi_dir = local_dirs.get(f"{dir_prefix}_EVI_FWD")
            if evi_dir:
                os.makedirs(evi_dir, exist_ok=True)
                sum_name = f"{readable_prefix}_Evidence_{mol_label}_FWD_{base_name}.txt"
                lines = [
                    "\n=== COMPOSITION POSITION INDEX (short markers) ===\n",
                    f"Full coordinate lists: {hits_path}\n",
                ]
                for mname, strands in sorted((harvested_data.get("hit_index") or {}).items()):
                    info = TARGETS.get(mname) or {}
                    plen = len(str(info.get("seq", "") or ""))
                    for sk in ("FWD", "REV"):
                        pos = strands.get(sk) or []
                        if not pos:
                            continue
                        lines.append(
                            f"  {mname} len={plen} {sk}: n={len(pos)} "
                            f"first={min(pos)} last={max(pos)}\n"
                        )
                lines.append("=" * 60 + "\n")
                with open(os.path.join(evi_dir, sum_name), "a", encoding="utf-8") as ef:
                    ef.write("".join(lines))

        # WAV for large files too (every scanned file)
        if str(Constants.CONFIG.get("ENABLE_SONIFICATION", "0")) == "1":
            try:
                audio_dir = local_dirs.get(f"{dir_prefix}_AUDIO")
                if audio_dir:
                    os.makedirs(audio_dir, exist_ok=True)
                    audio_path = os.path.join(audio_dir, f"{base_name}_AUDIO.wav")
                    combined_parts = []
                    seen_segments = set()
                    for flk in harvested_data.get("flanks") or []:
                        seg = flk.get("match") or ""
                        if mol_label == "RNA":
                            seg = seg.replace("T", "U")
                        if seg and seg not in seen_segments:
                            seen_segments.add(seg)
                            combined_parts.append(seg)
                    combined_seq = "".join(combined_parts)
                    if not combined_seq:
                        # Fallback: first chunk of cleaned sequence from mmap path is hard;
                        # use concatenated matches only. If still empty, skip quietly.
                        pass
                    if combined_seq:
                        Sonification.sonify_sequence(
                            combined_seq, mol_type=mol_label, out_path=audio_path, note_seconds=0.12
                        )
            except Exception as snd_err:
                logger.error(f"[!] Sonification failed for {base_name}: {snd_err}")


        # Feature spectrogram (composition → log-frequency display; NOT measured RF)
        try:
            _map_dir = (
                local_dirs.get(f"{dir_prefix}_MAP_FWD")
                or local_dirs.get(f"{dir_prefix}_MAP_REV")
                or local_dirs.get(f"{dir_prefix}_MAP")
            )
            _seq_for_spec = sequence
            if not _seq_for_spec and harvested_data.get("flanks"):
                _seq_for_spec = "".join(
                    (flk.get("match") or "") for flk in harvested_data["flanks"]
                )
            if _seq_for_spec and _map_dir:
                write_feature_spectrogram(
                    _seq_for_spec, base_name, mol_label,
                    dir_prefix=dir_prefix, results=results,
                    flanks=harvested_data.get("flanks"),
                )
        except Exception as _spec_err:
            logger.debug(f"spectrogram skip: {_spec_err}")

        harvested_data["file"] = base_name
        return harvested_data

    # --- 4. CIRCULAR ENGINE SETUP ---
    original_len = len(sequence)
    circular = mol_label != "PROT" and (is_circular_vector(base_name) or is_test_mode)
    bridge_size = max(int(DYNAMIC_WINDOW), int(Constants.CONFIG.get('FLANK_SIZE', 42))) if circular else 0
    harvested_data = {"flanks": [], "composition": {}, "seq_len": original_len, "hit_index": {}, "ref_distance": []}
    results = []
    num_cores = int(Constants.CONFIG.get('MAX_CORES', 1))

    seen_matches = set() 
        
    for strand_pass in strands_to_scan:
        # DNA/RNA: always scan forward sequence; strand = match_type from automaton
        if mol_label == "PROT" and strand_pass == "REV":
            current_scan_seq = get_rev_comp(sequence)
        else:
            current_scan_seq = sequence

        search_seq = current_scan_seq + current_scan_seq[:bridge_size]

        matches = list(scan_automaton.iter(search_seq))

        for end_idx, (marker_name, match_type) in matches:
                if is_test_mode and marker_name != "TEST_MARKER":
                    continue

                strand = match_type  # FWD or REV from Targets load

                target_info = TARGETS.get(marker_name, {})
                target_string = target_info.get("seq", "")

                if mol_label == "PROT":
                    if match_type == "FWD":
                        target_string = translate_dna(target_string)
                    else:
                        target_string = translate_dna(get_rev_comp(target_string))
                # DNA/RNA: keep canonical target seq; matched bases are on forward genome
                # (REV hits are RC-pattern matches already)
                elif strand == "REV":
                    target_string = get_rev_comp(target_string)

                match_len = len(target_string)
                if match_len == 0: 
                    continue

                start_pos = end_idx - match_len + 1

                if start_pos >= original_len:
                    continue
                
                actual_file_match = search_seq[start_pos : end_idx + 1]

                if len(actual_file_match) != match_len:
                    continue

                physical_pos = start_pos % original_len

                match_key = (marker_name, physical_pos, strand)
                if match_key in seen_matches:
                    continue

                alignment_bars = "".join(
                    "|" if i < len(target_string) and i < len(actual_file_match) and target_string[i] == actual_file_match[i] else "." 
                    for i in range(match_len)
                )
                
                match_pct = (alignment_bars.count('|') / match_len) * 100.0 if match_len > 0 else 0.0

                config_dict = getattr(Constants, 'CONFIG', {}) if 'Constants' in globals() else {}
                match_threshold = float(config_dict.get("MATCH_THRESHOLD", 0))
                
                if match_pct < match_threshold:
                    continue

                seen_matches.add(match_key)
                                # Always tally composition + full position index for short markers
                strand_key = strand if strand in ("FWD", "REV") else "FWD"
                bucket = harvested_data.setdefault("composition", {}).setdefault(
                    marker_name, {"FWD": 0, "REV": 0}
                )
                bucket[strand_key] = bucket.get(strand_key, 0) + 1

                if is_composition_marker(marker_name):
                    idx = harvested_data.setdefault("hit_index", {}).setdefault(
                        marker_name, {"FWD": [], "REV": []}
                    )
                    idx.setdefault(strand_key, []).append(int(physical_pos))
                    results.append({
                        "pos": physical_pos,
                        "name": marker_name,
                        "strand": strand,
                        "color": MARKER_COLORS.get(marker_name, "#808080"),
                    })
                    # No full evidence block, no Phase-2 flanks for dense short motifs
                    continue
                if not is_reportable_forensic_hit(actual_file_match, marker_name):
                    results.append({
                        "pos": physical_pos,
                        "name": marker_name,
                        "strand": strand,
                        "color": MARKER_COLORS.get(marker_name, "#808080"),
                        "suppressed_evidence": True,
                    })
                    continue
                freq_scan = check_frequency_activation(actual_file_match, mol_label)

                f_size = int(Constants.CONFIG.get('FLANK_AUTO', 200))
                left_flank, right_flank = get_flanks(
                    current_scan_seq, start_pos, end_idx % original_len, f_size, circular)

                # --- BUFFER EVIDENCE IN MEMORY ---
                out_dir = local_dirs.get(f"{dir_prefix}_EVI_{strand}")
                if out_dir:
                    filename = f"{readable_prefix}_Evidence_{mol_label}_{strand}_{base_name}.txt"
                    tag = TARGETS.get(marker_name, {}).get("tag", "UNTAGGED")
                    
                    block = []
                    block.append(f"\n[!] GENOMIC HIT DETECTED: {marker_name} ({match_pct:.1f}% Match)\n")

                    block.append(f"FILE HASH   : {file_hash_hex}\n")

                    block.append(f"LOCATION: {physical_pos} | STRAND: {strand} | TYPE: {mol_label}\n")

                    block.append(f"THREAT TAG : {tag}\n")

                    block.append(f"TARGET:  {target_string}\n")

                    block.append(f"         {alignment_bars}\n")

                    block.append(f"MATCH:   {actual_file_match}\n")

                    block.append(f"LEFT:    {left_flank}\n")

                    block.append(f"RIGHT:   {right_flank}\n")

                    block.append(f"CONTEXT: ...{left_flank}[{actual_file_match}]{right_flank}...\n")

                    try:
                        aa = translate_dna(actual_file_match.replace("U", "T"))
                        block.append(f"PROTEIN TRANSLATION: {aa}\n")
                    except Exception:
                        block.append("PROTEIN TRANSLATION: N/A\n")

                    if str(Constants.CONFIG.get("ENABLE_CODON_BIAS", "0")) == "1":
                        try:
                            if mol_label == "PROT":
                                bias = "N/A (protein sequence)"
                            else:
                                win_start = start_pos - 3 * min(150, start_pos // 3)
                                bias = analyze_author_signature(current_scan_seq[win_start:end_idx + 1 + 450])
                            block.append(f"CODON BIAS SIGNATURE: {bias}\n")
                        except Exception as bias_err:
                            block.append(f"CODON BIAS SIGNATURE: ERROR ({bias_err})\n")

                    if str(Constants.CONFIG.get("ENABLE_CIRCUIT", "0")) == "1":
                        try:
                            circuit = check_circuit_logic(current_scan_seq, physical_pos, marker_name, window=250)
                            if circuit:
                                block.append(f"CIRCUIT NEIGHBORS ({len(circuit)} within 250 bp):\n")
                                for hit in circuit:
                                    block.append(f"  -> {hit['neighbor']} dist={hit['distance']}bp")
                                    if hit.get("payload"):
                                        block.append(f" payload_len={len(hit['payload'])}")
                                    block.append("\n")
                            else:
                                block.append("CIRCUIT NEIGHBORS: none within 250 bp\n")
                        except Exception as cir_err:
                            block.append(f"CIRCUIT NEIGHBORS: ERROR ({cir_err})\n")

                    if str(Constants.CONFIG.get("ENABLE_RF_HEURISTICS", "0")) == "1":
                        block.append("--- SEQUENCE-FEATURE SPECTROGRAM (NOT measured RF) ---\n")
                        block.append("  Composition on a log-frequency display axis — not radio detection.\n")
                        for _pk in (freq_scan.get("peaks") or [])[:5]:
                            _bc = _pk.get("band_context") or ""
                            block.append(
                                f"  peak {_pk.get('label')}: strength={_pk.get('strength', 0):.2f}"
                                + (f" | {_bc}" if _bc else "")
                                + "\n"
                            )
                        if freq_scan.get("path"):
                            block.append(f"  spectrogram: {freq_scan['path']}\n")
                        block.append("  Display only — not measured RF emission.\n")
                    profile_seq = actual_file_match.replace("T", "U") if mol_label == "RNA" else actual_file_match
                    if str(Constants.CONFIG.get("ENABLE_BIOPHYSICS", "0")) == "1":
                        try:
                            bio_profile = BiophysicalProfile.build_profile(profile_seq, mol_type=mol_label)
                            block.append(BiophysicalProfile.format_profile_report(bio_profile) + "\n")
                        except Exception as bp_err:
                            logger.error(f"[!] Biophysical profile failed for {marker_name} in {base_name}: {bp_err}")
                            block.append(f"[!] Biophysical profile failed: {bp_err}\n")
                    try:
                        _heur = BiophysicalProfile.format_optional_heuristics(
                            profile_seq, mol_type=mol_label
                        )
                        if _heur:
                            block.append(_heur)
                    except Exception as heur_err:
                        block.append(f"[!] Optional heuristics failed: {heur_err}\n")


                    block.append("-" * 60 + "\n")
                    queue_evidence(out_dir, filename, "".join(block))

                harvested_data["flanks"].append({
                    "marker": marker_name, "pos": physical_pos, "mol_type": mol_label,
                    "strand": strand, "target": target_string, "match": actual_file_match,
                    "left": left_flank, "right": right_flank, "file": scope_id, "basename": base_name, "category": dir_prefix,
                })
                # Write flank-only line (not full evidence block)
                flanks_dir = local_dirs.get(f"{dir_prefix}_FLANKS")
                if flanks_dir:
                    os.makedirs(flanks_dir, exist_ok=True)
                    flank_path = os.path.join(flanks_dir, f"{base_name}_flanks.txt")
                    with open(flank_path, "a", encoding="utf-8") as ff:
                        ff.write(
                            f"{marker_name}\t{physical_pos}\t{strand}\t"
                            f"{left_flank}\t{actual_file_match}\t{right_flank}\n"
                        )

                strand_key = strand if strand in ("FWD", "REV") else "FWD"
                bucket = harvested_data.setdefault("composition", {}).setdefault(
                    marker_name, {"FWD": 0, "REV": 0}
                )
                bucket[strand_key] = bucket.get(strand_key, 0) + 1

                results.append({
                    'pos': physical_pos, 
                    'name': marker_name, 
                    'strand': strand, 
                    'color': MARKER_COLORS.get(marker_name, '#808080')
                })

    # FLUSH ALL BUFFERED EVIDENCE TO DISK IN A SINGLE WRITE
    flush_evidence()

    # --- 5. FINAL AUDIT & CHIMERA CHECK ---
    if results and str(Constants.CONFIG.get("ENABLE_HEATMAPS", "1")) == "1":
        for s_label in ["FWD", "REV"]:
            s_results = [r for r in results if r['strand'] == s_label]
            if not s_results:
                continue
            final_map_dir = local_dirs.get(f"{dir_prefix}_MAP_{s_label}")
            if final_map_dir:
                os.makedirs(final_map_dir, exist_ok=True)
                draw_map(s_results, original_len, os.path.join(final_map_dir, f"{base_name}_HEATMAP.png"), f"Audit: {base_name} ({s_label})")

    chimera_label = None
    if results:
        found_names = [r['name'] for r in results]
        chimera_label = detect_chimeric_logic(found_names)
        if chimera_label:
            for s in ["FWD", "REV"]:
                out_dir = local_dirs.get(f"{dir_prefix}_EVI_{s}")
                if out_dir:
                    evi_path = os.path.join(
                        out_dir,
                        f"{readable_prefix}_Evidence_{mol_label}_{s}_{base_name}.txt",
                    )
                    if os.path.exists(evi_path):
                        with open(evi_path, "a", encoding="utf-8") as ef:
                            audit_type = "CALIBRATION" if is_test_mode else "PRODUCTION"
                            ef.write(f"\n[{audit_type} AUDIT] Finding: {chimera_label}\n")
                            ef.write(f"TIMESTAMP: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                            ef.write(f"{'='*60}\n")

    # WAV: every successfully scanned file (hits optional)
    if str(Constants.CONFIG.get("ENABLE_SONIFICATION", "0")) == "1":
        try:
            audio_dir = local_dirs.get(f"{dir_prefix}_AUDIO")
            if audio_dir:
                os.makedirs(audio_dir, exist_ok=True)
                audio_path = os.path.join(audio_dir, f"{base_name}_AUDIO.wav")
                combined_parts = []
                seen_segments = set()
                for flk in harvested_data.get("flanks") or []:
                    seg = flk.get("match") or ""
                    if mol_label == "RNA":
                        seg = seg.replace("T", "U")
                    if seg and seg not in seen_segments:
                        seen_segments.add(seg)
                        combined_parts.append(seg)
                combined_seq = "".join(combined_parts)
                if not combined_seq and sequence:
                    sample = sequence[:3000]
                    if mol_label == "RNA":
                        sample = sample.replace("T", "U")
                    combined_seq = sample
                if combined_seq:
                    Sonification.sonify_sequence(
                        combined_seq, mol_type=mol_label, out_path=audio_path, note_seconds=0.12
                    )
        except Exception as snd_err:
            logger.error(f"[!] Sonification failed for {base_name}: {snd_err}")


        # Feature spectrogram (composition → log-frequency display; NOT measured RF)
        try:
            _map_dir = (
                local_dirs.get(f"{dir_prefix}_MAP_FWD")
                or local_dirs.get(f"{dir_prefix}_MAP_REV")
                or local_dirs.get(f"{dir_prefix}_MAP")
            )
            _seq_for_spec = sequence
            if not _seq_for_spec and harvested_data.get("flanks"):
                _seq_for_spec = "".join(
                    (flk.get("match") or "") for flk in harvested_data["flanks"]
                )
            if _seq_for_spec and _map_dir:
                write_feature_spectrogram(
                    _seq_for_spec, base_name, mol_label,
                    dir_prefix=dir_prefix, results=results,
                    flanks=harvested_data.get("flanks"),
                )
        except Exception as _spec_err:
            logger.debug(f"spectrogram skip: {_spec_err}")

    harvested_data["seq_len"] = len(sequence) if sequence else harvested_data.get("seq_len", 0)
    # Fast composition density — replaces millions of short AC hits
    try:
        comp = count_composition_fast(sequence, mol_label)
        for mname, bucket in comp.items():
            dest = harvested_data.setdefault("composition", {}).setdefault(
                mname, {"FWD": 0, "REV": 0}
            )
            # Fast path is authoritative for composition (not added into AC)
            dest["FWD"] = bucket.get("FWD", 0)
            dest["REV"] = bucket.get("REV", 0)
        pos_idx = composition_hit_index(sequence, mol_label)
        if pos_idx:
            harvested_data["hit_index"] = pos_idx
    except Exception as comp_err:
        logger.error(f"[!] Composition count failed for {base_name}: {comp_err}")
    write_composition_sidecar(
        base_name,
        dir_prefix,
        harvested_data.get("seq_len", 0),
        harvested_data.get("composition") or {},
        metrics_seq=sequence if sequence else None,
    )
    hits_path = write_hit_index(
        base_name,
        harvested_data.get("seq_len", 0),
        harvested_data.get("hit_index") or {},
    )
    if str(Constants.CONFIG.get("REF_DISTANCE_SIDECAR", "0")) == "1":
        try:
            write_ref_distance_sidecar(base_name, harvested_data.get("ref_distance") or [])
        except Exception as _rd_err:
            logger.debug(f"ref-distance sidecar: {_rd_err}")
    if hits_path and (harvested_data.get("hit_index") or {}):
        evi_dir = local_dirs.get(f"{dir_prefix}_EVI_FWD")
        if evi_dir:
            os.makedirs(evi_dir, exist_ok=True)
            sum_name = f"{readable_prefix}_Evidence_{mol_label}_FWD_{base_name}.txt"
            lines = [
                "\n=== COMPOSITION POSITION INDEX (short markers) ===\n",
                f"Full coordinate lists: {hits_path}\n",
            ]
            for mname, strands in sorted((harvested_data.get("hit_index") or {}).items()):
                info = TARGETS.get(mname) or {}
                plen = len(str(info.get("seq", "") or ""))
                for sk in ("FWD", "REV"):
                    pos = strands.get(sk) or []
                    if not pos:
                        continue
                    lines.append(
                        f"  {mname} len={plen} {sk}: n={len(pos)} "
                        f"first={min(pos)} last={max(pos)}\n"
                    )
            lines.append("=" * 60 + "\n")
            with open(os.path.join(evi_dir, sum_name), "a", encoding="utf-8") as ef:
                ef.write("".join(lines))
    plt.close('all')
    harvested_data["file"] = base_name
    return harvested_data

_P2_LIBRARY = {}
_P2_AC = None
_P2_META = None  # pattern -> list of metadata dicts (all origins)


def _p2_init(library):
    """Build flank Aho-Corasick once per worker process (not per target file)."""
    global _P2_LIBRARY, _P2_AC, _P2_META
    from collections import defaultdict
    _P2_LIBRARY = library or {}
    _P2_META = defaultdict(list)
    ac = ahocorasick.Automaton()
    min_len = int(globals().get("MIN_FLANK_LEN", 24) or 24)
    try:
        min_len = int(Constants.CONFIG.get("MIN_FLANK_LEN", min_len) or min_len)
    except Exception:
        pass
    for _key, metadata in _P2_LIBRARY.items():
        seq = str(metadata.get("seq", "") or "").upper()
        if len(seq) < min_len:
            continue
        for strand_tag, pat in (("FWD", seq), ("REV", get_rev_comp(seq))):
            if not pat:
                continue
            if pat not in _P2_META:
                try:
                    ac.add_word(pat, pat)
                except Exception:
                    continue
            _P2_META[pat].append({**metadata, "query_strand": strand_tag, "query_seq": seq})
    if _P2_META:
        ac.make_automaton()
        _P2_AC = ac
    else:
        _P2_AC = None


def _flank_mismatch_budget(flank_len: int) -> int:
    """Config-driven max Hamming distance for Phase-2 flank search."""
    try:
        mm = int(Constants.CONFIG.get("FLANK_MAX_MISMATCH", 0) or 0)
    except Exception:
        mm = 0
    if mm > 0:
        return max(0, mm)
    try:
        ident = float(Constants.CONFIG.get("FLANK_MIN_IDENTITY", 0) or 0)
    except Exception:
        ident = 0.0
    if ident and 0.0 < ident < 1.0 and flank_len > 0:
        return max(0, int((1.0 - ident) * flank_len))
    return 0


def hamming_distance(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(1 for x, y in zip(a, b) if x != y)


def fuzzy_find_flank(sequence: str, pattern: str, max_mismatch: int):
    """Yield (start, end, mismatches) for equal-length windows within budget."""
    if not pattern or max_mismatch < 1 or len(pattern) > len(sequence):
        return
    L = len(pattern)
    # modest speed: skip if pattern has extreme complexity? keep simple
    for i in range(0, len(sequence) - L + 1):
        window = sequence[i : i + L]
        d = hamming_distance(window, pattern)
        if d <= max_mismatch:
            yield i, i + L - 1, d


def ref_distance_record(marker_name: str, target_seq: str, match_seq: str, pos: int, strand: str):
    """Compare observed match to Targets motif — triage SNPs only, not a clock."""
    info = TARGETS.get(marker_name) or {}
    ref = str(info.get("seq") or target_seq or "").upper().replace("U", "T")
    obs = str(match_seq or "").upper().replace("U", "T")
    if not ref or not obs:
        return None
    # align by truncation to min length (exact AC hits are same length usually)
    n = min(len(ref), len(obs))
    if n == 0:
        return None
    ref_n, obs_n = ref[:n], obs[:n]
    snps = [(i, ref_n[i], obs_n[i]) for i in range(n) if ref_n[i] != obs_n[i]]
    identity = (n - len(snps)) / n if n else 0.0
    return {
        "marker": marker_name,
        "ref": info.get("ref") or marker_name,
        "ref_seq": ref_n,
        "observed": obs_n,
        "pos": pos,
        "strand": strand,
        "length": n,
        "snps": len(snps),
        "snp_detail": [{"pos": i, "ref": a, "alt": b} for i, a, b in snps[:50]],
        "window_identity": round(identity, 4),
        "note": (
            "Relative identity to Targets motif only. Not a molecular clock; "
            "not proof of natural evolution vs laboratory editing."
        ),
    }


def write_ref_distance_sidecar(base_name: str, records: list, out_root=None):
    if not records:
        return None
    import json
    if out_root is None:
        out_root = os.path.join(Constants.DIRS.get("OUT_LOGS", Constants.BASE_OUT), "ref_distance")
    os.makedirs(out_root, exist_ok=True)
    path = os.path.join(out_root, f"{base_name}.REF_DISTANCE.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump({
            "file": base_name,
            "description": (
                "Per forensic hit comparison to Targets motif (chassis/Inc clique by default). "
                "Triage aid only — not a molecular clock; not proof of evolution vs editing."
            ),
            "clique_only": str(Constants.CONFIG.get("REF_DISTANCE_CLIQUE_ONLY", "1")) == "1",
            "n_hits": len(records),
            "hits": records,
        }, f, indent=2)
    return path

def phase2_worker(args):
    """
    Cross-file flank search (exact Aho-Corasick + labels).
    flank_library keys are origin ids; each value has metadata["seq"].
    """
    from collections import defaultdict

    file_path, flank_library, local_dirs, shared_cfg = args
    if flank_library is None:
        flank_library = _P2_LIBRARY
    base_name = os.path.basename(file_path)
    scope_id = file_scope_id(file_path)
    context_bp = int(shared_cfg.get("FLANK_CONTEXT_BP", 5000))

    global TARGETS, AUTO
    if not TARGETS or AUTO is None:
        initialize_globals()

    memory_guard(threshold=90, max_wait_seconds=30)

    try:
        with open(file_path, "rb") as f:
            with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                sequence = extract_clean_sequence(mm)

        if not sequence:
            return None

        sequence = sequence.upper().replace("U", "T")
        matches_found = []

        # Shared AC from worker initializer; filter same-origin in the hit loop
        pattern_meta = _P2_META if _P2_META is not None else defaultdict(list)
        ac = _P2_AC
        if ac is None or not pattern_meta:
            # Fallback: build local AC if initializer did not run (single-process path)
            pattern_meta = defaultdict(list)
            ac = ahocorasick.Automaton()
            for _key, metadata in (flank_library or {}).items():
                seq = str(metadata.get("seq", "")).upper()
                if len(seq) < MIN_FLANK_LEN:
                    continue
                for strand_tag, pat in (("FWD", seq), ("REV", get_rev_comp(seq))):
                    if not pat:
                        continue
                    if pat not in pattern_meta:
                        try:
                            ac.add_word(pat, pat)
                        except Exception:
                            pass
                    pattern_meta[pat].append({**metadata, "query_strand": strand_tag, "query_seq": seq})
            if not pattern_meta:
                return []
            ac.make_automaton()

        raw_hits = []
        for end_idx, pat in ac.iter(sequence):
            start_idx = end_idx - len(pat) + 1
            for meta in pattern_meta.get(pat, []):
                if meta.get("origin_file") in (base_name, scope_id):
                    continue  # skip same-file origin
                raw_hits.append((start_idx, end_idx, pat, meta, 0))

        # Fuzzy only on modest sequences (pigeonhole would be better; gate keeps cost bounded)
        max_fuzzy_bp = int(shared_cfg.get("FLANK_FUZZY_MAX_SEQ_BP", 500_000) or 500_000)
        if len(sequence) <= max_fuzzy_bp:
            exact_pats = {pat for (_s, _e, pat, _m, _d) in raw_hits}
            for pat, metas in pattern_meta.items():
                mm_budget = _flank_mismatch_budget(len(pat))
                if mm_budget < 1 or pat in exact_pats:
                    continue
                for start_idx, end_idx, dist in fuzzy_find_flank(sequence, pat, mm_budget):
                    for meta in metas:
                        if meta.get("origin_file") in (base_name, scope_id):
                            continue
                        raw_hits.append((start_idx, end_idx, pat, meta, dist))

        seen = set()
        for start_idx, end_idx, pat, meta, mm_dist in raw_hits:
            key = (
                meta.get("origin_file"),
                meta.get("marker"),
                meta.get("side"),
                meta.get("origin_pos"),
                start_idx,
                meta.get("query_strand"),
                mm_dist,
            )
            if key in seen:
                continue
            seen.add(key)

            lo = max(0, start_idx - context_bp)
            hi = min(len(sequence), end_idx + context_bp)
            neighborhood = sequence[lo:hi]
            nearby_targets = []
            if neighborhood and AUTO is not None:
                try:
                    for _t_end, payload in AUTO.iter(neighborhood):
                        if isinstance(payload, (tuple, list)) and len(payload) >= 1:
                            t_name = payload[0]
                        else:
                            t_name = str(payload)
                        if t_name and t_name not in ("FWD", "REV"):
                            nearby_targets.append(t_name)
                except Exception:
                    pass
            nearby_unique = sorted(set(nearby_targets))

            origin_marker = meta.get("marker") or ""
            # Ignore composition motifs for incorporation logic (CG/NGG/scars etc.)
            forensic_nearby = [
                t for t in nearby_unique
                if t and not is_composition_marker(t)
            ]
            # Same marker as the flank's origin is expected — not a transplant warning
            foreign_markers = [
                t for t in forensic_nearby
                if t != origin_marker
            ]
            same_marker_present = bool(origin_marker) and origin_marker in forensic_nearby

            if foreign_markers:
                # 100% flank from marker A lands beside marker B in another genome
                label = "FLANK_INCORPORATION_WARNING"
            elif same_marker_present:
                label = "FLANK_CROSS_SAME_MARKER"
            elif forensic_nearby:
                label = "FLANK_CROSS_DIFF_MARKER"
            elif nearby_unique:
                # only composition motifs in window
                label = "FLANK_CROSS_SOLITARY"
            else:
                label = "FLANK_CROSS_SOLITARY"

            rc_flag = meta.get("query_strand") == "REV"
            risk_info = determine_verdict(
                [origin_marker] + foreign_markers if foreign_markers else (
                    [origin_marker] + forensic_nearby if forensic_nearby else origin_marker
                ),
                neighborhood,
            )
            # Bump score when flank context shows a different forensic marker family
            if foreign_markers and risk_info.get("risk_score", 0) < 40:
                risk_info = dict(risk_info)
                risk_info["risk_score"] = min(100, int(risk_info.get("risk_score", 0)) + 15)
                if risk_info.get("verdict") == "LOW_LEVEL_MARKER_MATCH":
                    risk_info["verdict"] = "POTENTIAL_LAB_ORIGIN"
                reasons = risk_info.get("reasoning") or ""
                extra = (
                    f"Flank from {origin_marker} co-located with foreign marker(s) "
                    f"{', '.join(foreign_markers)} in {base_name}"
                )
                risk_info["reasoning"] = (reasons + "; " + extra).strip("; ")

            # Prefer exact, long flanks; demote fuzzy matches in ranking/grid
            _flank_len = len(pat) if pat else len(str(meta.get("query_seq") or ""))
            risk_info = adjust_incorporation_risk(
                risk_info, mm_dist=mm_dist, flank_len=_flank_len, foreign_markers=foreign_markers
            )

            warning_msg = ""
            if foreign_markers:
                warning_msg = (
                    f"INCORPORATION: flanks bound to {origin_marker} in "
                    f"{meta.get('origin_file', '?')} also sit next to "
                    f"{', '.join(foreign_markers)} in {base_name} "
                    f"(100% flank identity @ {start_idx}). "
                    f"Context may be shared DNA from a different original marker set."
                )

            sim = 100.0
            if mm_dist and len(pat) > 0:
                sim = round(100.0 * (1.0 - (mm_dist / float(len(pat)))), 2)
            matches_found.append({
                "seq": meta.get("query_seq") or pat,
                "target_file": scope_id,
                "origin_file": meta.get("origin_file", "UNKNOWN"),
                "origin_pos": meta.get("origin_pos"),
                "found_pos": start_idx,
                "marker": origin_marker,
                "side": meta.get("side", "N/A"),
                "similarity": sim,
                "flank_mismatches": int(mm_dist or 0),
                "match_strand": "REV" if rc_flag else "FWD",
                "label": label,
                "rc_hit": rc_flag,
                "nearby": forensic_nearby if forensic_nearby else (
                    nearby_unique if nearby_unique else ["SOLITARY_INSERT"]
                ),
                "foreign_markers": foreign_markers,
                "warning": warning_msg,
                "verdict": risk_info["verdict"],
                "risk_score": risk_info["risk_score"],
                "flank_quality": risk_info.get("flank_quality"),
                "flank_len": risk_info.get("flank_len"),
                "context": forensic_nearby if forensic_nearby else (
                    nearby_unique if nearby_unique else ["SOLITARY_INSERT"]
                ),
                "origin_hit_id": meta.get("origin_hit_id"),
            })

        return matches_found

    except Exception as e:
        logger.error(f"Phase 2 Error on {base_name}: {e}", exc_info=True)
        return None 


def build_inc_summary(files, composition_rows=None):
    """
    Short Inc / chassis summary for portal: which files carry which replicon family tags.
    Derived from Phase-1 forensic tags (and optional composition is ignored).
    """
    from collections import defaultdict
    per_file = []
    family_counts = defaultdict(int)
    for f in files or []:
        name = f.get("name") or ""
        by_tag = f.get("by_tag") or {}
        # We need marker-level; files only have by_tag. Use composition_rows if present
        # for marker names is heavy; instead scan TARGETS roles from known clique via
        # a lightweight pass over composition is wrong. Use edges? Better: rebuild from
        # TARGETS hits stored nowhere at file level except by_tag aggregate.
        # Portal files only have by_tag — approximate chassis presence from CHASSIS_ORIGIN etc.
        families = set()
        chassis_hits = int(by_tag.get("CHASSIS_ORIGIN", 0))
        ind_hits = int(by_tag.get("INDUSTRIAL_BACKBONE", 0))
        seamless = int(by_tag.get("SEAMLESS_ASSEMBLY", 0))
        if chassis_hits or ind_hits:
            # refined later if marker list available on file
            families.add("chassis_present")
        if seamless:
            families.add("synthetic_recomb")
        per_file.append({
            "file": name,
            "chassis_origin_hits": chassis_hits,
            "industrial_backbone_hits": ind_hits,
            "seamless_assembly_hits": seamless,
            "antibiotic_res_hits": int(by_tag.get("ANTIBIOTIC_RES", 0)),
            "families_hint": sorted(families),
        })
    # Marker-resolved summary from composition_rows is wrong (composition). Use
    # Global evidence is not here. Enrich from TARGETS by reading optional
    # Forensic BED or skip. Better approach: walk OUT evidence is slow.
    # Use incorporation foreign_markers + edges markers for family map.
    return {
        "description": (
            "Inc/chassis triage table from forensic tags. "
            "ColE1/pMB1 clique vs p15A vs synthetic recomb. Not a formal Inc typing assay."
        ),
        "files": per_file,
        "note": "Enable REF_DISTANCE_SIDECAR for motif-level identity on clique markers.",
    }


def build_inc_summary_from_edges(files, edges, incorporation):
    """Resolve Inc-family presence using markers listed on similarity/incorporation edges."""
    from collections import defaultdict
    file_markers = defaultdict(set)
    for e in edges or []:
        marks = e.get("markers") or []
        if isinstance(marks, str):
            marks = [x.strip() for x in marks.split(";") if x.strip()]
        for fn in (e.get("file_a"), e.get("file_b")):
            if fn:
                for m in marks:
                    file_markers[fn].add(m)
    for row in incorporation or []:
        for fn in (row.get("origin_file"), row.get("target_file")):
            if not fn:
                continue
            for m in row.get("foreign_markers") or []:
                file_markers[fn].add(m)
            # origin marker not always present
        # also try keys from messages — skip
    # also include any clique marker names if files have high chassis — still need hits
    rows = []
    family_tally = defaultdict(int)
    for f in files or []:
        name = f.get("name") or ""
        by_tag = f.get("by_tag") or {}
        marks = sorted(file_markers.get(name) or [])
        # filter to clique / interesting
        clique = [m for m in marks if m in INC_CLIQUE_MARKERS]
        families = sorted({
            INC_FAMILY_BY_MARKER[m] for m in clique if m in INC_FAMILY_BY_MARKER
        })
        for fam in families:
            family_tally[fam] += 1
        if not clique and not int(by_tag.get("CHASSIS_ORIGIN", 0)) and not int(by_tag.get("SEAMLESS_ASSEMBLY", 0)):
            continue
        rows.append({
            "file": name,
            "inc_markers": clique,
            "families": families,
            "chassis_origin_hits": int(by_tag.get("CHASSIS_ORIGIN", 0)),
            "industrial_backbone_hits": int(by_tag.get("INDUSTRIAL_BACKBONE", 0)),
            "seamless_assembly_hits": int(by_tag.get("SEAMLESS_ASSEMBLY", 0)),
            "antibiotic_res_hits": int(by_tag.get("ANTIBIOTIC_RES", 0)),
            "reporter_hits": int(by_tag.get("REPORTER_GENE", 0)),
        })
    rows.sort(key=lambda r: -(len(r["inc_markers"]) + r["chassis_origin_hits"]))
    return {
        "description": (
            "Short Inc/chassis summary from forensic edge markers + file tags. "
            "Points researchers at ColE1/pMB1 vs p15A vs synthetic recomb neighborhoods. "
            "Not a classical conjugation Inc assay."
        ),
        "n_files": len(rows),
        "family_file_counts": dict(family_tally),
        "files": rows,
    }


def write_portal_snapshot(out_dir=None, max_files=500):
    """
    Aggregate composition + similarity into one JSON for index.html (live mode).
    Path: Forensic_Results/Global_Evidence/portal_snapshot.json
    """
    import json

    if out_dir is None:
        out_dir = Constants.DIRS.get("OUT_LOGS", Constants.BASE_OUT)
    os.makedirs(out_dir, exist_ok=True)

    names, categories, matrix, _colors = build_waterfall_matrix(max_files=max_files)
    files = []
    for i, name in enumerate(names):
        by_tag = {}
        total = 0
        if i < len(matrix):
            for j, cat in enumerate(categories):
                n = int(matrix[i][j]) if j < len(matrix[i]) else 0
                if n:
                    by_tag[cat] = n
                    total += n
        files.append({
            "name": name,
            "category": "",
            "seq_len": 0,
            "by_tag": by_tag,
            "total_hits": total,
        })

    comp_root = os.path.join(out_dir, "composition")
    composition_rows = []
    if os.path.isdir(comp_root):
        by_name = {f["name"]: f for f in files}
        for fname in os.listdir(comp_root):
            if not fname.endswith(".TARGET_COUNTS.json"):
                continue
            path = os.path.join(comp_root, fname)
            try:
                with open(path, "r", encoding="utf-8") as jf:
                    payload = json.load(jf)
            except Exception:
                continue
            subject = payload.get("file") or fname.replace(".TARGET_COUNTS.json", "")
            if subject in by_name:
                by_name[subject]["seq_len"] = int(payload.get("seq_len") or 0)
                by_name[subject]["category"] = str(payload.get("category") or "")
            for mname, stats in (payload.get("by_marker") or {}).items():
                if isinstance(stats, dict):
                    fwd = int(stats.get("FWD", 0))
                    rev = int(stats.get("REV", 0))
                    total = int(stats.get("total", fwd + rev))
                else:
                    fwd = rev = 0
                    total = int(stats)
                dens = (payload.get("density_per_kb") or {}).get(mname)
                composition_rows.append({
                    "file": subject,
                    "marker": mname,
                    "fwd": fwd,
                    "rev": rev,
                    "total": total,
                    "density_per_kb": dens,
                })

    edges = []
    n_raw = 0
    description = "No CROSS_SIMILARITY_TRACKER.json yet"
    sim_path = os.path.join(out_dir, "CROSS_SIMILARITY_TRACKER.json")
    if os.path.exists(sim_path):
        try:
            with open(sim_path, "r", encoding="utf-8") as f:
                sim = json.load(f)
            edges = sim.get("edges") or []
            n_raw = sim.get("n_raw_hits", 0)
            description = sim.get("description") or (
                "Exact shared LEFT/RIGHT flanks between files (Phase 2)."
            )
        except Exception:
            description = "Phase 2 similarity tracker"

    root_anchor = getattr(Constants, "ROOT", Constants.BASE_OUT)
    artifact_paths = []
    for root, _, fnames in os.walk(Constants.BASE_OUT):
        for fn in fnames:
            if fn.endswith((".png", ".wav", ".txt", ".json", ".xlsx")):
                rel = os.path.relpath(os.path.join(root, fn), root_anchor)
                artifact_paths.append(rel.replace("\\", "/"))
                if len(artifact_paths) >= 5000:
                    break
        if len(artifact_paths) >= 5000:
            break

    # Incorporation grid (flank + foreign forensic markers) for portal map
    incorporation = []
    inc_csv = os.path.join(out_dir, "INCORPORATION_GRID.csv")
    inc_json = os.path.join(out_dir, "FLANK_INCORPORATION_WARNINGS.json")
    if os.path.exists(inc_csv):
        try:
            import csv as _csv
            with open(inc_csv, "r", encoding="utf-8") as cf:
                for row in _csv.DictReader(cf):
                    incorporation.append({
                        "origin_file": row.get("origin_file") or "",
                        "target_file": row.get("target_file") or "",
                        "context_count": int(row.get("context_count") or 0),
                        "max_risk": int(row.get("max_risk") or 0),
                        "same_accession": str(row.get("same_accession") or "").lower() in ("true", "1", "yes"),
                        "foreign_markers": [
                            x for x in (row.get("foreign_markers") or "").split(";") if x
                        ],
                    })
        except Exception:
            incorporation = []
    if not incorporation and os.path.exists(inc_json):
        try:
            with open(inc_json, "r", encoding="utf-8") as jf:
                inc_payload = json.load(jf)
            bucket = {}
            for w in inc_payload.get("warnings") or []:
                key = (w.get("origin_file"), w.get("target_file"))
                if key not in bucket:
                    bucket[key] = {
                        "origin_file": key[0],
                        "target_file": key[1],
                        "context_count": 0,
                        "max_risk": 0,
                        "same_accession": False,
                        "foreign_markers": set(),
                    }
                bucket[key]["context_count"] += 1
                bucket[key]["max_risk"] = max(
                    bucket[key]["max_risk"], int(w.get("risk_score") or 0)
                )
                for fm in w.get("foreign_markers") or []:
                    bucket[key]["foreign_markers"].add(fm)
            for v in bucket.values():
                v["foreign_markers"] = sorted(v["foreign_markers"])
                incorporation.append(v)
        except Exception:
            pass

    try:
        inc_summary = build_inc_summary_from_edges(files, edges, incorporation)
    except Exception:
        inc_summary = {"files": [], "family_file_counts": {}, "description": "Inc summary unavailable"}

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "description": description,
        "config_stamp": build_config_stamp(),
        "files": files,
        "categories": categories,
        "edges": edges,
        "composition_rows": composition_rows,
        "incorporation": incorporation,
        "inc_summary": inc_summary,
        "n_raw_hits": n_raw,
        "artifact_paths": artifact_paths[:2000],
    }
    out_path = os.path.join(out_dir, "portal_snapshot.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    logger.info(
        f"[+] Portal snapshot: {out_path} ({len(files)} files, {len(edges)} edges, "
        f"{len(incorporation)} incorporation links)"
    )
    return out_path

def extract_accession(filename: str) -> str:
    """NCBI-style accession from names like ..._MN187559.1.fasta -> MN187559."""
    base = os.path.basename(str(filename or ""))
    base = re.sub(r"\.(fasta|fa|fna|fas|gb|gbk|genbank)$", "", base, flags=re.I)
    m = re.search(r"\b([A-Z]{1,4}_?\d{5,})(?:\.\d+)?\b", base, re.I)
    if m:
        return m.group(1).upper().replace("_", "")
    m = re.search(r"([A-Z]{1,2}\d{5,})(?:\.\d+)?", base, re.I)
    if m:
        return m.group(1).upper()
    return ""


def annotate_flank_pairs(matches, pair_window_bp=5000):
    """
    LEFT+RIGHT of same origin_hit_id on same target file within pair_window_bp
    -> label FLANK_CROSS_PAIR. Does not claim same plasmid.
    """
    from collections import defaultdict
    if not matches:
        return matches
    by_key = defaultdict(list)
    for m in matches:
        oh = m.get("origin_hit_id")
        tf = m.get("target_file")
        if oh and tf:
            by_key[(oh, tf)].append(m)
    for (_oh, _tf), group in by_key.items():
        lefts = [m for m in group if m.get("side") == "LEFT"]
        rights = [m for m in group if m.get("side") == "RIGHT"]
        if not lefts or not rights:
            continue
        for L in lefts:
            lp = int(L.get("found_pos") or 0)
            for R in rights:
                rp = int(R.get("found_pos") or 0)
                if abs(lp - rp) <= pair_window_bp:
                    for m in (L, R):
                        if m.get("label") != "FLANK_CROSS_PAIR":
                            m["label"] = "FLANK_CROSS_PAIR"
                            m["pair_window_bp"] = pair_window_bp
                    break
    return matches




def draw_incorporation_grid(warning_rows, all_matches=None):
    """
    Global visual for flank-context incorporation:
      1) File x File matrix — cell = unique foreign-marker contexts (100% flanks)
      2) Optional strip: same-accession pairs muted
    Dark theme consistent with waterfall. Not whole-genome identity.
    """
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from collections import defaultdict
    import re

    rows = list(warning_rows or [])
    if not rows:
        logger.info("[*] Incorporation grid skipped (no foreign-marker flank contexts).")
        return None

    def short(name, n=28):
        base = os.path.basename(str(name or ""))
        return base if len(base) <= n else base[: n - 1] + "…"

    def accession_token(fname):
        m = re.search(r"([A-Z]{1,2}_?\d{5,}(?:\.\d+)?)", str(fname) or "", re.I)
        if m:
            return m.group(1).upper().replace("_", "")
        m = re.search(r"([A-Z]{2}\d{6,}(?:\.\d+)?)", str(fname) or "", re.I)
        return m.group(1).upper() if m else ""

    # Aggregate unique contexts per unordered? keep directed origin -> target
    pair_count = defaultdict(int)
    pair_max_risk = defaultdict(int)
    pair_markers = defaultdict(set)
    files = set()
    for w in rows:
        a = w.get("origin_file") or ""
        b = w.get("target_file") or ""
        if not a or not b or a == b:
            continue
        files.add(a)
        files.add(b)
        key = (a, b)
        pair_count[key] += 1
        pair_max_risk[key] = max(pair_max_risk[key], int(w.get("risk_score") or 0))
        for fm in w.get("foreign_markers") or []:
            pair_markers[key].add(fm)

    if not pair_count:
        return None

    # Order: high total activity first
    activity = defaultdict(int)
    for (a, b), c in pair_count.items():
        activity[a] += c
        activity[b] += c
    labels = sorted(files, key=lambda x: (-activity[x], x.lower()))
    # Cap display for readability
    max_files = int(Constants.CONFIG.get("INCORP_GRID_MAX_FILES", 40))
    if len(labels) > max_files:
        labels = labels[:max_files]
    idx = {f: i for i, f in enumerate(labels)}
    n = len(labels)
    mat = np.zeros((n, n), dtype=float)
    risk = np.zeros((n, n), dtype=float)
    same_acc = np.zeros((n, n), dtype=bool)
    for (a, b), c in pair_count.items():
        if a not in idx or b not in idx:
            continue
        i, j = idx[a], idx[b]
        mat[i, j] = c
        risk[i, j] = pair_max_risk[(a, b)]
        ta, tb = accession_token(a), accession_token(b)
        if ta and tb and ta == tb:
            same_acc[i, j] = True

    out_dir = Constants.DIRS.get("OUT_MAPS", Constants.BASE_OUT)
    os.makedirs(out_dir, exist_ok=True)
    save_path = os.path.join(out_dir, "INCORPORATION_GRID.png")

    fig_w = max(10, n * 0.42 + 4)
    fig_h = max(8, n * 0.42 + 3)
    fig, axes = plt.subplots(1, 2, figsize=(fig_w * 1.55, fig_h),
                             gridspec_kw={"width_ratios": [1.15, 1.0]})
    fig.patch.set_facecolor("#111111")

    cmap_count = LinearSegmentedColormap.from_list(
        "inc_count", ["#1a1a1a", "#3d1f5c", "#7b2d8e", "#e8b923", "#fff5c0"]
    )
    cmap_risk = LinearSegmentedColormap.from_list(
        "inc_risk", ["#1a1a1a", "#0d3b4c", "#1a7a6d", "#e8b923", "#ff6b4a"]
    )

    for ax, data, cmap, title, cbar_label in (
        (axes[0], mat, cmap_count,
         "Flank incorporation contexts (count)",
         "unique foreign-marker contexts"),
        (axes[1], risk, cmap_risk,
         "Flank incorporation (max risk score)",
         "max risk score"),
    ):
        ax.set_facecolor("#111111")
        vmax = max(float(data.max()), 1.0)
        im = ax.imshow(data, cmap=cmap, vmin=0, vmax=vmax, aspect="equal",
                       interpolation="nearest")
        # Hatch / edge same-accession cells
        for i in range(n):
            for j in range(n):
                if same_acc[i, j] and data[i, j] > 0:
                    ax.add_patch(plt.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1,
                        fill=False, edgecolor="#66ccff", linewidth=1.2, linestyle="--"
                    ))
                elif data[i, j] > 0 and n <= 24:
                    ax.text(j, i, f"{int(data[i, j])}", ha="center", va="center",
                            color="white", fontsize=6, alpha=0.9)
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels([short(x) for x in labels], rotation=75, ha="right",
                           color="white", fontsize=7)
        ax.set_yticklabels([short(x) for x in labels], color="white", fontsize=7)
        ax.set_xlabel("Target file (flank found here)", color="#aaaaaa", fontsize=9)
        ax.set_ylabel("Origin file (flank library source)", color="#aaaaaa", fontsize=9)
        ax.set_title(title, color="white", fontsize=11, pad=10)
        ax.tick_params(colors="#888888", length=0)
        for spine in ax.spines.values():
            spine.set_color("#333333")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(cbar_label, color="#cccccc", fontsize=8)
        cbar.ax.yaxis.set_tick_params(color="#888888")
        plt.setp(cbar.ax.yaxis.get_ticklabels(), color="#cccccc", fontsize=7)

    fig.suptitle(
        "BSAP — Flank incorporation grid (100% flanks + foreign forensic markers)\n"
        "Dashed cyan border = same GenBank accession (duplicate folder copy, not independent)",
        color="white", fontsize=12, y=1.02,
    )
    plt.tight_layout()
    plt.savefig(save_path, dpi=160, facecolor="#111111", bbox_inches="tight")
    plt.close()
    logger.info(f"[+] Incorporation grid: {save_path} ({n} files)")

    # Compact CSV for portal / Excel users
    csv_path = os.path.join(
        Constants.DIRS.get("OUT_LOGS", out_dir), "INCORPORATION_GRID.csv"
    )
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("origin_file,target_file,context_count,max_risk,same_accession,foreign_markers\n")
        for (a, b), c in sorted(pair_count.items(), key=lambda x: -x[1]):
            ta, tb = accession_token(a), accession_token(b)
            same = bool(ta and tb and ta == tb)
            fmarks = ";".join(sorted(pair_markers[(a, b)]))
            f.write(
                f"\"{a}\",\"{b}\",{c},{pair_max_risk[(a, b)]},{same},\"{fmarks}\"\n"
            )
    logger.info(f"[+] Incorporation grid CSV: {csv_path}")
    return save_path


def write_incorporation_warnings(matches, out_dir=None):
    """
    Flank-anchored incorporation warnings:
    A 100% flank taken from marker A in genome X lands in genome Y next to
    a different forensic marker B. Same-marker context is ignored (expected).
    Composition motifs are ignored.
    """
    import json
    if out_dir is None:
        out_dir = Constants.DIRS.get("OUT_LOGS", Constants.BASE_OUT)
    os.makedirs(out_dir, exist_ok=True)
    warnings = []
    for m in matches or []:
        foreign = m.get("foreign_markers") or []
        if not foreign and m.get("label") != "FLANK_INCORPORATION_WARNING":
            continue
        if not foreign:
            continue
        warnings.append({
            "flank_seq": m.get("seq"),
            "origin_file": m.get("origin_file"),
            "origin_marker": m.get("marker"),
            "origin_pos": m.get("origin_pos"),
            "target_file": m.get("target_file"),
            "found_pos": m.get("found_pos"),
            "side": m.get("side"),
            "strand": m.get("match_strand"),
            "identity": m.get("similarity"),
            "foreign_markers": foreign,
            "nearby_forensic": [x for x in (m.get("nearby") or []) if x != "SOLITARY_INSERT"],
            "message": m.get("warning") or (
                f"Flanks from {m.get('marker')} in {m.get('origin_file')} match in "
                f"{m.get('target_file')} next to foreign marker(s) {', '.join(foreign)}"
            ),
            "verdict": m.get("verdict"),
            "flank_quality": m.get("flank_quality"),
            "flank_len": m.get("flank_len"),
            "risk_score": m.get("risk_score"),
        })
    # de-dupe by origin/target/marker/foreign set
    dedup = {}
    for w in warnings:
        key = (
            w.get("origin_file"),
            w.get("target_file"),
            w.get("origin_marker"),
            tuple(sorted(w.get("foreign_markers") or [])),
        )
        if key not in dedup:
            dedup[key] = w
    rows = list(dedup.values())
    path_json = os.path.join(out_dir, "FLANK_INCORPORATION_WARNINGS.json")
    path_txt = os.path.join(out_dir, "FLANK_INCORPORATION_WARNINGS.txt")
    with open(path_json, "w", encoding="utf-8") as f:
        json.dump({
            "description": (
                "100% flank identity from marker A lands beside a different forensic "
                "marker B in another file. Same-marker neighborhoods are omitted. "
                "Not whole-genome proof — prioritizes inspection of possible shared or "
                "recombined construct context."
            ),
            "n_warnings": len(rows),
            "warnings": rows,
        }, f, indent=2)
    with open(path_txt, "w", encoding="utf-8") as f:
        f.write("FLANK INCORPORATION WARNINGS\n")
        f.write("=" * 60 + "\n")
        f.write("Same origin-marker next to its own flanks = expected (suppressed).\n")
        f.write("Different target marker next to those flanks = listed below.\n\n")
        if not rows:
            f.write("(none)\n")
        for i, w in enumerate(rows, 1):
            f.write(f"[{i}] {w.get('message')}\n")
            f.write(
                f"    ORIGIN: {w.get('origin_file')} marker={w.get('origin_marker')} "
                f"@ {w.get('origin_pos')}\n"
            )
            f.write(
                f"    FOUND:  {w.get('target_file')} @ {w.get('found_pos')} "
                f"({w.get('side')}/{w.get('strand')})\n"
            )
            f.write(f"    FOREIGN MARKERS: {', '.join(w.get('foreign_markers') or [])}\n")
            f.write(f"    VERDICT: {w.get('verdict')} ({w.get('risk_score')}/100)\n")
            f.write("-" * 40 + "\n")
    logger.info(
        f"[+] Incorporation warnings: {path_json} ({len(rows)} unique contexts)"
    )
    try:
        draw_incorporation_grid(rows, matches)
    except Exception as grid_err:
        logger.error(f"[!] Incorporation grid failed: {grid_err}")
    return path_json

def write_triage_summary(out_dir, matches=None):
    """Ranked where-to-look report (not a diagnosis). TRIAGE_SUMMARY.json

    Primary score blends:
      - long forensic hit totals (motif len >= FORENSIC_EVIDENCE_MIN_LEN)
      - unique long-marker tags (diversity of lab signal)
      - Phase-2 flank support (files that appear in cross-match edges)
    Short forensic hits stay in short_forensic_hits_excluded and do not inflate rank.
    """
    from collections import defaultdict
    os.makedirs(out_dir, exist_ok=True)
    matches = matches or []
    comp_root = os.path.join(out_dir, "composition")
    min_len = forensic_evidence_min_len()
    stamp = build_config_stamp()

    file_lab = defaultdict(int)       # long motif totals
    file_short = defaultdict(int)
    file_tags = defaultdict(lambda: defaultdict(int))
    file_marker_set = defaultdict(set)

    if os.path.isdir(comp_root):
        for fname in os.listdir(comp_root):
            if not fname.endswith(".TARGET_COUNTS.json") and not fname.endswith(".composition.json"):
                continue
            path = os.path.join(comp_root, fname)
            try:
                with open(path, encoding="utf-8") as jf:
                    payload = json.load(jf)
            except Exception:
                continue
            subject = payload.get("file") or fname
            # normalize subject to basename for matching phase2
            subject_base = os.path.basename(str(subject))
            for mname, stats in (payload.get("by_marker") or {}).items():
                total = int(stats.get("total", 0)) if isinstance(stats, dict) else int(stats or 0)
                if total <= 0:
                    continue
                info = TARGETS.get(mname) or {}
                tag = info.get("tag", "UNKNOWN")
                role = str(info.get("role", "")).lower()
                if role == "composition" or is_composition_marker(mname):
                    continue
                mlen = len(str(info.get("seq", "") or ""))
                if mlen < min_len:
                    file_short[subject_base] += total
                    continue
                file_tags[subject_base][tag] += total
                file_lab[subject_base] += total
                file_marker_set[subject_base].add(mname)

    # Phase-2 flank support: count edges / shared_flank_hits per basename
    flank_support = defaultdict(int)
    flank_edge_count = defaultdict(int)
    sim_path = os.path.join(out_dir, "CROSS_SIMILARITY_TRACKER.json")
    if os.path.exists(sim_path):
        try:
            with open(sim_path, encoding="utf-8") as sf:
                sim = json.load(sf)
            for e in (sim.get("edges") or []):
                a = os.path.basename(str(e.get("file_a") or ""))
                b = os.path.basename(str(e.get("file_b") or ""))
                shared = int(e.get("shared_flank_hits") or e.get("similarity_score") or 0)
                if a:
                    flank_support[a] += shared
                    flank_edge_count[a] += 1
                if b:
                    flank_support[b] += shared
                    flank_edge_count[b] += 1
        except Exception:
            pass
    # also credit raw phase2 match list if present
    for m in matches:
        bn = os.path.basename(str(m.get("file") or m.get("basename") or m.get("target_file") or ""))
        if bn:
            flank_support[bn] += 1

    def triage_score(f):
        hits = int(file_lab.get(f, 0))
        ntags = len(file_tags.get(f) or {})
        fl = int(flank_support.get(f, 0))
        # hits dominate weakly; tag diversity and flank reuse boost "start here"
        return hits + (5 * ntags) + min(fl, 500)

    all_subjects = set(file_lab) | set(file_short) | set(flank_support)
    top_files = []
    for f in all_subjects:
        hits = int(file_lab.get(f, 0))
        short = int(file_short.get(f, 0))
        ntags = len(file_tags.get(f) or {})
        fl = int(flank_support.get(f, 0))
        edges = int(flank_edge_count.get(f, 0))
        score = triage_score(f) if hits or fl else 0
        row = {
            "file": f,
            "triage_score": score,
            "lab_marker_hits": hits,
            "unique_long_tags": ntags,
            "phase2_flank_support": fl,
            "phase2_edge_count": edges,
            "short_forensic_hits_excluded": short,
            "triage_min_len": min_len,
            "long_markers": sorted(file_marker_set.get(f) or []),
        }
        if hits == 0 and short > 0 and fl == 0:
            row["note"] = "short motifs only — excluded from primary rank"
            row["triage_score"] = 0
        top_files.append(row)

    top_files.sort(
        key=lambda r: (-int(r.get("triage_score") or 0), -int(r.get("lab_marker_hits") or 0), r["file"])
    )
    # keep ranked with score>0 first, then a few short-only notes
    primary = [r for r in top_files if int(r.get("triage_score") or 0) > 0][:50]
    short_only = [r for r in top_files if int(r.get("triage_score") or 0) == 0 and int(r.get("short_forensic_hits_excluded") or 0) > 0][:15]
    top_files = primary + short_only

    multi = []
    viral_tags = {"VIRAL_VECTOR", "VIRAL_ORIGIN", "VIRAL_EXPORT", "VACCINE_VECTOR"}
    bact_tags = {"CHASSIS_ORIGIN", "INDUSTRIAL_BACKBONE", "ANTIBIOTIC_RES", "LAB_MARKER"}
    tag_tags = {"PURIFICATION_TAG", "REPORTER_GENE"}
    for subject, tags in file_tags.items():
        has_v = any(tags.get(t, 0) > 0 for t in viral_tags)
        has_b = any(tags.get(t, 0) > 0 for t in bact_tags)
        has_t = any(tags.get(t, 0) > 0 for t in tag_tags)
        domains = sum([has_v, has_b, has_t])
        if domains >= 2:
            multi.append({
                "file": subject,
                "domains": domains,
                "viral": has_v,
                "bacterial_or_lab": has_b,
                "tag_or_reporter": has_t,
                "by_tag": dict(tags),
                "triage_score": triage_score(subject),
            })
    multi.sort(key=lambda r: (-r["domains"], -r.get("triage_score", 0), r["file"]))

    top_pairs = []
    if os.path.exists(sim_path):
        try:
            with open(sim_path, encoding="utf-8") as sf:
                sim = json.load(sf)
            edges = sorted(
                sim.get("edges") or [],
                key=lambda e: -int(e.get("similarity_score") or 0),
            )[:40]
            for e in edges:
                top_pairs.append({
                    "file_a": e.get("file_a"),
                    "file_b": e.get("file_b"),
                    "similarity_score": e.get("similarity_score"),
                    "shared_flank_hits": e.get("shared_flank_hits"),
                    "pair_support_count": e.get("pair_support_count"),
                    "same_accession": e.get("same_accession", False),
                    "markers": e.get("markers") or [],
                })
        except Exception:
            pass

    payload = {
        "description": (
            "Triage ranking only — points researchers at files/pairs to inspect. "
            "Not a medical or legal determination of origin or disease. "
            f"triage_score = long_hits(>={min_len}bp) + 5*unique_long_tags + "
            "capped Phase-2 flank support. Short motifs listed but score 0."
        ),
        "config_stamp": stamp,
        "triage_min_len": min_len,
        "scoring": {
            "primary": "triage_score = lab_marker_hits + 5*unique_long_tags + min(phase2_flank_support, 500)",
            "lab_marker_hits": "non-composition forensic totals with len(seq) >= triage_min_len",
            "excluded": "short_forensic_hits_excluded = same class with len(seq) < triage_min_len",
            "config_key": "FORENSIC_EVIDENCE_MIN_LEN",
        },
        "top_files_by_lab_markers": top_files,
        "top_file_pairs_by_flank_score": top_pairs,
        "multi_domain_files": multi[:40],
        "n_phase2_hits": len(matches),
    }
    out_path = os.path.join(out_dir, "TRIAGE_SUMMARY.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    # also write stamp sidecar next to triage for easy packaging
    try:
        with open(os.path.join(out_dir, "RUN_CONFIG_STAMP.json"), "w", encoding="utf-8") as sf:
            json.dump(stamp, sf, indent=2)
    except Exception:
        pass
    logger.info(f"[+] Triage summary: {out_path}")
    return out_path


def write_coordinate_export(out_dir=None):
    """BED tracks for long forensic hits from Evidence TXT (composition skipped)."""
    if str(Constants.CONFIG.get("ENABLE_COORDINATE_EXPORT", "1")) != "1":
        return
    if out_dir is None:
        out_dir = Constants.DIRS.get("OUT_LOGS", Constants.BASE_OUT)
    bed_root = os.path.join(out_dir, "coordinates")
    os.makedirs(bed_root, exist_ok=True)
    min_len = int(Constants.CONFIG.get("MIN_EVIDENCE_LEN", Constants.CONFIG.get("FORENSIC_EVIDENCE_MIN_LEN", 24)))
    from collections import defaultdict
    by_file = defaultdict(list)
    for key, path in Constants.DIRS.items():
        m = re.fullmatch(r"OUT_(CLN|IND|REF)_EVI_(FWD|REV)", key)
        if not m or not os.path.isdir(path):
            continue
        strand = m.group(2)
        strand_bed = "+" if strand == "FWD" else "-"
        for root, _, files in os.walk(path):
            for fname in files:
                if "_Evidence_" not in fname or not fname.endswith(".txt"):
                    continue
                after = fname.split("_Evidence_", 1)[-1]
                parts = after.split("_", 2)
                subject = parts[2][:-4] if len(parts) >= 3 and parts[2].endswith(".txt") else fname
                fpath = os.path.join(root, fname)
                try:
                    with open(fpath, encoding="utf-8", errors="replace") as ef:
                        content = ef.read()
                except Exception:
                    continue
                for block in content.split("-" * 60):
                    tm = re.search(r"GENOMIC HIT DETECTED:\s+([A-Za-z0-9_-]+)", block)
                    pm = re.search(r"LOCATION:\s+(\d+)", block)
                    if not tm or not pm:
                        continue
                    mname = tm.group(1).strip()
                    if is_composition_marker(mname):
                        continue
                    info = TARGETS.get(mname) or {}
                    mlen = len(str(info.get("seq", "") or "")) or 1
                    if mlen < min_len:
                        continue
                    start = int(pm.group(1))
                    end = start + mlen
                    by_file[subject].append(
                        f"{subject}\t{start}\t{end}\t{mname}\t0\t{strand_bed}\n"
                    )
    n_files = 0
    for subject, lines in by_file.items():
        if not lines:
            continue
        safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in subject)
        bed_path = os.path.join(bed_root, f"{safe}.bed")
        with open(bed_path, "w", encoding="utf-8") as bf:
            bf.write("# BED: long forensic hits only (composition motifs excluded)\n")
            bf.writelines(sorted(set(lines), key=lambda x: int(x.split("\t")[1])))
        n_files += 1
    if n_files:
        logger.info(f"[+] Coordinate BED export: {n_files} files under {bed_root}")

def write_cross_similarity_tracker(matches, out_dir):
    """
    Aggregates Phase 2 flank hits into a file-pair similarity graph.
    pair_support = both LEFT+RIGHT of same origin hit in same target.
    same_accession = shared NCBI accession token (duplicate sequence warning).
    """
    from collections import defaultdict

    if not matches:
        return

    edge_count = defaultdict(int)
    edge_pair_support = defaultdict(int)
    edge_labels = defaultdict(lambda: defaultdict(int))
    edge_markers = defaultdict(set)
    sides_seen = defaultdict(set)

    for m in matches:
        origin = m.get("origin_file") or "UNKNOWN"
        target = m.get("target_file") or "UNKNOWN"
        if origin == target:
            continue
        key = tuple(sorted((origin, target)))
        edge_count[key] += 1
        edge_labels[key][m.get("label", "FLANK_CROSS")] += 1
        if m.get("marker"):
            edge_markers[key].add(m["marker"])
        oh = m.get("origin_hit_id")
        if oh:
            sides_seen[(oh, target)].add(m.get("side"))

    for (oh, target), sides in sides_seen.items():
        if "LEFT" in sides and "RIGHT" in sides:
            origin = str(oh).split("|", 1)[0]
            key = tuple(sorted((origin, target)))
            edge_pair_support[key] += 1

    rows = []
    for key, n in sorted(edge_count.items(), key=lambda kv: (-kv[1], kv[0])):
        a, b = key
        pair_n = edge_pair_support.get(key, 0)
        acc_a, acc_b = extract_accession(a), extract_accession(b)
        same_acc = bool(acc_a and acc_b and acc_a == acc_b)
        score = n + 8 * pair_n
        top_label = max(edge_labels[key], key=edge_labels[key].get) if edge_labels[key] else ""
        rows.append({
            "file_a": a,
            "file_b": b,
            "shared_flank_hits": n,
            "pair_support_count": pair_n,
            "similarity_score": score,
            "labels": dict(edge_labels[key]),
            "markers": sorted(edge_markers[key]),
            "same_accession": same_acc,
            "accession": acc_a if same_acc else "",
            "top_label": top_label,
            "note": (
                "DUPLICATE_SEQUENCE_WARNING: same accession in both filenames — "
                "not independent biological cross-link"
                if same_acc else ""
            ),
        })

    os.makedirs(out_dir, exist_ok=True)
    json_path = os.path.join(out_dir, "CROSS_SIMILARITY_TRACKER.json")
    csv_path = os.path.join(out_dir, "CROSS_SIMILARITY_TRACKER.csv")

    payload = {
        "description": (
            "Cross-file flank reuse graph from Phase 2 exact Aho-Corasick hits. "
            "Not whole-genome identity. same_accession=true means shared accession "
            "token (e.g. CLN/IND copy of same GenBank record) — treat as duplicate, "
            "not independent evidence."
        ),
        "n_edges": len(rows),
        "n_raw_hits": len(matches),
        "edges": rows,
    }
    with open(json_path, "w", encoding="utf-8") as jf:
        json.dump(payload, jf, indent=2)

    with open(csv_path, "w", newline="", encoding="utf-8") as cf:
        w = csv.writer(cf)
        w.writerow([
            "file_a", "file_b", "shared_flank_hits", "pair_support_count",
            "similarity_score", "markers", "top_label", "same_accession", "accession", "note",
        ])
        for r in rows:
            w.writerow([
                r["file_a"], r["file_b"], r["shared_flank_hits"], r["pair_support_count"],
                r["similarity_score"], ";".join(r["markers"]), r.get("top_label", ""),
                r.get("same_accession", False), r.get("accession", ""), r.get("note", ""),
            ])

    logger.info(f"[+] Cross-similarity tracker: {json_path} ({len(rows)} file pairs)")

def draw_flank_cross_heatmaps(all_cross):
    """
    Phase 2 only: for each file WHERE a foreign LEFT/RIGHT flank was found,
    plot those hit positions. Label = origin marker + side + short origin file.
    """
    if not all_cross:
        logger.info("[*] Flank-cross heatmaps: no Phase 2 hits.")
        return

    from collections import defaultdict

    # target_file (where flank was FOUND) -> list of draw_map rows
    by_found = defaultdict(list)
    for m in all_cross:
        found = m.get("target_file")
        pos = m.get("found_pos")
        if not found or pos is None:
            continue
        origin = os.path.basename(str(m.get("origin_file") or "?"))
        # short label so the bottom strip stays readable
        if len(origin) > 28:
            origin = origin[:25] + "..."
        label = f"{m.get('marker', '?')}|{m.get('side', '?')}|{origin}"
        strand = m.get("match_strand") or "FWD"
        color = (TARGETS.get(m.get("marker") or "") or {}).get("color", "#00CED1")
        by_found[found].append({
            "pos": int(pos),
            "name": label,
            "strand": strand,
            "color": color,
        })

    # map basename -> category key for OUT_*_MAP_FLANKS
    file_to_prefix = {}
    for key, folder in (("CLN", "OUT_CLN"), ("IND", "OUT_IND"), ("REF", "OUT_REF")):
        d = Constants.DIRS.get(key)
        if not d or not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            file_to_prefix[f] = folder

    for found_name, hits in by_found.items():
        # scope_id may be "OUT_CLN/foo.fasta" — split category from basename
        raw = str(found_name or "")
        if "/" in raw or "\\" in raw:
            parts = raw.replace("\\", "/").split("/")
            cat_hint = parts[0] if parts[0].startswith("OUT_") else None
            bare = parts[-1]
        else:
            cat_hint = None
            bare = raw

        prefix = cat_hint or file_to_prefix.get(bare) or file_to_prefix.get(raw)
        if not prefix:
            map_dir = os.path.join(Constants.DIRS.get("OUT_MAPS", Constants.BASE_OUT), "FLANKS")
        else:
            map_dir = Constants.DIRS.get(f"{prefix}_MAP_FLANKS")
            if not map_dir:
                map_dir = os.path.join(Constants.DIRS.get("OUT_MAPS", Constants.BASE_OUT), "FLANKS")
        os.makedirs(map_dir, exist_ok=True)

        max_pos = max(h["pos"] for h in hits) + 500
        seq_len = max_pos
        comp_dir = os.path.join(Constants.DIRS.get("OUT_LOGS", Constants.BASE_OUT), "composition")
        safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in bare)
        for suffix in (f"{safe}.TARGET_COUNTS.json", f"{bare}.TARGET_COUNTS.json"):
            jp = os.path.join(comp_dir, suffix)
            if os.path.isfile(jp):
                try:
                    with open(jp, encoding="utf-8") as jf:
                        seq_len = max(int(json.load(jf).get("seq_len") or 0), max_pos)
                except Exception:
                    pass
                break

        # Never put path separators in the PNG filename
        out_png = os.path.join(map_dir, f"{safe}_FLANKS_HEATMAP.png")
        draw_map(
            hits,
            seq_len,
            out_png,
            f"Flank cross-hits on {bare} (from other files)",
        )
        logger.info(f"[+] Flank-cross heatmap: {out_png} ({len(hits)} hits)")

def compare_flanks(phase1_harvest):
    if not phase1_harvest:
        logger.warning("Phase 2: No harvest data received. Skipping.")
        return
        
    logger.info(f"[*] Analyzing {len(phase1_harvest)} hits for orientation and synthetic tags...")
    
    if float(Constants.CONFIG.get('MATCH_THRESHOLD', 100)) < 100:
        logger.warning("MATCH_THRESHOLD < 100 has no effect in Phase 2: flank cross-matching is exact-only.")

    # --- 2. PRESERVED CROSS-MATCH LOGIC ---
    setup_directories()
    # Written into Global_Evidence so it surfaces in the UI's audit log list
    report_path = os.path.join(Constants.DIRS["OUT_LOGS"], "CROSS_MATCH_LOG.txt")
    if os.path.exists(report_path): os.remove(report_path)

    flank_library = {}
    hit_id = 0
    for result in phase1_harvest:
        source_file = result.get("file", "UNKNOWN")
        for flk in result.get("flanks", []):
            marker_name = flk.get("marker")
            origin_pos = flk.get("pos")
            source_file = flk.get("file") or result.get("file") or "UNKNOWN"
            hit_id += 1
            origin_hit_id = f"{source_file}|{marker_name}|{origin_pos}|{hit_id}"
            for side, seq in (("LEFT", flk.get("left", "")), ("RIGHT", flk.get("right", ""))):
                seq = str(seq).upper().replace("U", "T")
                if len(seq) < MIN_FLANK_LEN or is_low_complexity(seq):
                    continue
                lib_key = f"{origin_hit_id}|{side}"
                flank_library[lib_key] = {
                    "seq": seq,
                    "marker": marker_name,
                    "side": side,
                    "origin_file": source_file,
                    "origin_pos": origin_pos,
                    "origin_hit_id": origin_hit_id,
                }

    if not flank_library: return

    # Recursive Search in simplified folders
    all_files = list_input_files()

    num_cores = int(Constants.CONFIG.get('MAX_CORES', 2))
    task_queue = [(f, None, Constants.DIRS, Constants.CONFIG) for f in all_files]

    all_cross = []
    with tqdm(
        total=len(all_files),
        desc="[CROSS-AUDIT]",
        leave=True,
        dynamic_ncols=True,
        mininterval=0.3,
        colour=get_random_hex(),
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}]",
    ) as pbar:
        with mp.Pool(processes=num_cores, initializer=_p2_init, initargs=(flank_library,)) as pool:
            for worker_result in pool.imap_unordered(phase2_worker, task_queue):
                if worker_result and isinstance(worker_result, list):
                    all_cross.extend(worker_result)
                    with open(report_path, "a") as master_log:
                        for m in worker_result:
                            warn = m.get("warning") or ""
                            foreign = m.get("foreign_markers") or []
                            master_log.write(
                                f"[!] CROSS-LINK: {m['seq']}\n"
                                f"LABEL: {m.get('label', 'FLANK_CROSS')}{' | RC' if m.get('rc_hit') else ''}\n"
                                f"MARKER: {m['marker']} | SIDE: {m['side']} | STRAND: {m['match_strand']}\n"
                                f"IDENTITY: {m['similarity']:.1f}%\n"
                                f"FOUND: {m['target_file']} @ {m.get('found_pos', '?')}\n"
                                f"ORIGIN: {m['origin_file']} @ {m.get('origin_pos', '?')}\n"
                                f"VERDICT: {m['verdict']} (Risk Score: {m['risk_score']}/100)\n"
                                f"NEARBY CONTEXT: {', '.join(m['context'])}\n"
                            )
                            if foreign:
                                master_log.write(
                                    f"FOREIGN MARKERS NEXT TO FLANK: {', '.join(foreign)}\n"
                                )
                            if warn:
                                master_log.write(f"[!] {warn}\n")
                            master_log.write(f"{'='*40}\n")
                pbar.update(1)

    pair_window = int(Constants.CONFIG.get("PAIR_WINDOW_BP", 5000))
    annotate_flank_pairs(all_cross, pair_window_bp=pair_window)
    write_cross_similarity_tracker(all_cross, Constants.DIRS["OUT_LOGS"])
    try:
        write_incorporation_warnings(all_cross, Constants.DIRS["OUT_LOGS"])
    except Exception as inc_err:
        logger.error(f"[!] Incorporation warnings failed: {inc_err}")
    try:
        write_triage_summary(Constants.DIRS["OUT_LOGS"], all_cross)
    except Exception as triage_err:
        logger.error(f"[!] Triage summary failed: {triage_err}")
    try:
        write_coordinate_export(Constants.DIRS["OUT_LOGS"])
    except Exception as bed_err:
        logger.error(f"[!] Coordinate export failed: {bed_err}")
    try:
        draw_flank_cross_heatmaps(all_cross)
    except Exception as e:
        logger.error(f"[!] Flank-cross heatmaps failed: {e}")

def determine_verdict(marker_context, sequence_data=None):
    """
    Universal Forensic Verdict Engine.
    Evaluates genetic risk by weighing cumulative markers and codon bias signatures.
    """
    # 1. Normalize input (Handles single marker name or list of neighboring markers)
    if isinstance(marker_context, list):
        log_content = "|".join(marker_context)
    else:
        log_content = str(marker_context)
    
    score = 0
    reasons = []

    # --- 1. HIGH THREAT / GAIN-OF-FUNCTION (GOF) ---
    if "Furin_Cleavage" in log_content: 
        score += 60
        reasons.append("High-consequence cleavage site (GOF Signal)")
    if "NLS_Signal" in log_content:     
        score += 40
        reasons.append("Nuclear Localization Signal (Viral/Human delivery)")
    if "CRISPR_Target" in log_content:  
        score += 35
        reasons.append("Engineered CRISPR-Cas targeting site")
    
    # --- 2. VIRAL VECTORS & MOBILIZATION ---
    viral_markers = ["Lenti_LTR", "RRE_Element", "WPRE_Element", "SV40_Origin", "SV40_Promoter"]
    found_viral = [m for m in viral_markers if m in log_content]
    if found_viral:
        score += (len(found_viral) * 20)
        reasons.append(f"Viral mobilization machinery: {', '.join(found_viral)}")

    # --- 3. ASSEMBLY SCARS & LINKERS (THE "SMOKING GUNS") ---
    scars = ["BsaI_Scar", "BsmBI_Scar", "NotI_Site", "BamHI_Site"]
    found_scars = [s for s in scars if s in log_content]
    if found_scars:
        score += (len(found_scars) * 15)
        reasons.append(f"Restriction/Cloning scars: {', '.join(found_scars)}")

    if "Gibson_Overlap" in log_content: score += 20
    if "attB1_Gateway" in log_content:  score += 15
    if "27nt_Bridge" in log_content or "Engineered_Linker" in log_content:
        score += 15
        reasons.append("Synthetic linker/bridge detected")

    # --- 4. CHASSIS & SELECTION MARKERS ---
    industrial = ["pET_Origin", "pUC_Origin", "KanR_Marker", "AmpR_Marker", "T7_Promoter"]
    found_ind = [m for m in industrial if m in log_content]
    if found_ind:
        score += (len(found_ind) * 12)
        reasons.append(f"Industrial/Lab chassis components: {', '.join(found_ind)}")

    # --- 5. CODON BIAS CROSS-REFERENCE ---
    # If we have sequence data, we verify if the host match supports the markers
    if sequence_data and score > 0:
        bias_verdict = analyze_author_signature(sequence_data)
        if "INDUSTRIAL_OPTIMIZATION" in bias_verdict:
            score += 20
            reasons.append("Host codon usage is consistent with an industrial chassis (supporting only, not proof)")

    # --- 6. CHIMERIC COMPLEXITY BONUS ---
    unique_domains = set()
    if any(x in log_content for x in ["Lenti", "SV40", "WPRE"]): unique_domains.add("VIRAL")
    if any(x in log_content for x in ["pET", "pUC", "KanR", "AmpR"]): unique_domains.add("BACTERIAL")
    if any(x in log_content for x in ["NLS", "Furin"]): unique_domains.add("MAMMALIAN")
    
    if len(unique_domains) >= 2:
        score += 30 
        reasons.append(f"Interspecies Chimerism detected: {list(unique_domains)}")

    # --- FINAL CATEGORIZATION ---
    if score >= 90:
        label = "CRITICAL_RISK_GOF_UPGRADE"
    elif score >= 60:
        label = "HIGH_CONFIDENCE_SYNTHETIC_CONSTRUCT"
    elif score >= 30:
        label = "POTENTIAL_LAB_ORIGIN"
    elif score > 0:
        label = "LOW_LEVEL_MARKER_MATCH"
    else:
        label = "NO_SIGNIFICANT_ENGINEERING_DETECTED"

    return {
        "verdict": label,
        "risk_score": min(score, 100), # Cap score at 100 for reporting
        "reasoning": "; ".join(reasons) if reasons else "Natural genomic baseline."
    }

# ---------------------------------------------------------------------------
# WATERFALL / MULTI-FILE COLOUR-BAND VIEW
# ---------------------------------------------------------------------------
CATEGORY_COLORS = {
    "RESTRICTION_SCAR": "#1E90FF",
    "SEAMLESS_ASSEMBLY": "#800080",
    "GATEWAY_CLONING": "#9932CC",
    "SYNTHETIC_LINKER": "#DA70D6",
    "MOBILE_ELEMENT": "#BA55D3",
    "LAB_MARKER": "#228B22",
    "VACCINE_VECTOR": "#006400",
    "WRC_SENSITIVE": "#FF0000",
    "GOF_UPGRADE": "#8B0000",
    "GENETIC_MOD": "#FF4500",
    "VIRAL_VECTOR": "#B22222",
    "REPORTER_GENE": "#FF69B4",
    "PURIFICATION_TAG": "#C71585",
    "INDUSTRIAL_BACKBONE": "#FF4500",
    "CHASSIS_ORIGIN": "#FF8C00",
    "VIRAL_ORIGIN": "#DAA520",
    "ANTIBIOTIC_RES": "#696969",
    "RNA_STABILITY": "#008080",
    "TRANSLATION_INIT": "#20B2AA",
    "VIRAL_EXPORT": "#48D1CC",
    "RNA_SCAFFOLD": "#00CED1",
    "TEST_VALIDATION": "#FFD700",
    "HUMAN_EDIT": "#FF6347",
    "EPIGENETIC_SITE": "#9370DB",
    "UNKNOWN": "#808080",
}


def build_waterfall_matrix(evidence_root=None, max_files=500, include_comp=None):
    """rows = files, columns = marker-tag categories, values = hit counts.
    Forensic hits from evidence logs; composition totals from *.composition.json.
    """
    if evidence_root is None:
        evidence_root = Constants.BASE_OUT
    from collections import defaultdict, OrderedDict
    file_hits = OrderedDict()

    def ensure_row(subject):
        if subject not in file_hits:
            if len(file_hits) >= max_files:
                return False
            file_hits[subject] = defaultdict(int)
        return True

    # Every hit (forensic and composition) is already in the per-file sidecar by tag.
    # --- B) Composition JSON sidecars (accurate totals, no spam) ---
    comp_roots = []
    for _r in (
        os.path.join(Constants.DIRS.get("OUT_LOGS", evidence_root), "composition"),
        os.path.join(evidence_root, "Global_Evidence", "composition"),
        os.path.join(evidence_root, "composition"),
    ):
        _r = os.path.abspath(_r)
        if _r not in comp_roots:
            comp_roots.append(_r)
    for comp_root in comp_roots:
        if not os.path.isdir(comp_root):
            continue
        for fname in os.listdir(comp_root):
            if not (fname.endswith(".composition.json") or fname.endswith(".TARGET_COUNTS.json")):
                continue
            path = os.path.join(comp_root, fname)
            try:
                with open(path) as jf:
                    payload = json.load(jf)
                subject = (
                    payload.get("file")
                    or fname.replace(".composition.json", "")
                    .replace(".TARGET_COUNTS.json", "")
                )
                if not ensure_row(subject):
                    continue
                if include_comp is None:
                    include_comp = str(Constants.CONFIG.get("WATERFALL_INCLUDE_COMPOSITION", "0")) == "1"
                # Prefer per-marker tallies so composition can be stripped from the decision matrix
                by_marker = payload.get("by_marker") or {}
                if by_marker:
                    for mname, stats in by_marker.items():
                        if isinstance(stats, dict):
                            total = int(stats.get("total", 0) or 0)
                            is_comp = stats.get("composition")
                            if is_comp is None:
                                is_comp = is_composition_marker(mname)
                            tag = stats.get("tag") or (TARGETS.get(mname) or {}).get("tag", "UNKNOWN")
                        else:
                            total = int(stats or 0)
                            is_comp = is_composition_marker(mname)
                            tag = (TARGETS.get(mname) or {}).get("tag", "UNKNOWN")
                        if not include_comp and is_comp:
                            continue
                        if total <= 0:
                            continue
                        file_hits[subject][tag] += total
                elif include_comp:
                    for tag, total in (payload.get("by_tag") or {}).items():
                        file_hits[subject][tag] += int(total)
            except Exception:
                continue

    categories = sorted({t for hits in file_hits.values() for t in hits})
    if not categories:
        categories = list(CATEGORY_COLORS.keys())
    matrix = []
    names = list(file_hits.keys())
    for subject in names:
        matrix.append([file_hits[subject].get(cat, 0) for cat in categories])
    colors = [CATEGORY_COLORS.get(c, "#808080") for c in categories]
    return names, categories, matrix, colors


def draw_waterfall(file_names, categories, matrix, colors, save_path, title="Forensic Waterfall"):
    if not file_names or not categories:
        return
    n_files = len(file_names)
    n_cats = len(categories)
    row_h = max(0.12, min(0.35, 40.0 / max(n_files, 1)))
    fig_h = max(4, n_files * row_h + 1.5)
    fig_w = max(10, n_cats * 0.9 + 4)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_facecolor("#111111")
    fig.patch.set_facecolor("#111111")

    arr = np.array(matrix, dtype=float)
    max_c = arr.max() if arr.size and arr.max() > 0 else 1.0

    for i, row in enumerate(matrix):
        for j, count in enumerate(row):
            if count <= 0:
                continue
            alpha = 0.25 + 0.75 * (count / max_c)
            ax.add_patch(plt.Rectangle(
                (j, n_files - i - 1), 1, 1,
                facecolor=colors[j], alpha=alpha,
                edgecolor="#222222", linewidth=0.3,
            ))

    ax.set_xlim(0, n_cats)
    ax.set_ylim(0, n_files)
    ax.set_xticks([x + 0.5 for x in range(n_cats)])
    ax.set_xticklabels(categories, rotation=55, ha="right", color="white", fontsize=7)
    step = max(1, n_files // 40)
    yticks = list(range(0, n_files, step))
    ax.set_yticks([n_files - i - 0.5 for i in yticks])
    ax.set_yticklabels([file_names[i][:40] for i in yticks], color="white", fontsize=6)
    ax.set_title(title, color="white", fontsize=12, pad=8)
    ax.tick_params(colors="gray")
    for spine in ax.spines.values():
        spine.set_color("#333333")
    plt.tight_layout()
    os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
    plt.savefig(save_path, dpi=150, facecolor="#111111", bbox_inches="tight")
    plt.close()
    logger.info(f"[+] Waterfall written: {save_path} ({n_files} files x {n_cats} categories)")


def clean_previous_run():
    """Outputs are written in append mode, so last run's files must be moved away first
    (otherwise every re-run duplicates hits). They are archived, not deleted, to keep
    chain of custody; set ARCHIVE_PREVIOUS_RUNS=0 to delete instead."""
    keep = str(Constants.CONFIG.get("ARCHIVE_PREVIOUS_RUNS", "1")) == "1"
    stale = [d for k, d in Constants.DIRS.items()
             if k.startswith("OUT_") and k not in ("OUT_LOGS", "OUT_MAPS")]
    stale += [
        os.path.join(Constants.DIRS["OUT_LOGS"], "composition"),
        os.path.join(Constants.DIRS["OUT_LOGS"], "hits"),
        os.path.join(Constants.DIRS["OUT_LOGS"], "CROSS_MATCH_LOG.txt"),
        os.path.join(Constants.DIRS["OUT_LOGS"], "CROSS_SIMILARITY_TRACKER.json"),
        os.path.join(Constants.DIRS["OUT_LOGS"], "CROSS_SIMILARITY_TRACKER.csv"),
        os.path.join(Constants.DIRS["OUT_LOGS"], "FLANK_INCORPORATION_WARNINGS.json"),
        os.path.join(Constants.DIRS["OUT_LOGS"], "FLANK_INCORPORATION_WARNINGS.txt"),
        os.path.join(Constants.DIRS["OUT_LOGS"], "INCORPORATION_GRID.csv"),
        os.path.join(Constants.DIRS["OUT_LOGS"], "TRIAGE_SUMMARY.json"),
        os.path.join(Constants.DIRS["OUT_LOGS"], "coordinates"),
        os.path.join(Constants.DIRS["OUT_LOGS"], "ref_distance"),
        os.path.join(Constants.DIRS["OUT_MAPS"], "WATERFALL_SUMMARY.png"),
        os.path.join(Constants.DIRS["OUT_MAPS"], "INCORPORATION_GRID.png"),
        os.path.join(Constants.DIRS["OUT_LOGS"], "portal_snapshot.json"),
        Constants.DIRS["REPORT"],
    ]
    stamp = time.strftime("%Y%m%d_%H%M%S")
    for path in stale:
        if os.path.isdir(path):
            if not any(files for _, _, files in os.walk(path)):
                continue
        elif not os.path.isfile(path):
            continue
        if keep:
            arch_root = getattr(Constants, "AUDIT_LIBRARY", None) or os.path.join(
                getattr(Constants, "ROOT", os.path.dirname(Constants.BASE_OUT)), "Audit_Library")
            dest = os.path.join(arch_root, "_working_archive", stamp, os.path.relpath(path, Constants.BASE_OUT))
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            shutil.move(path, dest)
        elif os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
        else:
            os.remove(path)


def export_run_bundle():
    """
    Copy the publishable snapshot set into <ROOT>/Audit_Library/<stamp>/
    (sibling of Forensic_Results). Deleting Forensic_Results does not wipe history.
    Zip any stamp folder + index.html for a static portal view.
    """
    import shutil
    stamp = time.strftime("%Y%m%d_%H%M%S")
    runs_root = getattr(Constants, "AUDIT_LIBRARY", None) or os.path.join(
        getattr(Constants, "ROOT", os.path.dirname(Constants.BASE_OUT)), "Audit_Library"
    )
    dest = os.path.join(runs_root, stamp)
    os.makedirs(dest, exist_ok=True)

    logs = Constants.DIRS.get("OUT_LOGS", os.path.join(Constants.BASE_OUT, "Global_Evidence"))
    maps = Constants.DIRS.get("OUT_MAPS", os.path.join(Constants.BASE_OUT, "Global_Heatmaps"))

    copies = [
        (os.path.join(logs, "portal_snapshot.json"), "portal_snapshot.json"),
        (os.path.join(logs, "TRIAGE_SUMMARY.json"), "TRIAGE_SUMMARY.json"),
        (os.path.join(logs, "RUN_CONFIG_STAMP.json"), "RUN_CONFIG_STAMP.json"),
        (os.path.join(logs, "CROSS_SIMILARITY_TRACKER.json"), "CROSS_SIMILARITY_TRACKER.json"),
        (os.path.join(logs, "CROSS_SIMILARITY_TRACKER.csv"), "CROSS_SIMILARITY_TRACKER.csv"),
        (os.path.join(logs, "INCORPORATION_GRID.csv"), "INCORPORATION_GRID.csv"),
        (os.path.join(logs, "FLANK_INCORPORATION_WARNINGS.json"), "FLANK_INCORPORATION_WARNINGS.json"),
        (os.path.join(logs, "FLANK_INCORPORATION_WARNINGS.txt"), "FLANK_INCORPORATION_WARNINGS.txt"),
        (os.path.join(maps, "WATERFALL_SUMMARY.png"), "WATERFALL_SUMMARY.png"),
        (os.path.join(maps, "INCORPORATION_GRID.png"), "INCORPORATION_GRID.png"),
        (Constants.DIRS.get("REPORT") or os.path.join(Constants.BASE_OUT, "FINAL_GENOMIC_REPORT.xlsx"),
         "FINAL_GENOMIC_REPORT.xlsx"),
    ]
    written = []
    for src, name in copies:
        if src and os.path.isfile(src):
            shutil.copy2(src, os.path.join(dest, name))
            written.append(name)

    # Manifest for index.html run picker
    manifest_path = os.path.join(runs_root, "index.json")
    entries = []
    if os.path.isfile(manifest_path):
        try:
            with open(manifest_path, encoding="utf-8") as mf:
                entries = json.load(mf).get("runs") or []
        except Exception:
            entries = []
    entries = [e for e in entries if e.get("id") != stamp]
    try:
        _cfg_stamp = build_config_stamp()
    except Exception:
        _cfg_stamp = {}
    # always write stamp into package even if missing from logs
    try:
        with open(os.path.join(dest, "RUN_CONFIG_STAMP.json"), "w", encoding="utf-8") as sf:
            json.dump(_cfg_stamp, sf, indent=2)
        if "RUN_CONFIG_STAMP.json" not in written:
            written.append("RUN_CONFIG_STAMP.json")
    except Exception:
        pass
    entries.insert(0, {
        "id": stamp,
        "path": f"Audit_Library/{stamp}",
        "snapshot": f"Audit_Library/{stamp}/portal_snapshot.json",
        "files": written,
        "created": stamp,
        "config_stamp": _cfg_stamp,
    })
    # keep last 30
    entries = entries[:30]
    with open(manifest_path, "w", encoding="utf-8") as mf:
        json.dump({"runs": entries, "latest": stamp}, mf, indent=2)

    # Also mirror latest snapshot next to index for simple hosting
    latest_snap = os.path.join(logs, "portal_snapshot.json")
    if os.path.isfile(latest_snap):
        try:
            root_copy = os.path.join(Constants.ROOT if hasattr(Constants, "ROOT") else os.path.dirname(Constants.BASE_OUT), "portal_snapshot.json")
            # prefer project root = parent of Forensic_Results
            root_copy = os.path.join(os.path.dirname(Constants.BASE_OUT), "portal_snapshot.json")
            shutil.copy2(latest_snap, root_copy)
        except Exception:
            pass

    
    # Maintain Audit_Library/latest/ for index.html without query string
    try:
        latest_dir = os.path.join(runs_root, "latest")
        if os.path.isdir(latest_dir):
            shutil.rmtree(latest_dir, ignore_errors=True)
        shutil.copytree(dest, latest_dir)
    except Exception as _lat_err:
        logger.debug(f"latest package copy skipped: {_lat_err}")

    logger.info(f"[+] Audit package: {dest} ({len(written)} files) · manifest {manifest_path}")
    print(f"[+] Run snapshot folder: {dest}")
    return dest



def run_audit():
    global TARGETS, MARKER_COLORS, AUTO, DIRS

    # --- STEP 1: ENVIRONMENT, FRESH OUTPUT FOLDERS, TARGETS ---
    ensure_environment()
    DIRS = Constants.DIRS
    clean_previous_run()
    setup_directories()

    if not initialize_globals() or not TARGETS:
        raise RuntimeError("Targets.json is missing, corrupted, or has no usable markers.")

    # --- STEP 2: LEGEND (after the folders exist) ---
    print("[*] Generating Legend PNG...")
    check_and_generate_legend(MARKER_COLORS)

    if str(Constants.CONFIG.get("SYNC_ON_START", "0")) == "1":
        sync_global_threats()
    disk_guard(min_gb=1.0)
    
    # --- STEP 4: FILE COLLECTION ---
    all_files = list_input_files()
    
    if not all_files: 
        print("[!] No files found to audit. Please check your source directories.")
        return

    # --- STEP 5: PHASE 1 - SCANNING & HEATMAP GENERATION ---
    config = Constants.CONFIG
    num_cores = int(config.get('MAX_CORES', 2))
    
    # Pack the task queue for the workers
    task_queue = [(f, i, DIRS) for i, f in enumerate(all_files)]
    phase1_harvest = []
    
    num_cores = int(Constants.CONFIG.get("MAX_CORES", 1))

    with tqdm(
        total=len(all_files),
        desc="[MASTER AUDIT]",
        leave=True,
        dynamic_ncols=True,
        mininterval=0.3,
        colour=get_random_hex(),
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}]",
    ) as pbar:
        pool = None
        try:
            # maxtasksperchild recycles workers so one fat file cannot pin RAM forever
            pool = mp.Pool(processes=num_cores, maxtasksperchild=4)
            for worker_output in pool.imap_unordered(index_worker, task_queue):
                if worker_output:
                    file_name = worker_output.get("file", "Unknown")
                    if not worker_output.get("skipped"):
                        tqdm.write(f"[+] Audit Complete for {file_name}")
                    if worker_output.get("flanks"):
                        phase1_harvest.append(worker_output)
                pbar.update(1)
        except MemoryError:
            tqdm.write("[!] MemoryError in Phase 1 — stopping remaining workers. "
                       "Lower MAX_CORES, raise COMPOSITION_AC_MIN_LEN, set COMPOSITION_FAST_ONLY=1, "
                       "or reduce LARGE_FILE_MAX_HITS.")
            if pool is not None:
                try:
                    pool.terminate()
                    pool.join()
                except Exception:
                    pass
                pool = None
            raise
        except KeyboardInterrupt:
            tqdm.write("[!] Interrupted — terminating worker pool…")
            if pool is not None:
                try:
                    pool.terminate()
                    pool.join()
                except Exception:
                    pass
                pool = None
            raise
        finally:
            if pool is not None:
                try:
                    pool.close()
                    pool.join()
                except Exception:
                    try:
                        pool.terminate()
                        pool.join()
                    except Exception:
                        pass

    # --- STEP 6: PHASE 2 - CROSS-MATCHING / FLANK ANALYSIS ---
    # This is where the engine looks for synthetic adapter scars around the hits
    if phase1_harvest:
        print("\n" + "="*60)
        print("[*] Phase 1 Complete. Initiating Phase 2: Flank Analysis...")
        # 
        compare_flanks(phase1_harvest)
    else:
        print("[!] Phase 1 found zero markers. Skipping Phase 2.")

    # --- STEP 7: PHASE 3 - FORENSIC REPORT COMPILATION ---
    print("\n" + "="*60)
    print("[*] Initiating Phase 3: Generating Final Forensic Report...")
    try:
        # Re-import to ensure we have the latest Compiler state
        importlib.reload(Compiler)
        Compiler.run_compiler()
    except Exception as e:
        print(f"[!] Error during report generation: {e}")

    # --- WATERFALL SUMMARY (multi-file colour-band view) ---
    try:
        print("[*] Building multi-file waterfall view...")
        names, cats, matrix, colors = build_waterfall_matrix(max_files=500)
        if names:
            wf_path = os.path.join(
                Constants.DIRS.get("OUT_MAPS", Constants.BASE_OUT), "WATERFALL_SUMMARY.png"
            )
            draw_waterfall(
                names, cats, matrix, colors, wf_path,
                title=f"Forensic Waterfall — {len(names)} files x {len(cats)} marker categories",
            )
            print(f"[+] Waterfall: {wf_path}")
        else:
            print("[!] Waterfall skipped (no evidence hits yet).")
    except Exception as wf_err:
        print(f"[!] Waterfall generation failed: {wf_err}")

    try:
        write_portal_snapshot()
    except Exception as snap_err:
        print(f"[!] Portal snapshot failed: {snap_err}")

    try:
        write_triage_summary(Constants.DIRS["OUT_LOGS"], None)  # keep file-tag ranking; matches optional
    except Exception as triage_err:
        print(f"[!] Triage summary failed: {triage_err}")
    try:
        write_coordinate_export(Constants.DIRS["OUT_LOGS"])
    except Exception as bed_err:
        print(f"[!] Coordinate export failed: {bed_err}")

    try:
        export_run_bundle()
    except Exception as bundle_err:
        print(f"[!] Run bundle export failed: {bundle_err}")

    # --- STEP 8: FINAL GLOBAL MEMORY CLEARING ---
    print("\n[*] Finalizing Audit: Releasing System Resources...")    
    AUTO = None
    TARGETS = {}
    MARKER_COLORS = {}    
    plt.close('all')
    gc.collect()    
    print("[+] System Clean. All forensic memory purged.")
    print("="*60 + "\n")

if __name__ == "__main__":
    ensure_environment()
    current_os = platform.system()   
    try:
        if current_os == "Windows":
            mp.set_start_method('spawn', force=True)
            logger.info("[*] OS Detected: Windows. Using 'spawn' mode.")
        else:
            mp.set_start_method('fork', force=True)
            logger.info(f"[*] OS Detected: {current_os}. Using 'fork' mode.")
    except RuntimeWarning:
        pass
    run_audit()