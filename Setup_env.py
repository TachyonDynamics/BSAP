import os
import json
import Constants
import importlib.util

def setup_forensic_environment():
    # 1. Load Configuration
    folders = list(Constants.DIRS.values())
    is_test_mode = str(Constants.CONFIG.get('TEST_MODE', '0')) == '1'
    
    print("=" * 60)
    print(f"[*] Initializing Environment (TEST_MODE={is_test_mode})")
    print("=" * 60)
    
    # 2. Build Directory Scaffolding
    created = 0
    for f in folders:
        if isinstance(f, str) and not f.endswith(".xlsx"):
            if not os.path.exists(f):
                os.makedirs(f, exist_ok=True)
                created += 1
            else:
                os.makedirs(f, exist_ok=True)
    print(f"[+] Directory tree ready ({created} new paths created).")

    # 3. The Logic Gate
    if is_test_mode:
        print("\n[!] TEST_MODE ACTIVE: Generating Calibration Matrix from Targets.json...")
        generate_test_files()
    else:
        print("\n[*] PRODUCTION MODE: Triggering Master Database Sync (DL.py)...")
        try:
            import DL
            DL.sync_all()
        except ImportError as exc:
            print(f"[!] Could not import DL.py ({exc}). Skipping production download.")

def generate_test_files():
    """
    STRESS TEST GENERATION:
    Creates 2x Forward and 2x Reverse hits per file to verify configuration.
    """
    try:
        import Alpha4
    except ImportError:
        print("[!] Critical: Alpha4.py not found.")
        return

    try:
        with open(Constants.TARGETS_PATH, "r", encoding="utf-8") as f:
            targets_data = json.load(f)
        target_dna = targets_data["FORENSIC_MARKERS"]["TEST_MARKER"]["seq"].upper()
    except Exception as e:
        print(f"[!] Error: {e}")
        return

    # 1. Prepare Payloads
    prot_fwd = Alpha4.translate_dna(target_dna)

    # CORRECT REV PROTEIN: translate the reverse-complement DNA, matching
    # exactly what Alpha4.initialize_globals() registers as the REV peptide.
    # (A naive string reversal of prot_fwd is NOT the same peptide and would
    # never be found by the scanner.)
    rc_target = Alpha4.get_rev_comp(target_dna)
    prot_rev = Alpha4.translate_dna(rc_target)

    dna_end, dna_start = target_dna[7:], target_dna[:7]
    prot_end, prot_start = prot_fwd[2:], prot_fwd[:2]
    n_gap = "N" * 50
    
    # 2. Build Content (2x FWD, 2x REV, plus the Split marker at ends)
    # Layout: [Split Suffix] - [REV] - [FWD] - [REV] - [FWD] - [Split Prefix]
    dna_data = (f">DNA_Test\n{dna_end}{n_gap}{rc_target}{n_gap}{target_dna}"
                f"{n_gap}{rc_target}{n_gap}{target_dna}{n_gap}{dna_start}\n")
    
    rna_data = ">RNA_Test\n" + dna_data.split("\n", 1)[1].replace("T", "U")
    
    prot_data = (f">Prot_Test\n{prot_end}{n_gap[:10]}{prot_rev}{n_gap[:10]}{prot_fwd}"
                 f"{n_gap[:10]}{prot_rev}{n_gap[:10]}{prot_fwd}{n_gap[:10]}{prot_start}\n")

    payloads = [("DNA_Test.fasta", dna_data), ("RNA_Test.fasta", rna_data), ("Prot_Test.fasta", prot_data)]

    # 3. Write to Sector Folders
    for key in ["CLN", "IND", "REF"]:
        dest = Constants.DIRS.get(key)
        if not dest: continue
        label = os.path.basename(dest).split('_')[0]
        for filename, content in payloads:
            path = os.path.join(dest, f"{label}_{filename}")
            with open(path, "w") as f: f.write(content)
            print(f"    [+] Created: {path}")

if __name__ == "__main__":
    setup_forensic_environment()