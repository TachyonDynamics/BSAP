#!/usr/bin/env python3
"""
Optional baseline FASTA fetcher (NCBI GenBank via E-utilities).

- Does NOT download on import or on a silent main() run: always asks Y/N.
- DATABASE is a small public teaching set (short genomes / plasmids only).
- Each entry is saved as {Name}_{Accession}.fasta under Reference_Archive,
  Industrial_Sources, or Clinical_Samples.
- Skip logic is per DATABASE entry (that exact file). Other FASTA files the user
  already placed in those folders are left alone; Alpha4 scans any file present,
  not only this list.
"""

import os
import random
import shutil
import sys
import time
import urllib.request

from tqdm import tqdm


def get_random_hex():
    """Same style as Alpha4 master/cross audit bars."""
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))


# --- TARGET DIRECTORIES (prefer absolute Constants.DIRS; else next to this file) ---
def _resolve_dirs():
    root = os.path.dirname(os.path.abspath(__file__))
    try:
        import Constants
        d = getattr(Constants, "DIRS", {}) or {}
        out = {}
        for k, rel in (
            ("REF", "Reference_Archive"),
            ("IND", "Industrial_Sources"),
            ("CLN", "Clinical_Samples"),
        ):
            pth = d.get(k)
            out[k] = pth if pth else os.path.join(root, rel)
        return out
    except Exception:
        return {
            "REF": os.path.join(root, "Reference_Archive"),
            "IND": os.path.join(root, "Industrial_Sources"),
            "CLN": os.path.join(root, "Clinical_Samples"),
        }

DIRS = _resolve_dirs()

# --- PUBLIC BASELINE SET (small files only; one useful example per niche) ---
# All accessions are standard NCBI public records. Entries are offered only
# after the user answers Y at the prompt (or passes --yes).
DATABASE = {
    "REF": {
        # Classic small phage (~5.4 kb) — positive-control genome
        "PhiX174_Phage": "NC_001422.1",
        # Enterobacteria phage lambda (~48.5 kb) — historic molecular-biology reference
        "Lambda_Phage": "NC_001416.1",
        # Small DNA virus (~5.2 kb)
        "SV40_Virus": "NC_001669.1",
        # Human mitochondrion (~16.6 kb) — animal organelle baseline
        "Human_Mitochondrion": "NC_012920.1",
        # Yeast mitochondrion (~85.8 kb) — fungal organelle baseline
        "Yeast_Mitochondrion": "NC_001224.1",
        # Escherichia coli K-12 MG1655 is huge; use a small natural plasmid instead
        # ColE1-type natural plasmid ancestor family (pBR/pUC derive from this world)
        "ColE1_Plasmid": "J01566.1",
    },
    "IND": {
        # Workhorse cloning plasmid (~2.7 kb)
        "pUC19_Plasmid": "L09137.1",
        # Older cloning / antibiotic-marker plasmid (~4.4 kb)
        "pBR322_Plasmid": "J01749.1",
        # AAV2 complete genome (~4.7 kb) — gene-therapy vector backbone family
        "AAV2_Vector": "AF043303.1",
        # T7 phage (~39 kb) — industrial expression / polymerase system
        "T7_Phage": "NC_001604.1",
    },
    "CLN": {
        # HIV-1 HXB2 reference (~9.7 kb) — LTR / lentiviral motif baseline
        "HIV1_HXB2": "K03455.1",
        # Hepatitis B ayw (~3.2 kb) — compact clinical DNA virus
        "HBV_ayw": "NC_003977.2",
        # HPV-16 reference (~7.9 kb) — DNA virus clinical baseline
        "HPV16_Reference": "NC_001526.4",
    },
}


def expected_path(name, accession, folder):
    return os.path.join(folder, f"{name}_{accession}.fasta")


def fetch_fasta(name, accession, folder, retries=3):
    """
    Returns "kept", "new", or "failed".

    "kept" only if THIS list entry's file already exists and looks non-empty.
    Other FASTAs in the same folder do not satisfy this entry (and are never deleted).
    Writes to a .part file and renames only after the payload starts with '>'.
    """
    file_path = expected_path(name, accession, folder)
    if os.path.exists(file_path) and os.path.getsize(file_path) > 100:
        return "kept"

    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={accession}&rettype=fasta"
    )
    part = file_path + ".part"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as response, open(part, "wb") as out:
                shutil.copyfileobj(response, out)
            with open(part, "rb") as f:
                if f.read(1) != b">":
                    raise ValueError("response is not FASTA")
            os.replace(part, file_path)
            time.sleep(0.4)  # NCBI: <=3 requests/s without an API key
            return "new"
        except Exception as e:
            err = str(e)
            tqdm.write(f"[!] FAIL {name} | {accession} | {err}")
            if any(code in err for code in ("400", "404", "HTTP Error 400", "HTTP Error 404")):
                break
            time.sleep(2 ** attempt)
    if os.path.exists(part):
        try:
            os.remove(part)
        except OSError:
            pass
    return "failed"


def _count_existing_any(folder):
    """How many sequence-like files the user already has in this folder (any name)."""
    if not os.path.isdir(folder):
        return 0
    n = 0
    for fn in os.listdir(folder):
        if fn.startswith("."):
            continue
        low = fn.lower()
        if low.endswith((".fasta", ".fa", ".fna", ".ffn", ".faa")):
            n += 1
    return n


def _prompt_yes_no(message, default_no=True):
    """Skip interactive prompt under Streamlit / pipes / non-TTY (default = No)."""
    if not sys.stdin or not getattr(sys.stdin, "isatty", lambda: False)():
        print(f"{message} [auto: {'N' if default_no else 'Y'} — non-interactive]")
        return not default_no
    hint = " [y/N] " if default_no else " [Y/n] "
    try:
        ans = input(message + hint).strip().lower()
    except EOFError:
        return not default_no
    if not ans:
        return not default_no
    return ans in ("y", "yes")


def _baseline_missing():
    """List (category, name, accession) for DATABASE entries not yet on disk."""
    missing = []
    for category, files in DATABASE.items():
        folder = DIRS[category]
        for name, acc in files.items():
            p = expected_path(name, acc, folder)
            if not (os.path.exists(p) and os.path.getsize(p) > 100):
                missing.append((category, name, acc))
    return missing


def sync_all(assume_yes=False):
    print("--- BSAP baseline FASTA sync (NCBI) ---")
    print("Public teaching set only (short plasmids / viruses).")
    print("Your own files in Reference_Archive / Industrial_Sources /")
    print("Clinical_Samples are never deleted or overwritten by name mismatch.\n")

    for d in DIRS.values():
        os.makedirs(d, exist_ok=True)

    total_list = 0
    total_present = 0
    for category, folder in DIRS.items():
        n_any = _count_existing_any(folder)
        entries = DATABASE.get(category, {})
        n_list = sum(
            1
            for name, acc in entries.items()
            if os.path.exists(expected_path(name, acc, folder))
            and os.path.getsize(expected_path(name, acc, folder)) > 100
        )
        total_list += len(entries)
        total_present += n_list
        print(
            f"  {category}: {n_any} FASTA(s) already on disk, "
            f"{n_list}/{len(entries)} baseline list present"
        )

    missing = _baseline_missing()
    if not missing:
        print(
            f"\n[+] All {total_present}/{total_list} baseline FASTAs already present — "
            "nothing to download."
        )
        print("[*] Alpha4 will scan every FASTA in those folders, not only this list.")
        return

    print(f"\n[*] Missing {len(missing)} baseline file(s):")
    for category, name, acc in missing:
        print(f"      - {category}/{name}_{acc}.fasta")

    if not assume_yes:
        if not _prompt_yes_no(
            "\nDownload missing baseline FASTAs from NCBI now?",
            default_no=True,
        ):
            print(
                "[*] Skipped downloads. Drop your own FASTA files into the "
                "three folders and run the audit."
            )
            return

    for category, files in DATABASE.items():
        folder = DIRS[category]
        print(f"\n[*] Syncing baseline list -> {category} ({folder})")
        tally = {"kept": 0, "new": 0, "failed": 0}

        for name, acc in tqdm(
            files.items(),
            colour=get_random_hex(),
            desc=f"Syncing {category}",
            leave=True,
            dynamic_ncols=True,
            bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}]",
        ):
            tally[fetch_fasta(name, acc, folder)] += 1

        print(
            f"[+] {category} Summary: {tally['kept']} kept, "
            f"{tally['new']} new, {tally['failed']} failed."
        )

    print(
        "\n[*] Done. Alpha4 will scan every FASTA in those folders, "
        "not only this list."
    )


if __name__ == "__main__":
    # Optional: python DL.py --yes  for non-interactive CI (still explicit)
    assume = any(a in ("--yes", "-y") for a in sys.argv[1:])
    sync_all(assume_yes=assume)
