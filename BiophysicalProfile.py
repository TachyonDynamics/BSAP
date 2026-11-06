"""
BiophysicalProfile.py
======================
Computes a real, physically-grounded profile for a matched sequence, using
Physics.py's existing functions. Physics.py itself is left untouched; this
module calls its functions correctly and adds the one correction that
matters for a biological (not vacuum) context: dielectric screening by
water.

WHY THE WATER-DIELECTRIC CORRECTION MATTERS:
UniversalPhysics.calculate_coulomb_force() computes electrostatic force as
if the two charges were in a vacuum. In an aqueous/cellular environment,
water's dielectric constant (already defined in Physics.py as
WATER_DIELECTRIC = 80.1, but never applied anywhere) screens electrostatic
interactions by roughly that same factor. This module reports BOTH the raw
vacuum figure Physics.py's function produces and the water-screened
estimate, side by side, so nothing is silently wrong or silently hidden.

WHAT THE PHOTON-ENERGY COMPARISON IS FOR:
get_photon_energy() and get_thermal_velocity() are real, correctly-applied
formulas (E=hf in eV; kT/h as the thermal-photon crossover frequency). Used
on a sequence's estimated vibrational frequency, they answer two genuinely
meaningful physics questions: (1) is this frequency anywhere near strong
enough to break a covalent bond (~eV range), and (2) does this frequency
sit above or below the ambient thermal noise floor at body temperature.
Both answers are real math, reported honestly -- including when the honest
answer is "no, nowhere close."
WHY MASS IS NOT COMPUTED VIA BiochemicalEngine.get_molecular_mass() DIRECTLY:
That function auto-detects molecule type per-sequence ("all ATCG -> DNA,
else -> PROTEIN"). Real sequencing data routinely contains IUPAC ambiguity
codes (N, R, Y, S, W, K, M, B, D, H, V -- this same codebase's own
get_rev_comp() explicitly handles all of these). A single 'N' in an
otherwise normal DNA sequence causes that auto-detection to misclassify the
whole sequence as a protein and silently return a mass computed from amino
acid residue formulas instead of nucleotide ones -- confirmed to be off by
more than 3x on a real test case, with no error or warning. Since this
module already knows the correct mol_type from the caller (Alpha4.py knows
exactly what it scanned), _typed_molecular_mass() below computes mass using
that known type for every residue, skipping (and honestly counting) any
residue it can't resolve rather than silently guessing a different
molecule type.
"""

import sys
import os
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Physics import UniversalPhysics, BiochemicalEngine, UniversalFrequencyMapper, AcousticForensics
from Sonification import get_residue_vibrational_hz

WATER_DA = 18.01528

def _typed_molecular_mass(sequence, mol_type):
    """Mass using the CALLER-specified mol_type for every residue (see module docstring).
    Returns (mass_da, unresolved_count); unresolved residues (N, '*', ...) are counted, not guessed."""
    total = 0.0
    unresolved = 0
    for residue, n in Counter(sequence.upper()).items():
        if mol_type in ("DNA", "RNA"):
            formula = BiochemicalEngine.NUCLEOTIDE_MAP.get(mol_type, {}).get(residue)
        else:
            formula = BiochemicalEngine.PROTEIN_MAP.get(residue, {}).get("formula")
        if not formula:
            unresolved += n
            continue
        total += BiochemicalEngine.get_atomic_mass(formula) * n
    if total:
        total += WATER_DA      # residues are chain units (monomer minus H2O); one water closes the ends
    return round(total, 4), unresolved

def build_profile(sequence, mol_type="DNA", ph=7.4, temp_k=310.15):
    sequence = sequence.upper()
    profile = {
        "sequence": sequence,
        "length": len(sequence),
        "mol_type": mol_type,
    }

    # --- Molecular mass (type-safe -- see module docstring for why this
    # doesn't call BiochemicalEngine.get_molecular_mass() directly) ---
    mass, unresolved = _typed_molecular_mass(sequence, mol_type)
    profile["molecular_mass_da"] = mass
    profile["unresolved_residue_count"] = unresolved

    # --- Net charge (protein only -- Henderson-Hasselbalch, real biochemistry) ---
    if mol_type == "PROT":
        profile["net_charge_at_ph"] = UniversalPhysics.calculate_net_charge(sequence, ph=ph)
    else:
        profile["net_charge_at_ph"] = None  # not meaningful for nucleic acid backbone here

    # --- Electrostatics ---
    # Protein: Coulomb force between the chain halves, separated by the contour length (vacuum vs water-screened).
    # DNA/RNA: that estimate (q ~ L, r ~ L) cancels to a constant, so report total backbone charge
    # (1 negative charge per phosphate) and charge density instead.
    if mol_type == "PROT":
        raw_force_n = UniversalPhysics.calculate_coulomb_force(profile["net_charge_at_ph"], len(sequence))
        profile["coulomb_force_vacuum_newtons"] = raw_force_n
        profile["coulomb_force_water_screened_newtons"] = raw_force_n / UniversalPhysics.WATER_DIELECTRIC
        profile["total_charge_e"] = None
        profile["charge_density_e_per_nm"] = None
    else:
        profile["coulomb_force_vacuum_newtons"] = float("nan")
        profile["coulomb_force_water_screened_newtons"] = float("nan")
        profile["total_charge_e"] = -len(sequence)
        profile["charge_density_e_per_nm"] = -1 / 0.34 if sequence else 0.0

    # --- Vibrational estimate + honest physical context for the FIRST residue ---
    # (see Sonification.py docstring for why per-residue mass is used, not whole-sequence mass)
    if sequence:
        first_residue, vib_hz = sequence[0], None
        for _res in sequence:
            vib_hz = get_residue_vibrational_hz(_res, mol_type)
            if vib_hz:
                first_residue = _res
                break
        profile["sample_residue"] = first_residue
        profile["sample_vibrational_hz"] = vib_hz

        if vib_hz:
            photon_info = UniversalPhysics.get_photon_energy(vib_hz)
            profile["sample_photon_energy_ev"] = photon_info["eV"]
            profile["sample_bond_breaking_capable"] = photon_info["threat"] == "IONIZING"

            thermal_info = UniversalPhysics.get_thermal_velocity(vib_hz, temp_k=temp_k)
            profile["thermal_noise_floor_hz"] = thermal_info["limit_hz"]
            profile["above_thermal_noise_floor"] = thermal_info["detectable"]

    return profile


def format_profile_report(profile):
    """Human-readable text block for inclusion in an evidence file."""
    lines = []
    lines.append("--- BIOPHYSICAL PROFILE (estimates; see Sonification.py / BiophysicalProfile.py docstrings for methodology and caveats) ---")
    lines.append(f"Sequence: {profile['sequence']} ({profile['mol_type']}, {profile['length']}bp/aa)")
    lines.append(f"Molecular mass: {profile['molecular_mass_da']} Da")
    if profile.get("unresolved_residue_count", 0) > 0:
        lines.append(f"  (NOTE: {profile['unresolved_residue_count']} residue(s) could not be resolved -- e.g. an ambiguity code like 'N', or a stop codon '*' -- and were excluded from mass, not guessed)")
    if profile["net_charge_at_ph"] is not None:
        lines.append(f"Net charge (Henderson-Hasselbalch, pH 7.4): {profile['net_charge_at_ph']}")
    force = profile["coulomb_force_vacuum_newtons"]
    if force == force:   # not NaN (protein)
        lines.append(f"Electrostatic force -- vacuum estimate: {force:.3e} N")
        lines.append(f"Electrostatic force -- water-screened estimate: {profile['coulomb_force_water_screened_newtons']:.3e} N")
    else:
        lines.append(f"Backbone charge: {profile['total_charge_e']} e ({profile['charge_density_e_per_nm']:.2f} e/nm)")
    if profile.get("sample_vibrational_hz"):
        lines.append(
            f"Sample residue '{profile['sample_residue']}' single-oscillator vibrational estimate: "
            f"{profile['sample_vibrational_hz']/1e12:.3f} THz "
            f"(published DNA THz spectroscopy range: 0.03-10 THz)"
        )
        lines.append(
            f"  -> Photon-equivalent energy: {profile['sample_photon_energy_ev']} eV "
            f"(bond-breaking/ionizing threshold ~10 eV -> "
            f"{'CAPABLE' if profile['sample_bond_breaking_capable'] else 'NOT capable of bond-breaking at this frequency'})"
        )
        lines.append(
            f"  -> Body-temp thermal noise floor: {profile['thermal_noise_floor_hz']/1e12:.3f} THz -> "
            f"{'above' if profile['above_thermal_noise_floor'] else 'below'} thermal noise floor"
        )
    return "\n".join(lines)


def format_optional_heuristics(sequence, mol_type="DNA", vib_hz=None, mass_da=None):
    """
    Optional Physics.py helpers gated by ENABLE_* config keys.
    All blocks are labeled HEURISTIC — illustrative models, not measurements
    and not proof of lab origin, emission, or absorption in vivo.
    """
    import Constants
    seq = (sequence or "").upper()
    lines = []
    cfg = getattr(Constants, "CONFIG", {}) or {}

    def on(key):
        return str(cfg.get(key, "0")).strip() == "1"

    if mass_da is None and seq:
        try:
            mass_da, _ = _typed_molecular_mass(seq, mol_type)
        except Exception:
            mass_da = None
    if vib_hz is None and seq:
        try:
            for res in seq:
                vib_hz = get_residue_vibrational_hz(res, mol_type)
                if vib_hz:
                    break
        except Exception:
            vib_hz = None

    if on("ENABLE_ZETA_HEURISTICS"):
        lines.append("--- HEURISTIC: STATIC ZETA / BACKBONE CHARGE (not measured zeta-potential) ---")
        try:
            z = BiochemicalEngine.get_static_zeta_potential(ph=7.4)
            lines.append(f"  Model output @ pH 7.4: {z}")
            lines.append("  HEURISTIC only — coarse slogan, not an electrophoretic measurement.")
        except Exception as e:
            lines.append(f"  unavailable: {e}")
        lines.append("")

    if on("ENABLE_DIELECTRIC_HEURISTICS"):
        lines.append("--- HEURISTIC: DIELECTRIC HEATING MODEL (not measured RF absorption) ---")
        try:
            if mass_da and vib_hz:
                heat = BiochemicalEngine.calculate_dielectric_heating(
                    mass_da, float(vib_hz), exposure_seconds=60
                )
                lines.append(
                    f"  Input mass: {mass_da} Da | stand-in freq: {float(vib_hz):.3e} Hz | exposure: 60 s"
                )
                lines.append(f"  Model delta-T: {heat}")
                lines.append("  HEURISTIC only — does not prove heating, emission, or biological effect.")
            else:
                lines.append("  skipped (need mass + vibrational estimate)")
        except Exception as e:
            lines.append(f"  unavailable: {e}")
        lines.append("")

    if on("ENABLE_CODON_FOOTPRINT_HEURISTICS") and mol_type in ("DNA", "RNA"):
        lines.append("--- HEURISTIC: CODON FOOTPRINT SCAN (BiochemicalEngine.scan_sequence) ---")
        try:
            dna = seq.replace("U", "T")
            report, footprint, aa_counts = BiochemicalEngine.scan_sequence(dna)
            lines.append(f"  Codons scored: {len(report)} | AA counts: {dict(list(aa_counts.items())[:12])}")
            lines.append(f"  Atom footprint (summed DNA triplets): {footprint}")
            for row in (report or [])[:5]:
                lines.append(
                    f"    {row.get('codon')} -> {row.get('aa')} "
                    f"prot_mass={row.get('prot_mass')} dna_mass={row.get('dna_mass')}"
                )
            lines.append("  HEURISTIC only — oligo composition, not expression proof.")
        except Exception as e:
            lines.append(f"  unavailable: {e}")
        lines.append("")

    if on("ENABLE_MOLECULE_TYPE_HEURISTICS"):
        lines.append("--- HEURISTIC: MOLECULE-TYPE GUESS (alphabet heuristic) ---")
        try:
            guess = BiochemicalEngine.identify_molecule_type(seq)
            lines.append(f"  Engine guess: {guess} | declared mol_type: {mol_type}")
            lines.append("  HEURISTIC only — alphabet statistics, not lab assay.")
        except Exception as e:
            lines.append(f"  unavailable: {e}")
        lines.append("")

    if on("ENABLE_FREQ_CONTEXT_HEURISTICS") and vib_hz:
        lines.append("--- HEURISTIC: FREQUENCY CONTEXT (display taxonomy, not RF detection) ---")
        try:
            ctx = UniversalFrequencyMapper.get_forensic_context(float(vib_hz), medium="LAB")
            lines.append(f"  analyze_signal_physics: {ctx.get('physics')}")
            lines.append(f"  optical context: {ctx.get('optical')}")
            lines.append(f"  sonar-band tag: {ctx.get('sonar')}")
            lines.append(f"  medium multiplier (LAB): {UniversalFrequencyMapper.get_medium_multiplier('LAB')}")
            lines.append(
                f"  acoustic transposition helper: "
                f"{UniversalFrequencyMapper.get_acoustic_transposition(float(vib_hz))}"
            )
            lines.append("  HEURISTIC / display only — not measured radio, optical, or sonar emission.")
        except Exception as e:
            lines.append(f"  unavailable: {e}")
        lines.append("")

    if on("ENABLE_ACOUSTIC_HEURISTICS"):
        lines.append("--- HEURISTIC: ACOUSTIC ABSORPTION MODEL (material toy model) ---")
        try:
            thickness_m = max(len(seq), 1) * 0.34e-9
            density = 1000.0
            bulk_mod = 2.2e9
            peak = AcousticForensics.calculate_absorption_peak(thickness_m, density, bulk_mod)
            coef = AcousticForensics.get_absorption_coefficient(0.5, 1.0)
            lines.append(f"  Toy thickness from length: {thickness_m:.3e} m ({len(seq)} bp x 0.34 nm)")
            lines.append(f"  Model absorption-peak Hz: {peak}")
            lines.append(f"  Model absorption coefficient (R=0.5 demo): {coef}")
            lines.append("  HEURISTIC only — not ultrasound data, not proof of a cavity or implant.")
        except Exception as e:
            lines.append(f"  unavailable: {e}")
        lines.append("")

    if on("ENABLE_ALT_MASS_HEURISTICS") and seq and mol_type in ("DNA", "RNA", "PROT"):
        lines.append("--- HEURISTIC: ALTERNATE get_molecular_mass PATH ---")
        try:
            alt = BiochemicalEngine.get_molecular_mass(seq)
            lines.append(f"  BiochemicalEngine.get_molecular_mass: {alt}")
            lines.append(f"  Typed path mass (primary): {mass_da}")
            lines.append("  HEURISTIC comparison — primary path remains typed residue mass.")
        except Exception as e:
            lines.append(f"  unavailable: {e}")
        lines.append("")

    if not lines:
        return ""
    return "\n".join(lines).rstrip() + "\n"
