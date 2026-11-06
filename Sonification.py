"""
Sonification.py
================
Turns a nucleotide/protein sequence into real, playable audio.

WHAT THIS IS: a "sonification" — a well-established technique (used for LIGO
gravitational wave data, radio-astronomy datasets, accessibility tools, etc.)
that maps a physical quantity onto an audible tone so a human can listen to
data instead of only reading it. It does NOT claim that DNA emits sound, or
that these tones have any biological effect on the sequence being played.

WHERE THE FREQUENCY PER BASE COMES FROM (and what it does/doesn't mean):
Physics.py's UniversalPhysics.get_vibrational_resonance() models a single
mass on a spring (f = 1/(2*pi) * sqrt(k / m)) using the DNA_BOND_K constant.
Fed the mass of a whole sequence, that model is meaningless (a real DNA
molecule's mass has no simple relationship to a single spring constant).
Fed the mass of a single nucleotide residue (~305-345 Da), however, this
same unmodified formula lands at ~4.7-5.0 THz per base -- which sits
squarely inside the real, independently measured range for DNA collective
vibrational modes reported by terahertz time-domain spectroscopy in the
literature (published measurements: ~0.03-10 THz; see e.g. Fischer et al.,
"Far-infrared vibrational modes of DNA components studied by terahertz
time-domain spectroscopy"; Xie et al., "THz spectra and corresponding
vibrational modes of DNA base pair cocrystals and polynucleotides").

This is still a simplified single-oscillator approximation, not a full
normal-mode / molecular-dynamics simulation -- real measured spectra are
collective, many-atom motions and the literature notes sequence-dependent
differences concentrate mostly below ~1 THz, not spread evenly across the
whole per-base range this toy model produces. Treat the per-base frequency
here as an order-of-magnitude-correct, illustrative estimate, not a
laboratory measurement of your specific sequence.

The 20Hz-20kHz audio tone assigned to each base used to come from Physics.py's
get_acoustic_transposition() (a fixed octave-halving of the real frequency).
That preserves the real *ratio* between residues -- but since all four DNA
bases' real vibrational estimates sit within ~13% of each other (they're all
similar-mass nucleotides), a fixed-ratio octave shift keeps them clustered
within a couple Hz of each other no matter how far down you shift them --
audibly indistinguishable. That's not a bug, it's what the real mass
similarity predicts; it's just not a *useful* sonification.

The fix used here is the standard sonification remedy: map each residue's
real vibrational frequency PROPORTIONALLY across a wide, clearly audible
pitch range (220Hz-1760Hz, three octaves), based on where it falls between
the lightest and heaviest residue in its own alphabet. This still preserves
the real physical ordering (heavier residue -> lower real frequency -> lower
pitch) -- it just stretches the real, tiny differences out to where a human
ear can actually hear them, exactly like any other data sonification does.

A second, independent axis of real distinction is layered on top for timbre:
purines (A/G -- real two-ring structures) get a genuine second harmonic
mixed in with the fundamental tone; pyrimidines (T/C/U -- real single-ring
structures) stay a pure sine tone. For protein sequences, the same idea is
applied to real aromatic side chains (F/W/Y/H) vs. everything else. Neither
of these is an invented category -- purine/pyrimidine ring count and amino
acid aromaticity are standard, textbook structural classifications.
"""

import functools
import math
import wave
import sys
import os

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Physics import UniversalPhysics, BiochemicalEngine, UniversalFrequencyMapper

AUDIBLE_MIN_HZ = 220.0   # A3
AUDIBLE_MAX_HZ = 1760.0  # A6 -- three octaves of clearly perceptible range
MAX_TONES = 3000         # about 7.5 minutes of audio; longer inputs are truncated
PURINES = {"A", "G"}          # real two-ring nucleobases
PYRIMIDINES = {"T", "C", "U"}  # real single-ring nucleobases
AROMATIC_AA = {"F", "W", "Y", "H"}  # real aromatic-side-chain amino acids

_range_cache = {}

@functools.lru_cache(maxsize=None)

def get_residue_vibrational_hz(residue, mol_type="DNA"):
    """
    Real per-residue vibrational frequency estimate (Hz), using Physics.py's
    unmodified get_vibrational_resonance() fed the ACTUAL mass of a single
    nucleotide or amino acid residue (not the whole sequence -- see module
    docstring for why that distinction matters).
    """
    residue = residue.upper()

    if mol_type in ("DNA", "RNA"):
        formula = BiochemicalEngine.NUCLEOTIDE_MAP.get(mol_type, {}).get(residue)
        if not formula:
            return None
    else:  # PROT
        formula = BiochemicalEngine.PROTEIN_MAP.get(residue, {}).get("formula")
        if not formula:
            return None

    mass_daltons = BiochemicalEngine.get_atomic_mass(formula)
    if mass_daltons <= 0:
        return None

    return UniversalPhysics.get_vibrational_resonance(mass_daltons)


def _residue_alphabet(mol_type):
    if mol_type in ("DNA", "RNA"):
        return list(BiochemicalEngine.NUCLEOTIDE_MAP.get(mol_type, {}).keys())
    return list(BiochemicalEngine.PROTEIN_MAP.keys())


def _real_hz_range(mol_type):
    """Real min/max vibrational estimate across this molecule type's own
    alphabet, computed from Physics.py's own formulas -- not hardcoded."""
    if mol_type in _range_cache:
        return _range_cache[mol_type]
    values = [get_residue_vibrational_hz(r, mol_type) for r in _residue_alphabet(mol_type)]
    values = [v for v in values if v]
    result = (min(values), max(values)) if values else (0, 0)
    _range_cache[mol_type] = result
    return result


def map_to_audible_range(real_hz, mol_type):
    """Proportionally maps a real vibrational frequency onto AUDIBLE_MIN_HZ..
    AUDIBLE_MAX_HZ based on where it falls between the lightest and heaviest
    residue of its own alphabet. See module docstring for why this replaces
    a fixed-ratio octave shift."""
    lo, hi = _real_hz_range(mol_type)
    if hi <= lo:
        return (AUDIBLE_MIN_HZ + AUDIBLE_MAX_HZ) / 2
    fraction = (real_hz - lo) / (hi - lo)
    fraction = max(0.0, min(1.0, fraction))
    return AUDIBLE_MIN_HZ + fraction * (AUDIBLE_MAX_HZ - AUDIBLE_MIN_HZ)


def has_ring_harmonic(residue, mol_type):
    """Real structural classification used for timbre: two-ring purines /
    aromatic amino acids get a second harmonic; everything else is a pure tone."""
    if mol_type in ("DNA", "RNA"):
        return residue in PURINES
    return residue in AROMATIC_AA


def sequence_to_tones(sequence, mol_type="DNA"):
    """
    Maps every residue to (residue, real_vibrational_hz, audible_hz,
    has_harmonic) using the physically-grounded per-residue frequency,
    proportionally mapped into the audible range (see module docstring),
    with a real structural flag for whether a second harmonic should be
    mixed in for timbre.
    """
    tones = []
    for residue in sequence.upper():
        real_hz = get_residue_vibrational_hz(residue, mol_type)
        if real_hz is None:
            continue
        audible_hz = map_to_audible_range(real_hz, mol_type)
        harmonic = has_ring_harmonic(residue, mol_type)
        tones.append((residue, real_hz, audible_hz, harmonic))
    return tones


def synthesize_wave(tones, note_seconds=0.15, sample_rate=44100, amplitude=0.4):
    """Renders (residue, real_hz, audible_hz, has_harmonic) tuples into PCM samples (float array in -1..1):
    one note per residue. Purines/aromatics mix in a second harmonic (2x, 35% amplitude); everything
    else is a pure sine. A 10 ms linear fade in/out per note avoids clicks."""
    n = int(sample_rate * note_seconds)
    t = np.arange(n) / sample_rate
    fade = max(1, int(sample_rate * 0.01))
    envelope = np.ones(n)
    envelope[:fade] = np.linspace(0, 1, fade, endpoint=False)
    envelope[-fade:] = np.linspace(1, 0, fade)
    harmonic_ratio = 0.35

    cache, notes = {}, []
    for _residue, _real_hz, audible_hz, harmonic in tones:
        if audible_hz <= 0:
            continue
        key = (round(audible_hz, 3), harmonic)
        if key not in cache:
            wave_ = np.sin(2 * np.pi * audible_hz * t)
            if harmonic:
                wave_ = (wave_ + harmonic_ratio * np.sin(4 * np.pi * audible_hz * t)) / (1 + harmonic_ratio)
            cache[key] = amplitude * wave_ * envelope
        notes.append(cache[key])
    return np.concatenate(notes) if notes else np.zeros(0)

def export_wav(samples, path, sample_rate=44100):
    """Writes real PCM samples to a standard, playable 16-bit mono .wav file."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    pcm = (np.clip(np.asarray(samples), -1.0, 1.0) * 32767).astype("<i2")
    with wave.open(path, "wb") as wav_file:
        wav_file.setnchannels(1)
        wav_file.setsampwidth(2)  # 16-bit
        wav_file.setframerate(sample_rate)
        wav_file.writeframes(pcm.tobytes())




# =============================================================================
# FIXED-DURATION DENSE SONIFICATION
# =============================================================================
# Any length sequence maps to a short clip (default ~2.5 s). Local windows
# drive noise floor, crackle rate, and brief tones — a compressed overview,
# not one-note-per-base (impossible for multi-GB genomes).

DEFAULT_TARGET_SECONDS = 2.5
DEFAULT_MAX_SECONDS = 6.0
N_WINDOWS = 512  # resolution across the clip


def _window_features(sequence, mol_type, n_windows=N_WINDOWS):
    """Per-window GC/purine/entropy-like scores in [0,1]."""
    seq = sequence.upper()
    L = len(seq)
    if L == 0:
        return []
    n_windows = max(8, min(n_windows, L))
    win = max(1, L // n_windows)
    feats = []
    for i in range(n_windows):
        a = i * win
        b = L if i == n_windows - 1 else min(L, (i + 1) * win)
        chunk = seq[a:b]
        if not chunk:
            feats.append({"gc": 0.0, "pur": 0.0, "run": 0.0, "len": 0})
            continue
        if mol_type == "PROT":
            charged = sum(1 for c in chunk if c in "KRHDE") / len(chunk)
            aromatic = sum(1 for c in chunk if c in "FWYH") / len(chunk)
            feats.append({"gc": charged, "pur": aromatic, "run": aromatic, "len": len(chunk)})
        else:
            gc = sum(1 for c in chunk if c in "GC") / len(chunk)
            pur = sum(1 for c in chunk if c in "AG") / len(chunk)
            # longest same-base run fraction
            best = 1
            cur = 1
            for j in range(1, len(chunk)):
                if chunk[j] == chunk[j - 1]:
                    cur += 1
                    best = max(best, cur)
                else:
                    cur = 1
            run = min(1.0, best / max(8, len(chunk) * 0.05))
            feats.append({"gc": gc, "pur": pur, "run": run, "len": len(chunk)})
    return feats


def synthesize_dense(sequence, mol_type="DNA", target_seconds=None, sample_rate=44100):
    """
    Compress any sequence into a short noisy overview clip.
    ~1 MB DNA → fraction of a second to ~1 s; multi-GB → still a few seconds.
    """
    if target_seconds is None:
        target_seconds = DEFAULT_TARGET_SECONDS
    # Slight length scaling: tiny files shorter, huge files cap
    L = len(sequence)
    if L < 500:
        target_seconds = min(target_seconds, 0.4)
    elif L < 50_000:
        target_seconds = min(target_seconds, 1.2)
    elif L < 5_000_000:
        target_seconds = min(target_seconds, 2.5)
    else:
        target_seconds = min(max(target_seconds, 2.5), DEFAULT_MAX_SECONDS)

    n_samples = max(1, int(sample_rate * target_seconds))
    feats = _window_features(sequence, mol_type, n_windows=min(N_WINDOWS, max(16, n_samples // 64)))
    if not feats:
        return np.zeros(n_samples, dtype=np.float64)

    # Upsample feature tracks to sample rate
    idx = np.linspace(0, len(feats) - 1, n_samples)
    gc = np.array([feats[int(i)].get("gc", 0) for i in idx], dtype=np.float64)
    pur = np.array([feats[int(i)].get("pur", 0) for i in idx], dtype=np.float64)
    run = np.array([feats[int(i)].get("run", 0) for i in idx], dtype=np.float64)

    rng = np.random.default_rng(abs(hash(sequence[:64])) % (2**32))
    # Broadband static modulated by GC
    noise = rng.standard_normal(n_samples) * (0.08 + 0.35 * gc)
    # Occasional crackles (impulses) where runs/complexity spike
    crackle_prob = 0.002 + 0.02 * run + 0.01 * (gc - 0.5) ** 2
    crackles = (rng.random(n_samples) < crackle_prob).astype(np.float64)
    crackles *= rng.uniform(0.3, 1.0, n_samples)
    # Brief tonal blips from purine density (audible carrier)
    t = np.arange(n_samples) / sample_rate
    carrier_hz = 180 + 900 * pur  # wanders with purine content
    # Integrate phase for continuous frequency modulation
    phase = 2 * np.pi * np.cumsum(carrier_hz) / sample_rate
    tone = 0.12 * np.sin(phase) * (0.3 + 0.7 * pur)
    # Sparse high ticks
    ticks = (rng.random(n_samples) < (0.0005 + 0.008 * gc * pur)).astype(np.float64)
    ticks *= rng.choice([-1.0, 1.0], size=n_samples) * 0.5

    samples = noise + crackles * 0.6 + tone + ticks
    # Soft limit
    peak = np.max(np.abs(samples)) or 1.0
    samples = 0.85 * samples / peak
    return samples


def sonify_sequence(
    sequence,
    mol_type="DNA",
    out_path=None,
    note_seconds=0.15,
    sample_rate=44100,
    mode="dense",
    target_seconds=None,
):
    """
    mode=\"dense\" (default): fixed short duration overview for any length file.
    mode=\"legacy\": one tone per residue (truncated), useful only for short motifs.
    """
    if mode == "legacy":
        tones = sequence_to_tones(sequence, mol_type)[:MAX_TONES]
        if out_path:
            samples = synthesize_wave(tones, note_seconds=note_seconds, sample_rate=sample_rate)
            export_wav(samples, out_path, sample_rate=sample_rate)
        return tones

    samples = synthesize_dense(
        sequence, mol_type=mol_type, target_seconds=target_seconds, sample_rate=sample_rate
    )
    if out_path:
        export_wav(samples, out_path, sample_rate=sample_rate)
    return {
        "mode": "dense",
        "n_samples": int(len(samples)),
        "duration_s": round(len(samples) / sample_rate, 3),
        "seq_len": len(sequence),
    }


# =============================================================================
# SEQUENCE-FEATURE SPECTROGRAM (NOT measured RF)
# =============================================================================
# Maps local sequence composition onto a log-frequency axis from ~1 Hz to
# optical (~3e14 Hz) for visualisation only. Peaks are where local GC / purine
# / run-length features concentrate after a deliberate band mapping — they are
# NOT radio detections of 2.4 GHz Wi-Fi or cellular bands.

def feature_spectrogram(
    sequence,
    mol_type="DNA",
    out_png=None,
    n_time=512,
    n_freq=160,
    f_min=1.0,
    f_max=6.0e14,
    marker_positions=None,
    title=None,
):
    """
    Heatmap-style sequence-feature spectrogram.

    X = genomic position, Y = log frequency *display* axis (1 Hz → ~600 THz).
    Color = local composition energy mapped onto author-chosen bands.

    This is NOT measured electromagnetic absorption or emission.
    marker_positions: optional list of {pos, name, color} for a thin track.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    seq = (sequence or "").upper()
    disclaimer = (
        "Sequence-feature spectrogram — NOT measured RF/EM emission. "
        "Y-axis is a log frequency DISPLAY axis; intensity comes from local "
        "GC / purine / run structure (and optional marker proximity), not a spectrum analyzer."
    )
    if not seq:
        return {"peaks": [], "path": None, "disclaimer": disclaimer}

    # Cap work for multi-GB: even subsample keeps relative structure
    max_bp = 2_000_000
    if len(seq) > max_bp:
        step = len(seq) // max_bp
        seq = seq[::step]

    L = len(seq)
    n_time = max(64, min(n_time, L))
    feats = _window_features(seq, mol_type, n_windows=n_time)
    freqs = np.logspace(np.log10(f_min), np.log10(f_max), n_freq)
    logf = np.log10(freqs)
    grid = np.zeros((n_freq, n_time), dtype=np.float64)

    # Display-band kernels (author-chosen channels, not physics measurements)
    for ti, f in enumerate(feats):
        gc, pur, run = f["gc"], f["pur"], f["run"]
        centers = [
            (3.0, 0.35 * run),                 # kHz-ish weight from runs
            (6.0, 0.55 * run),                 # MHz-ish
            (9.38, 0.90 * gc),                 # ~2.4e9 display channel × GC
            (9.70, 0.55 * gc),                 # ~5e9 display channel × GC
            (11.0, 0.40 * (gc * pur)),         # tens of GHz-ish
            (12.5, 0.85 * pur),                # THz-ish × purines
            (13.5, 0.50 * ((gc + pur) / 2)),   # far-IR-ish
            (14.5, 0.45 * ((gc + pur) / 2)),   # near optical edge
            (14.78, 0.35 * pur),               # ~600 THz display channel
        ]
        for c, w in centers:
            if w <= 0:
                continue
            grid[:, ti] += w * np.exp(-0.5 * ((logf - c) / 0.38) ** 2)

    # Optional boost near marker positions (highlight forensic hits on the same canvas)
    if marker_positions and L > 0:
        for m in marker_positions:
            try:
                pos = float(m.get("pos", 0))
            except (TypeError, ValueError):
                continue
            ti = int(max(0, min(n_time - 1, (pos / max(L, 1)) * (n_time - 1))))
            grid[:, ti] += 0.25

    profile = grid.mean(axis=1)
    peaks = []
    if profile.max() > 0:
        for i in range(2, len(profile) - 2):
            if profile[i] >= profile[i - 1] and profile[i] >= profile[i + 1] and profile[i] > 0.12 * profile.max():
                peaks.append({
                    "hz": float(freqs[i]),
                    "strength": float(profile[i] / profile.max()),
                    "label": _freq_label(freqs[i]),
                })
        peaks.sort(key=lambda p: -p["strength"])
        peaks = peaks[:10]
        # Band labels from Physics.UniversalFrequencyMapper (display taxonomy only)
        for p in peaks:
            try:
                ctx = UniversalFrequencyMapper.get_forensic_context(float(p["hz"]), medium="LAB")
                phys = ctx.get("physics") or ()
                band = phys[2] if len(phys) >= 3 else ""
                note = phys[3] if len(phys) >= 4 else ""
                opt = ctx.get("optical") or {}
                extra = []
                if band:
                    extra.append(str(band))
                if note:
                    extra.append(str(note))
                if opt.get("color"):
                    extra.append(f"optical~{opt.get('color')}")
                if opt.get("impact"):
                    extra.append(str(opt.get("impact")))
                if extra:
                    p["band_context"] = " · ".join(extra)
                    p["label"] = f"{p.get('label', '')} [{extra[0]}]"
            except Exception:
                p.setdefault("band_context", "")

    # Position-axis energy (for top density panel)
    pos_energy = grid.mean(axis=0)
    if pos_energy.max() > 0:
        pos_energy = pos_energy / pos_energy.max()

    path = None
    if out_png:
        os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
        fig = plt.figure(figsize=(14, 7.2), facecolor="#0b1220")
        gs = GridSpec(3, 1, height_ratios=[1.4, 3.2, 0.55], hspace=0.12)
        ax_top = fig.add_subplot(gs[0])
        ax_mid = fig.add_subplot(gs[1], sharex=ax_top)
        ax_bot = fig.add_subplot(gs[2], sharex=ax_top)

        for ax in (ax_top, ax_mid, ax_bot):
            ax.set_facecolor("#0b1220")
            ax.tick_params(colors="#94a3b8")
            for sp in ax.spines.values():
                sp.set_color("#334155")

        x = np.linspace(0, 1, n_time)
        ax_top.plot(x, pos_energy, color="#4AF626", lw=1.3, alpha=0.9)
        ax_top.fill_between(x, 0, pos_energy, color="#4AF626", alpha=0.15)
        ax_top.set_ylim(0, 1.15)
        ax_top.set_ylabel("Feature\nenergy", color="#94a3b8", fontsize=8)
        ax_top.set_title(
            title or "Sequence-feature spectrogram (NOT RF measurement)",
            color="#e2e8f0",
            fontsize=12,
            pad=8,
        )
        ax_top.set_xticks([])

        im = ax_mid.imshow(
            grid,
            aspect="auto",
            origin="lower",
            cmap="magma",
            extent=[0, 1, 0, 1],
            interpolation="bilinear",
        )
        yticks = np.linspace(0, 1, 7)
        ylabels = [
            _freq_label(10 ** (np.log10(f_min) + y * (np.log10(f_max) - np.log10(f_min))))
            for y in yticks
        ]
        ax_mid.set_yticks(yticks)
        ax_mid.set_yticklabels(ylabels, fontsize=7, color="#94a3b8")
        ax_mid.set_ylabel("Display frequency (log)", color="#94a3b8", fontsize=8)
        ax_mid.set_xticks([])

        # Marker track
        ax_bot.set_ylim(0, 1)
        ax_bot.set_yticks([])
        if marker_positions and L > 0:
            for m in marker_positions:
                try:
                    pos = float(m.get("pos", 0))
                except (TypeError, ValueError):
                    continue
                xf = pos / max(L, 1)
                if 0 <= xf <= 1:
                    ax_bot.axvline(xf, color=m.get("color") or "#22d3ee", alpha=0.7, lw=0.8)
        ax_bot.set_xlabel("Sequence position (start → end)", color="#94a3b8", fontsize=9)
        ax_bot.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax_bot.set_xticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=8, color="#94a3b8")
        ax_bot.text(0.01, 0.35, "markers", color="#64748b", fontsize=7, transform=ax_bot.transAxes)

        fig.text(0.01, 0.005, disclaimer, fontsize=6.5, color="#64748b", wrap=True)
        plt.savefig(out_png, dpi=130, facecolor=fig.get_facecolor(), bbox_inches="tight")
        plt.close(fig)
        path = out_png

    return {
        "peaks": peaks,
        "path": path,
        "disclaimer": disclaimer,
        "grid_max": float(grid.max()) if grid.size else 0.0,
        "n_windows": n_time,
        "seq_len_used": L,
    }



def _freq_label(hz):
    hz = float(hz)
    if hz >= 1e14:
        return f"{hz/1e14:.1f}e14 Hz (opt.)"
    if hz >= 1e12:
        return f"{hz/1e12:.2f} THz"
    if hz >= 1e9:
        return f"{hz/1e9:.2f} GHz"
    if hz >= 1e6:
        return f"{hz/1e6:.2f} MHz"
    if hz >= 1e3:
        return f"{hz/1e3:.1f} kHz"
    return f"{hz:.1f} Hz"
