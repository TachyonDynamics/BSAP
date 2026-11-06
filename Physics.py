import math

class UniversalPhysics:
    # --- CONSTANTS ---
    C = 299792458            
    H = 6.62607015e-34       
    K_B = 1.380649e-23       
    N_A = 6.02214076e23      
    E_Q = 1.60217663e-19     
    DNA_BOND_K = 500.0       
    SPECIFIC_HEAT_DNA = 1.45 
    WATER_DIELECTRIC = 80.1  

    @staticmethod
    def get_photon_energy(hz):
        ev = (UniversalPhysics.H * hz) / UniversalPhysics.E_Q
        return {"eV": round(ev, 4), "threat": "IONIZING" if ev > 10 else "NON-IONIZING"}

    @staticmethod
    def get_vibrational_resonance(mass_daltons):
        m_kg = mass_daltons * 1.660539e-27
        res_hz = (1 / (2 * math.pi)) * math.sqrt(UniversalPhysics.DNA_BOND_K / m_kg)
        return round(res_hz, 2)
        
    @staticmethod
    def get_thermal_velocity(hz, temp_k=310.15):
        """
        Calculates if the signal frequency exceeds the 'Thermal Noise' floor.
        If hz < thermal_noise, the signal is lost in body heat.
        """
        thermal_noise_limit = (UniversalPhysics.K_B * temp_k) / UniversalPhysics.H
        is_detectable = hz > thermal_noise_limit
        return {"limit_hz": round(thermal_noise_limit, 2), "detectable": is_detectable}

    @staticmethod
    def calculate_coulomb_force(net_charge, length_bp):
        """
        Applies Coulomb's Law to calculate the electrostatic force 
        across the DNA fragment.
        """
        if length_bp == 0: return 0.0
        
        # DNA length is approx 0.34 nm per base pair
        r_meters = length_bp * 0.34e-9
        
        # Elementary charge conversion
        q_total = net_charge * UniversalPhysics.E_Q
        
        # Coulomb's constant (k_e) = 8.987e9 N*m^2/C^2
        k_e = 8.9875517923e9
        
        # Assuming the charge is split across the dipole poles
        q1 = q_total / 2
        q2 = -(q_total / 2)
        
        force_newtons = k_e * abs(q1 * q2) / (r_meters**2)
        return force_newtons

    @staticmethod
    def calculate_net_charge(sequence, ph=7.4):
        charge = (1 / (1 + 10**(ph - BiochemicalEngine.AA_PKA['N_term']))) - \
                 (1 / (1 + 10**(BiochemicalEngine.AA_PKA['C_term'] - ph)))
        for aa in sequence:
            if aa in ['R', 'K', 'H']:
                charge += 1 / (1 + 10**(ph - BiochemicalEngine.AA_PKA[aa]))
            elif aa in ['D', 'E', 'C', 'Y']:
                charge -= 1 / (1 + 10**(BiochemicalEngine.AA_PKA[aa] - ph))
        return round(charge, 4)
        
class UniversalFrequencyMapper:  # Used by Sonification feature-spectrogram peak labels (display taxonomy)
    SPECTRUM_MAP = {
        "EXTREMELY_LOW_GEOLOGICAL": {"range": (0.0001, 3.0), "markers": {0.0: "DC_BIOPOTENTIAL_STATIC", 1.5: "DELTA_DEEP_SLEEP"}},
        "BRAIN_CONSCIOUSNESS": {"range": (3.0, 45.0), "markers": {5.5: "THETA_HYPNAGOGIC", 7.83: "SCHUMANN_EARTH", 10.5: "ALPHA_RELAXED", 20.0: "BETA_ACTIVE", 40.0: "GAMMA_COHERENCE"}},
        "POWER_GRID_INDUCTION": {"range": (45.0, 100.0),"markers": {50.0: "EURO_GRID_HUM", 60.0: "US_GRID_HUM"}},
        "ACOUSTIC_ULTRASONIC": {"range": (100.0, 20000.0), "markers": {432.0: "VERDI_TUNING", 528.0: "SOLFEGGIO_REPAIR"}},
        "RADIO_COMM_INFRA": {"range": (20000.0, 3e11), "markers": {915e6: "SMART_METER_RF", 2.44e9: "WIFI_BT_EMF", 12e9: "STARLINK_SPACE", 60e9: "OXYGEN_RESONANCE"}},
        "OPTICAL_BIO_GROWTH": {"range": (3e11, 7.5e14), "markers": {4.3e14: "697nm_DEEP_RED", 5.2e14: "576nm_GREEN_REFL", 6.8e14: "440nm_DEEP_BLUE", 7.1e14: "422nm_VIOLET"}},
        "ULTRAVIOLET_DNA_DAMAGE": {"range": (7.5e14, 2.4e15), "markers": {1.18e15: "UVC_254nm_PHOTOLESION"}},
        "IONIZING_RADIATION": {"range": (2.4e15, 1e22), "markers": {2.4e15: "IONIZATION_THRESHOLD_10eV", 3e19: "COSMIC_RAY_STRIKE"}}
    }

    # --- OPERATIONAL DOMAINS ---
    RADIO_DOMAIN_MAP = {
        "SEA_MARITIME": {"range": (2e4, 5e5), "markers": {24000.0: "NAA_CUTLER_SUB_COMM", 100000.0: "LORAN_C_NAV"}},
        "LAND_AM_HAM":  {"range": (5e5, 3e7), "markers": {4625000.0: "UVB_76_BUZZER", 14e6: "HAM_20M_BAND", 27e6: "CB_RADIO_CITIZENS"}},
        "AIR_AVIONICS": {"range": (1.08e8, 1.37e8), "markers": {121.5e6: "AIRCRAFT_EMERGENCY", 109e6: "ADS_B_TRANSPONDER"}},
        "LAND_MOBILE":  {"range": (3e8, 6e9), "markers": {915e6: "INDUSTRIAL_ISM_BAND", 2.44e9: "WIFI_BT_CONSUMER", 5.8e9: "WIFI_5G_UPPER"}},
        "SPACE_SAT":    {"range": (6e9, 1e11), "markers": {12e9: "STARLINK_S_BAND", 24e9: "POLICE_RADAR_K", 60e9: "5G_OXYGEN_ABSORB"}}
    }

    # --- HIGH-RES OPTICAL MAP ---
    OPTICAL_MAP = {
        "FAR_RED": {"nm": 750, "hz": 4e14, "tag": "FAR_RED_GROWTH"},
        "RED_PEAK": {"nm": 660, "hz": 4.5e14, "tag": "CHLOROPHYLL_A_MAX"},
        "ORANGE": {"nm": 610, "hz": 4.9e14, "tag": "BIOLUM_ORANGE"},
        "YELLOW": {"nm": 580, "hz": 5.1e14, "tag": "REFLECTANCE_PEAK"},
        "GREEN": {"nm": 530, "hz": 5.6e14, "tag": "GREEN_GAP_REFLECT"},
        "BLUE_PEAK": {"nm": 450, "hz": 6.6e14, "tag": "CIRCADIAN_TRIGGER"},
        "VIOLET": {"nm": 400, "hz": 7.5e14, "tag": "DNA_STRESS_THRESHOLD"},
        "UV_C": {"nm": 254, "hz": 1.1e15, "tag": "DNA_SHATTER_UV_C"}
    }

    # --- SUB-SURFACE ACOUSTICS ---
    SONAR_MAP = {
        "LF_SONAR": {"range": (100.0, 1000.0), "markers": {450.0: "SURTASS_LFA_PING"}},
        "MF_SONAR": {"range": (1000.0, 10000.0), "markers": {3500.0: "AN_SQS_53C_TACTICAL"}},
        "HF_SONAR": {"range": (20000.0, 500000.0), "markers": {200000.0: "COMMERCIAL_ECHOSOUNDER"}}
    }

    @staticmethod
    def analyze_signal_physics(hz):
        abs_hz = abs(hz)
        direction = "FORWARD" if hz >= 0 else "REVERSE"
        if 45.0 <= abs_hz <= 65.0:
            return direction, abs_hz, "AC (GRID_INDUCTION)", "Power Grid Interference"
        return direction, abs_hz, "AC (OSCILLATORY)", "Active Signal"

    @staticmethod
    def get_medium_multiplier(environment="LAB"):
        """Adjusts resonance for environmental density."""
        multipliers = {"VACUUM": 1.0, "AIR": 0.9997, "WATER": 0.75, "BONE": 0.45, "LAB": 0.85}
        return multipliers.get(environment, 1.0)

    @staticmethod
    def get_optical_context(hz):
        """Translates Hz to Nanometers for light-based forensics."""
        if hz < 3e11: return None
        c = 299792458
        nm = (c / hz) * 1e9
        for label, data in UniversalFrequencyMapper.OPTICAL_MAP.items():
            if abs(nm - data['nm']) / data['nm'] < 0.05:
                return {"wavelength": f"{round(nm, 2)}nm", "color": label, "impact": data['tag']}
        return {"wavelength": f"{round(nm, 2)}nm", "color": "OUT_OF_BAND"}

    @staticmethod
    def detect_sonar_interference(hz, medium="WATER"):
        """Checks for acoustic sonar pings in fluid mediums."""
        if medium != "WATER": return None
        for band, data in UniversalFrequencyMapper.SONAR_MAP.items():
            if data['range'][0] <= hz <= data['range'][1]:
                for p_hz, name in data['markers'].items():
                    if abs(hz - p_hz) / p_hz < 0.02:
                        return f"MARITIME_ALERT: {name} Sonar"
        return None

    @staticmethod
    def get_forensic_context(hz, medium="LAB"):
        """The Master Classifier: Flowing Hz through all sub-maps."""
        abs_hz = abs(hz)
        results = {"physics": UniversalFrequencyMapper.analyze_signal_physics(hz)}
        
        # Check Optical
        opt = UniversalFrequencyMapper.get_optical_context(abs_hz)
        if opt: results["optical"] = opt
        
        # Check Sonar
        snr = UniversalFrequencyMapper.detect_sonar_interference(abs_hz, medium)
        if snr: results["sonar"] = snr
        
        return results

    @staticmethod
    def get_acoustic_transposition(hz):
        """
        Transposes biological resonances into the audible spectrum.
        Ensures results are mathematically harmonized within 20Hz - 20kHz.
        """
        if hz <= 0: return 0
        
        # Initial 40-octave drop
        audio_hz = hz / (2**40)
        
        # Automatic Octave Sifting (Audible Normalization)
        # Shift UP if too low
        while audio_hz > 0 and audio_hz < 20:
            audio_hz *= 2
        
        # Shift DOWN if too high
        while audio_hz > 20000:
            audio_hz /= 2
            
        return round(audio_hz, 2)

class BiochemicalEngine:
    # --- ATOMIC WEIGHTS ---
    ELEMENTS = {'H': 1.00784, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'P': 30.9738, 'S': 32.065}

    # --- DNA & RNA ATOMIC MAPS ---
    NUCLEOTIDE_MAP = {
        # Residue in a chain = nucleotide monophosphate minus H2O (same convention as PROTEIN_MAP).
        'DNA': {
            'A': {'C': 10, 'H': 12, 'N': 5, 'O': 5, 'P': 1},   # dAMP residue 313.21 Da
            'T': {'C': 10, 'H': 13, 'N': 2, 'O': 7, 'P': 1},   # dTMP residue 304.20 Da
            'C': {'C': 9,  'H': 12, 'N': 3, 'O': 6, 'P': 1},   # dCMP residue 289.18 Da
            'G': {'C': 10, 'H': 12, 'N': 5, 'O': 6, 'P': 1}    # dGMP residue 329.21 Da
        },
        'RNA': {
            'A': {'C': 10, 'H': 12, 'N': 5, 'O': 6, 'P': 1},   # AMP residue 329.21 Da
            'U': {'C': 9,  'H': 11, 'N': 2, 'O': 8, 'P': 1},   # UMP residue 306.17 Da
            'C': {'C': 9,  'H': 12, 'N': 3, 'O': 7, 'P': 1},   # CMP residue 305.18 Da
            'G': {'C': 10, 'H': 12, 'N': 5, 'O': 7, 'P': 1}    # GMP residue 345.21 Da
        }
    }

    # Protein Residue Formulas (Formulas minus H2O for peptide bonding)
    PROTEIN_MAP = {
        'A': {'name': 'Alanine',       'formula': {'C': 3, 'H': 5, 'N': 1, 'O': 1}},
        'R': {'name': 'Arginine',      'formula': {'C': 6, 'H': 12, 'N': 4, 'O': 1}},
        'N': {'name': 'Asparagine',    'formula': {'C': 4, 'H': 6, 'N': 2, 'O': 2}},
        'D': {'name': 'Aspartic Acid', 'formula': {'C': 4, 'H': 5, 'N': 1, 'O': 3}},
        'C': {'name': 'Cysteine',      'formula': {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'S': 1}},
        'E': {'name': 'Glutamic Acid', 'formula': {'C': 5, 'H': 7, 'N': 1, 'O': 3}},
        'Q': {'name': 'Glutamine',     'formula': {'C': 5, 'H': 8, 'N': 2, 'O': 2}},
        'G': {'name': 'Glycine',       'formula': {'C': 2, 'H': 3, 'N': 1, 'O': 1}},
        'H': {'name': 'Histidine',     'formula': {'C': 6, 'H': 7, 'N': 3, 'O': 1}},
        'I': {'name': 'Isoleucine',    'formula': {'C': 6, 'H': 11, 'N': 1, 'O': 1}},
        'L': {'name': 'Leucine',       'formula': {'C': 6, 'H': 11, 'N': 1, 'O': 1}},
        'K': {'name': 'Lysine',        'formula': {'C': 6, 'H': 12, 'N': 2, 'O': 1}},
        'M': {'name': 'Methionine',    'formula': {'C': 5, 'H': 9, 'N': 1, 'O': 1, 'S': 1}},
        'F': {'name': 'Phenylalanine', 'formula': {'C': 9, 'H': 9, 'N': 1, 'O': 1}},
        'P': {'name': 'Proline',       'formula': {'C': 5, 'H': 7, 'N': 1, 'O': 1}},
        'S': {'name': 'Serine',        'formula': {'C': 3, 'H': 5, 'N': 1, 'O': 2}},
        'T': {'name': 'Threonine',     'formula': {'C': 4, 'H': 7, 'N': 1, 'O': 2}},
        'W': {'name': 'Tryptophan',    'formula': {'C': 11, 'H': 10, 'N': 2, 'O': 1}},
        'Y': {'name': 'Tyrosine',      'formula': {'C': 9, 'H': 9, 'N': 1, 'O': 2}},
        'V': {'name': 'Valine',        'formula': {'C': 5, 'H': 9, 'N': 1, 'O': 1}},
        '_': {'name': 'STOP',          'formula': {'C': 0, 'H': 0, 'N': 0, 'O': 0}}
    }

    CODON_TO_AA = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    # --- AMINO ACID pKa VALUES ---
    AA_PKA = {
        'R': 12.48, 'K': 10.53, 'H': 6.0,
        'D': 3.65,  'E': 4.25,  'C': 8.33, 'Y': 10.07,
        'N_term': 9.6, 'C_term': 2.34
    }

    # --- UNIFIED BIOME DATA ---
    # format: { "ph": value, "conductivity": S/m, "permittivity": epsilon_r }
    ORGAN_BIOMES = {
        "SKIN":        {"ph": 5.5, "conductivity": 0.0002, "permittivity": 2.0},
        "BLOOD":       {"ph": 7.4, "conductivity": 1.50,   "permittivity": 80.0},
        "BONE":        {"ph": 7.4, "conductivity": 0.02,   "permittivity": 5.0},
        "STOMACH":     {"ph": 2.0, "conductivity": 0.50,   "permittivity": 75.0},
        "LIVER":       {"ph": 7.2, "conductivity": 0.15,   "permittivity": 40.0},
        "AGAR_DISH":   {"ph": 7.0, "conductivity": 2.00,   "permittivity": 78.0},
        "FERMENT_VAT": {"ph": 3.5, "conductivity": 4.50,   "permittivity": 65.0},
        "DIRT_SOIL":   {"ph": 6.5, "conductivity": 0.01,   "permittivity": 4.0},
        "AIR_AEROSOL": {"ph": 7.0, "conductivity": 0.00,   "permittivity": 1.0}
    }

    # --- CODON SPECIAL DATA ---
    CODON_PHYSICS = {
        'NNN': {'mass': 110.0, 'resonance': 'UNKNOWN'},
        'MASKED': {'mass': 0.0, 'resonance': 'INERT'}
    }

    @staticmethod
    def get_molecular_mass(sequence):
        """Unified mass calculator for DNA, RNA, and Protein."""
        seq = sequence.upper()
        # Logic: If contains 'U', it's RNA. If only ATCG, it's DNA. Else, Protein.
        m_type = "RNA" if 'U' in seq else ("DNA" if all(b in 'ATCG' for b in seq) else "PROT")
        
        total_mass = 0
        if m_type in ['DNA', 'RNA']:
            for base in seq:
                # Uses the Nucleotide map
                formula = BiochemicalEngine.NUCLEOTIDE_MAP[m_type].get(base, {})
                total_mass += sum(BiochemicalEngine.ELEMENTS[a] * c for a, c in formula.items())
        else:
            for aa in seq:
                # Uses the Protein map
                formula = BiochemicalEngine.PROTEIN_MAP.get(aa, {}).get('formula', {})
                total_mass += sum(BiochemicalEngine.ELEMENTS[a] * c for a, c in formula.items())
        if seq:
            total_mass += 18.01528   # terminal H2O (chain = sum of residues + one water)        
        return round(total_mass, 4)

    @staticmethod
    def calculate_dielectric_heating(mass_daltons, freq_hz, exposure_seconds=60):
        m_g = mass_daltons * 1.66e-24
        absorption = (freq_hz * 1e-12) * 0.05 
        delta_t = (absorption * exposure_seconds) / (m_g * UniversalPhysics.SPECIFIC_HEAT_DNA)
        return round(delta_t, 6)

    @staticmethod
    def identify_molecule_type(sequence):
        """Auto-detects if sequence is DNA, RNA, or Protein."""
        seq = sequence.upper()
        if 'U' in seq: return "RNA"
        if all(base in 'ATCG' for base in seq): return "DNA"
        return "PROTEIN_OR_PEPTIDE"

    @staticmethod
    def get_atomic_mass(formula):
        """Calculates mass in Daltons for any given atomic formula."""
        return round(sum(BiochemicalEngine.ELEMENTS[atom] * count for atom, count in formula.items()), 4)

    @staticmethod
    def scan_sequence(dna_string):
        """The Master Scan: Returns atomic footprint and AA arrangement."""
        report = []
        total_footprint = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0, 'S': 0}
        aa_counts = {}

        for i in range(0, len(dna_string), 3):
            codon = dna_string[i:i+3].upper()
            if len(codon) < 3: break
            
            aa_char = BiochemicalEngine.CODON_TO_AA.get(codon, '_')
            aa_info = BiochemicalEngine.PROTEIN_MAP.get(aa_char) or {
                "name": "UNKNOWN", "formula": {}
            }
            
            # DNA Triplet Footprint
            dna_atoms = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'P': 0, 'S': 0}
            for base in codon:
                base_atoms = BiochemicalEngine.NUCLEOTIDE_MAP['DNA'].get(base, {})
                for atom, count in base_atoms.items():
                    dna_atoms[atom] += count
                    total_footprint[atom] += count

            aa_counts[aa_info['name']] = aa_counts.get(aa_info['name'], 0) + 1
            
            report.append({
                'codon': codon,
                'aa': aa_info['name'],
                'prot_mass': BiochemicalEngine.get_atomic_mass(aa_info['formula']),
                'dna_mass': BiochemicalEngine.get_atomic_mass(dna_atoms)
            })
        
        return report, total_footprint, aa_counts

    @staticmethod
    def get_static_zeta_potential(ph=7.4):
        """
        Calculates the expected electrostatic charge of the DNA backbone.
        Static logic: DNA is negatively charged due to the phosphate groups.
        """
        # DNA phosphate backbone has ~1 negative charge per base
        # At clinical pH (7.4), this is stable. 
        # Low pH (Stomach) neutralizes this, causing "clumping".
        if ph < 3.0: return "NEUTRALIZED_CLUMPING_RISK"
        return "NEGATIVE_CHARGE_ELECTROPHORETIC"

class AcousticForensics:  # Optional HEURISTIC via ENABLE_ACOUSTIC_HEURISTICS
    @staticmethod
    def calculate_absorption_peak(thickness_m, material_density, bulk_modulus):
        """
        Predicts the 'Null Signal' where the material absorbs sound energy.
        """
        try:
            speed_of_sound = math.sqrt(bulk_modulus / material_density)
            fundamental_hz = speed_of_sound / (2 * thickness_m)
            return round(fundamental_hz, 2)
        except ZeroDivisionError:
            return 0.0

    @staticmethod
    def get_absorption_coefficient(reflected_p, incident_p):
        """
        The math for your Phase-Inversion idea. 
        Calculates how much energy the 'Resonating Cavity' ate.
        """
        if incident_p == 0: return 0.0
        return round(1 - (reflected_p / incident_p)**2, 4)