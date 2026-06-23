"""
common.py  --  Shared constants and helpers for xs.py and xs_arithmetic.py
=========================================================================

Material constants, atomic masses, BISO geometry, blanket material builders,
cross-section library helpers, and flux-spectrum loaders used by both scripts.
"""

import os
import sys
import csv
import numpy as np

import openmc
import openmc.data


# =============================================================================
# MATERIAL CONSTANTS  (from Python/parameters.py and Python/utilities.py)
# =============================================================================

TARGET_TEMP_K  = 900                # match the model temperature

# --- FLiBe breeder ----------------------------------------------------------
DENSITY_FLIBE = 1.9505   # [g/cm3]  2(LiF)-BeF2 EOS at 900 K, 101 kPa
ENRICH_FLIBE  = 7.50     # [at%]    Li-6 enrichment
DENSITY_UF4   = 6.88     # [g/cm3]  UF4 molar-volume density in FLiBe at 900 K
DENSITY_ThF4  = 6.61     # [g/cm3]  ThF4 density in FLiBe (Cantor 1966/1973)

# --- DCLL (Pb-17Li) blanket -------------------------------------------------
DENSITY_DCLL  = 9.40     # [g/cm3]  Pb-17Li
ENRICH_DCLL   = 90.00    # [at%]    Li-6 enrichment in the Pb-17Li
DCLL_VF_FS_NOM = 0.019   # ferritic steel (F82H) volume fraction of blanket
DCLL_VF_LL_NOM = 0.808   # lead-lithium (breeder) volume fraction of blanket
DCLL_VF_SI_NOM = 0.076   # silicon carbide volume fraction of blanket
DCLL_VF_HE_NOM = 0.097   # helium volume fraction of blanket
DCLL_BLANKET_NDENS = 0.03541604638  # [atom/b-cm]

# --- HCPB (He-cooled pebble bed, Li4SiO4 + Be) blanket ---------------------
DENSITY_LI4SIO4 = 2.42     # [g/cm3]  Li4SiO4 ceramic
DENSITY_BE      = 1.85     # [g/cm3]  Beryllium metal
ENRICH_HCPB     = 60.00    # [at%]    Li-6 enrichment in Li4SiO4 ceramic
HCPB_VF_LI_NOM = 0.1304   # Li4SiO4 (ceramic)
HCPB_VF_BE_NOM = 0.3790   # Be (metal)
HCPB_VF_EU_NOM = 0.1176   # Eurofer
HCPB_VF_HE_NOM = 0.3730   # He-4 (gas)

# --- fertile materials ------------------------------------------------------
ENRICH_U     =  0.71     # [wt%]   U-235 in the UO2 BISO kernel
DENSITY_UO2  = 10.50     # [g/cm3]
DENSITY_ThO2 = 10.00     # [g/cm3]
DENSITY_SIC  =  3.20     # [g/cm3]

# --- BISO particle geometry -------------------------------------------------
BISO_KERNEL_RADIUS = 0.04   # [cm]
BISO_RADIUS        = 0.05   # [cm]
BISO_VOLUME          = (4/3) * 3.14159265359 * BISO_RADIUS**3
KERNEL_VOLUME        = (4/3) * 3.14159265359 * BISO_KERNEL_RADIUS**3
BISO_KERNEL_VOL_FRAC = KERNEL_VOLUME / BISO_VOLUME   # = 0.512
BISO_COAT_VOL_FRAC   = 1.0 - BISO_KERNEL_VOL_FRAC    # = 0.488

# --- atomic masses (from atom.kaeri.re.kr) ----------------------------------
AMU_F19   =  18.9984
AMU_O     =  15.999
AMU_U     = 238.02891          # natural-U average mass
AMU_U238  = 238.05078826
AMU_Th232 = 232.0381
AMU_UF4   = AMU_U + 4 * AMU_F19
AMU_ThF4  = AMU_Th232 + 4 * AMU_F19
AMU_UO2   = AMU_U + 2 * AMU_O
AMU_ThO2  = AMU_Th232 + 2 * AMU_O
XF4_X_RATIO  = AMU_UF4 / AMU_U238    # ~1.31914
XF4_Th_RATIO = AMU_ThF4 / AMU_Th232  # ~1.32691


# =============================================================================
# FERTILE-LOADING HELPERS  (from Python/utilities.py)
# =============================================================================

def fertile_kgm3_to_biso_per_cc(fertile_kgm3, fertile_isotope='U238'):
    """kg/m3 of fertile isotope -> number of BISO particles per cm3 of breeder."""
    if fertile_isotope == 'U238':
        return (fertile_kgm3 * AMU_UO2 / AMU_U238
                / KERNEL_VOLUME / DENSITY_UO2 / 100**3 * 1000)
    elif fertile_isotope == 'Th232':
        return (fertile_kgm3 * AMU_ThO2 / AMU_Th232
                / KERNEL_VOLUME / DENSITY_ThO2 / 100**3 * 1000)
    sys.exit(f"Unknown fertile_isotope {fertile_isotope!r}")


def calc_biso_breeder_vol_fracs(fertile_kgm3, fertile_isotope='U238'):
    """Volume fractions of BISO and breeder, relative to the nominal breeder volume.

    Returns (vf_biso, vf_breeder, biso_per_cc).
    """
    biso_per_cc     = fertile_kgm3_to_biso_per_cc(fertile_kgm3, fertile_isotope)
    vf_biso_breeder = biso_per_cc * BISO_VOLUME

    if vf_biso_breeder > 1.0:
        sys.exit(f"Loading {fertile_kgm3} kg/m3 exceeds what fits in the breeder "
                 f"volume (vf_biso_breeder = {vf_biso_breeder:.3f} > 1).")

    vf_biso    = vf_biso_breeder / (vf_biso_breeder + 1)
    vf_breeder = 1 / (vf_biso_breeder + 1)
    return vf_biso, vf_breeder, biso_per_cc


# =============================================================================
# MATERIAL CONSTRUCTION
# =============================================================================

def build_flibe_blanket(fertile_kgm3, fertile_isotope='U238'):
    """FLiBe + actinide-fluoride blanket.

    fertile_isotope='U238'  -> UF4 dissolved in 2(LiF)-BeF2 (mirrors flibe.py)
    fertile_isotope='Th232' -> ThF4 dissolved in 2(LiF)-BeF2

    Returns (material, vf_xf4).
    """
    breeder = openmc.Material(name='breeder', temperature=TARGET_TEMP_K)
    breeder.set_density('g/cm3', DENSITY_FLIBE)
    breeder.add_elements_from_formula('F4Li2Be', 'ao',
                                      enrichment_target='Li6',
                                      enrichment_type='ao',
                                      enrichment=ENRICH_FLIBE)

    if fertile_kgm3 <= 0.0:                       # pure FLiBe baseline
        breeder.name = 'blanket'
        return breeder, 0.0

    fertile = openmc.Material(name='fertile', temperature=TARGET_TEMP_K)
    if fertile_isotope == 'U238':
        fertile.add_elements_from_formula('UF4', 'ao')   # natural U (U-234/235/238)
        fertile.set_density('g/cm3', DENSITY_UF4)
        rho_xf4  = DENSITY_UF4
        xf4_ratio = XF4_X_RATIO
    elif fertile_isotope == 'Th232':
        fertile.add_elements_from_formula('ThF4', 'ao')  # natural Th (~100% Th-232)
        fertile.set_density('g/cm3', DENSITY_ThF4)
        rho_xf4  = DENSITY_ThF4
        xf4_ratio = XF4_Th_RATIO
    else:
        sys.exit(f"Unknown fertile_isotope {fertile_isotope!r}")

    vf_xf4   = fertile_kgm3 * xf4_ratio / (rho_xf4 * 1000.0)  # 1 g/cm3 = 1000 kg/m3
    vf_flibe = 1.0 - vf_xf4
    if vf_xf4 >= 1.0:
        sys.exit(f"Loading {fertile_kgm3} kg/m3 gives vf_xf4={vf_xf4:.3f} >= 1 -- "
                 f"more fertile fluoride than fits in the volume.")

    blanket = openmc.Material.mix_materials([breeder, fertile],
                                            [vf_flibe, vf_xf4], 'vo')
    blanket.name = 'blanket'
    blanket.temperature = TARGET_TEMP_K
    return blanket, vf_xf4


def build_dcll_blanket(fertile_kgm3, fertile_isotope='U238'):
    """DCLL Pb-17Li + BISO blanket (mirrors Python/dcll.py).

    Returns (material, vf_biso_bl), where vf_biso_bl is the BISO volume
    fraction of the whole blanket.
    """
    # --- Pb-17Li breeder ----------------------------------------------------
    pbli = openmc.Material(name='breeder', temperature=TARGET_TEMP_K)
    pbli.set_density('g/cm3', DENSITY_DCLL)
    pbli.add_element('Pb', 0.83, percent_type='ao')
    pbli.add_element('Li', 0.17, percent_type='ao',
                     enrichment_target='Li6', enrichment_type='ao',
                     enrichment=ENRICH_DCLL)

    # --- F82H steel ---------------------------------------------------------
    f82h = openmc.Material(name='f82h', temperature=TARGET_TEMP_K)
    f82h.add_element('Fe', 89.3686, percent_type='wo')
    f82h.add_element('C',   0.1000, percent_type='wo')
    f82h.add_element('Si',  0.1000, percent_type='wo')
    f82h.add_element('Mn',  0.1300, percent_type='wo')
    f82h.add_element('Cr',  8.1600, percent_type='wo')
    f82h.add_element('W',   1.9400, percent_type='wo')
    f82h.add_element('V',   0.2000, percent_type='wo')
    f82h.add_element('N',   0.0014, percent_type='wo')
    f82h.set_density('g/cm3', 7.78)

    # --- SiC and He ---------------------------------------------------------
    sic = openmc.Material(name='SiC', temperature=TARGET_TEMP_K)
    sic.add_elements_from_formula('SiC')
    sic.set_density('g/cm3', DENSITY_SIC)

    he = openmc.Material(name='helium', temperature=TARGET_TEMP_K)
    he.set_density('g/cm3', 0.004)
    he.add_element('He', 1)

    # --- BISO fuel kernel + SiC coating -------------------------------------
    if fertile_isotope == 'U238':
        kernel = openmc.Material(name='UO2', temperature=TARGET_TEMP_K)
        kernel.add_elements_from_formula('UO2', enrichment=ENRICH_U)
        kernel.set_density('g/cm3', DENSITY_UO2)
    elif fertile_isotope == 'Th232':
        kernel = openmc.Material(name='ThO2', temperature=TARGET_TEMP_K)
        kernel.add_elements_from_formula('ThO2')
        kernel.set_density('g/cm3', DENSITY_ThO2)
    else:
        sys.exit(f"Unknown fertile_isotope {fertile_isotope!r}")

    biso = openmc.Material.mix_materials([kernel, sic],
                                         [BISO_KERNEL_VOL_FRAC, BISO_COAT_VOL_FRAC],
                                         'vo')

    # --- volume fractions of each component in the blanket ------------------
    vf_biso_br, vf_pbli_br, _ = calc_biso_breeder_vol_fracs(fertile_kgm3,
                                                            fertile_isotope)
    vf_biso_bl = vf_biso_br * DCLL_VF_LL_NOM
    vf_pbli_bl = vf_pbli_br * DCLL_VF_LL_NOM

    blanket = openmc.Material.mix_materials(
        [biso, pbli, f82h, sic, he],
        [vf_biso_bl, vf_pbli_bl, DCLL_VF_FS_NOM, DCLL_VF_SI_NOM, DCLL_VF_HE_NOM],
        'vo')
    blanket.set_density('atom/b-cm', DCLL_BLANKET_NDENS)
    blanket.temperature = TARGET_TEMP_K
    blanket.name = 'blanket'
    return blanket, vf_biso_bl


def build_hcpb_blanket(fertile_kgm3, fertile_isotope='U238'):
    """HCPB Li4SiO4+Be + BISO blanket (mirrors Python/hcpb.py).

    Returns (material, vf_biso_bl), where vf_biso_bl is the BISO volume
    fraction of the whole blanket.
    """
    # --- Li4SiO4 ceramic breeder ------------------------------------------------
    li4sio4 = openmc.Material(name='Li4SiO4', temperature=TARGET_TEMP_K)
    li4sio4.set_density('g/cm3', DENSITY_LI4SIO4)
    li4sio4.add_elements_from_formula('Li4SiO4',
                                       enrichment_target='Li6',
                                       enrichment_type='ao',
                                       enrichment=ENRICH_HCPB)

    # --- Beryllium neutron multiplier -------------------------------------------
    be = openmc.Material(name='Beryllium', temperature=TARGET_TEMP_K)
    be.set_density('g/cm3', DENSITY_BE)
    be.add_element('Be', 1, percent_type='wo')

    # --- Eurofer steel ----------------------------------------------------------
    eurofer = openmc.Material(name='Eurofer', temperature=TARGET_TEMP_K)
    eurofer.set_density('g/cm3', 7.8)
    eurofer.add_element('Fe', 89.36, percent_type='wo')
    eurofer.add_element('C',   0.11, percent_type='wo')
    eurofer.add_element('Cr',  9.00, percent_type='wo')
    eurofer.add_element('W',   1.10, percent_type='wo')
    eurofer.add_element('Mn',  0.40, percent_type='wo')
    eurofer.add_element('N',   0.03, percent_type='wo')

    # --- Helium coolant ---------------------------------------------------------
    he = openmc.Material(name='helium', temperature=TARGET_TEMP_K)
    he.set_density('atom/b-cm', 0.00049800000)
    he.add_element('He', 1)

    # --- SiC coating ------------------------------------------------------------
    sic = openmc.Material(name='SiC', temperature=TARGET_TEMP_K)
    sic.add_elements_from_formula('SiC')
    sic.set_density('g/cm3', DENSITY_SIC)

    # --- BISO fuel kernel -------------------------------------------------------
    if fertile_isotope == 'U238':
        kernel = openmc.Material(name='UO2', temperature=TARGET_TEMP_K)
        kernel.add_elements_from_formula('UO2', enrichment=ENRICH_U)
        kernel.set_density('g/cm3', DENSITY_UO2)
    elif fertile_isotope == 'Th232':
        kernel = openmc.Material(name='ThO2', temperature=TARGET_TEMP_K)
        kernel.add_elements_from_formula('ThO2')
        kernel.set_density('g/cm3', DENSITY_ThO2)
    else:
        sys.exit(f"Unknown fertile_isotope {fertile_isotope!r}")

    biso = openmc.Material.mix_materials([kernel, sic],
                                          [BISO_KERNEL_VOL_FRAC, BISO_COAT_VOL_FRAC],
                                          'vo')

    # --- volume fractions (mirrors hcpb.py) ------------------------------------
    vf_biso_br, vf_libe_br, _ = calc_biso_breeder_vol_fracs(fertile_kgm3,
                                                              fertile_isotope)
    breeder_nom = HCPB_VF_LI_NOM + HCPB_VF_BE_NOM   # 0.5094

    vf_biso_bl = vf_biso_br * breeder_nom
    vf_libe_bl = vf_libe_br * breeder_nom
    vf_li_bl   = vf_libe_bl * HCPB_VF_LI_NOM / breeder_nom
    vf_be_bl   = vf_libe_bl * HCPB_VF_BE_NOM / breeder_nom

    blanket = openmc.Material.mix_materials(
        [biso, li4sio4, be, eurofer, he],
        [vf_biso_bl, vf_li_bl, vf_be_bl, HCPB_VF_EU_NOM, HCPB_VF_HE_NOM],
        'vo')
    blanket.temperature = TARGET_TEMP_K
    blanket.name = 'blanket'
    return blanket, vf_biso_bl


def build_blanket(blanket, fertile_kgm3, fertile_isotope='U238'):
    """Dispatch to the right builder.  Returns (material, vol_frac_of_fuel)."""
    if blanket == 'FLiBe':
        return build_flibe_blanket(fertile_kgm3, fertile_isotope)
    if blanket == 'HCPB':
        return build_hcpb_blanket(fertile_kgm3, fertile_isotope)
    if blanket == 'DCLL':
        return build_dcll_blanket(fertile_kgm3, fertile_isotope)
    sys.exit(f"Unknown blanket {blanket!r} (expected 'FLiBe', 'HCPB', or 'DCLL').")


def atom_densities(mat):
    """{nuclide: atom density [atom/b-cm]} -- robust to OpenMC version differences."""
    raw = mat.get_nuclide_atom_densities()
    out = {}
    for k, v in raw.items():
        out[k] = float(v[1]) if isinstance(v, (tuple, list)) else float(v)
    return out


# =============================================================================
# CROSS-SECTION DATA HELPERS
# =============================================================================

def open_library():
    xs_xml = (openmc.config.get('cross_sections')
              if hasattr(openmc, 'config') else None) or os.environ.get('OPENMC_CROSS_SECTIONS')
    if not xs_xml or not os.path.isfile(str(xs_xml)):
        sys.exit("Could not find cross_sections.xml. Set OPENMC_CROSS_SECTIONS "
                 "or openmc.config['cross_sections'] to your ENDF/B-VIII.0 library.")
    return openmc.data.DataLibrary.from_xml(xs_xml)


def find_path(lib, name):
    """Path to the neutron HDF5 file for `name`, or None."""
    for entry in lib.libraries:
        if entry.get('type') == 'neutron' and name in entry.get('materials', []):
            return entry['path']
    return None


def pick_temp(nuc, target_k=TARGET_TEMP_K):
    """Temperature string ('900K') closest to target among those available."""
    temps = list(nuc.temperatures)
    if not temps:
        raise RuntimeError(f"{nuc.name}: no temperatures in data file")
    want = f"{int(target_k)}K"
    if want in temps:
        return want
    nearest = min(temps, key=lambda t: abs(int(t.rstrip('K')) - target_k))
    print(f"  ! {nuc.name}: {want} unavailable, using {nearest}")
    return nearest


# =============================================================================
# FLUX-SPECTRUM HELPERS
# =============================================================================

FLUX_DIR = './Data'
FLUX_FILES = {
    'FLiBe': 'FLiBe_900K_Li07.5_U238_flux.csv',
    'HCPB':  'HCPB_900K_Li60.0_U238_flux.csv',
    'DCLL':  'DCLL_900K_Li90.0_U238_flux.csv',
}


def load_flux_spectra(flux_dir=FLUX_DIR, flux_files=FLUX_FILES):
    """Read the flux CSVs.  Returns {blanket: {loading: (E_mid, phi)}} or {}
    for blankets whose file is missing / None.  Only non-zero flux bins are
    kept (zero-flux bins cannot contribute to the weighted average)."""
    spectra = {}
    for bname, fname in flux_files.items():
        if fname is None:
            continue
        path = os.path.join(flux_dir, fname)
        if not os.path.isfile(path):
            print(f"  ! flux file not found: {path}  -- falling back to "
                  f"point evaluation for {bname}")
            continue
        rows_by_loading = {}
        with open(path) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                ld = float(row['fertile_kg/m3'])
                if ld not in rows_by_loading:
                    rows_by_loading[ld] = ([], [])
                rows_by_loading[ld][0].append(float(row['energy mid [eV]']))
                rows_by_loading[ld][1].append(float(row['mean']))
        blk = {}
        for ld, (emids, phis) in rows_by_loading.items():
            E = np.array(emids)
            P = np.array(phis)
            order = np.argsort(E)
            E, P = E[order], P[order]
            keep = P > 0
            blk[ld] = (E[keep], P[keep])
        spectra[bname] = blk
    return spectra


def match_loading(available, target, tol=1.0):
    """Closest loading in *available* within *tol* of *target*, or None."""
    best = min(available, key=lambda x: abs(x - target))
    return best if abs(best - target) <= tol else None
