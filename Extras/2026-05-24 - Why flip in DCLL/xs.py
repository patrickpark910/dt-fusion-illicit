#!/usr/bin/env python3
"""
xs.py  --  Elemental macroscopic cross sections of the FLiBe, HCPB, and DCLL blankets
=============================================================================

Plots the *macroscopic* cross section  Sigma(E)  [1/cm]  for the reactions
listed below, so you can rank which reactions matter most in each breeding
blanket.  Sigma for a single nuclide is  N_i * sigma_i(E); for an ELEMENT it is

        Sigma_element(E) = sum_i  N_i * sigma_i(E)        (sum over isotopes i)

i.e. the per-isotope macroscopic contributions are added (microscopic cross
sections of different isotopes cannot be added directly -- only the N_i*sigma_i
products can).  Uranium is reported as the element "U" (U-235 + U-238) and lead
as the element "Pb" (Pb-204/206/207/208); lithium is kept split into Li-6 / Li-7
because the (n,t) vs (n,Xt) tritium channels are the whole point of the study.

Each blanket material is built EXACTLY the way it is in the full D-shaped models:
  * FLiBe (Python/flibe.py): UF4 is volume-mixed into 2(LiF)-BeF2 with the rule
        vf_xf4 = fertile_kgm3 * (AMU_UF4/AMU_U238) / (rho_UF4 * 1000)
  * HCPB  (Python/hcpb.py): a 5-component volume mix of UO2/SiC BISO particles,
        Li4SiO4 ceramic, Be metal, Eurofer steel, and He, with the BISO loading
        set by calc_biso_breeder_vol_fracs(...) and the breeder fraction
        (Li4SiO4 + Be) replacing the nominal HCPB_VF_LI_NOM + HCPB_VF_BE_NOM.
  * DCLL  (Python/dcll.py): a 5-component volume mix of UO2/SiC BISO particles,
        Pb-17Li breeder, F82H steel, SiC and He, with the BISO loading set by
        calc_biso_breeder_vol_fracs(...) and the total atom density pinned to the
        fixed value from the MCNP model (set_density('atom/b-cm', ...)).
The per-nuclide atom densities are taken from the mixed material, so Sigma
reflects both the cross section AND the actual loading-dependent concentration.

Output: a single 3 x 2 grid (macro_xs_grid.png/.pdf).
        rows    = blankets   (FLiBe, HCPB, DCLL)
        columns = loadings   (FERTILE_CASES; 0.1 and 1000 kg/m3 by default)
A ranked table of Sigma at the source energy is also printed for each panel.

        A second 3 x 2 grid (macro_xs_ratio_grid.png/.pdf) shows the
        energy-dependent ratio  Sigma(n,elastic) / Sigma(n,gamma)  of the
        fertile element, with both loading cases overlaid on each panel:
            rows    = blankets   (FLiBe, HCPB, DCLL)
            columns = fertile    (U, Th)
        The U columns sum U-235 + U-238; the Th columns use Th-232 alone.
        FLiBe-Th uses ThF4 dissolved in 2(LiF)-BeF2 (analogous to UF4);
        HCPB-Th uses ThO2 BISO particles (as in hcpb.py);
        DCLL-Th uses ThO2 BISO particles (as in dcll.py).

This file is STANDALONE: every constant and helper it needs is copied in here,
so it does not import anything from the Python/ package.

Requirements: openmc (with openmc.data) + an ENDF/B-VIII.0 HDF5 library that
includes 900 K data.  Point OPENMC_CROSS_SECTIONS (env var) or
openmc.config['cross_sections'] at your cross_sections.xml -- the same library
the rest of your model uses.

    python xs.py

"""

import os
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import openmc
import openmc.data


# =============================================================================
# CONFIGURATION
# =============================================================================

BLANKETS       = ['FLiBe', 'HCPB', 'DCLL']  # rows of the grid
FERTILE_CASES  = [0.1, 1000.0]      # kg/m3 of fertile -> columns of the grid.
                                    # Add 0.0 for a pure-breeder column.
FERTILE_ISOTOPE = 'U238'            # 'U238' or 'Th232'
                                    #   FLiBe always uses natural U via UF4;
                                    #   this switch selects the DCLL BISO kernel
                                    #   (UO2 vs ThO2) and the fertile element.
TARGET_TEMP_K  = 900                # match the model temperature

E_MIN, E_MAX   = 1e0, 2e7           # eV  (model tally range is 1e-3..20e6 eV)
SIGMA_MIN      = 1e-6               # 1/cm  -- floor; below this is negligible
SIGMA_MAX      = 1e1               # 1/cm

E_MARK         = 14.0e6            # eV  -- the 14 MeV DT source energy (marked)
RANK_ENERGY    = 14.0e6            # eV  -- energy at which the ranking is taken

OUTDIR         = "."               # where to drop the PNG/PDF
OUTNAME_BASE   = "macro_xs"        # file becomes macro_xs_grid.{png,pdf}

# --- CSV export -------------------------------------------------------------
WRITE_CSV      = True              # also dump the plotted spectra to CSV
CSV_OUTDIR     = OUTDIR            # where the CSVs go (one per blanket+loading)
CSV_E_RANGE    = (E_MIN, E_MAX)    # energy window written to CSV, matching the
                                   # plotted x-axis.  Set to None to export each
                                   # curve's FULL native energy grid instead.
CSV_FLOAT_FMT  = "%.9g"            # number format for energies and sigmas
CSV_DOWNSAMPLE = 100               # keep ~1 of every N points of the shared
                                   # energy grid, so the CSV is ~1/N the size of
                                   # the full-resolution dump.  The two endpoints
                                   # are always kept and the native grid's
                                   # relative point density is preserved (the
                                   # resonance region stays the most sampled).
                                   # Set to 1 to write every grid point.

# --- Flux-spectrum CSVs (for flux-weighted ranking) -------------------------
# When present the printed ranking uses <Sigma>_phi instead of a point value
# at RANK_ENERGY.  Set a path to None to fall back to point evaluation.
FLUX_DIR = './Data'
FLUX_FILES = {
    'FLiBe': 'FLiBe_900K_Li07.5_U238_flux.csv',
    'HCPB':  'HCPB_900K_Li60.0_U238_flux.csv',
    'DCLL':  'DCLL_900K_Li90.0_U238_flux.csv',
}


# =============================================================================
# CONSTANTS  (copied from Python/parameters.py and Python/utilities.py)
# =============================================================================

# --- FLiBe breeder ----------------------------------------------------------
DENSITY_FLIBE = 1.9505   # [g/cm3]  2(LiF)-BeF2 EOS at 900 K, 101 kPa
ENRICH_FLIBE  = 7.50     # [at%]    Li-6 enrichment
DENSITY_UF4   = 6.88     # [g/cm3]  UF4 molar-volume density in FLiBe at 900 K
DENSITY_ThF4  = 6.32     # [g/cm3]  ThF4 solid density (900 K FLiBe value TBD)

# --- DCLL (Pb-17Li) blanket -------------------------------------------------
DENSITY_DCLL  = 9.40     # [g/cm3]  Pb-17Li
ENRICH_DCLL   = 90.00    # [at%]    Li-6 enrichment in the Pb-17Li
DCLL_VF_FS_NOM = 0.019   # ferritic steel (F82H) volume fraction of blanket
DCLL_VF_LL_NOM = 0.808   # lead-lithium (breeder) volume fraction of blanket
DCLL_VF_SI_NOM = 0.076   # silicon carbide volume fraction of blanket
DCLL_VF_HE_NOM = 0.097   # helium volume fraction of blanket
# Fixed total atom density of the homogenised DCLL blanket, from Glaser et al.
# (2025) MCNP.  dcll.py pins the mix to this value, so the per-nuclide atom
# densities come out as (atom fraction from the volume mix) * this total.
DCLL_BLANKET_NDENS = 0.03541604638  # [atom/b-cm]

# --- HCPB (He-cooled pebble bed, Li4SiO4 + Be) blanket ---------------------
DENSITY_LI4SIO4 = 2.42     # [g/cm3]  Li4SiO4 ceramic
DENSITY_BE      = 1.85     # [g/cm3]  Beryllium metal
ENRICH_HCPB     = 60.00    # [at%]    Li-6 enrichment in Li4SiO4 ceramic
# Volume fractions of material in breeding regions (Lu et al. 2017, Tb.2)
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

# --- atomic masses (utilities.py, from atom.kaeri.re.kr) ---------------------
AMU_F19   =  18.9984
AMU_O     =  15.999
AMU_U     = 238.02891          # natural-U average mass (used to build UF4)
AMU_U238  = 238.05078826
AMU_Th232 = 232.0381
AMU_UF4   = AMU_U + 4 * AMU_F19      # = 314.02251 g/mol
AMU_ThF4  = AMU_Th232 + 4 * AMU_F19  # = 307.9317  g/mol
AMU_UO2   = AMU_U + 2 * AMU_O        # = 270.02691 g/mol
AMU_ThO2  = AMU_Th232 + 2 * AMU_O    # = 264.0361  g/mol

# UF4 -> U-238 mass ratio used to convert kg(U-238)/m3 into a UF4 volume fraction
XF4_X_RATIO = AMU_UF4 / AMU_U238  # ~1.31914
# ThF4 -> Th-232 mass ratio (analogous)
XF4_Th_RATIO = AMU_ThF4 / AMU_Th232  # ~1.32691


# =============================================================================
# FERTILE-LOADING HELPERS  (copied verbatim from Python/utilities.py)
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
# REACTION REGISTRY
# =============================================================================
#
# For every nuclide, a list of (label, kind, mt):
#   kind == 'mt'   -> single ENDF MT number `mt`
#   kind == 'inel' -> total inelastic, built by SUMMING MT 51..91
#                     (OpenMC mis-handles MT=4 directly -- see reactor.py note --
#                      so we reconstruct it from the discrete-level + continuum
#                      channels, the same way the model's tallies do)
#   kind == 'abs'  -> absorption: tries MT 27, then 101, then the sum of the
#                     disappearance partials (102..117)
#
# Nuclides absent from a given blanket are skipped automatically (atom density 0),
# so the same registry serves all three blankets:
#   FLiBe  ->  F-19, Li-6, Li-7, Be-9, U-235, U-238
#   HCPB   ->  Li-6, Li-7, Be-9, O, Si, C, Fe, U-235, U-238  (+ Eurofer trace)
#   DCLL   ->  Li-6, Li-7, Pb, C, Si, Fe, O, U-235, U-238  (+ F82H trace)
# Isotopes within an element are summed into a single elemental Sigma at plot
# time: U = U-235+U-238, Pb = Pb-204/206/207/208, Fe = Fe-54/56/57/58,
# O = O-16/17/18, Si = Si-28/29/30, C = C-12/13  (see GROUPS below).

MT_ELASTIC, MT_2N, MT_FISSION, MT_GAMMA = 2, 16, 18, 102
MT_NT      = 105      # (n,t)  -- the clean two-body Li-6(n,t)alpha
MT_XT      = 205      # (n,Xt) -- total tritium production (Li-7 breakup channel)

REACTIONS = {
    'F19':   [('(n,elastic)',    'mt',  MT_ELASTIC),
              ('(n,absorption)', 'abs', None)],

    'Li6':   [('(n,elastic)',    'mt',   MT_ELASTIC),
              ('(n,inelastic)',  'inel', None),
              ('(n,t)',          'mt',   MT_NT)],

    'Li7':   [('(n,elastic)',    'mt',   MT_ELASTIC),
              ('(n,inelastic)',  'inel', None),
              ('(n,Xt)',         'mt',   MT_XT)],

    'Be9':   [('(n,elastic)',    'mt',  MT_ELASTIC),
              ('(n,absorption)', 'abs', None)],

    'U235':  [('(n,elastic)',    'mt',   MT_ELASTIC),
              ('(n,inelastic)',  'inel', None),
              ('(n,2n)',         'mt',   MT_2N),
              ('(n,gamma)',      'mt',   MT_GAMMA),
              ('(n,fission)',    'mt',   MT_FISSION)],

    'U238':  [('(n,elastic)',    'mt',   MT_ELASTIC),
              ('(n,inelastic)',  'inel', None),
              ('(n,2n)',         'mt',   MT_2N),
              ('(n,gamma)',      'mt',   MT_GAMMA),
              ('(n,fission)',    'mt',   MT_FISSION)],

    'Th232': [('(n,elastic)',    'mt',  MT_ELASTIC),
              ('(n,inelastic)',  'inel', None),
              ('(n,2n)',         'mt',  MT_2N),
              ('(n,gamma)',      'mt',  MT_GAMMA),
              ('(n,fission)',    'mt',  MT_FISSION)],

    'Pb204': [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,inelastic)', 'inel', None),
              ('(n,2n)', 'mt', MT_2N), ('(n,absorption)', 'abs', None)],
    'Pb206': [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,inelastic)', 'inel', None),
              ('(n,2n)', 'mt', MT_2N), ('(n,absorption)', 'abs', None)],
    'Pb207': [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,inelastic)', 'inel', None),
              ('(n,2n)', 'mt', MT_2N), ('(n,absorption)', 'abs', None)],
    'Pb208': [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,inelastic)', 'inel', None),
              ('(n,2n)', 'mt', MT_2N), ('(n,absorption)', 'abs', None)],

    # --- structural / ceramic nuclides (>= 1 % mod power in HCPB or DCLL) ---
    'Fe54':  [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,inelastic)', 'inel', None),
              ('(n,absorption)', 'abs', None)],
    'Fe56':  [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,inelastic)', 'inel', None),
              ('(n,absorption)', 'abs', None)],
    'Fe57':  [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,inelastic)', 'inel', None),
              ('(n,absorption)', 'abs', None)],
    'Fe58':  [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,inelastic)', 'inel', None),
              ('(n,absorption)', 'abs', None)],

    'O16':   [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,absorption)', 'abs', None)],
    'O17':   [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,absorption)', 'abs', None)],
    'O18':   [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,absorption)', 'abs', None)],

    'Si28':  [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,absorption)', 'abs', None)],
    'Si29':  [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,absorption)', 'abs', None)],
    'Si30':  [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,absorption)', 'abs', None)],

    'C12':   [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,absorption)', 'abs', None)],
    'C13':   [('(n,elastic)', 'mt', MT_ELASTIC), ('(n,absorption)', 'abs', None)],
}


# =============================================================================
# ELEMENTAL GROUPING
# =============================================================================
# Isotopes summed into a single element line.  Anything not listed here is its
# own species (F-19, Li-6, Li-7, Be-9, Th-232).  The isotopes inside one group
# share an identical reaction list in REACTIONS, so the group inherits it.

GROUPS = {
    'U':  ['U235', 'U238'],
    'Pb': ['Pb204', 'Pb206', 'Pb207', 'Pb208'],
    'Fe': ['Fe54', 'Fe56', 'Fe57', 'Fe58'],
    'O':  ['O16', 'O17', 'O18'],
    'Si': ['Si28', 'Si29', 'Si30'],
    'C':  ['C12', 'C13'],
}
NUCLIDE_TO_SPECIES = {nuc: el for el, members in GROUPS.items() for nuc in members}


def species_of(nuclide):
    return NUCLIDE_TO_SPECIES.get(nuclide, nuclide)


def members_of(species):
    return GROUPS.get(species, [species])


def species_reactions(species):
    """(label, kind, mt) list for a species; grouped members share one set."""
    for m in members_of(species):
        if m in REACTIONS:
            return REACTIONS[m]
    return []


# =============================================================================
# PLOT STYLING
# =============================================================================
#  colour  -> species (element or kept-split nuclide)
#  style   -> reaction type    (so e.g. all elastic curves are solid)

SPECIES_COLOR = {
    'F19':  '#2ca02c',   # green
    'Li6':  '#1f77b4',   # blue
    'Li7':  '#17becf',   # cyan
    'Be9':  '#ff7f0e',   # orange
    'Th232':'#9467bd',   # purple
    'U':    '#d62728',   # red    (U-235 + U-238)
    'Pb':   '#7f7f7f',   # grey   (Pb-204/206/207/208)
    'Fe':   '#8c564b',   # brown  (Fe-54/56/57/58)
    'O':    '#e377c2',   # pink   (O-16/17/18)
    'Si':   '#bcbd22',   # olive  (Si-28/29/30)
    'C':    '#393b79',   # indigo (C-12/13)
}
SPECIES_PRETTY = {
    'F19':'F-19', 'Li6':'Li-6', 'Li7':'Li-7', 'Be9':'Be-9', 'Th232':'Th-232',
    'U':'U', 'Pb':'Pb', 'Fe':'Fe', 'O':'O', 'Si':'Si', 'C':'C',
}
# legend / draw order
SPECIES_ORDER = ['F19', 'Li6', 'Li7', 'Be9', 'C', 'O', 'Si', 'Fe',
                 'Th232', 'U', 'Pb']

# per-nuclide pretty names (used only for data-file warnings during caching)
PRETTY = {'F19':'F-19', 'Li6':'Li-6', 'Li7':'Li-7', 'Be9':'Be-9',
          'U235':'U-235', 'U238':'U-238', 'Th232':'Th-232',
          'Pb204':'Pb-204', 'Pb206':'Pb-206', 'Pb207':'Pb-207', 'Pb208':'Pb-208',
          'Fe54':'Fe-54', 'Fe56':'Fe-56', 'Fe57':'Fe-57', 'Fe58':'Fe-58',
          'O16':'O-16', 'O17':'O-17', 'O18':'O-18',
          'Si28':'Si-28', 'Si29':'Si-29', 'Si30':'Si-30',
          'C12':'C-12', 'C13':'C-13'}

# reaction label -> linestyle
REACTION_STYLE = {
    '(n,elastic)':    '-',
    '(n,inelastic)':  '--',
    '(n,2n)':         ':',
    '(n,gamma)':      '-.',
    '(n,fission)':    (0, (5, 1, 1, 1, 1, 1)),  # dash-dot-dot
    '(n,t)':          (0, (1, 1)),              # densely dotted  (tritium)
    '(n,Xt)':         (0, (1, 1)),              # densely dotted  (tritium)
    '(n,absorption)': (0, (6, 2)),              # long dash
}
# how the reaction legend is ordered + labelled ((n,t)/(n,Xt) share one entry)
RXN_LEGEND_ORDER = ['(n,elastic)', '(n,inelastic)', '(n,2n)', '(n,gamma)',
                    '(n,fission)', '(n,t)', '(n,absorption)']
RXN_LEGEND_LABEL = {'(n,t)': '(n,t) / (n,Xt)'}

# light, clean rcParams (no external font dependency -> stays standalone)
plt.rcParams.update({
    'font.size': 8, 'axes.labelsize': 9, 'axes.titlesize': 9,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7.5,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'axes.linewidth': 0.7, 'mathtext.default': 'regular', 'pdf.fonttype': 42,
})


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
    # pin the total atom density to the model's fixed value (see note above);
    # the per-nuclide atom fractions from the mix are preserved and rescaled.
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
    """Path to the neutron HDF5 file for `name`, or None. (Version-stable lookup.)"""
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


def _xs_at(nuc, T, mt, energies):
    """sigma(E) for a single MT on `energies`, or None if absent."""
    rx = nuc.reactions.get(mt)
    if rx is None or T not in rx.xs:
        return None
    return rx.xs[T](energies)


def micro_xs(nuc, T, kind, mt, energies):
    """Microscopic sigma(E) [barn] for one (kind, mt) on the nuclide energy grid."""
    if kind == 'mt':
        return _xs_at(nuc, T, mt, energies)

    if kind == 'inel':                                  # sum MT 51..91
        sig, found = np.zeros_like(energies), False
        for m in range(51, 92):
            s = _xs_at(nuc, T, m, energies)
            if s is not None:
                sig = sig + s
                found = True
        return sig if found else None

    if kind == 'abs':                                   # MT 27 -> 101 -> partials
        for m in (27, 101):
            s = _xs_at(nuc, T, m, energies)
            if s is not None:
                return s
        sig, found = np.zeros_like(energies), False
        for m in (102, 103, 104, 105, 106, 107, 108, 109,
                  111, 112, 113, 114, 115, 116, 117):
            s = _xs_at(nuc, T, m, energies)
            if s is not None:
                sig = sig + s
                found = True
        return sig if found else None

    raise ValueError(kind)


def ensure_cached(names, lib, cache):
    """Read microscopic XS once per nuclide into `cache[name] = (E, {label: sigma}, T)`."""
    for name in names:
        if name in cache:
            continue
        path = find_path(lib, name)
        if path is None:
            print(f"  ! no data file for {PRETTY.get(name, name)} in the library -- skipped")
            cache[name] = None
            continue
        nuc = openmc.data.IncidentNeutron.from_hdf5(path)
        T   = pick_temp(nuc, TARGET_TEMP_K)
        E   = nuc.energy[T]
        sig = {}
        for label, kind, mt in REACTIONS[name]:
            s = micro_xs(nuc, T, kind, mt, E)
            if s is None:
                print(f"  ! {PRETTY[name]} {label}: not in data file -- skipped")
                continue
            sig[label] = np.asarray(s)
        cache[name] = (E, sig, T)


# =============================================================================
# ELEMENTAL MACROSCOPIC CROSS SECTION
# =============================================================================

def elemental_sigma(species, label, dens, cache):
    """Macroscopic Sigma(E) [1/cm] for one reaction, summed over the element's
    isotopes:  Sigma = sum_i N_i * sigma_i(E).

    Each isotope's sigma is evaluated on its own (cached) energy grid; the sum is
    formed on the union of those grids (linear interpolation between an isotope's
    native points -- matching the lin-lin tabulation of the broadened XS).
    For a single-isotope species this is exactly N * sigma on the native grid.
    Returns (energy_grid, Sigma) or (None, None) if nothing contributes.
    """
    members = [m for m in members_of(species)
               if cache.get(m) and label in cache[m][1] and dens.get(m, 0.0) > 0.0]
    if not members:
        return None, None

    if len(members) == 1:                       # fast path: no interpolation
        E, sig, _ = cache[members[0]]
        return E, dens[members[0]] * sig[label]

    Egrid = np.unique(np.concatenate([cache[m][0] for m in members]))
    Sigma = np.zeros_like(Egrid)
    for m in members:
        Em, sig, _ = cache[m]
        Sigma = Sigma + dens[m] * np.interp(Egrid, Em, sig[label])
    return Egrid, Sigma


def present_species(blanket, densities):
    """Ordered species (elements / split nuclides) present in `blanket`."""
    present_nuc = {n for n in REACTIONS
                   if any(densities[(blanket, c)].get(n, 0.0) > 0.0
                          for c in FERTILE_CASES)}
    return [sp for sp in SPECIES_ORDER
            if any(m in present_nuc for m in members_of(sp))]


def fertile_species(blanket):
    """Element tracked in a panel annotation."""
    if blanket == 'FLiBe':
        return 'U'
    return 'U' if FERTILE_ISOTOPE == 'U238' else 'Th232'


# =============================================================================
# CSV EXPORT  (one file per blanket + loading, holding every plotted curve)
# =============================================================================

def write_spectra_csvs(densities, blk_species, cache,
                       outdir=CSV_OUTDIR, base=OUTNAME_BASE):
    """Write one CSV per (blanket, loading case) with every per-isotope Sigma(E)
    contribution behind the plotted curves, in WIDE form -- one row per energy,
    one column per isotope+reaction:

        energy_eV, <isotope><reaction>, ...

    e.g. columns  F19(n,elastic), Li6(n,t), U235(n,fission), U238(n,fission),
    Pb208(n,2n), ...   The grouped element curves in the figure (U = U235+U238,
    Pb = Pb204/206/207/208) are NOT summed here: each isotope's macroscopic
    contribution N_i*sigma_i(E) [1/cm] is written as its own column, so an
    element curve from the plot is just the row-wise sum of its isotope columns.
    No SIGMA_MIN floor is applied to the values.

    All columns in a panel are placed on a shared energy grid = the sorted union
    of the isotopes' native grids, and each is then linearly interpolated onto it
    (the same lin-lin interpolation elemental_sigma uses), so the values match
    the broadened cross sections point for point.  Where an isotope has no native
    data at a grid energy the cell is left blank rather than extrapolated.  The
    grid is restricted to CSV_E_RANGE (the plotted x-limits) unless that is None,
    then thinned to ~1/CSV_DOWNSAMPLE of its points to keep the file small.
    The blanket and loading are encoded in each file name rather than repeated as
    constant columns; add them back if you need to concatenate the four files
    into one combined table.

    Returns the list of file paths written.
    """
    os.makedirs(outdir, exist_ok=True)
    written = []

    def fmt_col(arr):
        """Vectorised '%g' formatting; NaNs become blank cells."""
        s = np.char.mod(CSV_FLOAT_FMT, arr)
        s[np.isnan(arr)] = ''
        return s

    for blanket in BLANKETS:
        for case in FERTILE_CASES:
            dens  = densities[(blanket, case)]
            fname = os.path.join(outdir,
                                 f"{base}_{blanket}_{case}kgm3.csv")

            # 1. gather the panel's columns (header -> (E, Sigma)) in figure
            #    order, but kept PER ISOTOPE: the U235/U238 and Pb204/206/207/208
            #    contributions are written separately instead of being summed into
            #    the "U" / "Pb" element curves shown in the plot.  elemental_sigma
            #    on a single nuclide returns that isotope's N*sigma on its own
            #    native grid (the no-interpolation fast path).
            curves = []
            for sp in blk_species[blanket]:
                for nuc in members_of(sp):
                    if nuc not in REACTIONS:            # e.g. natural-U U234 -> skip
                        continue
                    for label, kind, mt in REACTIONS[nuc]:
                        E, Sigma = elemental_sigma(nuc, label, dens, cache)
                        if E is None:                   # isotope absent / no data
                            continue
                        curves.append((f"{nuc}{label}", E, Sigma))

            if not curves:                              # emit header-only stub
                with open(fname, 'w', newline='') as fh:
                    csv.writer(fh).writerow(['energy_eV'])
                written.append(fname)
                print(f"Saved: {fname}  (no curves in panel)")
                continue

            # 2. shared energy grid = union of every curve's native grid, clipped
            #    to the plotted window.
            E_union = np.unique(np.concatenate([E for _, E, _ in curves]))
            if CSV_E_RANGE is not None:
                lo, hi  = CSV_E_RANGE
                E_union = E_union[(E_union >= lo) & (E_union <= hi)]

            # 2b. thin the grid to ~1/CSV_DOWNSAMPLE of its points so the file is
            #     ~1/CSV_DOWNSAMPLE the size.  Keep evenly spaced indices (the two
            #     endpoints included), which preserves the native grid's relative
            #     density -- the resonance region, where points cluster, stays the
            #     most finely sampled.
            if CSV_DOWNSAMPLE and CSV_DOWNSAMPLE > 1 and E_union.size > 2:
                n_keep = max(2, round(E_union.size / CSV_DOWNSAMPLE))
                if n_keep < E_union.size:
                    idx = np.unique(np.linspace(0, E_union.size - 1, n_keep)
                                    .round().astype(int))
                    E_union = E_union[idx]

            # 3. interpolate + format each curve onto that grid (blank outside the
            #    curve's own native range so nothing is silently extrapolated).
            names    = [name for name, _, _ in curves]
            col_strs = []
            for _, E, Sigma in curves:
                S = np.interp(E_union, E, Sigma)
                S[(E_union < E[0]) | (E_union > E[-1])] = np.nan
                col_strs.append(fmt_col(S))

            # 4. write the wide CSV.
            table = np.column_stack([fmt_col(E_union)] + col_strs)
            with open(fname, 'w', newline='') as fh:
                writer = csv.writer(fh)
                writer.writerow(['energy_eV'] + names)
                writer.writerows(table.tolist())
            written.append(fname)
            print(f"Saved: {fname}  ({len(E_union)} rows x {len(names)} curves)")
    return written


# =============================================================================
# FLUX-SPECTRUM HELPERS  (shared by the ranking and any future flux-averaged
# quantities -- mirrors the approach in xs_arithmetic.py)
# =============================================================================

def _load_flux_spectra():
    """Read the flux CSVs.  Returns {blanket: {loading: (E_mid, phi)}} or {}
    for blankets whose file is missing / None.  Only non-zero flux bins are
    kept (zero-flux bins cannot contribute to the weighted average)."""
    spectra = {}
    for bname, fname in FLUX_FILES.items():
        if fname is None:
            continue
        path = os.path.join(FLUX_DIR, fname)
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


def _match_loading(available, target, tol=1.0):
    """Closest loading in *available* within *tol* of *target*, or None."""
    best = min(available, key=lambda x: abs(x - target))
    return best if abs(best - target) <= tol else None


def _flux_weighted_sigma(E_xs, Sigma_xs, phi_E, phi_w, phi_sum):
    """<Sigma>_phi = sum(Sigma(E_i) * phi_i) / sum(phi_i).

    Interpolates the macroscopic XS onto the flux-bin midpoints.
    """
    Sigma_at_bins = np.interp(phi_E, E_xs, Sigma_xs)
    return float(np.dot(Sigma_at_bins, phi_w) / phi_sum)


# =============================================================================
# 3 x 2 GRID  (rows = blankets, cols = loadings)
# =============================================================================

def make_grid(lib, cache):
    n_rows, n_cols = len(BLANKETS), len(FERTILE_CASES)

    # --- 1. build every material and grab atom densities --------------------
    densities, vf_info = {}, {}
    for blanket in BLANKETS:
        for case in FERTILE_CASES:
            mat, vf = build_blanket(blanket, case, FERTILE_ISOTOPE)
            densities[(blanket, case)] = atom_densities(mat)
            vf_info[(blanket, case)]   = vf

    # --- 2. read microscopic cross sections once per needed nuclide ---------
    needed = [n for n in REACTIONS
              if any(densities[k].get(n, 0.0) > 0.0 for k in densities)]
    ensure_cached(needed, lib, cache)

    # species present per blanket + the union for the shared legend
    blk_species = {b: present_species(b, densities) for b in BLANKETS}
    legend_species = [sp for sp in SPECIES_ORDER
                      if any(sp in blk_species[b] for b in BLANKETS)]

    # --- 3. plot ------------------------------------------------------------
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(6.1 * n_cols, 5.4 * n_rows),
                             sharex=True, sharey=True, squeeze=False)

    styles_seen = set()
    for i, blanket in enumerate(BLANKETS):
        fert = fertile_species(blanket)
        for j, case in enumerate(FERTILE_CASES):
            ax   = axes[i][j]
            dens = densities[(blanket, case)]

            for sp in blk_species[blanket]:
                for label, kind, mt in species_reactions(sp):
                    E, Sigma = elemental_sigma(sp, label, dens, cache)
                    if E is None:
                        continue
                    ax.loglog(E, Sigma, color=SPECIES_COLOR[sp],
                              linestyle=REACTION_STYLE[label], lw=1.3,
                              alpha=0.9, rasterized=True)
                    styles_seen.add('(n,t)' if label in ('(n,t)', '(n,Xt)') else label)

            # 14 MeV source marker
            ax.axvline(E_MARK, color='0.55', lw=0.8, ls=(0, (6, 3)), zorder=0)
            ax.text(E_MARK, SIGMA_MAX, ' 14 MeV', color='0.45', fontsize=6.5,
                    ha='left', va='top', rotation=90)

            # panel title: blanket + loading + fuel annotation
            if case > 0:
                nfert = sum(dens.get(m, 0.0) for m in members_of(fert))
                if blanket == 'FLiBe':
                    sub = (f"UF$_4$ = {vf_info[(blanket, case)]*100:.3g} vol%,  "
                           f"N({SPECIES_PRETTY[fert]}) = {nfert:.2e} /b-cm")
                else:
                    sub = (f"BISO = {vf_info[(blanket, case)]*100:.3g} vol%,  "
                           f"N({SPECIES_PRETTY[fert]}) = {nfert:.2e} /b-cm")
            else:
                sub = "pure breeder"
            ax.set_title(f"{blanket}  --  {case:g} kg/m$^3$\n{sub}")

            ax.set_xlim(E_MIN, E_MAX)
            ax.set_ylim(SIGMA_MIN, SIGMA_MAX)
            ax.grid(True, which='major', alpha=0.30, lw=0.5)
            ax.grid(True, which='minor', alpha=0.12, lw=0.3)
            if i == n_rows - 1:
                ax.set_xlabel("Neutron energy [eV]")
            if j == 0:
                ax.set_ylabel(r"$\Sigma = N\,\sigma$  [cm$^{-1}$]")

    # --- two-part legend (colour=species, style=reaction) -------------------
    sp_handles = [Line2D([0], [0], color=SPECIES_COLOR[sp], lw=2.4)
                  for sp in legend_species]
    sp_labels  = [SPECIES_PRETTY[sp] for sp in legend_species]
    rxn_order  = [r for r in RXN_LEGEND_ORDER if r in styles_seen]
    rxn_handles = [Line2D([0], [0], color='0.25', lw=1.5,
                          linestyle=REACTION_STYLE[r]) for r in rxn_order]
    rxn_labels  = [RXN_LEGEND_LABEL.get(r, r) for r in rxn_order]

    leg1 = fig.legend(sp_handles, sp_labels, title="Element / nuclide (colour)",
                      loc='lower center', bbox_to_anchor=(0.27, 0.0),
                      ncol=min(6, len(legend_species)), frameon=False)
    fig.legend(rxn_handles, rxn_labels, title="Reaction (line style)",
               loc='lower center', bbox_to_anchor=(0.74, 0.0),
               ncol=4, frameon=False)
    fig.add_artist(leg1)

    fig.suptitle("Elemental macroscopic cross sections -- FLiBe, HCPB, DCLL breeding "
                 f"blankets ({TARGET_TEMP_K} K, ENDF/B-VIII.0)", fontsize=11)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.95, bottom=0.11,
                        wspace=0.06, hspace=0.32)

    os.makedirs(OUTDIR, exist_ok=True)
    stem = os.path.join(OUTDIR, f"{OUTNAME_BASE}_grid")
    fig.savefig(stem + ".png", dpi=200)
    fig.savefig(stem + ".pdf", dpi=200)
    print(f"Saved: {stem}.png\n       {stem}.pdf")
    plt.close(fig)

    # --- 3b. CSV dump of the plotted spectra (one file per blanket+loading) --
    if WRITE_CSV:
        write_spectra_csvs(densities, blk_species, cache)

    # --- 4. printed ranking (flux-weighted when spectra are available) --------
    flux = _load_flux_spectra()

    for blanket in BLANKETS:
        blk_flux = flux.get(blanket, {})
        for case in FERTILE_CASES:
            dens = densities[(blanket, case)]

            # resolve the flux spectrum for this (blanket, loading)
            phi_E = phi_w = phi_sum = None
            method_tag = f"at {RANK_ENERGY/1e6:.0f} MeV"
            if blk_flux:
                matched = _match_loading(blk_flux.keys(), case)
                if matched is not None:
                    phi_E, phi_w = blk_flux[matched]
                    phi_sum = phi_w.sum()
                    if phi_sum > 0:
                        method_tag = f"flux-weighted (matched {matched:g} kg/m\u00b3)"
                    else:
                        phi_E = phi_w = phi_sum = None

            rows = []
            for sp in blk_species[blanket]:
                for label, kind, mt in species_reactions(sp):
                    E, Sigma = elemental_sigma(sp, label, dens, cache)
                    if E is None:
                        continue
                    if phi_E is not None:
                        val = _flux_weighted_sigma(E, Sigma, phi_E, phi_w,
                                                   phi_sum)
                    else:
                        val = float(np.interp(RANK_ENERGY, E, Sigma))
                    rows.append((val, f"{SPECIES_PRETTY[sp]:<5} {label}"))
            rows.sort(reverse=True)
            print(f"\n=== {blanket}: ranking by \u03a3,  {method_tag} "
                  f"-- {case:g} kg/m\u00b3 ===")
            for val, lab in rows:
                print(f"   {lab:<22} {val:.3e} /cm")


# =============================================================================
# SELF-SHIELDING GRID  (3 x 2: host Sigma_el, U Sigma, flux ratio)
# =============================================================================

def _total_macro_xs(dens, cache):
    """Approximate total macroscopic XS by summing all cached reaction channels
    for every nuclide present in `dens`.  Returns (E, Sigma_t) or (None, None).
    """
    nuclides = [n for n in REACTIONS
                if dens.get(n, 0.0) > 0.0 and cache.get(n) is not None]
    if not nuclides:
        return None, None

    Egrid = np.unique(np.concatenate([cache[n][0] for n in nuclides]))
    Sigma_t = np.zeros_like(Egrid)

    for n in nuclides:
        E_n, sig_dict, _ = cache[n]
        N = dens[n]
        for s in sig_dict.values():
            Sigma_t += N * np.interp(Egrid, E_n, s)

    return Egrid, Sigma_t


def make_shielding_grid(lib, cache):
    """3x2 grid: host Sigma_el, U Sigma_el, U Sigma_gamma on primary axis;
    narrow-resonance flux ratio phi(rho)/phi_0 underlaid on secondary axis.

    Rows: FLiBe, HCPB, DCLL.  Cols: 0.1 and 1000 kg/m3.
    phi(rho)/phi_0 = Sigma_t(0) / Sigma_t(rho)  (NR approximation).
    """
    n_rows, n_cols = len(BLANKETS), len(FERTILE_CASES)
    FERTILE_NUC = {'U235', 'U238', 'Th232'}

    # --- 1. build materials at each loading + zero baseline --------------------
    densities = {}
    for blanket in BLANKETS:
        mat0, _ = build_blanket(blanket, 0.0, FERTILE_ISOTOPE)
        densities[(blanket, 0.0)] = atom_densities(mat0)
        for case in FERTILE_CASES:
            mat, _ = build_blanket(blanket, case, FERTILE_ISOTOPE)
            densities[(blanket, case)] = atom_densities(mat)

    # --- 2. ensure micro XS are cached ----------------------------------------
    needed = [n for n in REACTIONS
              if any(densities[k].get(n, 0.0) > 0.0 for k in densities)]
    ensure_cached(needed, lib, cache)

    # --- 3. plot ---------------------------------------------------------------
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(6.5 * n_cols, 4.8 * n_rows),
                              squeeze=False)

    for i, blanket in enumerate(BLANKETS):
        for j, case in enumerate(FERTILE_CASES):
            ax  = axes[i][j]
            dens  = densities[(blanket, case)]
            dens0 = densities[(blanket, 0.0)]

            # --- secondary axis (underlaid): flux ratio -----------------------
            ax2 = ax.twinx()
            # put primary on top with transparent background so ax2 shows through
            ax.set_zorder(ax2.get_zorder() + 1)
            ax.patch.set_visible(False)

            Et,  Sig_t  = _total_macro_xs(dens,  cache)
            Et0, Sig_t0 = _total_macro_xs(dens0, cache)

            if Et is not None and Et0 is not None:
                E_common = np.unique(np.concatenate([Et, Et0]))
                mask = (E_common >= E_MIN) & (E_common <= E_MAX)
                Ep = E_common[mask]

                St  = np.interp(Ep, Et,  Sig_t)
                St0 = np.interp(Ep, Et0, Sig_t0)

                with np.errstate(divide='ignore', invalid='ignore'):
                    phi_ratio = np.where(St > 0, St0 / St, 1.0)

                ax2.fill_between(Ep, phi_ratio, 1.0,
                                  where=(phi_ratio < 1.0),
                                  color='#FFD700', alpha=0.30,
                                  interpolate=True, rasterized=True)
                ax2.plot(Ep, phi_ratio,
                         color='#DAA520', lw=0.5, alpha=0.5,
                         rasterized=True)

            ax2.set_xscale('log')
            ax2.set_ylim(0, 1.15)
            if j == n_cols - 1:
                ax2.set_ylabel(r'$\phi(\rho)\,/\,\phi_0$',
                               color='#B8860B', fontsize=8)
                ax2.tick_params(axis='y', colors='#B8860B', labelsize=7)
            else:
                ax2.set_yticks([])

            # --- primary axis: cross sections ---------------------------------

            # Host Sigma_elastic (all non-fertile nuclides)
            host_nuc = [n for n in REACTIONS
                        if dens.get(n, 0.0) > 0.0
                        and cache.get(n) is not None
                        and n not in FERTILE_NUC]
            if host_nuc:
                Eh = np.unique(np.concatenate(
                    [cache[n][0] for n in host_nuc]))
                Sig_h = np.zeros_like(Eh)
                for n in host_nuc:
                    E_n, sig_d, _ = cache[n]
                    if '(n,elastic)' in sig_d:
                        Sig_h += dens[n] * np.interp(
                            Eh, E_n, sig_d['(n,elastic)'])
                mh = (Eh >= E_MIN) & (Eh <= E_MAX)
                ax.loglog(Eh[mh], Sig_h[mh],
                          color='#7f7f7f', ls='-', lw=1.6, alpha=0.85,
                          rasterized=True)

            # U Sigma_elastic
            E_uel, S_uel = elemental_sigma('U', '(n,elastic)', dens, cache)
            if E_uel is not None:
                mu = (E_uel >= E_MIN) & (E_uel <= E_MAX)
                ax.loglog(E_uel[mu], S_uel[mu],
                          color='#1f77b4', ls='-', lw=1.3, alpha=0.9,
                          rasterized=True)

            # U Sigma_gamma
            E_ung, S_ung = elemental_sigma('U', '(n,gamma)', dens, cache)
            if E_ung is not None:
                mg = (E_ung >= E_MIN) & (E_ung <= E_MAX)
                ax.loglog(E_ung[mg], S_ung[mg],
                          color='#d62728', ls='-.', lw=1.3, alpha=0.9,
                          rasterized=True)

            # 14 MeV marker
            ax.axvline(E_MARK, color='0.55', lw=0.8,
                       ls=(0, (6, 3)), zorder=1)
            ax.text(E_MARK, SIGMA_MAX, ' 14 MeV', color='0.45',
                    fontsize=6.5, ha='left', va='top', rotation=90)

            ax.set_title(f"{blanket}  --  {case:g} kg/m$^3$")
            ax.set_xlim(E_MIN, E_MAX)
            ax.set_ylim(SIGMA_MIN, SIGMA_MAX)
            ax.grid(True, which='major', alpha=0.30, lw=0.5)
            ax.grid(True, which='minor', alpha=0.12, lw=0.3)

            if i == n_rows - 1:
                ax.set_xlabel("Neutron energy [eV]")
            if j == 0:
                ax.set_ylabel(r"$\Sigma$  [cm$^{-1}$]")

    # --- shared legend --------------------------------------------------------
    handles = [
        Line2D([0], [0], color='#7f7f7f', lw=2.0),
        Line2D([0], [0], color='#1f77b4', lw=1.5),
        Line2D([0], [0], color='#d62728', lw=1.5, ls='-.'),
        Patch(facecolor='#FFD700', alpha=0.35, edgecolor='#DAA520', lw=0.8),
    ]
    labels = [
        r'Host $\Sigma_{\mathrm{el}}$',
        r'U $\Sigma_{\mathrm{el}}$',
        r'U $\Sigma_{\gamma}$',
        r'$\phi(\rho)/\phi_0$  (NR approx.)',
    ]
    fig.legend(handles, labels,
               loc='lower center', bbox_to_anchor=(0.5, 0.0),
               ncol=4, frameon=False, fontsize=8)

    fig.suptitle(
        r"Self-shielding: host $\Sigma_{\mathrm{el}}$, U cross sections,"
        r" and NR flux ratio $\phi(\rho)/\phi_0$"
        f"  ({TARGET_TEMP_K} K, ENDF/B-VIII.0)", fontsize=11)
    fig.subplots_adjust(left=0.07, right=0.92, top=0.95, bottom=0.08,
                        wspace=0.18, hspace=0.28)

    os.makedirs(OUTDIR, exist_ok=True)
    stem = os.path.join(OUTDIR, f"{OUTNAME_BASE}_shielding_grid")
    fig.savefig(stem + ".png", dpi=200)
    fig.savefig(stem + ".pdf", dpi=200)
    print(f"Saved: {stem}.png\n       {stem}.pdf")
    plt.close(fig)


# =============================================================================
# ELASTIC / CAPTURE RATIO  (3 x 2 grid: rows = blankets, cols = fertile)
# =============================================================================

def make_ratio_grid(lib, cache):
    """6-panel plot of  Sigma(n,elastic) / Sigma(n,gamma)  for the fertile
    element, with both loading cases overlaid on each panel.

    Layout:
        (0,0) FLiBe-U   Sigma_U(n,el)  / Sigma_U(n,gamma)   0.1 & 1000 kg/m3
        (0,1) FLiBe-Th  Sigma_Th(n,el) / Sigma_Th(n,gamma)  0.1 & 1000 kg/m3
        (1,0) HCPB-U    Sigma_U(n,el)  / Sigma_U(n,gamma)   0.1 & 1000 kg/m3
        (1,1) HCPB-Th   Sigma_Th(n,el) / Sigma_Th(n,gamma)  0.1 & 1000 kg/m3
        (2,0) DCLL-U    Sigma_U(n,el)  / Sigma_U(n,gamma)   0.1 & 1000 kg/m3
        (2,1) DCLL-Th   Sigma_Th(n,el) / Sigma_Th(n,gamma)  0.1 & 1000 kg/m3

    The elemental macroscopic XS for each channel is summed over isotopes
    (U-235 + U-238 for U; Th-232 alone for Th) before dividing.

    Output: macro_xs_ratio_grid.{png,pdf}
    """
    FERTILE_ISOTOPES = ['U238', 'Th232']
    # species key used by elemental_sigma (grouped element or single nuclide)
    FERT_SPECIES = {'U238': 'U', 'Th232': 'Th232'}
    FERT_COL_LABEL = {'U238': 'U', 'Th232': 'Th'}
    FERT_SIGMA_TEX = {
        'U238':  (r'$\Sigma_{\mathrm{U}}$',  r'$\Sigma_{\mathrm{U}}$(n,el) / '
                  r'$\Sigma_{\mathrm{U}}$(n,$\gamma$)'),
        'Th232': (r'$\Sigma_{\mathrm{Th}}$', r'$\Sigma_{\mathrm{Th}}$(n,el) / '
                  r'$\Sigma_{\mathrm{Th}}$(n,$\gamma$)'),
    }

    # line style / colour per loading case
    CASE_STYLE = {
        0.1:    dict(color='#1f77b4', ls='-',  lw=1.3),   # blue, solid
        1000.0: dict(color='#d62728', ls='--', lw=1.3),   # red,  dashed
    }

    n_rows, n_cols = len(BLANKETS), len(FERTILE_ISOTOPES)

    # --- 1. build every material (blanket x fertile x loading) ----------------
    densities, vf_info = {}, {}
    for blanket in BLANKETS:
        for fiso in FERTILE_ISOTOPES:
            for case in FERTILE_CASES:
                mat, vf = build_blanket(blanket, case, fiso)
                densities[(blanket, fiso, case)] = atom_densities(mat)
                vf_info[(blanket, fiso, case)]   = vf

    # --- 2. ensure micro XS are cached (idempotent if make_grid ran first) ----
    needed = [n for n in REACTIONS
              if any(densities[k].get(n, 0.0) > 0.0 for k in densities)]
    ensure_cached(needed, lib, cache)

    # --- 3. plot --------------------------------------------------------------
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(6.1 * n_cols, 5.4 * n_rows),
                             sharex=True, sharey=True, squeeze=False)

    for i, blanket in enumerate(BLANKETS):
        for j, fiso in enumerate(FERTILE_ISOTOPES):
            ax = axes[i][j]
            sp = FERT_SPECIES[fiso]

            for case in FERTILE_CASES:
                dens = densities[(blanket, fiso, case)]

                E_el, Sig_el = elemental_sigma(sp, '(n,elastic)', dens, cache)
                E_ng, Sig_ng = elemental_sigma(sp, '(n,gamma)',   dens, cache)

                if E_el is None or E_ng is None:
                    continue

                # common energy grid (union of the two native grids)
                E_common = np.unique(np.concatenate([E_el, E_ng]))
                Sig_el_c = np.interp(E_common, E_el, Sig_el)
                Sig_ng_c = np.interp(E_common, E_ng, Sig_ng)

                # quotient, suppressing division warnings for zero capture
                with np.errstate(divide='ignore', invalid='ignore'):
                    ratio = Sig_el_c / Sig_ng_c
                ratio[Sig_ng_c <= 0.0] = np.nan

                # restrict to the plotted energy window
                mask = (E_common >= E_MIN) & (E_common <= E_MAX)
                sty = CASE_STYLE[case]
                ax.loglog(E_common[mask], ratio[mask],
                          color=sty['color'], ls=sty['ls'], lw=sty['lw'],
                          alpha=0.9, label=f"{case:g} kg/m$^3$",
                          rasterized=True)

            # horizontal guide at ratio = 1
            ax.axhline(1.0, color='0.55', lw=0.7, ls='--', zorder=0)

            # 14 MeV source marker (data-x, axes-y coordinates)
            ax.axvline(E_MARK, color='0.55', lw=0.8, ls=(0, (6, 3)), zorder=0)
            ax.text(E_MARK, 0.98, ' 14 MeV', color='0.45', fontsize=6.5,
                    ha='left', va='top', rotation=90,
                    transform=ax.get_xaxis_transform())

            ax.set_title(f"{blanket} -- {FERT_COL_LABEL[fiso]}")

            ax.set_xlim(E_MIN, E_MAX)
            ax.grid(True, which='major', alpha=0.30, lw=0.5)
            ax.grid(True, which='minor', alpha=0.12, lw=0.3)
            ax.legend(loc='best', fontsize=7, framealpha=0.7)
            if i == n_rows - 1:
                ax.set_xlabel("Neutron energy [eV]")
            if j == 0:
                ax.set_ylabel(
                    r"$\Sigma$(n,el) / $\Sigma$(n,$\gamma$)")

    fig.suptitle(
        "Fertile-element elastic / capture ratio  "
        r"$\Sigma$(n,el) / $\Sigma$(n,$\gamma$)"
        f"  ({TARGET_TEMP_K} K, ENDF/B-VIII.0)", fontsize=11)
    fig.subplots_adjust(left=0.08, right=0.985, top=0.95, bottom=0.06,
                        wspace=0.08, hspace=0.32)

    os.makedirs(OUTDIR, exist_ok=True)
    stem = os.path.join(OUTDIR, f"{OUTNAME_BASE}_ratio_grid")
    fig.savefig(stem + ".png", dpi=200)
    fig.savefig(stem + ".pdf", dpi=200)
    print(f"Saved: {stem}.png\n       {stem}.pdf")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    lib   = open_library()
    cache = {}   # sigma(E) per nuclide is read from HDF5 only once
    make_grid(lib, cache)
    make_shielding_grid(lib, cache)
    make_ratio_grid(lib, cache)


if __name__ == "__main__":
    main()