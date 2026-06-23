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

Constants and material builders live in common.py (imported below); this file
handles the cross-section computation and plotting.

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

from common import *


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
# CROSS-SECTION DATA HELPERS
# =============================================================================



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
    flux = load_flux_spectra()

    for blanket in BLANKETS:
        blk_flux = flux.get(blanket, {})
        for case in FERTILE_CASES:
            dens = densities[(blanket, case)]

            # resolve the flux spectrum for this (blanket, loading)
            phi_E = phi_w = phi_sum = None
            method_tag = f"at {RANK_ENERGY/1e6:.0f} MeV"
            if blk_flux:
                matched = match_loading(blk_flux.keys(), case)
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