"""
plot_spectra_2x2.py

2x2 grid of cumulative normalized energy spectra, one panel per blanket/case:
  (0,0) HCPB Model A   (0,1) HCPB Model C
  (1,0) DCLL Model A   (1,1) DCLL Model C

Each panel overlays five spectra (both 0.1 and 1000 kg/m³ plotted together):
  1. Flux
  2. U-238 (n,gamma)
  3. Li TBR  [Li6 + Li7  (n,Xt)]
  4. U-235 + U-238 fission
  5. U-235 + U-238 (n,2n)

Cell selection:
  HCPB  ->  cell 23  (outboard breeding channel)
  DCLL  ->  cell 34  (first outboard breeding channel)

Usage:
    python plot_spectra_2x2.py
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from prism_utilities import *

# ──────────────────────────────────────────────────────────────────────────────
# FONTS & RCPARAMS  (same style as plot_hist.py)
# ──────────────────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────────────────
# USER CONFIGURATION
# ──────────────────────────────────────────────────────────────────────────────

ISOTOPE = 'U238'
FERTILE_LOADINGS = [0.1, 999.99]
OPENMC_BASE = './OpenMC'
LINEWIDTH = 0.8

BLANKETS = {
    'DCLL': {'prefix': 'dcll', 'enrich': '90.0'},
    'HCPB': {'prefix': 'hcpb', 'enrich': '60.0'},
}

# Cell IDs to filter on for each blanket type
#   HCPB: cell 23 = outboard breeding channel
#   DCLL: cell 34 = first outboard breeding channel
TARGET_CELL_ID = {
    'HCPB': 23,
    'DCLL': 34,
}

# Tally names
TALLY_FLUX    = 'flux spectrum'
TALLY_FERTILE = 'Fertile rxn rates spectrum'
TALLY_LI      = 'Li rxn rates spectrum'

# Cross-section underlay paths
XS_LI6_PATH       = './Figures/XSPlot/Li6nt.txt'
XS_LI7_PATH       = './Figures/XSPlot/Li7nt.txt'
XS_U238_FPR_PATH  = './Figures/XSPlot/U238gamma.txt'
XS_U238_FIS_PATH  = './Figures/XSPlot/U238fis.txt'

# ──────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

def read_xs_txt(filepath):
    """Read two-column XS text file -> (energy, log10_xs)."""
    energy, xs = [], []
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if i < 2:
                continue
            parts = line.split()
            if len(parts) >= 2:
                e, x = float(parts[0]), float(parts[1])
                if x > 0:
                    energy.append(e)
                    xs.append(np.log10(x))
    return np.array(energy), np.array(xs)


def shift_to_01(arr, lo=0.1, hi=0.9):
    """Min-max scale into [lo, hi]."""
    mn, mx = arr.min(), arr.max()
    if mx == mn:
        return np.full_like(arr, 0.5 * (lo + hi))
    return (arr - mn) * (hi - lo) / (mx - mn) + lo


def find_statepoint(blanket_key, case, fertile_kgm3):
    """Locate the latest statepoint file."""
    info = BLANKETS[blanket_key]
    fertile_str = f"{fertile_kgm3:06.2f}"
    pattern = os.path.join(
        OPENMC_BASE,
        f"{info['prefix']}_Li{info['enrich']}_wedge{case}_{ISOTOPE}_{fertile_str}kgm3_*",
        "statepoint.*.h5",
    )
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No statepoint for pattern:\n  {pattern}")
    return matches[-1]


def extract_spectrum_cell(sp_path, tally_name, score, nuclides, target_cell_id):
    """
    Extract energy-binned reaction rate from a single cell.

    Returns energy_lo, energy_hi, energy_mid, mean  (1-D arrays).
    """
    import openmc

    sp = openmc.StatePoint(sp_path)
    tally = sp.get_tally(name=tally_name)

    # Identify filters
    energy_filter = None
    cell_filter = None
    for f in tally.filters:
        if isinstance(f, openmc.EnergyFilter):
            energy_filter = f
        elif isinstance(f, openmc.CellFilter):
            cell_filter = f

    n_ebins = len(energy_filter.bins)

    # Find the index of target_cell_id inside the CellFilter
    cell_ids = [int(c) if isinstance(c, (int, np.integer)) else c.id for c in cell_filter.bins]
    if target_cell_id not in cell_ids:
        sp.close()
        raise ValueError(
            f"Cell {target_cell_id} not in CellFilter bins {cell_ids} "
            f"for tally '{tally_name}'"
        )
    cell_idx = cell_ids.index(target_cell_id)

    total = np.zeros(n_ebins)
    for nuc in nuclides:
        sl = tally.get_slice(scores=[score], nuclides=[nuc])
        vals = sl.mean.flatten()
        n_cells = len(vals) // n_ebins
        vals_2d = vals.reshape(n_cells, n_ebins)
        total += vals_2d[cell_idx, :]

    edges = energy_filter.bins
    energy_lo  = edges[:, 0]
    energy_hi  = edges[:, 1]
    energy_mid = np.sqrt(energy_lo * energy_hi)

    sp.close()
    return energy_lo, energy_hi, energy_mid, total


def extract_flux_cell(sp_path, target_cell_id):
    """Extract flux spectrum from a single cell."""
    import openmc

    sp = openmc.StatePoint(sp_path)
    tally = sp.get_tally(name=TALLY_FLUX)

    energy_filter = None
    cell_filter = None
    for f in tally.filters:
        if isinstance(f, openmc.EnergyFilter):
            energy_filter = f
        elif isinstance(f, openmc.CellFilter):
            cell_filter = f

    n_ebins = len(energy_filter.bins)

    cell_ids = [int(c) if isinstance(c, (int, np.integer)) else c.id for c in cell_filter.bins]
    cell_idx = cell_ids.index(target_cell_id)

    vals = tally.mean.flatten()
    n_cells = len(vals) // n_ebins
    vals_2d = vals.reshape(n_cells, n_ebins)
    flux = vals_2d[cell_idx, :]

    edges = energy_filter.bins
    energy_lo  = edges[:, 0]
    energy_hi  = edges[:, 1]
    energy_mid = np.sqrt(energy_lo * energy_hi)

    sp.close()
    return energy_lo, energy_hi, energy_mid, flux


def cum_norm_hist(ax, energy_mid, mean, bins, **kwargs):
    """Plot a single cumulative normalized step histogram."""
    total = mean.sum()
    if total == 0:
        return
    norm_mean = mean / total
    ax.hist(energy_mid, bins=bins, weights=norm_mean,
            cumulative=True, histtype='step', linewidth=LINEWIDTH, **kwargs)


# ──────────────────────────────────────────────────────────────────────────────
# SPECTRUM DEFINITIONS
# ──────────────────────────────────────────────────────────────────────────────
# Each entry: (label, tally_name, score, nuclides, color, linestyle)

SPECTRA = [
    ('Flux',               TALLY_FLUX,    'flux',      None,                '#888888', '-'),
    (r'U238 (n,$\gamma$)', TALLY_FERTILE, '(n,gamma)', ['U238'],            '#ff1f5b', '-'),
    ('(n,Xt)',             TALLY_LI,      '(n,Xt)',    ['Li6', 'Li7'],      '#0c9edd', '-'),
    ('U fis.',             TALLY_FERTILE, 'fission',   ['U235', 'U238'],    '#f48628', '-'),
    (r'U (n,2n)',          TALLY_FERTILE, '(n,2n)',    ['U235', 'U238'],    '#04cc6c', '--'),
]

# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():

    # ── Cross-section underlays ───────────────────────────────────────────
    # li6_energy,  li6_logxs  = read_xs_txt(XS_LI6_PATH)
    # li7_energy,  li7_logxs  = read_xs_txt(XS_LI7_PATH)
    # li6_shifted = shift_to_01(li6_logxs)
    # li7_shifted = shift_to_01(li7_logxs)
    # 
    # u238g_energy, u238g_logxs = read_xs_txt(XS_U238_FPR_PATH)
    # u238g_shifted = shift_to_01(u238g_logxs)
    # 
    # u238f_energy, u238f_logxs = read_xs_txt(XS_U238_FIS_PATH)
    # u238f_shifted = shift_to_01(u238f_logxs)

    # ── Panel layout ──────────────────────────────────────────────────────
    # Rows = blanket type, Cols = fertile loading
    # Each panel overlays Case A (homog., solid) vs Case C (lattice, dashed)
    # (row, col, blanket_key, loading, title)
    panels = [
        (0, 0, 'HCPB',   0.1, r'HCPB — 0.1 kg/m³ (cell 23)'),
        (0, 1, 'HCPB', 999.99, r'HCPB — 1000 kg/m³ (cell 23)'),
        (1, 0, 'DCLL',   0.1, r'DCLL — 0.1 kg/m³ (cell 34)'),
        (1, 1, 'DCLL', 999.99, r'DCLL — 1000 kg/m³ (cell 34)'),
    ]

    CASES = [
        ('A', '-',  'homog.'),
        ('C', '--', 'lattice'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharex=True, sharey=True)

    for row, col, blanket_key, loading, title in panels:
        ax = axes[row, col]
        cell_id = TARGET_CELL_ID[blanket_key]

        # ── Gray XS underlays ────────────────────────────────────────────
        # ax.plot(li6_energy,   li6_shifted,   lw=0.8, color='gray', alpha=0.20, label=r'Li-6 (n,t) xs')
        # ax.plot(li7_energy,   li7_shifted,   lw=0.8, color='gray', alpha=0.20, ls='--', label=r'Li-7 (n,t) xs')
        # ax.plot(u238g_energy, u238g_shifted, lw=0.8, color='gray', alpha=0.20, ls=':', label=r'U-238 (n,$\gamma$) xs')
        # ax.plot(u238f_energy, u238f_shifted, lw=0.8, color='gray', alpha=0.20, ls='-.', label=r'U-238 (n,f) xs')

        # ── Loop over cases (A=homog, C=lattice) ─────────────────────────
        bins = None
        for case, case_ls, case_tag in CASES:
            sp_path = find_statepoint(blanket_key, case, loading)
            loading_tag = f'{loading:.1f}' if loading < 1 else f'{loading:.0f}'
            print(f"  {blanket_key} case {case} | {loading_tag} kg/m³ → {sp_path}")

            for label, tally_name, score, nuclides, color, _unused_ls in SPECTRA:

                # Flux is special (no nuclide filter)
                if score == 'flux':
                    elo, ehi, emid, mean = extract_flux_cell(sp_path, cell_id)
                else:
                    elo, ehi, emid, mean = extract_spectrum_cell(
                        sp_path, tally_name, score, nuclides, cell_id
                    )

                if bins is None:
                    bins = np.sort(np.unique(np.concatenate([elo, ehi])))

                hist_label = f'{label} — {case_tag}'
                cum_norm_hist(ax, emid, mean, bins,
                              color=color, linestyle=case_ls,
                              label=hist_label)

        # ── Axis formatting ───────────────────────────────────────────────
        ax.set_xscale('log')
        ax.set_xlim(0.5e1, 1.5e7)
        ax.set_ylim(-0.03, 1.03)

        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

        leg = ax.legend(title=title, title_fontsize=12, fontsize=7,
                        fancybox=False, edgecolor='black', frameon=True,
                        framealpha=0.75, ncol=2, loc='upper left')
        leg.get_frame().set_linewidth(0.5)

    # ── Shared axis labels ────────────────────────────────────────────────
    for ax in axes[1, :]:
        ax.set_xlabel('Incident neutron energy [eV]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel('Cumulative fraction of reactions', fontsize=14)

    fig.tight_layout()

    # ── Save ──────────────────────────────────────────────────────────────
    os.makedirs('./Figures/PDF', exist_ok=True)
    os.makedirs('./Figures/PNG', exist_ok=True)
    plt.savefig(f'./Figures/PDF/fig_spectra_hist_{ISOTOPE}.pdf', bbox_inches='tight')
    plt.savefig(f'./Figures/PNG/fig_spectra_hist_{ISOTOPE}.png', bbox_inches='tight')
    print(f"\nSaved to ./Figures/PDF/ and ./Figures/PNG/")
    plt.show()


if __name__ == '__main__':
    main()