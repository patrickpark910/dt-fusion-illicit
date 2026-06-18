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

Cell selection (summed over all breeding channels):
  HCPB  ->  cells 13, 23  (inboard + outboard)
  DCLL  ->  cells 23, 25, 34, 36, 38  (2 inboard + 3 outboard)

Usage:
    python plot_spectra_2x2.py
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

from prism_utilities import *

plt.rcParams.update({
    'font.size':        8,
    'axes.titlesize':   9,
    'axes.labelsize':   8,
    'xtick.labelsize':  7,
    'ytick.labelsize':  7,
    'legend.fontsize':  7,
    'axes.linewidth':   0.6,
    'xtick.direction':  'in',
    'ytick.direction':  'in',
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size':  3.5,
    'ytick.major.size':  3.5,
    'xtick.minor.size':  2.0,
    'ytick.minor.size':  2.0,
    'grid.color':       '#DBDBDB',
})


# ──────────────────────────────────────────────────────────────────────────────
# USER CONFIGURATION
# ──────────────────────────────────────────────────────────────────────────────

ISOTOPES = ['U238', 'Th232']
OPENMC_BASE = './OpenMC'
LINEWIDTH = 0.8

BLANKETS = {
    'DCLL': {'prefix': 'dcll', 'enrich': '90.0'},
    'HCPB': {'prefix': 'hcpb', 'enrich': '60.0'},
}

# Breeding channel cell IDs to sum over
#   HCPB: cell 13 (inboard), cell 23 (outboard)
#   DCLL: cells 23, 25 (inboard), cells 34, 36, 38 (outboard)
BREEDING_CELL_IDS = {
    'HCPB': [13, 23],
    'DCLL': [23, 25, 34, 36, 38],
}

# Tally names
TALLY_FLUX    = 'flux spectrum'
TALLY_FERTILE = 'Fertile rxn rates spectrum'
TALLY_LI      = 'Li rxn rates spectrum'

# ──────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

def find_statepoint(blanket_key, case, fertile_kgm3, isotope):
    """Locate the latest statepoint file."""
    info = BLANKETS[blanket_key]
    fertile_str = f"{fertile_kgm3:06.2f}"
    pattern = os.path.join(
        OPENMC_BASE,
        f"{info['prefix']}_Li{info['enrich']}_wedge{case}_{isotope}_{fertile_str}kgm3_*",
        "statepoint.*.h5",
    )
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No statepoint for pattern:\n  {pattern}")
    return matches[-1]


def extract_spectrum_cells(sp_path, tally_name, score, nuclides, cell_ids):
    """
    Extract energy-binned reaction rate summed over specified breeding cells.

    Returns energy_lo, energy_hi, energy_mid, mean  (1-D arrays).
    """
    import openmc

    sp = openmc.StatePoint(sp_path)
    tally = sp.get_tally(name=tally_name)

    energy_filter = None
    cell_filter = None
    for f in tally.filters:
        if isinstance(f, openmc.EnergyFilter):
            energy_filter = f
        elif isinstance(f, openmc.CellFilter):
            cell_filter = f

    n_ebins = len(energy_filter.bins)

    all_cell_ids = [int(c) if isinstance(c, (int, np.integer)) else c.id
                    for c in cell_filter.bins]

    total = np.zeros(n_ebins)
    for nuc in nuclides:
        sl = tally.get_slice(scores=[score], nuclides=[nuc])
        vals = sl.mean.flatten()
        n_cells = len(vals) // n_ebins
        vals_2d = vals.reshape(n_cells, n_ebins)
        for cid in cell_ids:
            if cid in all_cell_ids:
                total += vals_2d[all_cell_ids.index(cid), :]

    edges = energy_filter.bins
    energy_lo  = edges[:, 0]
    energy_hi  = edges[:, 1]
    energy_mid = np.sqrt(energy_lo * energy_hi)

    sp.close()
    return energy_lo, energy_hi, energy_mid, total


def extract_flux_cells(sp_path, cell_ids):
    """Extract flux spectrum summed over specified breeding cells."""
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

    all_cell_ids = [int(c) if isinstance(c, (int, np.integer)) else c.id
                    for c in cell_filter.bins]

    vals = tally.mean.flatten()
    n_cells = len(vals) // n_ebins
    vals_2d = vals.reshape(n_cells, n_ebins)

    flux = np.zeros(n_ebins)
    for cid in cell_ids:
        if cid in all_cell_ids:
            flux += vals_2d[all_cell_ids.index(cid), :]

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

SPECTRA = {
    'U238': [
        ('Flux',               TALLY_FLUX,    'flux',      None,                '#888888', '-'),
        (r'U238 (n,$\gamma$)', TALLY_FERTILE, '(n,gamma)', ['U238'],            '#ff1f5b', '-'),
        ('(n,Xt)',             TALLY_LI,      '(n,Xt)',    ['Li6', 'Li7'],      '#0c9edd', '-'),
        ('U fis.',             TALLY_FERTILE, 'fission',   ['U235', 'U238'],    '#f48628', '-'),
        (r'U (n,2n)',          TALLY_FERTILE, '(n,2n)',    ['U235', 'U238'],    '#04cc6c', '--'),
    ],
    'Th232': [
        ('Flux',                 TALLY_FLUX,    'flux',      None,              '#888888', '-'),
        (r'Th232 (n,$\gamma$)',  TALLY_FERTILE, '(n,gamma)', ['Th232'],         '#ff1f5b', '-'),
        ('(n,Xt)',               TALLY_LI,      '(n,Xt)',    ['Li6', 'Li7'],    '#0c9edd', '-'),
        ('Th fis.',              TALLY_FERTILE, 'fission',   ['Th232'],         '#f48628', '-'),
        (r'Th (n,2n)',           TALLY_FERTILE, '(n,2n)',    ['Th232'],         '#04cc6c', '--'),
    ],
}

# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():

    # ── Panel layout ──────────────────────────────────────────────────────
    # Rows = blanket type, Cols = fertile loading
    # Each panel overlays Case A (homog., solid) vs Case C (lattice, dashed)
    # (row, col, blanket_key, loading)
    panels = [
        (0, 0, 'HCPB',   0.1),
        (0, 1, 'HCPB', 999.99),
        (1, 0, 'DCLL',   0.1),
        (1, 1, 'DCLL', 999.99),
    ]

    CASES = [
        ('A', '-',  'homog.'),
        ('C', (0, (8, 2)), 'lattice'),
    ]

    os.makedirs('./Figures/PDF', exist_ok=True)
    os.makedirs('./Figures/PNG', exist_ok=True)

    for isotope in ISOTOPES:
        print(f"\n{'='*60}\n  Isotope: {isotope}\n{'='*60}")

        fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.5), sharex=True, sharey=True)

        for row, col, blanket_key, loading in panels:
            ax = axes[row, col]
            cell_ids = BREEDING_CELL_IDS[blanket_key]

            bins = None
            for case, case_ls, case_tag in CASES:
                loading_tag = f'{loading:.1f}' if loading < 1 else f'{loading:.0f}'
                try:
                    sp_path = find_statepoint(blanket_key, case, loading, isotope)
                except FileNotFoundError:
                    print(f"  {blanket_key} case {case} | {loading_tag} kg/m³ → SKIPPED (no data)")
                    continue
                print(f"  {blanket_key} case {case} | {loading_tag} kg/m³ → {sp_path}")

                for label, tally_name, score, nuclides, color, _unused_ls in SPECTRA[isotope]:

                    if score == 'flux':
                        elo, ehi, emid, mean = extract_flux_cells(sp_path, cell_ids)
                    else:
                        elo, ehi, emid, mean = extract_spectrum_cells(
                            sp_path, tally_name, score, nuclides, cell_ids
                        )

                    if bins is None:
                        bins = np.sort(np.unique(np.concatenate([elo, ehi])))

                    cum_norm_hist(ax, emid, mean, bins,
                                  color=color, linestyle=case_ls)

            # ── Axis formatting ───────────────────────────────────────────
            ax.set_xscale('log')
            ax.set_xlim(0.5e1, 1.5e7)
            ax.set_ylim(-0.03, 1.03)

            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.yaxis.set_minor_locator(MultipleLocator(0.05))
            ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
            ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

            # ── Custom legend: reaction colors + linestyle keys ───────────
            handles = []
            for label, _, _, _, color, _ in SPECTRA[isotope]:
                handles.append(Line2D([], [], color=color, lw=LINEWIDTH, label=label))
            handles.append(Line2D([], [], color='black', lw=LINEWIDTH,
                                  linestyle='-', label='Homog.'))
            handles.append(Line2D([], [], color='black', lw=LINEWIDTH,
                                  linestyle=(0, (8, 2)), label='Lattice'))

            loading_str = f'{loading:.1f}' if loading < 1 else f'{round(loading)}'
            oxide = r'UO$_2$' if isotope == 'U238' else r'ThO$_2$'
            panel_title = f'{blanket_key}-{oxide}\n({loading_str} kg/m³)'
            leg = ax.legend(handles=handles, title=panel_title, title_fontsize=7,
                            fontsize=6, fancybox=False, edgecolor='black',
                            frameon=True, framealpha=0.75, loc='upper left')
            leg.get_frame().set_linewidth(0.5)

        # ── Shared axis labels ────────────────────────────────────────────
        for ax in axes[1, :]:
            ax.set_xlabel('Incident neutron energy [eV]')
        for ax in axes[:, 0]:
            ax.set_ylabel('Cumulative fraction of reactions')

        fig.tight_layout()

        # ── Save ──────────────────────────────────────────────────────────
        plt.savefig(f'./Figures/PDF/fig_wedge_spectra_hist_{isotope}.pdf', bbox_inches='tight')
        plt.savefig(f'./Figures/PNG/fig_wedge_spectra_hist_{isotope}.png', bbox_inches='tight')
        print(f"\n  Saved fig_wedge_spectra_hist_{isotope} to ./Figures/PDF/ and ./Figures/PNG/")
        plt.show()


if __name__ == '__main__':
    main()