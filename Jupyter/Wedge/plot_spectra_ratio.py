"""
plot_spectra_ratio_2x2.py

2x2 grid of C/A (lattice / homogeneous) ratios of reaction-rate spectra.

Layout (same as original):
  (0,0) HCPB — 0.1 kg/m³   (0,1) HCPB — 1000 kg/m³
  (1,0) DCLL — 0.1 kg/m³    (1,1) DCLL — 1000 kg/m³

Each panel overlays the C/A ratio for:
  1. Flux                        (optional, controlled by PLOT_FLUX)
  2. U-238 (n,gamma)
  3. Li TBR  [Li6 + Li7 (n,Xt)]
  4. U-235 + U-238 fission
  5. U-235 + U-238 (n,2n)

Cell selection:
  Summed over all cells in the tally.

Usage:
    python plot_spectra_ratio_2x2.py
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from prism_utilities import *


# ──────────────────────────────────────────────────────────────────────────────
# USER CONFIGURATION
# ──────────────────────────────────────────────────────────────────────────────

ISOTOPE = 'U238'
FERTILE_LOADINGS = [0.1, 999.99]
OPENMC_BASE = './OpenMC'
LINEWIDTH = 1.25
PLOT_FLUX = True           # Set to True to include flux ratio
N_COARSE_BINS = 100        # Number of log-uniform coarse energy bins, or None to skip rebinning

BLANKETS = {
    'DCLL': {'prefix': 'dcll', 'enrich': '90.0'},
    'HCPB': {'prefix': 'hcpb', 'enrich': '60.0'},
}

# Per-panel y-axis limits: (blanket_key, loading) -> (ymin, ymax)
# Set to None to let matplotlib auto-scale that panel.
LOG_Y = True              # Set to True for log-scale y-axis

YLIMS_LIN = {
    ('HCPB',   0.1):    (-0.05, 1.05),
    ('HCPB', 999.99):   (0.675, 1.425),
    ('DCLL',   0.1):    (0.50, 1.50),
    ('DCLL', 999.99):   (0.50, 1.50),
}

YLIMS_LOG = {
    ('HCPB',   0.1):    None,
    ('HCPB', 999.99):   None,
    ('DCLL',   0.1):    None,
    ('DCLL', 999.99):   None,
}

# Tally names
TALLY_FLUX    = 'flux spectrum'
TALLY_FERTILE = 'Fertile rxn rates spectrum'
TALLY_LI      = 'Li rxn rates spectrum'

# ──────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

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


def extract_spectrum(sp_path, tally_name, score, nuclides):
    """
    Extract energy-binned reaction rate summed over all cells and listed nuclides.

    Returns energy_lo, energy_hi, energy_mid, mean  (1-D arrays).
    """
    import openmc

    sp = openmc.StatePoint(sp_path)
    tally = sp.get_tally(name=tally_name)

    energy_filter = None
    for f in tally.filters:
        if isinstance(f, openmc.EnergyFilter):
            energy_filter = f
            break

    n_ebins = len(energy_filter.bins)

    total = np.zeros(n_ebins)
    for nuc in nuclides:
        sl = tally.get_slice(scores=[score], nuclides=[nuc])
        vals = sl.mean.flatten()
        n_cells = len(vals) // n_ebins
        total += vals.reshape(n_cells, n_ebins).sum(axis=0)

    edges = energy_filter.bins
    energy_lo  = edges[:, 0]
    energy_hi  = edges[:, 1]
    energy_mid = np.sqrt(energy_lo * energy_hi)

    sp.close()
    return energy_lo, energy_hi, energy_mid, total


def extract_flux(sp_path):
    """Extract flux spectrum summed over all cells."""
    import openmc

    sp = openmc.StatePoint(sp_path)
    tally = sp.get_tally(name=TALLY_FLUX)

    energy_filter = None
    for f in tally.filters:
        if isinstance(f, openmc.EnergyFilter):
            energy_filter = f
            break

    n_ebins = len(energy_filter.bins)

    sl = tally.get_slice(scores=['flux'])
    vals = sl.mean.flatten()
    n_cells = len(vals) // n_ebins
    flux = vals.reshape(n_cells, n_ebins).sum(axis=0)

    edges = energy_filter.bins
    energy_lo  = edges[:, 0]
    energy_hi  = edges[:, 1]
    energy_mid = np.sqrt(energy_lo * energy_hi)

    sp.close()
    return energy_lo, energy_hi, energy_mid, flux




def rebin_log_uniform(energy_lo, energy_hi, values, n_coarse):
    """
    Rebin fine-bin data into n_coarse log-uniform energy groups.

    Sums *values* (reaction rates) falling into each coarse bin.
    Fine bins are assigned to whichever coarse bin contains their midpoint.

    Returns coarse_mid, coarse_values  (1-D arrays of length n_coarse).
    """
    # Build coarse edges spanning the full fine-bin range
    e_min = energy_lo.min()
    e_max = energy_hi.max()
    coarse_edges = np.logspace(np.log10(e_min), np.log10(e_max), n_coarse + 1)

    # Midpoints of the fine bins (geometric mean, already log-centered)
    fine_mid = np.sqrt(energy_lo * energy_hi)

    # Assign each fine bin to a coarse bin via its midpoint
    # np.digitize returns 1-based index; clip to [0, n_coarse-1]
    idx = np.digitize(fine_mid, coarse_edges) - 1
    idx = np.clip(idx, 0, n_coarse - 1)

    coarse_values = np.zeros(n_coarse)
    np.add.at(coarse_values, idx, values)

    coarse_mid = np.sqrt(coarse_edges[:-1] * coarse_edges[1:])
    return coarse_mid, coarse_values


# ──────────────────────────────────────────────────────────────────────────────
# SPECTRUM DEFINITIONS
# ──────────────────────────────────────────────────────────────────────────────
# Each entry: (label, tally_name, score, nuclides, color, linestyle)

SPECTRA = [
    ('Flux',                TALLY_FLUX,    'flux',      None,             '#000000',  '-'),
    (r'U-238 (n,$\gamma$)', TALLY_FERTILE, '(n,gamma)', ['U238'],         '#ff1f5b', '-'),
    (r'Li (n,Xt)',          TALLY_LI,      '(n,Xt)',    ['Li6', 'Li7'],   '#0c9edd', '-'),
    (r'U (n,f)',            TALLY_FERTILE, 'fission',   ['U235', 'U238'], '#f48628', '-'),
    (r'U (n,2n)',           TALLY_FERTILE, '(n,2n)',    ['U235', 'U238'], '#04cc6c', '-'),
]

# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():

    panels = [
        (0, 0, 'HCPB',   0.1, r'HCPB-UO$_2$ (0.1 kg/m³)'),
        (0, 1, 'HCPB', 999.99, r'HCPB-UO$_2$ (1000 kg/m³)'),
        (1, 0, 'DCLL',   0.1, r'DCLL-UO$_2$ (0.1 kg/m³)'),
        (1, 1, 'DCLL', 999.99, r'DCLL-UO$_2$ (1000 kg/m³)'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharex=True, sharey=False)

    for row, col, blanket_key, loading, title in panels:
        ax = axes[row, col]

        loading_tag = f'{loading:.1f}' if loading < 1 else f'{loading:.0f}'

        # ── Get statepoints for both cases ────────────────────────────────
        sp_path_A = find_statepoint(blanket_key, 'A', loading)
        sp_path_C = find_statepoint(blanket_key, 'C', loading)
        print(f"  {blanket_key} {loading_tag} kg/m³")
        print(f"    Case A → {sp_path_A}")
        print(f"    Case C → {sp_path_C}")

        # ── Compute and plot C/A for each reaction ───────────────────────
        for label, tally_name, score, nuclides, color, ls in SPECTRA:

            if score == 'flux' and not PLOT_FLUX:
                continue

            # Extract spectra for Case A and Case C
            if score == 'flux':
                elo, ehi, emid, mean_A = extract_flux(sp_path_A)
                _,   _,   _,    mean_C = extract_flux(sp_path_C)
            else:
                elo, ehi, emid, mean_A = extract_spectrum(
                    sp_path_A, tally_name, score, nuclides
                )
                _, _, _, mean_C = extract_spectrum(
                    sp_path_C, tally_name, score, nuclides
                )

            # Rebin into coarse log-uniform groups (or use fine bins if None)
            if N_COARSE_BINS is not None:
                coarse_mid, coarse_A = rebin_log_uniform(elo, ehi, mean_A, N_COARSE_BINS)
                _,          coarse_C = rebin_log_uniform(elo, ehi, mean_C, N_COARSE_BINS)
            else:
                coarse_mid, coarse_A, coarse_C = emid, mean_A, mean_C

            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.where(coarse_A > 0, coarse_C / coarse_A, np.nan)

            ax.step(coarse_mid, ratio, where='mid', linewidth=LINEWIDTH,
                    color=color, linestyle=ls, label=label)

        # ── Reference line at ratio = 1 ──────────────────────────────────
        ax.axhline(1.0, color='black', linewidth=0.5, linestyle=':')

        # ── Axis formatting ───────────────────────────────────────────────
        ax.set_xscale('log')
        ax.set_xlim(0.5e1, 1.5e7)
        if LOG_Y:
            ax.set_yscale('log')
            ylim = YLIMS_LOG.get((blanket_key, loading))
        else:
            ylim = YLIMS_LIN.get((blanket_key, loading))
        if ylim is not None:
            ax.set_ylim(*ylim)

        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

        leg = ax.legend(title=title, title_fontsize=12, fontsize=9,
                        fancybox=False, edgecolor='black', frameon=True,
                        framealpha=0.75, loc='upper right')
        leg.get_frame().set_linewidth(0.5)

    # ── Shared axis labels ────────────────────────────────────────────────
    for ax in axes[1, :]:
        ax.set_xlabel('Incident neutron energy [eV]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel('C / A  (lattice / homogeneous)', fontsize=14)

    fig.tight_layout()

    # ── Save ──────────────────────────────────────────────────────────────
    os.makedirs('./Figures/PDF', exist_ok=True)
    os.makedirs('./Figures/PNG', exist_ok=True)
    plt.savefig(f'./Figures/PDF/fig_spectra_ratio_{ISOTOPE}.pdf', bbox_inches='tight')
    plt.savefig(f'./Figures/PNG/fig_spectra_ratio_{ISOTOPE}.png', bbox_inches='tight')
    print(f"\nSaved to ./Figures/PDF/ and ./Figures/PNG/")
    plt.show()


if __name__ == '__main__':
    main()