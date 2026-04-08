"""
plot_wedge_histograms.py

2x2 cumulative normalized energy histogram of (n,Xt) for TBR and (n,gamma) for FPR,
comparing Case A vs Case C at 0.1 and 1000 kg/m3, for DCLL and HCPB wedge models.

Gray cross-section underlays:
  - TBR panels: Li-6(n,t) + Li-7(n,t)
  - FPR panels: U-238(n,gamma) or Th-232(n,gamma)

Usage:
    python plot_wedge_histograms.py
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt

from prism_utilities import *


# ──────────────────────────────────────────────────────────────────────────────
# USER CONFIGURATION
# ──────────────────────────────────────────────────────────────────────────────

ISOTOPE = 'U238'          # 'U238' or 'Th232'
FERTILE_LOADINGS = [0.1, 999.99]   # kg/m3 values to compare
CASES = ['A', 'C']

# Base directory where OpenMC run folders live
OPENMC_BASE = './OpenMC'

# Naming convention from prism_dcll.py:
#   {prefix}_Li{enrich}_wedge{case}_{isotope}_{fertile:06.2f}kgm3_{nps_str}/
BLANKETS = {
    'DCLL': {'prefix': 'dcll', 'enrich': '90.0'},
    'HCPB': {'prefix': 'hcpb', 'enrich': '60.0'},
}

# Tally names (must match tallies() in your model scripts)
TALLY_TBR = 'Li rxn rates spectrum'       # scores: (n,Xt); nuclides: Li6, Li7
TALLY_FPR = 'Fertile rxn rates spectrum'   # scores: (n,gamma); nuclides: U238 or Th232
TALLY_FISSION = 'Fertile rxn rates spectrum'  # scores: fission; nuclides: U238 or Th232

# Cross-section text files for the gray underlay
XS_LI6_PATH   = './Figures/XSPlot/Li6nt.txt'
XS_LI7_PATH   = './Figures/XSPlot/Li7nt.txt'
XS_BE9_PATH   = './Figures/XSPlot/Be9n2n.txt'
XS_PB_PATH    = './Figures/XSPlot/Pbn2n.txt'
XS_U238_FPR_PATH  = './Figures/XSPlot/U238gamma.txt'
XS_TH232_FPR_PATH = './Figures/XSPlot/Th232gamma.txt'
XS_U238_FIS_PATH  = './Figures/XSPlot/U238fis.txt'
XS_TH232_FIS_PATH = './Figures/XSPlot/Th232fis.txt'


# ──────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

def read_xs_txt(filepath):
    """
    Read a two-column cross-section text file (Energy [eV], XS [b]).
    Skips the first two header lines.
    Returns (energy, log10_xs) arrays.  log10 is taken so the resonance
    structure compresses into [0,1] nicely when shifted.
    """
    energy, xs = [], []
    with open(filepath, 'r') as f:
        for i, line in enumerate(f):
            if i < 2:          # skip header lines
                continue
            parts = line.split()
            if len(parts) >= 2:
                e, x = float(parts[0]), float(parts[1])
                if x > 0:
                    energy.append(e)
                    xs.append(np.log10(x))
    return np.array(energy), np.array(xs)


def shift_to_01(arr, lo=0.1, hi=0.9):
    """Min-max scale *arr* into [lo, hi]."""
    mn, mx = arr.min(), arr.max()
    if mx == mn:
        return np.full_like(arr, 0.5 * (lo + hi))
    return (arr - mn) * (hi - lo) / (mx - mn) + lo


def find_statepoint(blanket_key, case, fertile_kgm3):
    """
    Locate the statepoint file for a given blanket/case/loading combination.
    """
    info = BLANKETS[blanket_key]
    fertile_str = f"{fertile_kgm3:06.2f}"
    pattern = os.path.join(
        OPENMC_BASE,
        f"{info['prefix']}_Li{info['enrich']}_wedge{case}_{ISOTOPE}_{fertile_str}kgm3_*",
        "statepoint.*.h5"
    )
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No statepoint found for pattern:\n  {pattern}")
    return matches[-1]


def extract_spectrum(sp_path, tally_name, score, nuclides):
    """
    Open a statepoint and extract the energy-binned reaction rate for the
    given tally, score, and nuclide(s).  Sums over all cells and all listed
    nuclides to produce a 1-D array vs energy bin.

    Returns:
        energy_lo, energy_hi, energy_mid, mean   (all 1-D arrays)
    """
    import openmc

    sp = openmc.StatePoint(sp_path)
    tally = sp.get_tally(name=tally_name)

    # Identify the energy filter
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
        vals_2d = vals.reshape(n_cells, n_ebins)
        total += vals_2d.sum(axis=0)

    edges = energy_filter.bins
    energy_lo  = edges[:, 0]
    energy_hi  = edges[:, 1]
    energy_mid = np.sqrt(energy_lo * energy_hi)

    sp.close()
    return energy_lo, energy_hi, energy_mid, total


def cum_norm_hist(ax, energy_mid, mean, bins, color, label):
    """Plot a single cumulative normalized step histogram."""
    norm_mean = mean / mean.sum()
    ax.hist(energy_mid, bins=bins, weights=norm_mean,
            cumulative=True, histtype='step', linewidth=1.2,
            color=color, label=label)


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():

    # ── Load cross-section underlays ──────────────────────────────────────
    li6_energy,  li6_logxs  = read_xs_txt(XS_LI6_PATH)
    li7_energy,  li7_logxs  = read_xs_txt(XS_LI7_PATH)
    li6_shifted = shift_to_01(li6_logxs)
    li7_shifted = shift_to_01(li7_logxs)

    if ISOTOPE == 'U238':
        fert_energy, fert_logxs = read_xs_txt(XS_U238_FPR_PATH)
        fert_label = r'U-238 (n,$\gamma$)'
        fis_energy, fis_logxs = read_xs_txt(XS_U238_FIS_PATH)
        fis_label = r'U-238 (n,f)'
    else:
        fert_energy, fert_logxs = read_xs_txt(XS_TH232_FPR_PATH)
        fert_label = r'Th-232 (n,$\gamma$)'
        fis_energy, fis_logxs = read_xs_txt(XS_TH232_FIS_PATH)
        fis_label = r'Th-232 (n,f)'
    fert_shifted = shift_to_01(fert_logxs)
    fis_shifted = shift_to_01(fis_logxs)


    # ── Line style map ────────────────────────────────────────────────────
    line_styles = {
        ('A',   0.1  ): {'color': '#ff1f5b', 'label': r'   0.1 kg$/$m³ - homog.'},
        ('A', 999.99 ): {'color': '#f48628', 'label': r'1000 kg$/$m³ - homog.'},
        ('C',   0.1  ): {'color': '#04cc6c', 'label': r'   0.1 kg$/$m³ - lattice'},
        ('C', 999.99 ): {'color': '#0c9edd', 'label': r'1000 kg$/$m³ - lattice'},
    }

    # ── Nuclide / score config ────────────────────────────────────────────
    if ISOTOPE == 'U238':
        fpr_nuclides = ['U238']
        iso_label    = r'UO$_2$'
    else:
        fpr_nuclides = ['Th232']
        iso_label    = r'ThO$_2$'

    tbr_nuclides = ['Li6', 'Li7']
    tbr_score    = '(n,Xt)'
    fpr_score    = '(n,gamma)'
    fission_score = 'fission'

    # ── Create 2×2 figure ─────────────────────────────────────────────────
    fig, axes = plt.subplots(3, 2, figsize=(16, 18), sharex=True, sharey=True)

    panel_config = [
        # row, col, blanket, qty,   tally,     score,     nuclides,     title
        (0, 0, 'HCPB', 'TBR', TALLY_TBR, tbr_score, tbr_nuclides, rf'HCPB-{iso_label} TBR'),
        (0, 1, 'DCLL', 'TBR', TALLY_TBR, tbr_score, tbr_nuclides, rf'DCLL-{iso_label} TBR'),
        (1, 0, 'HCPB', 'FPR', TALLY_FPR, fpr_score, fpr_nuclides, rf'HCPB-{iso_label} FPR'),
        (1, 1, 'DCLL', 'FPR', TALLY_FPR, fpr_score, fpr_nuclides, rf'DCLL-{iso_label} FPR'),
        (2, 0, 'HCPB', 'FIS', TALLY_FISSION, fission_score, fpr_nuclides, rf'HCPB-{iso_label} fis.'),
        (2, 1, 'DCLL', 'FIS', TALLY_FISSION, fission_score, fpr_nuclides, rf'DCLL-{iso_label} fis.'),
    ]

    for row, col, blanket_key, qty, tally_name, score, nuclides, title in panel_config:
        ax = axes[row, col]

        # ── Gray cross-section underlay (plotted first so histograms sit on top)
        if qty == 'TBR':
            ax.plot(li6_energy, li6_shifted, linewidth=1.0, color='gray',
                    alpha=0.25, label=r'Li-6 (n,t)')
            ax.plot(li7_energy, li7_shifted, linewidth=1.0, color='gray',
                    alpha=0.25, linestyle='--', label=r'Li-7 (n,t)')
        elif qty == 'FPR':
            ax.plot(fert_energy, fert_shifted, linewidth=1.0, color='gray',
                    alpha=0.25, label=fert_label)
        elif qty == 'FIS':
            ax.plot(fis_energy, fis_shifted, linewidth=1.0, color='gray',
                    alpha=0.25, label=fis_label)

        # ── Histograms ────────────────────────────────────────────────────
        bins = None
        for loading in FERTILE_LOADINGS:
            for case in CASES:
                sp_path = find_statepoint(blanket_key, case, loading)
                print(f"  Loading {blanket_key} {qty} | case {case} | "
                      f"{loading} kg/m³ → {sp_path}")

                energy_lo, energy_hi, energy_mid, mean = extract_spectrum(
                    sp_path, tally_name, score, nuclides
                )

                if bins is None:
                    bins = np.sort(np.unique(
                        np.concatenate([energy_lo, energy_hi])))

                style = line_styles[(case, loading)]
                cum_norm_hist(ax, energy_mid, mean, bins,
                              color=style['color'], label=style['label'])

        # ── Axis formatting ───────────────────────────────────────────────
        ax.set_xscale('log')
        ax.set_xlim(0.5e1, 1.5e7)
        ax.set_ylim(-0.03, 1.03)

        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

        leg_loc = 'lower right' if qty == 'FPR' else 'upper left'
        leg = ax.legend(title=title, title_fontsize=13, fontsize=9,
                        fancybox=False, edgecolor='black', frameon=True,
                        framealpha=0.75, ncol=1, loc=leg_loc)
        leg.get_frame().set_linewidth(0.5)

    # ── Shared axis labels ────────────────────────────────────────────────
    for ax in axes[2, :]:
        ax.set_xlabel('Incident neutron energy [eV]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel('Cumulative fraction of reactions', fontsize=14)

    fig.tight_layout()

    # ── Save ──────────────────────────────────────────────────────────────
    os.makedirs('./Figures/PDF', exist_ok=True)
    os.makedirs('./Figures/PNG', exist_ok=True)
    plt.savefig(f'./Figures/PDF/fig_hist_{ISOTOPE}.pdf', bbox_inches='tight')
    plt.savefig(f'./Figures/PNG/fig_hist_{ISOTOPE}.png', bbox_inches='tight')
    print("\nSaved to ./Figures/PDF/ and ./Figures/PNG/")
    plt.show()


if __name__ == '__main__':
    main()