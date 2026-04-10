"""
plot_hist_unnorm.py

4x2 cumulative (unnormalized) energy histogram of (n,Xt) for TBR, (n,gamma) for FPR,
fission for FIS, and total (n,2n) for neutron multiplication.

Layout:
    Columns : 0.1 kg/m³ (left)  |  1000 kg/m³ (right)
    Rows    : TBR  |  FPR  |  fission  |  (n,2n)
    Lines   : HCPB homog., HCPB lattice, DCLL homog., DCLL lattice

Gray cross-section underlays are drawn on a secondary y-axis so they
don't interfere with the absolute reaction-rate scale.

Usage:
    python plot_hist_unnorm.py
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

ISOTOPE = 'U238'          # 'U238' or 'Th232'
FERTILE_LOADINGS = [0.1, 999.99]   # kg/m3 values to compare (left col, right col)
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
TALLY_N2N = {
    'DCLL': 'Pb rxn rates spectrum',   # nuclides: Pb204, Pb206, Pb207, Pb208
    'HCPB': 'Be rxn rates spectrum',   # nuclides: Be9
}
N2N_NUCLIDES = {
    'DCLL': ['Pb204', 'Pb206', 'Pb207', 'Pb208'],
    'HCPB': ['Be9'],
}

# Cross-section text files for the gray underlay
XS_LI6_PATH   = './Figures/XSPlot/Li6nt.txt'
XS_LI7_PATH   = './Figures/XSPlot/Li7nt.txt'
XS_BE9_PATH   = './Figures/XSPlot/Be9n2n.txt'
XS_PB_PATH    = './Figures/XSPlot/Pbn2n.txt'
XS_U238_FPR_PATH  = './Figures/XSPlot/U238gamma.txt'
XS_TH232_FPR_PATH = './Figures/XSPlot/Th232gamma.txt'
XS_U235_FIS_PATH  = './Figures/XSPlot/U235fis.txt'
XS_U238_FIS_PATH  = './Figures/XSPlot/U238fis.txt'
XS_TH232_FIS_PATH = './Figures/XSPlot/Th232fis.txt'


# ──────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

def read_xs_txt(filepath):
    """
    Read a two-column cross-section text file (Energy [eV], XS [b]).
    Skips the first two header lines.
    Returns (energy, log10_xs) arrays.
    """
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


def shift_to_01(arr, lo=0.1, hi=0.9, global_min=None, global_max=None):
    """Min-max scale *arr* into [lo, hi].
    If global_min / global_max are supplied they override the array's own
    extremes, which lets two arrays be scaled on a common basis."""
    mn = arr.min() if global_min is None else global_min
    mx = arr.max() if global_max is None else global_max
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


def cum_unnorm_hist(ax, energy_mid, mean, bins, color, label):
    """
    Plot a single **unnormalized** cumulative step histogram.
    """
    ax.hist(energy_mid, bins=bins, weights=mean,
            cumulative=True, histtype='step', linewidth=1.2,
            color=color, label=label)


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():

    # ── Load cross-section underlays ──────────────────────────────────────
    li6_energy,  li6_logxs  = read_xs_txt(XS_LI6_PATH)
    li7_energy,  li7_logxs  = read_xs_txt(XS_LI7_PATH)
    # Scale Li-6 and Li-7 together so their relative magnitudes are preserved
    li_mn = min(li6_logxs.min(), li7_logxs.min())
    li_mx = max(li6_logxs.max(), li7_logxs.max())
    li6_shifted = shift_to_01(li6_logxs, global_min=li_mn, global_max=li_mx)
    li7_shifted = shift_to_01(li7_logxs, global_min=li_mn, global_max=li_mx)

    if ISOTOPE == 'U238':
        fert_energy, fert_logxs = read_xs_txt(XS_U238_FPR_PATH)
        fert_label = r'U-238(n,$\gamma$)'
        fis_energy, fis_logxs = read_xs_txt(XS_U238_FIS_PATH)
        fis_label = r'U-238(n,fis)'
        fis235_energy, fis235_logxs = read_xs_txt(XS_U235_FIS_PATH)
        fis235_label = r'U-235(n,fis)'
    else:
        fert_energy, fert_logxs = read_xs_txt(XS_TH232_FPR_PATH)
        fert_label = r'Th-232(n,$\gamma$)'
        fis_energy, fis_logxs = read_xs_txt(XS_U238_FIS_PATH)
        fis_label = r'U-233(n,fis)'
        fis235_energy, fis235_logxs = read_xs_txt(XS_U235_FIS_PATH)
        fis235_label = r'U-235 (n,f)'
    fert_shifted = shift_to_01(fert_logxs)
    # Scale U-235 and U-238 fission together so relative magnitudes are preserved
    fis_mn = min(fis_logxs.min(), fis235_logxs.min())
    fis_mx = max(fis_logxs.max(), fis235_logxs.max())
    fis_shifted = shift_to_01(fis_logxs, global_min=fis_mn, global_max=fis_mx)
    fis235_shifted = shift_to_01(fis235_logxs, global_min=fis_mn, global_max=fis_mx)

    # Load Be-9 and Pb (n,2n) cross-section underlays
    be9_energy, be9_logxs = read_xs_txt(XS_BE9_PATH)
    pb_energy,  pb_logxs  = read_xs_txt(XS_PB_PATH)
    n2n_mn = min(be9_logxs.min(), pb_logxs.min())
    n2n_mx = max(be9_logxs.max(), pb_logxs.max())
    be9_shifted = shift_to_01(be9_logxs, global_min=n2n_mn, global_max=n2n_mx)
    pb_shifted  = shift_to_01(pb_logxs,  global_min=n2n_mn, global_max=n2n_mx)

    # ── Nuclide / score config ────────────────────────────────────────────
    if ISOTOPE == 'U238':
        fpr_nuclides = ['U238']
        fis_nuclides = ['U235', 'U238']
        XO2    = r'UO$_2$'
        X      = r'U'
        X_iso  = r'U-238'
    else:
        fpr_nuclides = ['Th232']
        fis_nuclides = ['Th232']
        XO2    = r'ThO$_2$'
        X      = r'Th'
        X_iso  = r'Th-232'

    # ── Line style map ────────────────────────────────────────────────────
    # Lines are now blanket × case (4 per panel)
    line_styles = {
        ('HCPB', 'A'): {'color': '#ff1f5b', 'label': rf'HCPB-{XO2} - homog.'},
        ('HCPB', 'C'): {'color': '#f48628', 'label': rf'HCPB-{XO2} - lattice'},
        ('DCLL', 'A'): {'color': '#04cc6c', 'label': rf'DCLL-{XO2} - homog.'},
        ('DCLL', 'C'): {'color': '#0c9edd', 'label': rf'DCLL-{XO2} - lattice'},
    }

    tbr_nuclides = ['Li6', 'Li7']
    tbr_score    = '(n,Xt)'
    fpr_score    = '(n,gamma)'
    fission_score = 'fission'
    n2n_score    = '(n,2n)'

    # ── Create 4×2 figure ─────────────────────────────────────────────────
    # Columns = fertile loading;  Rows = reaction type
    fig, axes = plt.subplots(4, 2, figsize=(16, 24), sharex=True, sharey=False)

    # column index → fertile loading
    col_loadings = {0: 0.1, 1: 999.99}
    col_labels   = {0: r'0.1 kg$/$m³', 1: r'1000 kg$/$m³'}

    # row config: (row, qty_key, tally_name, score, nuclides, title_suffix)
    # For N2N, tally_name and nuclides are resolved per-blanket below
    row_config = [
        (0, 'TBR', TALLY_TBR,     tbr_score,     tbr_nuclides, rf'Li(n,Xt)'),
        (1, 'FPR', TALLY_FPR,     fpr_score,     fpr_nuclides, rf'{X}(n,$\gamma$)'),
        (2, 'FIS', TALLY_FISSION, fission_score, fis_nuclides, rf'{X}(n,fis)'),
        (3, 'N2N', None,          n2n_score,     None,          r'Total (n,2n)'),
    ]

    for row, qty, tally_name, score, nuclides, title_suffix in row_config:
        for col, loading in col_loadings.items():
            ax = axes[row, col]
            title = rf'{title_suffix} — {col_labels[col]}'

            # ── Gray cross-section underlay on twin y-axis ────────────────
            ax2 = ax.twinx()
            ax2.set_ylim(0, 1)
            ax2.tick_params(right=False, labelright=False)

            if qty == 'TBR':
                ax2.plot(li6_energy, li6_shifted, linewidth=1.0, color='gray',
                         alpha=0.25, label=r'Li-6(n,t)')
                ax2.plot(li7_energy, li7_shifted, linewidth=1.0, color='gray',
                         alpha=0.25, linestyle='--', label=r'Li-7(n,Xt)')
            elif qty == 'FPR':
                ax2.plot(fert_energy, fert_shifted, linewidth=1.0, color='gray',
                         alpha=0.25, label=fert_label)
            elif qty == 'FIS':
                ax2.plot(fis235_energy, fis235_shifted, linewidth=1.0, color='gray',
                         alpha=0.25, linestyle='-', label=fis235_label)
                ax2.plot(fis_energy, fis_shifted, linewidth=1.0, color='gray',
                         alpha=0.25, linestyle='--', label=fis_label)
            elif qty == 'N2N':
                ax2.plot(pb_energy, pb_shifted, linewidth=1.0, color='gray',
                         alpha=0.25, linestyle='-', label=r'Pb(n,2n)')
                ax2.plot(be9_energy, be9_shifted, linewidth=1.0, color='gray',
                         alpha=0.25, linestyle='--', label=r'Be-9(n,2n)')


            ax.set_zorder(ax2.get_zorder() + 1)
            ax.patch.set_visible(False)

            # ── Histograms: loop over blanket × case ──────────────────────
            bins = None
            for blanket_key in ['HCPB', 'DCLL']:
                # Resolve per-blanket tally/nuclides for N2N row
                if qty == 'N2N':
                    tn = TALLY_N2N[blanket_key]
                    nucs = N2N_NUCLIDES[blanket_key]
                else:
                    tn = tally_name
                    nucs = nuclides

                for case in CASES:
                    sp_path = find_statepoint(blanket_key, case, loading)
                    print(f"  Loading {blanket_key} {qty} | case {case} | "
                          f"{loading} kg/m³ → {sp_path}")

                    energy_lo, energy_hi, energy_mid, mean = extract_spectrum(
                        sp_path, tn, score, nucs
                    )

                    if bins is None:
                        bins = np.sort(np.unique(
                            np.concatenate([energy_lo, energy_hi])))

                    style = line_styles[(blanket_key, case)]
                    cum_unnorm_hist(ax, energy_mid, mean, bins,
                                    color=style['color'], label=style['label'])

            # ── Axis formatting ───────────────────────────────────────────
            ax.set_xscale('log')
            ax.set_xlim(0.5e1, 5e7)

            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
            ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

            # ── Legend: merge handles from both axes ──────────────────────
            h1, l1 = ax.get_legend_handles_labels()
            h2, l2 = ax2.get_legend_handles_labels()
            leg_loc = 'lower right' if qty == 'FPR' else 'upper left'
            leg = ax.legend(h1 + h2, l1 + l2,
                            title=title, title_fontsize=13, fontsize=9,
                            fancybox=False, edgecolor='black', frameon=True,
                            framealpha=0.75, ncol=1, loc=leg_loc)
            leg.get_frame().set_linewidth(0.5)

    # ── Shared axis labels ────────────────────────────────────────────────
    for ax in axes[3, :]:
        ax.set_xlabel('Incident neutron energy [eV]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel(r'Cumulative reaction rate [rxn$/$src-n]', fontsize=14)

    fig.tight_layout()

    # ── Save ──────────────────────────────────────────────────────────────
    os.makedirs('./Figures/PDF', exist_ok=True)
    os.makedirs('./Figures/PNG', exist_ok=True)
    plt.savefig(f'./Figures/PDF/fig_hist_unnorm_{ISOTOPE}.pdf', bbox_inches='tight')
    plt.savefig(f'./Figures/PNG/fig_hist_unnorm_{ISOTOPE}.png', bbox_inches='tight')
    print("\nSaved to ./Figures/PDF/ and ./Figures/PNG/")
    plt.show()


if __name__ == '__main__':
    main()