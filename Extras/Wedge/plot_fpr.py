"""
Fissile production rate vs fertile loading for HCPB and DCLL — homogeneous (case A) vs lattice (case C).
Run from Extras/Wedge:  python plot_fpr.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import Akima1DInterpolator

fm.fontManager.addfont('./fonts/arial.ttf')
plt.rcParams['font.family'] = fm.FontProperties(fname='./fonts/arial.ttf').get_name()
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['pdf.fonttype'] = 42

LONG_DASH = (0, (8, 2))

# ── Physics constants ─────────────────────────────────────────────────────────
AVO              = 6.022141076e+23
AMU_PU239        = 239.0521634
AMU_U233         = 233.039635207
BLANKET_COVERAGE = 1.00
N_PER_MJ         = 1 / (17.6 * 1.602e-19)
NPS_FUS          = 1000 * N_PER_MJ

TOT_N_PER_YR = NPS_FUS * 3.156e+7
PU239_CONV   = (TOT_N_PER_YR / AVO * AMU_PU239 / 1e3) * BLANKET_COVERAGE
U233_CONV    = (TOT_N_PER_YR / AVO * AMU_U233  / 1e3) * BLANKET_COVERAGE

KEEP = [0.1, 100, 250, 500, 999.99]


def load_wedge(blanket, isotope):
    path = f'./Figures/extracted_rates_{blanket}_{isotope}.csv'
    df = pd.read_csv(path)
    df = df[df['loading_kg_m3'].apply(lambda x: any(np.isclose(x, k) for k in KEEP))]
    hom = df[df['case'] == 'A'].sort_values('loading_kg_m3')
    lat = df[df['case'] == 'C'].sort_values('loading_kg_m3')
    return hom, lat


if __name__ == '__main__':

    plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8,
                         'xtick.labelsize': 8, 'ytick.labelsize': 8,
                         'axes.linewidth': 0.25,
                         'xtick.direction': 'in', 'ytick.direction': 'in',
                         'xtick.major.width': 0.25, 'ytick.major.width': 0.25,
                         'xtick.minor.width': 0.25, 'ytick.minor.width': 0.25,
                         'grid.color': '#DBDBDB'})

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.0))
    x_fine = np.linspace(0.1, 1000, 1000)

    for blanket, ax, color, name in [('hcpb', axes[0], '#b41f24', 'HCPB'),
                                      ('dcll', axes[1], '#0047ba', 'DCLL')]:
        for isotope, ls_h, ls_l, label_iso, gamma_col, conv in [
                ('U238',  '-',        '-',        r'UO$_2$',  'u238_gamma',  PU239_CONV),
                ('Th232', LONG_DASH, LONG_DASH, r'ThO$_2$', 'th232_gamma', U233_CONV)]:
            hom, lat = load_wedge(blanket, isotope)
            xh, yh = hom['loading_kg_m3'].values, hom[gamma_col].values * conv
            xl, yl = lat['loading_kg_m3'].values, lat[gamma_col].values * conv

            if len(xh) >= 2:
                ax.plot(x_fine, Akima1DInterpolator(xh, yh, method='makima')(x_fine),
                        linestyle=ls_h, lw=0.75, color=color, label=f'{label_iso} homog.')
            if len(xl) >= 2:
                ax.plot(x_fine, Akima1DInterpolator(xl, yl, method='makima')(x_fine),
                        linestyle=ls_l, lw=0.75, color=color, alpha=0.5, label=f'{label_iso} lattice')

        ax.text(0.90, 0.42, name, transform=ax.transAxes, fontweight='bold',
                fontsize=10, ha='right', va='top')
        ax.set_ylabel(r'Initial fissile production rate [kg/yr]')

    # Shared formatting
    for ax in axes:
        ax.set_xlim(-25, 1025)
        ax.set_ylim(-27, 927)
        ax.set_xlabel(r'Fertile density [kg/m³]')
        ax.xaxis.set_major_locator(MultipleLocator(100))
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        ax.yaxis.set_major_locator(MultipleLocator(100))
        ax.yaxis.set_minor_locator(MultipleLocator(50))
        ax.tick_params(which='both', top=True, right=True)
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)
        ax.legend(loc='lower right', frameon=False, fontsize=8)

    plt.tight_layout(pad=0.3)
    fig.subplots_adjust(wspace=0.225)
    plt.savefig('./Figures/PDF/fig_wedge_fpr.pdf', bbox_inches='tight', pad_inches=0.02)
    plt.show()
