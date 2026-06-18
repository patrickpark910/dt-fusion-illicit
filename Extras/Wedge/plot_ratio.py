"""
Lattice/homogeneous ratio (wedge C / wedge A) vs fertile loading for HCPB and DCLL.
Left:  TBR ratio.   Right: fissile production rate (FPR) ratio.
Run from Extras/Wedge:  python plot_ratio.py
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

KEEP_TBR = [0.1, 100, 250, 500, 999.99]
KEEP_FPR = [0.1, 250, 500, 999.99]

# Marker shape/size per (blanket, isotope), matching the main plot.py schema.
MARKERS = {('hcpb', 'U238'):  ('s', 6),
           ('hcpb', 'Th232'): ('1', 10),
           ('dcll', 'U238'):  ('^', 6),
           ('dcll', 'Th232'): ('x', 8)}


def load_ratio(blanket, isotope, col, keep):
    """Return (loading, wedge-C / wedge-A) for the given value column."""
    path = f'./Figures/extracted_rates_{blanket}_{isotope}.csv'
    df = pd.read_csv(path)
    df = df[df['loading_kg_m3'].apply(lambda x: any(np.isclose(x, k) for k in keep))]
    hom = df[df['case'] == 'A'].sort_values('loading_kg_m3')
    lat = df[df['case'] == 'C'].sort_values('loading_kg_m3')
    merged = pd.merge(hom[['loading_kg_m3', col]], lat[['loading_kg_m3', col]],
                      on='loading_kg_m3', suffixes=('_A', '_C'))
    x = merged['loading_kg_m3'].values
    ratio = merged[f'{col}_C'].values / merged[f'{col}_A'].values
    return x, ratio


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

    # (axis, value column, y-axis label, mode, keep)
    #   mode 'dots'  -> markers only (4 shapes)
    #   mode 'lines' -> interpolated curves
    panels = [
        (axes[0], 'tbr',        'Lattice / homog. TBR', 'dots',  KEEP_TBR),
        (axes[1], 'u238_gamma', 'Lattice / homog. FPR', 'lines', KEEP_FPR),
    ]
    # FPR uses the breeder-specific (n,gamma) column.
    gamma_col = {'U238': 'u238_gamma', 'Th232': 'th232_gamma'}

    for ax, value_col, ylabel, mode, keep in panels:
        for blanket, color, name in [('hcpb', '#b41f24', 'HCPB'),
                                     ('dcll', '#0047ba', 'DCLL')]:
            for isotope, ls, label_iso in [('U238',  '-',        r'UO$_2$'),
                                           ('Th232', LONG_DASH, r'ThO$_2$')]:
                col = gamma_col[isotope] if value_col == 'u238_gamma' else value_col
                x, y = load_ratio(blanket, isotope, col, keep)
                label = f'{name}-{label_iso}'
                if mode == 'dots':
                    marker, ms = MARKERS[(blanket, isotope)]
                    ax.plot(x, y, linestyle='none', marker=marker, markersize=ms,
                            color=color, markerfacecolor='none', label=label)
                elif len(x) >= 2:
                    ax.plot(x_fine, Akima1DInterpolator(x, y, method='makima')(x_fine),
                            linestyle=ls, lw=0.75, color=color, label=label)
        ax.set_ylabel(ylabel)

    # Shared formatting
    for ax in axes:
        ax.set_xlim(-25, 1025)
        ax.set_xlabel(r'Fertile density [kg/m³]')
        ax.xaxis.set_major_locator(MultipleLocator(100))
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        ax.tick_params(which='both', top=True, right=True)
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

    # Left subplot: boxed legend on a translucent white background.
    leg = axes[0].legend(loc='lower right', frameon=True, fontsize=8,
                         facecolor='white', framealpha=0.5, edgecolor='black')
    leg.get_frame().set_linewidth(0.25)
    # Right subplot: plain (no box).
    axes[1].legend(loc='best', frameon=False, fontsize=8)

    # Left subplot (TBR ratio): y from 0.99 to 1.01 with a 3% buffer each end.
    _lo, _hi = 0.99, 1.01
    _buf = 0.03 * (_hi - _lo)
    axes[0].set_ylim(_lo - _buf, _hi + _buf)
    axes[0].yaxis.set_major_locator(MultipleLocator(0.002))
    axes[0].yaxis.set_minor_locator(MultipleLocator(0.001))

    # Right subplot (FPR ratio): y from 0.60 to 1.10 with a 3% buffer each end.
    _lo, _hi = 0.60, 1.10
    _buf = 0.03 * (_hi - _lo)
    axes[1].set_ylim(_lo - _buf, _hi + _buf)
    axes[1].yaxis.set_major_locator(MultipleLocator(0.05))
    axes[1].yaxis.set_minor_locator(MultipleLocator(0.025))

    plt.tight_layout(pad=0.3)
    fig.subplots_adjust(wspace=0.225)
    plt.savefig('./Figures/PDF/fig_wedge_ratio.pdf', bbox_inches='tight', pad_inches=0.02)
    plt.show()
