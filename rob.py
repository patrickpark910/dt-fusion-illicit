"""
ROB - Results Overview & Blanket assessment
TBR drop at 10%, heat load per neutron, and fission heat contribution.
Assesses FLiBe, DCLL, HCPB in same manner as plot.py.
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import Akima1DInterpolator

from Python.utilities import *
from Python.parameters import *

# Fission energy release [eV] - ~200 MeV per fission
MEV_PER_FISSION = 200e6  # eV


def read_and_sort(path: str):
    """Read CSV and sort by fertile_kg/m3."""
    if not os.path.isfile(path):
        return None
    df = pd.read_csv(path)
    if 'fertile_kg/m3' in df.columns:
        df = df.sort_values(by='fertile_kg/m3').reset_index(drop=True)
    return df


def find_loading_at_tbr_drop(df, tbr_drop_frac=0.10):
    """
    Find fertile loading [kg/m³] where TBR drops by tbr_drop_frac (e.g. 10% -> 0.90 of baseline).
    Returns (loading_kgm3, tbr_at_drop) or (None, None) if not found.
    """
    if df is None or len(df) < 2:
        return None, None
    tbr_0 = df[df['fertile_kg/m3'] == 0]['tbr'].values
    if len(tbr_0) == 0:
        tbr_0 = df.iloc[0]['tbr']
    else:
        tbr_0 = tbr_0[0]
    tbr_thresh = (1 - tbr_drop_frac) * tbr_0
    df_below = df[df['tbr'] <= tbr_thresh]
    if len(df_below) == 0:
        return None, tbr_thresh
    # Interpolate to find exact loading
    x = df['fertile_kg/m3'].values
    y = df['tbr'].values
    if np.any(y <= tbr_thresh):
        # Find first crossing
        idx = np.where(y <= tbr_thresh)[0][0]
        if idx == 0:
            return 0, y[0]
        # Linear interp between idx-1 and idx
        xa, xb = x[idx-1], x[idx]
        ya, yb = y[idx-1], y[idx]
        loading = xa + (tbr_thresh - ya) * (xb - xa) / (yb - ya) if yb != ya else xb
        return float(loading), tbr_thresh
    return None, tbr_thresh


def run_all():
    os.makedirs('./Figures/pdf', exist_ok=True)
    os.makedirs('./Figures/png', exist_ok=True)
    os.makedirs('./Figures/data', exist_ok=True)

    # Data config: (path, label, color, marker) - FLiBe, DCLL, HCPB
    configs = [
        ('./Figures/Data/FLiBe_900K_Li07.5_U238_rxns.csv',  r'FLiBe-UF$_4$ (Li7.5)',  '#66b420', 'o'),
        ('./Figures/Data/FLiBe_900K_Li07.5_Th232_rxns.csv', r'FLiBe-ThF$_4$ (Li7.5)', '#66b420', '+'),
        ('./Figures/Data/DCLL_900K_Li90.0_U238_rxns.csv',  r'DCLL-UO$_2$',   '#0047ba', '^'),
        ('./Figures/Data/DCLL_900K_Li90.0_Th232_rxns.csv',  r'DCLL-ThO$_2$',  '#0047ba', 'x'),
        ('./Figures/Data/HCPB_900K_Li60.0_U238_rxns.csv',  r'HCPB-UO$_2$',   '#b41f24', 's'),
        ('./Figures/Data/HCPB_900K_Li60.0_Th232_rxns.csv', r'HCPB-ThO$_2$',  '#b41f24', '1'),
    ]

    # Apply HCPB TBR correction if available (from plot.py)
    try:
        df_u = pd.DataFrame(HCPB_CONV_U_TBR, columns=['fertile_kgm3', 'ratio'])
        df_th = pd.DataFrame(HCPB_CONV_TH_TBR, columns=['fertile_kgm3', 'ratio'])
        def get_ratio(x, df_fertile):
            return np.interp(x, df_fertile['fertile_kgm3'], df_fertile['ratio'])
    except NameError:
        df_u = df_th = None
        get_ratio = lambda x, df: 1.0

    results = []
    datasets = []

    for path, label, color, marker in configs:
        df = read_and_sort(path)
        if df is None:
            print(f"{C.YELLOW}Skip {label}: no data at {path}{C.END}")
            continue
        if 'HCPB' in label and df_u is not None and df_th is not None:
            if 'Th' in label:
                df = df.copy()
                df['tbr'] = df['tbr'] * get_ratio(df['fertile_kg/m3'], df_th)
            else:
                df = df.copy()
                df['tbr'] = df['tbr'] * get_ratio(df['fertile_kg/m3'], df_u)

        loading_10, tbr_thresh = find_loading_at_tbr_drop(df, 0.10)
        tbr_0 = df[df['fertile_kg/m3'] == 0]['tbr'].values
        tbr_0 = tbr_0[0] if len(tbr_0) else np.nan
        results.append({
            'label': label,
            'color': color,
            'marker': marker,
            'tbr_0': tbr_0,
            'loading_10pct_kgm3': loading_10,
            'tbr_at_10pct_drop': tbr_thresh,
        })
        datasets.append((df, label, color, marker))

    # Print summary
    print(f"\n{C.CYAN}=== Blanket coverage: {BLANKET_COVERAGE*100:.0f}% ==={C.END}\n")
    print("Loading [kg/m³] at which TBR drops 10% (TBR = 0.9 × TBR₀):")
    for r in results:
        ld = r['loading_10pct_kgm3']
        s = f"  {r['label']}: " + (f"{ld:.1f} kg/m³" if ld is not None else "> max loading")
        print(s)
    print()

    # --- Plot 1: TBR vs loading (with 10% drop line per config) ---
    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    ax.axhspan(0, 1.00, color='#e0ded8')
    # Data curves (solid)
    for df, label, color, marker in datasets:
        if df is None:
            continue
        tbr_scaled = BLANKET_COVERAGE * df['tbr']
        plt.plot(df['fertile_kg/m3'], tbr_scaled, marker=marker, linestyle='-', color=color, label=label, markersize=5)
    # Colored dashed horizontal line at 10% TBR drop for each config (0.9 * TBR_0)
    for r in results:
        if np.isnan(r['tbr_0']) or r['tbr_at_10pct_drop'] is None:
            continue
        tbr_drop = BLANKET_COVERAGE * r['tbr_at_10pct_drop']
        ax.axhline(tbr_drop, color=r['color'], linestyle='--', alpha=0.85, linewidth=1.2)
        ax.plot([], [], color=r['color'], linestyle='--', linewidth=1.5, alpha=0.85,
                label='10% drop: ' + r['label'])
    plt.xlabel(r'Fertile isotope density in breeder [kg/m³]')
    plt.ylabel('Tritium breeding ratio')
    plt.xlim(-5, 1025)
    plt.ylim(0.75, 1.5)
    plt.legend(loc='best', fontsize=8)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('./Figures/pdf/fig_rob_tbr_drop.pdf', bbox_inches='tight')
    plt.savefig('./Figures/png/fig_rob_tbr_drop.png', bbox_inches='tight')
    print("Saved: fig_rob_tbr_drop")
    plt.close()

    # --- Plot 2: Heat per incident neutron ---
    has_heat = False
    has_fission = False
    for path, *_ in configs:
        d = read_and_sort(path)
        if d is not None:
            if 'heat_eV_src' in d.columns:
                has_heat = True
            if 'U238(n,fis)' in d.columns or 'Th232(n,fis)' in d.columns:
                has_fission = True

    if has_heat or has_fission:
        # Plot 2a: Heat per neutron [eV/src] - initial (0 loading) vs at each loading
        if has_heat:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))
        else:
            fig, ax2 = plt.subplots(1, 1, figsize=(8, 6))
            ax1 = None

        for path, label, color, marker in configs:
            df = read_and_sort(path)
            if df is None or len(df) == 0:
                continue
            fis_col = 'U238(n,fis)' if 'U238' in path else 'Th232(n,fis)'
            heat_col = 'heat_eV_src'
            if heat_col in df.columns and df[heat_col].notna().any():
                heat = df[heat_col].values
                heat_0 = df[df['fertile_kg/m3'] == 0][heat_col].values
                heat_0 = heat_0[0] if len(heat_0) else heat[0]
                factor = heat / heat_0
                if ax1 is not None:
                    ax1.plot(df['fertile_kg/m3'], heat / 1e6, marker=marker, linestyle='-', color=color, label=label)
                ax2.plot(df['fertile_kg/m3'], factor, marker=marker, linestyle='-', color=color, label=label)
            elif fis_col in df.columns:
                # Estimate: fission heat = N_fis * 200 MeV; baseline at 0 loading
                fis = df[fis_col].values
                heat_fis_eV = fis * MEV_PER_FISSION
                # Assume baseline heat ~3e7 eV/src (typical neutron heating), or use first point
                heat_0_est = 3e7  # placeholder
                factor = 1 + heat_fis_eV / heat_0_est
                ax2.plot(df['fertile_kg/m3'], factor, marker=marker, linestyle='-', color=color, label=label + ' (fission est.)')
        if ax1 is not None:
            ax1.set_xlabel(r'Fertile [kg/m³]')
            ax1.set_ylabel('Heat per neutron [MeV/src]')
            ax1.set_xlim(-5, 1025)
            ax1.legend()
            ax1.grid(True, alpha=0.3)
        ax2.set_xlabel(r'Fertile isotope density in breeder [kg/m³]')
        ax2.set_ylabel('Heat factor (vs 0 loading)')
        ax2.set_xlim(-5, 1025)
        ax2.axhline(1.0, color='gray', linestyle='--', alpha=0.5)
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('./Figures/pdf/fig_rob_heat_factor.pdf', bbox_inches='tight')
        plt.savefig('./Figures/png/fig_rob_heat_factor.png', bbox_inches='tight')
        print("Saved: fig_rob_heat_factor")
        plt.close()
    else:
        print(f"{C.YELLOW}Heat/fission data not in rxns.csv. Re-run main.py (extract_tallies + collate) to populate heat_eV_src and fission.{C.END}")

    # Export results table
    df_results = pd.DataFrame(results)
    df_results.to_csv('./Figures/data/rob_tbr_drop_summary.csv', index=False)
    print("Saved: rob_tbr_drop_summary.csv")
    print("\nDone.")


if __name__ == '__main__':
    run_all()
