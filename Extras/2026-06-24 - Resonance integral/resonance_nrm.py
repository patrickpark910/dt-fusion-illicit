#!/usr/bin/env python3

"""NRM-only plots: reads cached CSVs produced by resonance_full.py."""

from pathlib import Path
import argparse, os, sys, csv
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "Figures/Data"

os.chdir(ROOT)
sys.path.insert(0, str(ROOT))
from Python.parameters import *
OUTS = HERE / "Data"
FIGS = HERE / "Figures_NRM"

LOADINGS = [0.0, 0.01, 0.10, 0.50, 1, 10, 25, 50, 75, 100, 150, 250, 500, 750, 1000, 2000, 3000, 4000]
R2_MIN = 0.97

# OpenMC comparison data
CSV_RXNS_FLIBE_U  = DATA / "FLiBe_900K_Li07.5_U238_summary.csv"
CSV_RXNS_HCPB_U   = DATA / "HCPB_900K_Li60.0_U238_summary.csv"
CSV_RXNS_DCLL_U   = DATA / "DCLL_900K_Li90.0_U238_summary.csv"
CSV_RXNS_FLIBE_TH = DATA / "FLiBe_900K_Li07.5_Th232_summary.csv"
CSV_RXNS_HCPB_TH  = DATA / "HCPB_900K_Li60.0_Th232_summary.csv"
CSV_RXNS_DCLL_TH  = DATA / "DCLL_900K_Li90.0_Th232_summary.csv"

# Output figure paths
OUT_COMPARE_LINLOG = FIGS / "integral_vs_openmc_linlog"
OUT_COMPARE_LOGLOG = FIGS / "integral_vs_openmc_loglog"
OUT_COMPARE_LINEAR = FIGS / "integral_vs_openmc_linlin"
OUT_PSI_IEFF       = FIGS / "integral_psi_vs_ieff"

# Cached NRM results (produced by resonance_full.py)
OUT_CSV_FLIBE_U_NRM  = OUTS / "integral_nrm_flibe_uf4.csv"
OUT_CSV_DCLL_U_NRM   = OUTS / "integral_nrm_dcll_uo2.csv"
OUT_CSV_HCPB_U_NRM   = OUTS / "integral_nrm_hcpb_uo2.csv"
OUT_CSV_FLIBE_TH_NRM = OUTS / "integral_nrm_flibe_thf4.csv"
OUT_CSV_DCLL_TH_NRM  = OUTS / "integral_nrm_dcll_tho2.csv"
OUT_CSV_HCPB_TH_NRM = OUTS / "integral_nrm_hcpb_tho2.csv"


def read_cached_rows(csv_path):
    if not csv_path.exists():
        return None
    print(f"  cached: {csv_path}")
    with csv_path.open() as f:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]
    for r in rows:
        if "kgU238_per_m3_breeder" in r and "kg_fertile_per_m3_breeder" not in r:
            r["kg_fertile_per_m3_breeder"] = r.pop("kgU238_per_m3_breeder")
        if "N_U238_atom_per_b_cm" in r and "N_fertile_atom_per_b_cm" not in r:
            r["N_fertile_atom_per_b_cm"] = r.pop("N_U238_atom_per_b_cm")
    return rows


def read_openmc_generic(csv_path, ng_col='U238(n,g)'):
    with csv_path.open() as f:
        rows = list(csv.DictReader(f))
    return {
        "x": np.array([float(r["fertile_kg/m3"]) for r in rows]),
        "y": np.array([float(r[ng_col]) for r in rows]),
        "sd": np.array([float(r[f"{ng_col}_sd"]) for r in rows]),
    }


def plot_compare_all(cases, scale='linear', suffix=''):
    """Combined NRM and OpenMC for all blankets.

    scale: 'linear', 'linlog', or 'loglog'
    cases: list of (label, nrm_rows, x_key, openmc_csv_path, ng_col, color, linestyle, marker, fillstyle)
    """
    fig, ax = plt.subplots(figsize=(3.5, 3.0))

    if scale in ('linlog', 'loglog'):
        ax.set_xscale("log")
        if suffix:
            ax.set_xlim(10**(-1 - 0.03*5), 10**(4 + 0.03*5))
        else:
            ax.set_xlim(10**(-1 - 0.03*4), 10**(3 + 0.03*4))
        if scale == 'loglog':
            ax.set_yscale("log")
        out = OUT_COMPARE_LOGLOG if scale == 'loglog' else OUT_COMPARE_LINLOG
    else:
        max_loading = max(r[cases[0][2]] for c in cases for r in c[1])
        ax.set_xlim(-0.03 * max_loading, max_loading + 0.03 * max_loading)
        out = OUT_COMPARE_LINEAR
    out = Path(str(out) + suffix)

    def mask_for(x, y):
        return (x > 0) & (y > 0) if scale != 'linear' else np.ones_like(x, dtype=bool)

    dummy_cases = []
    for label, nrm_rows, x_key, openmc_csv, ng_col, color, linestyle, mkr, fillstyle in cases:
        omc = read_openmc_generic(openmc_csv, ng_col=ng_col)
        omc_mask = mask_for(omc["x"], omc["y"])
        ms, mew = (3.5, 0.5) if fillstyle == 'none' else (2, None)
        ax.plot(omc["x"][omc_mask], omc["y"][omc_mask],
                marker=mkr, color=color, lw=0.75, markersize=ms,
                markeredgewidth=mew, fillstyle=fillstyle, linestyle='none',
                label='_nolegend_')

        x_m = np.array([r[x_key] for r in nrm_rows])
        y_m = np.array([r["R_per_source_neutron"] for r in nrm_rows])
        mask = mask_for(x_m, y_m)
        ax.plot(x_m[mask], y_m[mask], ls=linestyle, color=color, lw=0.75,
                markersize=1, label='_nolegend_')

        dummy_cases.append((label, color, linestyle, mkr, fillstyle))

    ax.set_xlabel(r"Fertile isotope density [kg$/$m³]")
    ax.set_ylabel(r"Fertile isotope (n,$\gamma$) rate [rxns$/$src-n]")

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    if scale != 'loglog':
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        if scale == 'linear':
            ax.yaxis.set_minor_locator(MultipleLocator(0.25))
        else:
            ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    if scale == 'linear' and not suffix:
        ax.xaxis.set_minor_locator(MultipleLocator(50))
    ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
    ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

    # Dummy plots for legend — Th (left col) then U (right col)
    th = [d for d in dummy_cases if d[4] == 'none']
    u  = [d for d in dummy_cases if d[4] == 'full']
    for label, color, linestyle, mkr, fillstyle in th + u:
        ms, mew = (3.5, 0.5) if fillstyle == 'none' else (2, None)
        ax.plot([9e8, 9e9], [9e8, 9e9], marker=mkr, color=color, markersize=ms,
                markeredgewidth=mew, fillstyle=fillstyle, linestyle=linestyle,
                linewidth=0.75, label=label)

    if scale == 'loglog':
        ax.set_ylim(10**(-5 - 0.03*5), 10**(0 + 0.03*5))
    elif scale == 'linear' and suffix:
        ax.set_ylim(-0.03 * 2.0, 2.0 + 0.03 * 2.0)
    elif suffix:
        ax.set_ylim(-0.03, 1 + 0.03)
    else:
        ax.set_ylim(-0.03 * 0.5, 0.5 + 0.03 * 0.5)

    fig.tight_layout()
    legend_loc = 'lower right' if scale == 'loglog' else 'upper left'
    leg = ax.legend(loc=legend_loc, fancybox=False, edgecolor='none',
                    frameon=True, framealpha=0.0, ncol=2, fontsize=6.5)
    leg.get_frame().set_linewidth(0.5)

    fig.savefig(f'{out}.png', bbox_inches='tight', pad_inches=0.01, dpi=300)
    fig.savefig(f'{out}.pdf', bbox_inches='tight', pad_inches=0.01, dpi=300)

    print(f"wrote {out}")


def fit_power_law(x, y):
    """Least-squares fit of y = a * x**c via a linear fit in log-log space."""
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    c, log_a = np.polyfit(np.log(x[mask]), np.log(y[mask]), 1)
    return np.exp(log_a), c


def fit_power_law_progressive(x, y, r2_min=R2_MIN):
    """Fit I_eff = a * psi_0**c starting from the lowest psi_0 (most U-238,
    most self-shielded/"black") and expanding one point at a time toward
    higher psi_0 (more dilute), refitting each time. The window only ever
    grows to k points *after* the k-point fit's R^2 has already cleared
    r2_min -- so the returned fit's R^2 is guaranteed >= r2_min (asserted
    below), never the first window that dropped below it. This is an
    empirical boundary of the black-resonance (~sqrt(psi_0)) regime before
    the dilute/flat (I_infinity) regime starts pulling a single power-law
    exponent down toward 0.

    Returns (a, c, r2, n_used, psi_boundary) where n_used is how many
    (lowest-psi_0) points the returned fit used and psi_boundary is the
    largest psi_0 among them.
    """
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x, y = x[mask], y[mask]
    order = np.argsort(x)
    x, y = x[order], y[order]

    def r_squared(xk, yk, a, c):
        log_y = np.log(yk)
        yhat = np.log(a) + c * np.log(xk)
        ss_res = np.sum((log_y - yhat) ** 2)
        ss_tot = np.sum((log_y - log_y.mean()) ** 2)
        return 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

    # need >=2 points for any line, and R^2 is trivially 1.0 (not meaningful)
    # until there's a 3rd point to actually test the fit against
    n_used = min(2, len(x))
    a, c = fit_power_law(x[:n_used], y[:n_used])

    for k in range(3, len(x) + 1):
        a_k, c_k = fit_power_law(x[:k], y[:k])
        if r_squared(x[:k], y[:k], a_k, c_k) < r2_min:
            break
        a, c, n_used = a_k, c_k, k

    r2 = r_squared(x[:n_used], y[:n_used], a, c)
    assert n_used <= 2 or r2 >= r2_min, (
        f"fit_power_law_progressive invariant broken: returned R2={r2:.4f} "
        f"is below r2_min={r2_min} at n_used={n_used} -- a failing window got "
        f"accepted somewhere above; this should be unreachable")
    return a, c, r2, n_used, x[n_used - 1]


def plot_psi_vs_ieff(cases, r2_min=R2_MIN, suffix=''):
    """psi vs I_eff for all blankets, log-log.

    cases: list of (label, nrm_rows, color, marker, fillstyle, fit_ls)
    """
    fig, ax = plt.subplots(figsize=(3.5, 3.0))

    dummy_fits = []
    for label, nrm_rows, color, marker, fillstyle, fit_ls in cases:
        psi = np.array([r["sigma_bg"] for r in nrm_rows])
        i_eff = np.array([r["I_eff_b"] for r in nrm_rows])
        mask = np.isfinite(psi) & (psi > 0) & (i_eff > 0)

        ms, mew = (3.5, 0.5) if fillstyle == 'none' else (2, None)
        ax.plot(psi[mask], i_eff[mask], marker, color=color, markersize=ms,
                markeredgewidth=mew, fillstyle=fillstyle, linestyle="none",
                label='_nolegend_')

        a, c, r2, n_used, psi_boundary = fit_power_law_progressive(
            psi[mask], i_eff[mask], r2_min)
        n_total = int(mask.sum())
        print(f"{label}: I_eff = {a:.4g} * psi^{c:.4g}   "
              f"(R2={r2:.4f} using {n_used}/{n_total} points up to "
              f"psi={psi_boundary:.3g} b)")

        x_fit = np.geomspace(psi[mask].min(), psi_boundary, 100)
        ax.plot(x_fit, a * x_fit**c, ls=fit_ls, color=color, lw=0.75,
                label='_nolegend_')

        dummy_fits.append((label, a, c, color, marker, fillstyle, fit_ls))

    # Dummy plots for legend — Th (left col) then U (right col)
    th = [d for d in dummy_fits if d[5] == 'none']
    u  = [d for d in dummy_fits if d[5] == 'full']
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\psi_o$ [b]")
    ax.set_ylabel(r"I-eff [b]")
    ax.set_xlim(10**(1 - 0.03*6), 10**(8 + 0.03*7))
    ax.set_ylim(10**(-1 - 0.03*2), 10**(1 + 0.03*2))

    for label, a, c, color, mkr, fillstyle, fit_ls in th + u:
        ms, mew = (3.5, 0.5) if fillstyle == 'none' else (2, None)
        ax.plot([9e8, 9e9], [9e8, 9e9], marker=mkr, color=color, markersize=ms,
                markeredgewidth=mew, fillstyle=fillstyle, linestyle=fit_ls,
                linewidth=0.75,
                label=rf"{label}: ${a:.3g}\,\psi_o^{{{c:.3g}}}$")

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_major_formatter(lambda x, _: f'{x:g}')
    ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
    ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

    fig.tight_layout()
    leg = ax.legend(loc='lower right', fancybox=False, edgecolor='none',
                    frameon=True, framealpha=0.0, ncol=2, fontsize=6.5)
    leg.get_frame().set_linewidth(0.5)

    out = Path(str(OUT_PSI_IEFF) + suffix)
    fig.savefig(f'{out}.png', bbox_inches='tight', pad_inches=0.01, dpi=300)
    fig.savefig(f'{out}.pdf', bbox_inches='tight', pad_inches=0.01, dpi=300)
    print(f"wrote {out}")


def filter_rows(rows, loadings):
    loading_set = set(loadings)
    out = [r for r in rows if r["kg_fertile_per_m3_breeder"] in loading_set]
    missing = loading_set - {r["kg_fertile_per_m3_breeder"] for r in out}
    if missing:
        raise ValueError(f"loadings {sorted(missing)} not found in cached CSV")
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--full', action='store_true',
                        help='run/plot all loadings up to 4499 kg/m3')
    args = parser.parse_args()

    if args.full:
        loadings = LOADINGS
        suffix = '_full'
    else:
        loadings = [l for l in LOADINGS if l <= 1000]
        suffix = ''

    os.makedirs(FIGS, exist_ok=True)

    flibe_u_nrm  = read_cached_rows(OUT_CSV_FLIBE_U_NRM)
    hcpb_u_nrm   = read_cached_rows(OUT_CSV_HCPB_U_NRM)
    dcll_u_nrm   = read_cached_rows(OUT_CSV_DCLL_U_NRM)
    flibe_th_nrm = read_cached_rows(OUT_CSV_FLIBE_TH_NRM)
    hcpb_th_nrm  = read_cached_rows(OUT_CSV_HCPB_TH_NRM)
    dcll_th_nrm  = read_cached_rows(OUT_CSV_DCLL_TH_NRM)

    flibe_u_nrm  = filter_rows(flibe_u_nrm, loadings)
    hcpb_u_nrm   = filter_rows(hcpb_u_nrm, loadings)
    dcll_u_nrm   = filter_rows(dcll_u_nrm, loadings)
    flibe_th_nrm = filter_rows(flibe_th_nrm, loadings)
    hcpb_th_nrm  = filter_rows(hcpb_th_nrm, loadings)
    dcll_th_nrm  = filter_rows(dcll_th_nrm, loadings)

    x_key = "kg_fertile_per_m3_breeder"
    compare_cases = [
        (r"FLiBe-UF$_4$",  flibe_u_nrm,  x_key, CSV_RXNS_FLIBE_U,  'U238(n,g)',   '#66b420', '-',  'o', 'full'),
        (r"FLiBe-ThF$_4$", flibe_th_nrm, x_key, CSV_RXNS_FLIBE_TH, 'Th232(n,g)',  '#66b420', '--', 'o', 'none'),
        (r"HCPB-UO$_2$",   hcpb_u_nrm,   x_key, CSV_RXNS_HCPB_U,   'U238(n,g)',   '#b41f24', '-',  's', 'full'),
        (r"HCPB-ThO$_2$",  hcpb_th_nrm,  x_key, CSV_RXNS_HCPB_TH,  'Th232(n,g)',  '#b41f24', '--', 's', 'none'),
        (r"DCLL-UO$_2$",   dcll_u_nrm,   x_key, CSV_RXNS_DCLL_U,   'U238(n,g)',   '#0047ba', '-',  '^', 'full'),
        (r"DCLL-ThO$_2$",  dcll_th_nrm,  x_key, CSV_RXNS_DCLL_TH,  'Th232(n,g)',  '#0047ba', '--', '^', 'none'),
    ]
    plot_compare_all(compare_cases, scale='linlog', suffix=suffix)
    plot_compare_all(compare_cases, scale='loglog', suffix=suffix)
    plot_compare_all(compare_cases, scale='linear', suffix=suffix)

    psi_cases = [
        (r"FLiBe-UF$_4$",  flibe_u_nrm,  '#66b420', 'o', 'full',  '-'),
        (r"FLiBe-ThF$_4$", flibe_th_nrm, '#66b420', 'o', 'none',  '--'),
        (r"HCPB-UO$_2$",   hcpb_u_nrm,   '#b41f24', 's', 'full',  '-'),
        (r"HCPB-ThO$_2$",  hcpb_th_nrm,  '#b41f24', 's', 'none',  '--'),
        (r"DCLL-UO$_2$",   dcll_u_nrm,   '#0047ba', '^', 'full',  '-'),
        (r"DCLL-ThO$_2$",  dcll_th_nrm,  '#0047ba', '^', 'none',  '--'),
    ]
    plot_psi_vs_ieff(psi_cases, r2_min=R2_MIN, suffix=suffix)


if __name__ == "__main__":
    main()
