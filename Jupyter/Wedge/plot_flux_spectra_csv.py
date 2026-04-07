"""
Plot flux_spectra_AC_kernel_blanket_*.csv (no OpenMC — pandas + matplotlib only).

Collated CSVs should use the ``flux`` / ``flux spectrum`` tallies from ``prism_dcll.py``
breeding cells (23, 25, 34, 36, 38): sum OpenMC cell bins per energy bin for total blanket flux.

Usage:
  python plot_flux_spectra_csv.py
  python plot_flux_spectra_csv.py "C:/path/flux_spectra_AC_kernel_blanket_dcll_U238.csv"
  python plot_flux_spectra_csv.py "C:/path/data.csv" --shielding

With no path, the script looks for ``flux_spectra_AC_kernel_blanket_dcll_U238.csv`` in order:
current working directory, this script's folder, then ``~/Downloads/`` (Docker: use ``-w`` cwd or pass the file path).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DEFAULT_CSV_NAME = "flux_spectra_AC_kernel_blanket_dcll_U238.csv"


def _default_csv_candidates() -> list[Path]:
    """Resolve default CSV: cwd and script dir work in Docker; Downloads for local desktop."""
    here = Path(__file__).resolve().parent
    return [
        Path.cwd() / DEFAULT_CSV_NAME,
        here / DEFAULT_CSV_NAME,
        Path.home() / "Downloads" / DEFAULT_CSV_NAME,
    ]


def _resolve_csv_path(explicit: Path | None) -> tuple[Path, list[Path]]:
    """Return (path, tried) where path is the file to use or the first candidate for errors."""
    if explicit is not None:
        p = explicit.expanduser()
        return p, [p]
    tried: list[Path] = []
    for c in _default_csv_candidates():
        tried.append(c)
        if c.is_file():
            return c, tried
    return tried[0], tried


def _energy_col(df: pd.DataFrame) -> str:
    if "energy_mid_eV" in df.columns:
        return "energy_mid_eV"
    if "energy_mid [eV]" in df.columns:
        return "energy_mid [eV]"
    raise ValueError("CSV needs column energy_mid_eV or 'energy_mid [eV]'")


def load_flux_spectra_csv(csv_path: Path) -> tuple[pd.DataFrame, str]:
    """Read CSV and return cleaned frame + energy column name."""
    df = pd.read_csv(csv_path)
    emid = _energy_col(df)
    m = (df["flux_A"].fillna(0).abs() + df["flux_C"].fillna(0).abs()) > 0
    d = df.loc[m].copy()
    if d.empty:
        d = df.copy()
    return d, emid


def plot_self_shielding_diagnostic(
    d: pd.DataFrame,
    emid: str,
    out_path: Path,
    *,
    e_min_ev: float = 1.0,
    e_max_ev: float = 2e7,
    pct_ylim: tuple[float, float] | None = (-25.0, 25.0),
    deriv_ylim: tuple[float, float] | None = None,
) -> Path:
    """
    Explain C vs A bias in terms of *spectral* mismatch vs ~uniform systematic offset.

    - **Top:** percent difference ``100 * (C/A - 1)`` vs energy. A bias that is nearly
      independent of energy (flat band) is consistent with a **volume/path normalization**
      or geometry factor that is almost spectrally white; strong **variation** with
      energy—especially dips/peaks in epithermal U-238 resonances—signals **spatial
      self-shielding** in resolved TRISO (C) vs smeared cross sections in homogeneous (A).

    - **Bottom:** ``d(%Δ)/d(log10 E)``, i.e. sensitivity of the percent bias per decade of
      energy. Large |derivative| where cross sections vary rapidly (resonances) is the
      hallmark of **energy-dependent self-shielding**; near-zero derivative over a broad
      band supports a **roughly spectrally flat** ~few-% discrepancy.

    Requires columns: ``ratio_C_over_A`` (or flux_A/flux_C), ``loading_kg_m3``, energy column.
    """
    loadings = sorted(d["loading_kg_m3"].unique())

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

    ax_pct = axes[0]
    ax_d = axes[1]

    for L in loadings:
        sub = d[d["loading_kg_m3"] == L].sort_values(emid)
        e = sub[emid].to_numpy(dtype=float)
        if "ratio_C_over_A" in sub.columns:
            r = sub["ratio_C_over_A"].to_numpy(dtype=float)
        else:
            fa = sub["flux_A"].to_numpy(dtype=float)
            fc = sub["flux_C"].to_numpy(dtype=float)
            r = np.divide(fc, np.where(fa != 0, fa, np.nan))

        pct = 100.0 * (r - 1.0)
        log10e = np.log10(np.clip(e, 1e-300, None))
        # Derivative w.r.t. log10(E); NaNs handled by masking below
        dpct = np.gradient(pct, log10e, edge_order=2)

        lab = f"{L:g} kg/m³"
        ax_pct.semilogx(e, pct, "-", alpha=0.85, label=lab)
        ax_d.semilogx(e, dpct, "-", alpha=0.85, label=lab)

    ax_pct.axhline(0.0, color="k", ls="--", lw=0.8, alpha=0.45)
    ax_pct.set_ylabel(r"100 × (C/A − 1)   [%]")
    ax_pct.set_title(
        "Relative flux bias C vs A — flat offset ≈ bulk/systematic; "
        "structure vs E ≈ spectral self-shielding"
    )
    if pct_ylim is not None:
        ax_pct.set_ylim(pct_ylim)
    ax_pct.grid(True, which="both", ls=":", alpha=0.4)
    ax_pct.legend(fontsize=8, title="Loading")

    ax_d.set_ylabel("d(% bias) / d log10(E)   [% per decade]")
    ax_d.set_xlabel("Energy [eV]")
    ax_d.set_xlim(e_min_ev, e_max_ev)
    if deriv_ylim is not None:
        ax_d.set_ylim(deriv_ylim)
    ax_d.grid(True, which="both", ls=":", alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_standard_spectra(
    d: pd.DataFrame,
    emid: str,
    csv_path: Path,
    out_path: Path,
) -> Path:
    """Original two-panel: ratio C/A and log-log flux for max loading."""
    loadings = sorted(d["loading_kg_m3"].unique())

    e_min_ev, e_max_ev = 1.0, 2e7
    ratio_y_min, ratio_y_max = -0.05, 5.0
    flux_y_min, flux_y_max = 1e-5, 1e2

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

    ax0 = axes[0]
    for L in loadings:
        sub = d[d["loading_kg_m3"] == L]
        ax0.semilogx(sub[emid], sub["ratio_C_over_A"], "-", alpha=0.85, label=f"{L:g} kg/m³")
    ax0.axhline(1.0, color="k", ls="--", lw=0.8, alpha=0.5)
    ax0.set_ylabel("flux ratio C / A")
    ax0.set_title(csv_path.name)
    ax0.set_xlim(e_min_ev, e_max_ev)
    ax0.set_ylim(ratio_y_min, ratio_y_max)
    ax0.legend(fontsize=8, title="U-238 loading")
    ax0.grid(True, which="both", ls=":", alpha=0.4)

    L_one = max(loadings)
    s1 = d[d["loading_kg_m3"] == L_one]
    ax1 = axes[1]
    ax1.loglog(s1[emid], s1["flux_A"].clip(lower=1e-30), "-", label=f"A (homog.) L={L_one:g}")
    ax1.loglog(s1[emid], s1["flux_C"].clip(lower=1e-30), "--", label=f"C (TRISO) L={L_one:g}")
    ax1.set_xlim(e_min_ev, e_max_ev)
    ax1.set_ylim(flux_y_min, flux_y_max)
    ax1.set_xlabel("Energy [eV]")
    ax1.set_ylabel("Flux [n/cm²-s]")
    ax1.legend()
    ax1.grid(True, which="both", ls=":", alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot AC kernel blanket flux spectra CSV.")
    parser.add_argument(
        "csv_path",
        nargs="?",
        type=Path,
        default=None,
        help=f"Path to CSV (default: search {DEFAULT_CSV_NAME} in cwd, script dir, then ~/Downloads/)",
    )
    parser.add_argument(
        "--shielding",
        action="store_true",
        help="Also write *_self_shielding.png (%% bias vs E and d%%/d log10 E).",
    )
    args = parser.parse_args()

    csv_path, tried = _resolve_csv_path(args.csv_path)
    if not csv_path.is_file():
        print(f"File not found: {csv_path}")
        if len(tried) > 1:
            print("Tried:")
            for p in tried:
                print(f"  {p}")
        print("Pass the CSV path as the first argument, or copy the file into the working directory.")
        sys.exit(1)

    d, emid = load_flux_spectra_csv(csv_path)

    stem = csv_path.with_suffix("").name
    out_path = csv_path.parent / f"{stem}_plot.png"
    plot_standard_spectra(d, emid, csv_path, out_path)
    print(f"Saved figure: {out_path}")

    if args.shielding:
        shield_path = csv_path.parent / f"{stem}_self_shielding.png"
        plot_self_shielding_diagnostic(d, emid, shield_path)
        print(f"Saved figure: {shield_path}")


if __name__ == "__main__":
    main()
