from __future__ import annotations

import argparse
import os
import re
from collections import defaultdict
from typing import Optional, Tuple

import numpy as np
import openmc
import pandas as pd


def _li6_tag(enrichment_at: float) -> str:
    """Matches folder segment Li{04.1f} from prism_dcll / nat_prism_dcll."""
    return f"Li{enrichment_at:04.1f}"


def _matches_folder(name: str, breeder: str, isotope: str) -> bool:
    if isotope not in name:
        return False
    b = breeder.lower()
    if b == "hcpb":
        return name.startswith((f"{b}_wedgeA", f"{b}_wedgeC"))
    if b == "dcll":
        # Legacy: dcll_wedgeA_U238_...
        if name.startswith(("dcll_wedgeA", "dcll_wedgeC", "dcll_wedgeB")):
            return True
        # prism_dcll / nat_prism_dcll: wedgeA_dcll_Li07.5_U238_...
        if name.startswith(("wedgeA_dcll_", "wedgeB_dcll_", "wedgeC_dcll_")):
            return True
        return False
    return False


def _parse_case_loading(folder: str, breeder: str) -> Optional[Tuple[str, float]]:
    """Return (case letter, loading_kg_m3) or None."""
    # Loading: last token like 000.10kgm3 or 999.99kgm3
    m_load = re.search(r"(\d+\.\d+)kgm3", folder)
    if not m_load:
        return None
    loading = float(m_load.group(1))

    if breeder.lower() == "dcll" and folder.startswith("wedge") and "_dcll_" in folder:
        # wedgeA_dcll_Li07.5_U238_000.10kgm3
        m_case = re.match(r"wedge([ABC])_dcll_", folder)
        if m_case:
            return m_case.group(1), loading
    if folder.startswith("dcll_wedge"):
        # dcll_wedgeA_U238_000.10kgm3
        parts = folder.split("_")
        if len(parts) >= 2 and parts[1].startswith("wedge"):
            return parts[1][-1], loading
    if folder.startswith("hcpb_wedge"):
        parts = folder.split("_")
        if len(parts) >= 2:
            return parts[1][-1], loading
    return None


def _dcll_pair_key(folder: str) -> Optional[str]:
    """
    Suffix shared by wedge A and C for the same case (Li tag + isotope + loading).
    e.g. wedgeA_dcll_Li90.0_U238_000.10kgm3 -> Li90.0_U238_000.10kgm3
    """
    if re.match(r"wedge[ABC]_dcll_", folder):
        m = re.match(r"wedge[ABC]_dcll_(.+)", folder)
        return m.group(1) if m else None
    if folder.startswith("dcll_wedge"):
        m = re.match(r"dcll_wedge[ABC]_(.+)", folder)
        return m.group(1) if m else None
    return None


def _case_letter_from_folder(folder: str) -> Optional[str]:
    if folder.startswith("wedge") and "_dcll_" in folder:
        m = re.match(r"wedge([ABC])_dcll_", folder)
        return m.group(1) if m else None
    if folder.startswith("dcll_wedge"):
        parts = folder.split("_")
        if len(parts) >= 2 and parts[1].startswith("wedge"):
            return parts[1][-1]
    return None


def _latest_statepoint_path(folder_path: str) -> Optional[str]:
    sp_files = [f for f in os.listdir(folder_path) if f.startswith("statepoint.") and f.endswith(".h5")]
    if not sp_files:
        return None
    return os.path.join(folder_path, max(sp_files, key=lambda x: int(x.split(".")[1])))


def _summed_flux_spectrum_df(sp: openmc.StatePoint) -> pd.DataFrame:
    """Sum flux spectrum over all breeding blanket cells per energy bin (matches plot_flux_spectra_csv)."""
    flux_df = sp.get_tally(name="flux spectrum").get_pandas_dataframe()
    e_lo, e_hi = "energy low [eV]", "energy high [eV]"
    g = flux_df.groupby([e_lo, e_hi], as_index=False)["mean"].sum()
    g["energy_mid_eV"] = (g[e_lo] + g[e_hi]) / 2.0
    return g.rename(columns={"mean": "flux"})


def extract_dcll_flux_ac_spectra(
    isotope: str,
    base_dir: str = "./OpenMC",
    *,
    li6_enrichment_at: float | None = None,
    out_path: str | None = None,
) -> None:
    """
    Build ``flux_spectra_AC_kernel_blanket_dcll_<isotope>.csv`` for ``plot_flux_spectra_csv.py``.

    Pairs wedge A and wedge C runs with the same suffix (Li tag + isotope + loading), sums
    the ``flux spectrum`` tally over breeding cells for each case, and merges on energy bins.
    """
    if not os.path.isdir(base_dir):
        return

    folders = [f for f in os.listdir(base_dir) if _matches_folder(f, "dcll", isotope)]
    if li6_enrichment_at is not None:
        tag = _li6_tag(li6_enrichment_at)
        folders = [f for f in folders if tag in f]

    by_key: dict[str, dict[str, str]] = defaultdict(dict)
    for folder in folders:
        pk = _dcll_pair_key(folder)
        if pk is None:
            continue
        case = _case_letter_from_folder(folder)
        if case not in ("A", "C"):
            continue
        by_key[pk][case] = folder

    if out_path is None:
        if li6_enrichment_at is not None and abs(li6_enrichment_at - 7.5) < 0.01:
            out_path = f"flux_spectra_AC_kernel_blanket_dcll_nat_{isotope}.csv"
        else:
            out_path = f"flux_spectra_AC_kernel_blanket_dcll_{isotope}.csv"

    rows: list[pd.DataFrame] = []
    e_lo, e_hi = "energy low [eV]", "energy high [eV]"

    for pk in sorted(by_key.keys()):
        g = by_key[pk]
        if "A" not in g or "C" not in g:
            continue
        folder_a, folder_c = g["A"], g["C"]
        path_a = os.path.join(base_dir, folder_a)
        path_c = os.path.join(base_dir, folder_c)
        spa, spc = _latest_statepoint_path(path_a), _latest_statepoint_path(path_c)
        if not spa or not spc:
            print(f"Skip flux AC (no statepoint): {pk}")
            continue
        try:
            dfa = (
                _summed_flux_spectrum_df(openmc.StatePoint(spa))
                .drop(columns=["energy_mid_eV"], errors="ignore")
                .rename(columns={"flux": "flux_A"})
            )
            dfc = (
                _summed_flux_spectrum_df(openmc.StatePoint(spc))
                .drop(columns=["energy_mid_eV"], errors="ignore")
                .rename(columns={"flux": "flux_C"})
            )
        except Exception as e:
            print(f"Skip flux AC for {pk}: {e}")
            continue

        merged = pd.merge(dfa, dfc, on=[e_lo, e_hi], how="inner")
        merged["energy_mid_eV"] = (merged[e_lo] + merged[e_hi]) / 2.0
        if merged.empty:
            print(f"Skip flux AC (empty merge): {pk}")
            continue

        parsed = _parse_case_loading(folder_a, "dcll")
        loading = float(parsed[1]) if parsed else float("nan")
        merged["loading_kg_m3"] = loading
        merged["ratio_C_over_A"] = merged["flux_C"] / merged["flux_A"].replace(0, np.nan)
        merged["pair_key"] = pk
        merged["folder_A"] = folder_a
        merged["folder_C"] = folder_c
        rows.append(merged)
        print(f"Flux AC: {folder_a} + {folder_c}")

    if not rows:
        print(f"No A/C flux pairs collated for {isotope} (need wedge A and C with matching names).")
        return

    out = pd.concat(rows, ignore_index=True)
    cols = [
        e_lo,
        e_hi,
        "energy_mid_eV",
        "loading_kg_m3",
        "flux_A",
        "flux_C",
        "ratio_C_over_A",
        "pair_key",
        "folder_A",
        "folder_C",
    ]
    out = out[[c for c in cols if c in out.columns]]
    out.to_csv(out_path, index=False)
    print(f"Saved flux spectra A/C: {os.path.abspath(out_path)} ({len(out)} rows)")


def extract_rates(
    breeder: str,
    isotope: str,
    *,
    li6_enrichment_at: float | None = None,
    out_path: str | None = None,
    flux_ac: bool = True,
):
    """
    Extract TBR and fertile (n,gamma) from OpenMC statepoints under ./OpenMC.

    For ``breeder='dcll'``, also writes ``flux_spectra_AC_kernel_blanket_dcll_<isotope>.csv``
    (paired wedge A/C ``flux spectrum`` tallies) unless ``flux_ac=False``.

    Parameters
    ----------
    breeder
        ``hcpb`` or ``dcll``.
    isotope
        ``U238`` or ``Th232``.
    li6_enrichment_at
        If set (e.g. ``7.5`` or ``90.0``), only folders whose names include the
        segment ``Li{enrichment:04.1f}`` (e.g. ``Li07.5``, ``Li90.0``) are used.
        Use ``7.5`` for nat_prism_dcll outputs; ``90.0`` for prism_dcll.
    out_path
        CSV path. Default: ``extracted_rates_{breeder}_{isotope}.csv``, or
        ``extracted_rates_{breeder}_nat_{isotope}.csv`` when ``li6_enrichment_at=7.5``.
    flux_ac
        If True (default), collate A/C flux spectra for DCLL for ``plot_flux_spectra_csv.py``.
    """
    base_dir = "./OpenMC"
    results = []
    breeder = breeder.lower()

    if not os.path.exists(base_dir):
        print(f"Directory {base_dir} not found. Run from Jupyter/Wedge")
        return

    folders = [f for f in os.listdir(base_dir) if _matches_folder(f, breeder, isotope)]

    if li6_enrichment_at is not None:
        tag = _li6_tag(li6_enrichment_at)
        folders = [f for f in folders if tag in f]

    if out_path is None:
        if breeder == "dcll" and li6_enrichment_at is not None and abs(li6_enrichment_at - 7.5) < 0.01:
            out_path = f"extracted_rates_{breeder}_nat_{isotope}.csv"
        else:
            out_path = f"extracted_rates_{breeder}_{isotope}.csv"

    for folder in sorted(folders):
        path = os.path.join(base_dir, folder)
        if not os.path.isdir(path):
            continue

        parsed = _parse_case_loading(folder, breeder)
        if parsed is None:
            print(f"Skip (unparsed name): {folder}")
            continue
        case, loading = parsed

        sp_files = [f for f in os.listdir(path) if f.startswith("statepoint.") and f.endswith(".h5")]
        if not sp_files:
            print(f"No statepoint found in {folder}")
            continue

        sp_path = os.path.join(path, max(sp_files, key=lambda x: int(x.split(".")[1])))

        try:
            sp = openmc.StatePoint(sp_path)

            li_tally = sp.get_tally(name="Total (n,t) rxn rate")
            li_df = li_tally.get_pandas_dataframe()

            tbr_mean = float(li_df["mean"].sum())
            tbr_std = float(np.sqrt((li_df["std. dev."] ** 2).sum()))

            fertile_tally = sp.get_tally(name="Total fertile rxn rate")
            f_df = fertile_tally.get_pandas_dataframe()

            ng_mask = (f_df["nuclide"] == isotope) & (f_df["score"] == "(n,gamma)")

            ng_mean = float(f_df.loc[ng_mask, "mean"].sum())
            ng_std = float(np.sqrt((f_df.loc[ng_mask, "std. dev."] ** 2).sum()))

            row = {
                "folder": folder,
                "case": case,
                "loading_kg_m3": loading,
                "tbr": tbr_mean,
                "tbr_std": tbr_std,
                f"{isotope.lower()}_n_gamma": ng_mean,
                f"{isotope.lower()}_n_gamma_std": ng_std,
            }
            if li6_enrichment_at is not None:
                row["li6_enrichment_at"] = li6_enrichment_at
            results.append(row)

            print(f"Extracted {folder}")

        except Exception as e:
            print(f"Error processing {folder}: {e}")

    if breeder == "dcll" and flux_ac:
        extract_dcll_flux_ac_spectra(isotope, base_dir=base_dir, li6_enrichment_at=li6_enrichment_at)

    if results:
        df = pd.DataFrame(results)
        df.to_csv(out_path, index=False)
        print(f"\nSaved {len(results)} records to {os.path.abspath(out_path)}")
    else:
        print("No data extracted.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract TBR and fertile (n,gamma) from ./OpenMC statepoints.")
    parser.add_argument(
        "--nat-dcll",
        action="store_true",
        help="Natural-Li DCLL only (folder tag Li07.5); writes extracted_rates_dcll_nat_<isotope>.csv",
    )
    parser.add_argument(
        "--li6",
        type=float,
        default=None,
        metavar="AT_PERCENT",
        help="Filter folder names to this Li-6 at%% tag (e.g. 7.5 or 90). Overrides --nat-dcll if set.",
    )
    parser.add_argument(
        "--no-flux-ac",
        action="store_true",
        help="Do not write flux_spectra_AC_kernel_blanket_dcll_*.csv (paired A/C flux spectra).",
    )
    args = parser.parse_args()

    li_filter: Optional[float] = args.li6
    if args.nat_dcll and li_filter is None:
        li_filter = 7.5

    for b in ["hcpb", "dcll"]:
        for iso in ["U238", "Th232"]:
            if b == "dcll" and li_filter is not None:
                extract_rates(b, iso, li6_enrichment_at=li_filter, flux_ac=not args.no_flux_ac)
            else:
                extract_rates(b, iso, flux_ac=not args.no_flux_ac)
