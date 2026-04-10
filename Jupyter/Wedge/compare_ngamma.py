"""
compare_ngamma_C_vs_A.py

For each blanket type (HCPB, DCLL) and fertile loading, extract the
U-238(n,gamma) energy spectrum for cases A (homogeneous) and C (lattice),
then report which energy bins have a higher reaction rate in C than in A.

Usage:
    python compare_ngamma_C_vs_A.py
"""

import os
import glob
import numpy as np

# ── Configuration (mirrors plot_hist_unnorm.py) ──────────────────────────────

ISOTOPE = 'U238'
FERTILE_LOADINGS = [0.1, 999.99]
OPENMC_BASE = './OpenMC'

BLANKETS = {
    'DCLL': {'prefix': 'dcll', 'enrich': '90.0'},
    'HCPB': {'prefix': 'hcpb', 'enrich': '60.0'},
}

TALLY_NAME = 'Fertile rxn rates spectrum'
TALLY_FLUX = 'flux spectrum'
SCORE      = '(n,gamma)'
NUCLIDES   = ['U238']


# ── Helpers ──────────────────────────────────────────────────────────────────

def find_statepoint(blanket_key, case, fertile_kgm3):
    info = BLANKETS[blanket_key]
    fertile_str = f"{fertile_kgm3:06.2f}"
    pattern = os.path.join(
        OPENMC_BASE,
        f"{info['prefix']}_Li{info['enrich']}_wedge{case}_{ISOTOPE}_{fertile_str}kgm3_*",
        "statepoint.*.h5",
    )
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No statepoint found for:\n  {pattern}")
    return matches[-1]


def extract_spectrum(sp_path, tally_name, score, nuclides):
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
    """Extract the cell-summed flux spectrum (no nuclide filter)."""
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
    total = vals.reshape(n_cells, n_ebins).sum(axis=0)

    edges = energy_filter.bins
    energy_mid = np.sqrt(edges[:, 0] * edges[:, 1])

    sp.close()
    return energy_mid, total


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    os.makedirs('./Figures', exist_ok=True)
    outpath = f'./Figures/compare_ngamma_C_vs_A_{ISOTOPE}.txt'

    with open(outpath, 'w') as fout:
        def out(s=''):
            print(s)
            fout.write(s + '\n')

        for blanket in ['HCPB', 'DCLL']:
            for loading in FERTILE_LOADINGS:
                out(f"\n{'='*72}")
                out(f"  {blanket}  |  {loading} kg/m³")
                out(f"{'='*72}")

                sp_A = find_statepoint(blanket, 'A', loading)
                sp_C = find_statepoint(blanket, 'C', loading)

                elo, ehi, emid, rate_A = extract_spectrum(sp_A, TALLY_NAME, SCORE, NUCLIDES)
                _,   _,   _,    rate_C = extract_spectrum(sp_C, TALLY_NAME, SCORE, NUCLIDES)

                flux_emid_A, flux_A = extract_flux(sp_A)
                flux_emid_C, flux_C = extract_flux(sp_C)

                diff = rate_C - rate_A                    # positive ⟹ C > A
                mask = diff > 0                           # bins where lattice wins

                total_A = rate_A.sum()
                total_C = rate_C.sum()

                out(f"\n  Total (n,gamma) rate:  A = {total_A:.6e}   C = {total_C:.6e}"
                    f"   C/A = {total_C/total_A:.4f}")

                # Mode of flux spectrum
                flux_mode_A = flux_emid_A[np.argmax(flux_A)]
                flux_mode_C = flux_emid_C[np.argmax(flux_C)]
                out(f"  Flux mode energy:     A = {flux_mode_A:.4e} eV   "
                    f"C = {flux_mode_C:.4e} eV")

                # Mode of (n,gamma) rate spectrum
                ngamma_mode_A = emid[np.argmax(rate_A)]
                ngamma_mode_C = emid[np.argmax(rate_C)]
                out(f"  (n,gamma) mode energy: A = {ngamma_mode_A:.4e} eV   "
                    f"C = {ngamma_mode_C:.4e} eV")

                out(f"  Bins where C > A:  {mask.sum()} / {len(mask)}")

                if mask.any():
                    out(f"\n  {'Bin':>4s}   {'E_lo [eV]':>12s}   {'E_hi [eV]':>12s}"
                        f"   {'Rate A':>12s}   {'Rate C':>12s}   {'C - A':>12s}   {'C/A':>8s}")
                    out(f"  {'-'*80}")
                    for i in np.where(mask)[0]:
                        ratio = rate_C[i] / rate_A[i] if rate_A[i] > 0 else np.inf
                        out(f"  {i:4d}   {elo[i]:12.4e}   {ehi[i]:12.4e}"
                            f"   {rate_A[i]:12.4e}   {rate_C[i]:12.4e}"
                            f"   {diff[i]:12.4e}   {ratio:8.4f}")
                else:
                    out(f"\n  No bins where C > A.")

    print(f"\nSaved to {outpath}")


if __name__ == '__main__':
    main()