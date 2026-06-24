#!/usr/bin/env python3
"""
xs_actinides.py  --  Flux-averaged cross sections for U-235, U-238, Th-232
==========================================================================

Prints a table of flux-averaged microscopic σ [barn] and macroscopic Σ [1/cm]
for (n,gamma) and (n,fission) of U-235 and U-238/Th-232, across the three
blankets and both dopant types.

%change columns show the shift in flux-averaged σ relative to the lowest
loading, isolating the spectral self-shielding effect.

Uses the same material builders and flux spectra as xs.py / xs_arithmetic.py
(via common.py).

    python xs_actinides.py
"""

import os
import sys
import numpy as np

from common import *

# ── configuration ──────────────────────────────────────────────────────────

BLANKETS = ['FLiBe', 'HCPB', 'DCLL']
LOADINGS = [0.1, 150.0, 1000.0]

FLUX_DIR = './Data'
FLUX_FILES_U = {
    'FLiBe': 'FLiBe_900K_Li07.5_U238_flux.csv',
    'HCPB':  'HCPB_900K_Li60.0_U238_flux.csv',
    'DCLL':  'DCLL_900K_Li90.0_U238_flux.csv',
}
FLUX_FILES_Th = {
    'FLiBe': 'FLiBe_900K_Li07.5_Th232_flux.csv',
    'HCPB':  'HCPB_900K_Li60.0_Th232_flux.csv',
    'DCLL':  'DCLL_900K_Li90.0_Th232_flux.csv',
}

MT_GAMMA   = 102
MT_FISSION = 18

# ── intermediary-decay correction (Eq. 7) ─────────────────────────────────
HALF_LIFE_NP239 = 2.356    # [days]  U-238 → U-239 → *Np-239* → Pu-239
HALF_LIFE_PA233 = 26.98    # [days]  Th-232 → Th-233 → *Pa-233* → U-233
T_EVAL = 365.0             # [days]  evaluation period
FISSILE_PRODUCT = {'U238': 'Pu239', 'Th232': 'U233'}

# Average blanket scalar flux for in-situ destruction rate μ = <σ_a> · Φ
# Typical fusion blanket: ~1e13 (deep/outboard) to ~5e14 (near first wall)
PHI_ABS = 1e14                             # [n/cm²/s]  ← adjust to your geometry
SECONDS_PER_DAY = 86400.0


# ── helpers ────────────────────────────────────────────────────────────────

def load_nuclide(lib, name):
    """Return (IncidentNeutron, temperature_str) for *name*."""
    path = find_path(lib, name)
    if path is None:
        sys.exit(f"No data file for {name} in the library.")
    nuc = openmc.data.IncidentNeutron.from_hdf5(path)
    T = pick_temp(nuc, TARGET_TEMP_K)
    return nuc, T


def sigma_on_grid(nuc, T, mt, energies):
    """Microscopic sigma(E) [barn] for MT *mt* evaluated on *energies*."""
    rx = nuc.reactions.get(mt)
    if rx is None or T not in rx.xs:
        return None
    return rx.xs[T](energies)


def flux_avg(E_xs, xs_vals, phi_E, phi_w, phi_sum):
    """<sigma>_phi = sum(sigma(E_i) * phi_i) / sum(phi_i)."""
    xs_at_bins = np.interp(phi_E, E_xs, xs_vals)
    return float(np.dot(xs_at_bins, phi_w) / phi_sum)


def pct(val, ref):
    """% change string relative to *ref*, or 'ref' when they're the same."""
    if ref is None:
        return 'ref'
    return f"{(val / ref - 1) * 100:+.1f}%"


def fissile_correction(half_life_days, t_days=T_EVAL):
    """M(t)/(r·t) correction factor for intermediary-nuclide buildup (Eq. 7).

    Returns the fraction of naive production actually realised after *t_days*,
    accounting for the finite decay time of the longest-lived intermediary.
    """
    lam = np.log(2) / half_life_days
    return 1.0 - (1.0 - np.exp(-lam * t_days)) / (lam * t_days)


def fissile_inventory_factor(lam, mu, t):
    """Full correction factor N_fis(t)/(R·t) for BOTH intermediary decay AND
    in-situ fissile absorption.  Reduces to fissile_correction() when μ→0.

    lam : intermediary decay constant [1/time]
    mu  : fissile absorption rate = <σ_a>_fissile · Φ [1/time]
    t   : evaluation time [same time unit]
    """
    if mu * t < 1e-12:
        return 1.0 - (1.0 - np.exp(-lam * t)) / (lam * t)
    if abs(lam - mu) * t < 1e-8:
        elt = np.exp(-lam * t)
        return (1.0 - elt) / (lam * t) - elt
    return (
        (1.0 - np.exp(-mu * t)) / (mu * t)
        - (np.exp(-mu * t) - np.exp(-lam * t)) / ((lam - mu) * t)
    )


# ── main ───────────────────────────────────────────────────────────────────

def main():
    lib = open_library()

    nuc_data = {}
    for name in ['U235', 'U238', 'Th232', 'Pu239', 'U233']:
        nuc_data[name] = load_nuclide(lib, name)

    flux_u  = load_flux_spectra(flux_dir=FLUX_DIR, flux_files=FLUX_FILES_U)
    flux_th = load_flux_spectra(flux_dir=FLUX_DIR, flux_files=FLUX_FILES_Th)

    base_loading = min(LOADINGS)

    corr_rows = []

    # ── header ─────────────────────────────────────────────────────────────
    W = 10   # column width for data values
    P = 7    # column width for %change

    u5  = 'U-235'
    frt = 'U-238/Th-232'

    R = 10   # column width for R/src-n

    h1 = (f"{'':15} {'':>6}"
          f"  {'--- (n,gamma) ---':^{4*W+6}}"
          f"  {'--- (n,fission) ---':^{4*W+6}}"
          f"  {'--- R/src-n ---':^{4*R+6}}"
          f"  {'--- %Δσ ---':^{4*P+6}}")

    h2 = (f"{'Blanket-Dopant':15} {'kg/m³':>6}"
          f"  {'σγ_'+u5:>{W}}  {'σγ_fert':>{W}}"
          f"  {'Σγ_'+u5:>{W}}  {'Σγ_fert':>{W}}"
          f"  {'σf_'+u5:>{W}}  {'σf_fert':>{W}}"
          f"  {'Σf_'+u5:>{W}}  {'Σf_fert':>{W}}"
          f"  {'Rγ_'+u5:>{R}}  {'Rγ_fert':>{R}}"
          f"  {'Rf_'+u5:>{R}}  {'Rf_fert':>{R}}"
          f"  {'γ_'+u5:>{P}}  {'γ_fert':>{P}}"
          f"  {'f_'+u5:>{P}}  {'f_fert':>{P}}")

    units = (f"{'':15} {'':>6}"
             f"  {'[barn]':>{W}}  {'[barn]':>{W}}"
             f"  {'[1/cm]':>{W}}  {'[1/cm]':>{W}}"
             f"  {'[barn]':>{W}}  {'[barn]':>{W}}"
             f"  {'[1/cm]':>{W}}  {'[1/cm]':>{W}}"
             f"  {'[rxn]':>{R}}  {'[rxn]':>{R}}"
             f"  {'[rxn]':>{R}}  {'[rxn]':>{R}}"
             f"  {'':>{P}}  {'':>{P}}"
             f"  {'':>{P}}  {'':>{P}}")

    print(f"%change is on flux-averaged σ, relative to {base_loading:g} kg/m³")
    print(f"'fert' = U-238 when dopant is U, Th-232 when dopant is Th\n")
    print(h1)
    print(h2)
    print(units)
    print('─' * len(h2))

    # ── data rows ──────────────────────────────────────────────────────────
    for blanket in BLANKETS:
        for fiso, short in [('U238', 'U'), ('Th232', 'Th')]:
            flux_set = flux_u if fiso == 'U238' else flux_th
            blk_flux = flux_set.get(blanket, {})
            fert_nuc = fiso

            base = {}   # keyed by (nuclide, mt) -> base σ value

            for loading in LOADINGS:
                mat, _ = build_blanket(blanket, loading, fiso)
                dens = atom_densities(mat)

                phi_E = phi_w = phi_sum = None
                if blk_flux:
                    matched = match_loading(blk_flux.keys(), loading)
                    if matched is not None:
                        phi_E, phi_w = blk_flux[matched]
                        phi_sum = phi_w.sum()
                        if phi_sum <= 0:
                            phi_E = phi_w = phi_sum = None

                if phi_E is None:
                    continue

                has_u235 = (fiso == 'U238')
                blank = ' ' * W

                vals = {}  # (nuc_key, mt) -> (sig_micro, Sig_macro)
                nucs_to_eval = [(fert_nuc, MT_GAMMA), (fert_nuc, MT_FISSION)]
                if has_u235:
                    nucs_to_eval = [('U235', MT_GAMMA), ('U235', MT_FISSION)] + nucs_to_eval
                for nuc_key, mt in nucs_to_eval:
                    nuc, T = nuc_data[nuc_key]
                    E = nuc.energy[T]
                    sig = sigma_on_grid(nuc, T, mt, E)
                    if sig is None:
                        vals[(nuc_key, mt)] = (0.0, 0.0)
                        continue
                    s = flux_avg(E, sig, phi_E, phi_w, phi_sum)
                    N = dens.get(nuc_key, 0.0)
                    vals[(nuc_key, mt)] = (s, N * s)

                # store base-case sigmas for %change
                is_base = (loading == base_loading)
                if is_base:
                    for k, (s, _) in vals.items():
                        base[k] = s

                # unpack fertile (always present)
                sgf,  Sgf  = vals[(fert_nuc,  MT_GAMMA)]
                sff,  Sff  = vals[(fert_nuc,  MT_FISSION)]
                Rgf = Sgf * phi_sum
                Rff = Sff * phi_sum

                # %change strings
                if is_base:
                    pgf = pff = 'ref'
                else:
                    pgf = pct(sgf, base.get((fert_nuc, MT_GAMMA)))
                    pff = pct(sff, base.get((fert_nuc, MT_FISSION)))

                if has_u235:
                    sg5, Sg5 = vals[('U235', MT_GAMMA)]
                    sf5, Sf5 = vals[('U235', MT_FISSION)]
                    Rg5 = Sg5 * phi_sum
                    Rf5 = Sf5 * phi_sum
                    if is_base:
                        pg5 = pf5 = 'ref'
                    else:
                        pg5 = pct(sg5, base.get(('U235', MT_GAMMA)))
                        pf5 = pct(sf5, base.get(('U235', MT_FISSION)))
                    s_sg5 = f"{sg5:{W}.3e}"
                    s_Sg5 = f"{Sg5:{W}.3e}"
                    s_sf5 = f"{sf5:{W}.3e}"
                    s_Sf5 = f"{Sf5:{W}.3e}"
                    s_Rg5 = f"{Rg5:{W}.3e}"
                    s_Rf5 = f"{Rf5:{W}.3e}"
                else:
                    pg5 = pf5 = ''
                    s_sg5 = s_Sg5 = s_sf5 = s_Sf5 = blank
                    s_Rg5 = s_Rf5 = blank

                label = f"{blanket}-{short}"
                print(f"{label:15} {loading:>6g}"
                      f"  {s_sg5}  {sgf:{W}.3e}"
                      f"  {s_Sg5}  {Sgf:{W}.3e}"
                      f"  {s_sf5}  {sff:{W}.3e}"
                      f"  {s_Sf5}  {Sff:{W}.3e}"
                      f"  {s_Rg5}  {Rgf:{W}.3e}"
                      f"  {s_Rf5}  {Rff:{W}.3e}"
                      f"  {pg5:>{P}}  {pgf:>{P}}"
                      f"  {pf5:>{P}}  {pff:>{P}}")

                # Fissile-product absorption cross section (Pu-239 or U-233)
                fis_prod = FISSILE_PRODUCT[fiso]
                nuc_fp, T_fp = nuc_data[fis_prod]
                E_fp = nuc_fp.energy[T_fp]
                sig_f_fp = sigma_on_grid(nuc_fp, T_fp, MT_FISSION, E_fp)
                sig_g_fp = sigma_on_grid(nuc_fp, T_fp, MT_GAMMA, E_fp)
                sa_f = flux_avg(E_fp, sig_f_fp, phi_E, phi_w, phi_sum) if sig_f_fp is not None else 0.0
                sa_g = flux_avg(E_fp, sig_g_fp, phi_E, phi_w, phi_sum) if sig_g_fp is not None else 0.0

                mu_d = (sa_f + sa_g) * 1e-24 * PHI_ABS * SECONDS_PER_DAY
                hl = HALF_LIFE_NP239 if fiso == 'U238' else HALF_LIFE_PA233
                lam_d = np.log(2) / hl
                cf_eq7 = fissile_correction(hl)
                cf_full = fissile_inventory_factor(lam_d, mu_d, T_EVAL)

                corr_rows.append((label, loading, Rgf,
                                  sa_f, sa_g, mu_d, cf_eq7, cf_full))

        print()

    # ── corrected fissile production (decay + in-situ absorption) ───────
    print('=' * 90)
    print(f"Corrected fissile production, t = {T_EVAL:.0f} d, "
          f"Φ = {PHI_ABS:.1e} n/cm²/s")
    print()
    print("  dN_int/dt = R − λ·N_int            (intermediary buildup & decay)")
    print("  dN_fis/dt = λ·N_int − μ·N_fis      (fissile production − absorption)")
    print()
    print("  N_fis(t) = R/μ·(1−e^{−μt}) − R/(λ−μ)·(e^{−μt} − e^{−λt})")
    print("  Factor = N_fis(t)/(R·t)  →  corrected kg = r_naive [kg/yr] × factor")
    print()
    lam_np = np.log(2) / HALF_LIFE_NP239
    lam_pa = np.log(2) / HALF_LIFE_PA233
    print(f"  λ(Np-239) = {lam_np:.4f} /d  (t½ = {HALF_LIFE_NP239} d)   "
          f"λ(Pa-233) = {lam_pa:.4f} /d  (t½ = {HALF_LIFE_PA233} d)")
    print(f"  μ = <σ_a>_fissile · Φ   "
          f"(Pu-239 for U chain, U-233 for Th chain)")
    print()

    print(f"{'Case':15} {'kg/m³':>6}  {'σ_f(fis)':>9} {'σ_γ(fis)':>9} {'σ_a(fis)':>9}"
          f"  {'μ [1/d]':>10}  {'μt':>7}  {'Eq.7':>8}  {'Full':>8}")
    print(f"{'':15} {'':>6}  {'[barn]':>9} {'[barn]':>9} {'[barn]':>9}"
          f"  {'':>10}  {'':>7}  {'factor':>8}  {'factor':>8}")
    print('─' * 90)
    for label, ld, rn, sa_f, sa_g, mu_d, cf7, cff in corr_rows:
        sa = sa_f + sa_g
        mu_t = mu_d * T_EVAL
        print(f"{label:15} {ld:>6g}  {sa_f:>9.1f} {sa_g:>9.1f} {sa:>9.1f}"
              f"  {mu_d:>10.3e}  {mu_t:>7.4f}  {cf7:>8.6f}  {cff:>8.6f}")


if __name__ == '__main__':
    main()
