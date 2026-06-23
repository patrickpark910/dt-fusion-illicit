#!/usr/bin/env python3
"""
slowing_down.py  --  Elastic slowing-down parameters of blanket moderators
===========================================================================

Prints, for F-19, Be-9, Li-7 (from the FLiBe blanket), Be-9, Li-7 (from the HCPB
blanket), and elemental Pb (from the DCLL Pb-17Li blanket), the four moderation
figures of merit over the energy interval  E_0 -> E_final  (default 10 MeV -> 1 keV):

  Sigma_e   macroscopic elastic cross section            [1/cm]
            Sigma_e(E) = N * sigma_e(E);  for elemental Pb, summed over isotopes
            Sigma_e(E) = sum_i N_i * sigma_e,i(E).  Reported as the lethargy-
            averaged value over [E_final, E_0] (the average that enters the
            integrals below), with the endpoint values shown separately.

  xi        average logarithmic energy decrement per collision  ("e-folds/coll")
            xi_i = 1 + [a/(1-a)] ln a,   a = [(A-1)/(A+1)]^2.
            For elemental Pb: the elastic-scattering-weighted average of the
            isotopic xi_i.  (Energy-independent for elastic scattering.)
            Number of collisions to cross the interval: n = ln(E_0/E_final) / xi.

  tau_s     slowing-down time                              [s]
            tau_s = integral_{E_final}^{E_0}  dE / (E * xi * Sigma_e(E) * v(E)),
            v(E) = sqrt(2E/m_n).  Dominated by the low-energy end (slow v).

  l_s       slowing-down length                            [cm]
            l_s = sqrt(tau_age),  the Fermi age
            tau_age = integral_{E_final}^{E_0}  D(E)/(xi*Sigma_e(E)) * dE/E,
            D = 1/(3*Sigma_tr),  Sigma_tr = Sigma_e*(1 - mu0),  mu0 = 2/(3A).
            The rms crow-flight displacement to E_final is sqrt(6)*l_s.

All four are ELASTIC-only, classical-slowing-down quantities; see the caveats in
the printed header (inelastic scattering, anisotropy, absorption are ignored).

Constants and material builders live in common.py (imported below); this file
handles the slowing-down computation.

-----------------------------------------------------------------------------
DATA INPUTS YOU NEED TO PROVIDE
-----------------------------------------------------------------------------
1. An ENDF/B-VIII.0 HDF5 cross-section library that includes 900 K data -- the
   SAME library the rest of your model (and xs.py) uses.  Point it to either
   openmc.config['cross_sections'] or the OPENMC_CROSS_SECTIONS env var.
   It supplies the microscopic ELASTIC cross section sigma_e(E) = MT 2 for:
       F19, Be9, Li7, Pb204, Pb206, Pb207, Pb208.

2. Number densities N [atom/b-cm].  By default these are computed for you by
   building your FLiBe, HCPB, and DCLL materials (the builders below are copied
   verbatim from flibe.py / hcpb.py / dcll.py) at the loading FERTILE_KGM3 -- so
   you do NOT have to supply anything beyond the data library.  If you would rather
   hand the densities in directly (e.g. to decouple from OpenMC, or to use a
   specific tally-derived value), fill NDENS_OVERRIDE below.

That's it -- there are no other external inputs.  Compositions, atomic masses,
and the transport correction mu0 = 2/(3A) are all hard-coded / derived here.

    python slowing_down.py
"""

import os
import sys
import math
import numpy as np

import openmc
import openmc.data

from common import *


# =============================================================================
# CONFIGURATION
# =============================================================================

E0       = 10.0e6     # eV  -- start (top) energy of the slowing-down interval
E_FINAL  = 1.0e3      # eV  -- final (bottom) energy

# Loading at which the number densities are evaluated.  0.1 kg/m3 is the lowest
# nonzero loading in your sweep and is effectively the pure-breeder matrix
# (fuel is ~0.002 vol%), which is what you want to characterise a moderator.
FERTILE_KGM3   = 0.0
FERTILE_ISOTOPE = 'U238'   # only affects DCLL (UO2 vs ThO2 BISO); Pb is unaffected

# Targets: (display label, host blanket for N, [(nuclide, mass number A), ...]).
# F/Be/Li come from FLiBe; Be/Li from HCPB; elemental Pb (4 isotopes) from DCLL.
TARGETS = [
    ('F-19', 'FLiBe', [('F19', 19)]),
    ('Be-9', 'FLiBe', [('Be9',  9)]),
    ('Li-7', 'FLiBe', [('Li7',  7)]),
    ('Be-9', 'HCPB',  [('Be9',  9)]),
    ('Li-7', 'HCPB',  [('Li7',  7)]),
    ('Pb',   'DCLL',  [('Pb204', 204), ('Pb206', 206),
                        ('Pb207', 207), ('Pb208', 208)]),
]

# Optional manual number densities [atom/b-cm], keyed by nuclide name.  Anything
# listed here overrides the value taken from the built material.  Leave empty {}
# to use the material-derived densities.
NDENS_OVERRIDE = {}



# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

EV_J = 1.602176634e-19      # J per eV
M_N  = 1.67492749804e-27    # neutron mass [kg]

# np.trapz was renamed to np.trapezoid in NumPy 2.0 (trapz removed later).
# Use whichever this NumPy provides so the script runs on old and new versions.
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
if _trapz is None:                       # extremely old/odd NumPy: simple fallback
    def _trapz(y, x):
        y, x = np.asarray(y, float), np.asarray(x, float)
        return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))




# =============================================================================
# CROSS-SECTION DATA HELPERS
# =============================================================================



def load_elastic(names, lib):
    """name -> (sigma_callable, native_energy_grid[eV], temp_string), or None."""
    cache = {}
    for name in names:
        path = find_path(lib, name)
        if path is None:
            print(f"  ! no data file for {name} -- skipped")
            cache[name] = None
            continue
        nuc = openmc.data.IncidentNeutron.from_hdf5(path)
        T   = pick_temp(nuc, TARGET_TEMP_K)
        rx  = nuc.reactions.get(2)            # MT 2 = elastic
        if rx is None or T not in rx.xs:
            print(f"  ! {name}: no elastic (MT=2) at {T} -- skipped")
            cache[name] = None
            continue
        cache[name] = (rx.xs[T], np.asarray(nuc.energy[T]), T)
    return cache


# =============================================================================
# SLOWING-DOWN PHYSICS
# =============================================================================

def alpha_of(A):
    r = (A - 1.0) / (A + 1.0)
    return r * r


def xi_of(A):
    """Average logarithmic energy decrement per elastic collision."""
    a = alpha_of(A)
    if a <= 0.0:               # A == 1 limit
        return 1.0
    return 1.0 + a / (1.0 - a) * math.log(a)


def mu0_of(A):
    """Mean scattering cosine, s-wave (low-energy) elastic value."""
    return 2.0 / (3.0 * A)


def speed_cm_s(E_eV):
    """Neutron speed [cm/s] (classical) for energy E [eV]; array-safe."""
    E_J = np.asarray(E_eV, dtype=float) * EV_J
    return np.sqrt(2.0 * E_J / M_N) * 100.0


def moderation(isotopes, ndens, cache):
    """Compute the slowing-down parameters for one target (one or more isotopes).

    isotopes : list of (nuclide_name, mass_number_A)
    ndens    : {nuclide_name: number density [atom/b-cm]}
    cache    : output of load_elastic
    Returns a results dict, or None if no usable data.
    """
    members = [(n, A) for (n, A) in isotopes
               if cache.get(n) is not None and ndens.get(n, 0.0) > 0.0]
    if not members:
        return None

    # integration grid: union of the isotopes' native grids, clipped to range
    grids = [cache[n][1] for (n, _) in members]
    Eu = np.unique(np.concatenate(grids))
    Eu = Eu[(Eu >= E_FINAL) & (Eu <= E0)]
    Eu = np.unique(np.concatenate([[E_FINAL], Eu, [E0]]))

    # build Sigma_e(E), the xi-weighting numerator, and Sigma_tr(E)
    Sig   = np.zeros_like(Eu)     # Sigma_e            [1/cm]
    xinum = np.zeros_like(Eu)     # sum_i N_i sigma_i xi_i
    Str   = np.zeros_like(Eu)     # Sigma_tr           [1/cm]
    N_tot = 0.0
    for (n, A) in members:
        sig_fn, _, _ = cache[n]
        Ni = ndens[n]
        N_tot += Ni
        si = np.asarray(sig_fn(Eu), dtype=float)   # barns; N in atom/b-cm -> 1/cm
        contrib = Ni * si
        Sig   += contrib
        xinum += contrib * xi_of(A)
        Str   += contrib * (1.0 - mu0_of(A))

    if not np.all(Sig > 0.0):
        # elastic XS should be strictly positive across the range; guard anyway
        Sig = np.where(Sig > 0.0, Sig, np.finfo(float).tiny)

    xibar = xinum / Sig                            # effective xi(E)
    v     = speed_cm_s(Eu)                         # cm/s
    u     = math.log(E0 / E_FINAL)                 # total lethargy span

    # lethargy averages: <X> = integral X dlnE / u  = trapz(X/E, E) / u
    Sig_avg = _trapz(Sig / Eu, Eu) / u
    xi_avg  = _trapz(xibar / Eu, Eu) / u

    # slowing-down time:  tau_s = integral dE / (E * xi * Sigma_e * v)
    tau_s = _trapz(1.0 / (Eu * xibar * Sig * v), Eu)

    # Fermi age:  tau_age = integral D/(xi*Sigma_e) * dE/E ,  D = 1/(3 Sigma_tr)
    D   = 1.0 / (3.0 * Str)
    age = _trapz(D / (xibar * Sig * Eu), Eu)
    L_s = math.sqrt(age)

    return dict(
        N=N_tot,
        Sig_avg=Sig_avg,
        Sig_E0=float(np.interp(E0, Eu, Sig)),
        Sig_Ef=float(np.interp(E_FINAL, Eu, Sig)),
        xi=xi_avg,
        n_coll=u / xi_avg,
        tau_s=tau_s,
        age=age,
        L_s=L_s,
        lambda_s=1.0 / Sig_avg,
    )


# =============================================================================
# ELASTIC MODERATION POWER BREAKDOWN (console printout)
# =============================================================================

def _parse_nuclide(name):
    """'Li7' -> ('Li', 7), 'Pb208' -> ('Pb', 208), 'He4' -> ('He', 4)."""
    i = 0
    while i < len(name) and name[i].isalpha():
        i += 1
    elem = name[:i]
    A    = int(name[i:]) if i < len(name) else 0
    return elem, A




def print_moderation_breakdown(host_dens, lib):
    """Print ranked elastic Sigma and moderating-power breakdown for every
    element in each blanket material.

    When a flux-spectrum CSV is available for a blanket, the table reports
    the flux-weighted average

        <Sigma_el>_phi = N * sum[sigma_el(E_i) * phi_i] / sum[phi_i]

    where the sum runs over the tally energy bins.  xi is energy-independent
    for elastic scattering, so  <xi * Sigma_el>_phi = xi * <Sigma_el>_phi.

    When no flux file is found the function falls back to a point evaluation
    at E0 (the top of the slowing-down interval).

    Tables are printed at the baseline loading (FERTILE_KGM3, ~pure breeder)
    and at the high loading (1000 kg/m3) to show how fertile material shifts
    the moderation budget.  Both uranium (UF4/UO2) and thorium (ThF4/ThO2)
    fertile forms are shown.
    """
    LOADINGS = [FERTILE_KGM3, 1000.0]   # kg/m3
    FERTILE_ISOTOPES = ['U238', 'Th232']

    # fertile-compound display labels
    _FERT_LABEL = {
        ('FLiBe', 'U238'):  'UF4',   ('FLiBe', 'Th232'): 'ThF4',
        ('HCPB',  'U238'):  'UO2',   ('HCPB',  'Th232'): 'ThO2',
        ('DCLL',  'U238'):  'UO2',   ('DCLL',  'Th232'): 'ThO2',
    }

    # build materials at every (blanket, loading, fertile_isotope) combination
    all_dens = {}   # (blanket, loading, fiso) -> {nuclide: N}
    for fiso in FERTILE_ISOTOPES:
        for loading in LOADINGS:
            if loading == FERTILE_KGM3 and FERTILE_KGM3 <= 0.0:
                # At zero loading the fertile form doesn't matter; share one set
                if fiso == FERTILE_ISOTOPES[0]:
                    for bname in ('FLiBe', 'HCPB', 'DCLL'):
                        all_dens[(bname, loading, fiso)] = host_dens[bname]
                else:
                    for bname in ('FLiBe', 'HCPB', 'DCLL'):
                        all_dens[(bname, loading, fiso)] = host_dens[bname]
            else:
                all_dens[('FLiBe', loading, fiso)] = atom_densities(
                    build_flibe_blanket(loading, fiso)[0])
                all_dens[('HCPB',  loading, fiso)] = atom_densities(
                    build_hcpb_blanket(loading, fiso)[0])
                all_dens[('DCLL',  loading, fiso)] = atom_densities(
                    build_dcll_blanket(loading, fiso)[0])

    # load elastic XS for every nuclide across all cases
    all_nucs = sorted({n for dens in all_dens.values()
                       for n, N in dens.items() if N > 0})
    cache = load_elastic(all_nucs, lib)

    # load flux spectra (empty dict per blanket if file is missing)
    flux = load_flux_spectra()

    SIG  = "\u03a3"            # Σ
    XI   = "\u03be"            # ξ
    DOT  = "\u00b7"            # ·
    M3   = "m\u00b3"           # m³

    for bname in ('FLiBe', 'HCPB', 'DCLL'):
        for fiso in FERTILE_ISOTOPES:
            for loading in LOADINGS:
                # skip duplicate zero-loading tables for second fertile isotope
                if (loading == FERTILE_KGM3 and FERTILE_KGM3 <= 0.0
                        and fiso != FERTILE_ISOTOPES[0]):
                    continue

                dens = all_dens[(bname, loading, fiso)]
                flabel = _FERT_LABEL[(bname, fiso)]

                # --- resolve the flux spectrum for this case ------------------
                phi_E = phi_w = phi_sum = None
                weighting_label = f"at {E0 / 1e6:.0f} MeV"
                blk_flux = flux.get(bname, {})
                if blk_flux:
                    matched = match_loading(blk_flux.keys(), loading)
                    if matched is not None:
                        phi_E, phi_w = blk_flux[matched]
                        phi_sum = phi_w.sum()
                        if phi_sum > 0:
                            weighting_label = (f"flux-weighted  "
                                               f"(matched {matched:g} kg/{M3})")
                        else:
                            phi_E = phi_w = phi_sum = None   # fall back

                # --- accumulate per element -----------------------------------
                elem_Sig   = {}
                elem_xiSig = {}
                elem_N     = {}

                for nuc, N in dens.items():
                    if N <= 0 or cache.get(nuc) is None:
                        continue
                    sig_fn, _, _ = cache[nuc]

                    if phi_E is not None:
                        sig_at_bins = np.asarray(sig_fn(phi_E), dtype=float)
                        sigma_avg   = np.dot(sig_at_bins, phi_w) / phi_sum
                    else:
                        sigma_avg = float(sig_fn(E0))

                    Sigma_el = N * sigma_avg

                    elem, A = _parse_nuclide(nuc)
                    if A <= 0:
                        continue
                    xi = xi_of(A)

                    elem_Sig[elem]   = elem_Sig.get(elem, 0.0)   + Sigma_el
                    elem_xiSig[elem] = elem_xiSig.get(elem, 0.0) + xi * Sigma_el
                    elem_N[elem]     = elem_N.get(elem, 0.0)     + N

                if not elem_Sig:
                    continue

                total_Sig   = sum(elem_Sig.values())
                total_xiSig = sum(elem_xiSig.values())

                ranked = sorted(elem_Sig.keys(),
                                key=lambda e: elem_Sig[e], reverse=True)

                # --- header ---------------------------------------------------
                SIG_B = '\u03c3'   # σ  (lowercase for microscopic)
                W = 102
                print(f"\n{'=' * W}")
                print(f"  {bname} ({flabel}) -- {loading:g} kg/{M3} -- elastic "
                      f"{SIG} and moderation power,  {weighting_label}")
                print(f"{'=' * W}")
                print(f"  {'Element':<8}  {'N [/b-cm]':>12}  "
                      f"{'<' + SIG_B + '_el> [b]':>12}  "
                      f"{SIG + '_el [/cm]':>12}  {XI:>8}  "
                      f"{XI + DOT + SIG + '_el [/cm]':>14}  {'% mod power':>11}")
                print(f"  {'--------':<8}  {'------------':>12}  "
                      f"{'------------':>12}  "
                      f"{'------------':>12}  {'--------':>8}  "
                      f"{'--------------':>14}  {'-----------':>11}")

                for elem in ranked:
                    Sig   = elem_Sig[elem]
                    xiSig = elem_xiSig[elem]
                    N_e   = elem_N[elem]
                    sig_micro = Sig / N_e if N_e > 0 else 0.0
                    xi_eff = xiSig / Sig if Sig > 0 else 0.0
                    pct    = 100.0 * xiSig / total_xiSig if total_xiSig > 0 else 0.0

                    print(f"  {elem:<8}  {N_e:12.4e}  "
                          f"{sig_micro:12.4e}  "
                          f"{Sig:12.4e}  {xi_eff:8.5f}  "
                          f"{xiSig:14.4e}  {pct:10.4f}%")

                print(f"  {'--------':<8}  {'------------':>12}  "
                      f"{'------------':>12}  "
                      f"{'------------':>12}  {'--------':>8}  "
                      f"{'--------------':>14}  {'-----------':>11}")
                total_N = sum(elem_N.values())
                sig_micro_avg = total_Sig / total_N if total_N > 0 else 0.0
                xi_avg = total_xiSig / total_Sig if total_Sig > 0 else 0.0
                print(f"  {'Total':<8}  {total_N:12.4e}  "
                      f"{sig_micro_avg:12.4e}  "
                      f"{total_Sig:12.4e}  {xi_avg:8.5f}  "
                      f"{total_xiSig:14.4e}  {'100.0%':>11}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    lib = open_library()

    # number densities from the built materials (unless overridden)
    host_dens = {
        'FLiBe': atom_densities(build_flibe_blanket(FERTILE_KGM3, FERTILE_ISOTOPE)[0]),
        'HCPB':  atom_densities(build_hcpb_blanket(FERTILE_KGM3, FERTILE_ISOTOPE)[0]),
        'DCLL':  atom_densities(build_dcll_blanket(FERTILE_KGM3, FERTILE_ISOTOPE)[0]),
    }

    def ndens_for(blanket, isotopes):
        out = {}
        for (n, _) in isotopes:
            out[n] = NDENS_OVERRIDE.get(n, host_dens[blanket].get(n, 0.0))
        return out

    # load elastic cross sections once per needed nuclide
    names = sorted({n for (_, _, isos) in TARGETS for (n, _) in isos})
    cache = load_elastic(names, lib)

    # compute every target
    results = []
    for label, blanket, isos in TARGETS:
        res = moderation(isos, ndens_for(blanket, isos), cache)
        results.append((label, blanket, res))

    u = math.log(E0 / E_FINAL)

    # ---- header ------------------------------------------------------------
    print("\n" + "=" * 85)
    print("ELASTIC SLOWING-DOWN PARAMETERS")
    print("=" * 85)
    print(f"  Energy interval : {E0:.3e} eV  ->  {E_FINAL:.3e} eV")
    print(f"  Lethargy span u : ln(E0/E_final) = {u:.4f}")
    print(f"  Temperature     : {TARGET_TEMP_K} K   (ENDF/B-VIII.0, elastic = MT 2)")
    print(f"  Number densities: from the built materials at {FERTILE_KGM3:g} kg/m3"
          f" (~pure-breeder matrix)")
    for label, blanket, _ in results:
        print(f"      N({label:<5}) <- {blanket} blanket")
    print("\n  Definitions:")
    print("      Sigma_e = N*sigma_e  (elemental Pb: sum_i N_i*sigma_e,i),"
          " lethargy-averaged")
    print("      xi      = 1 + [a/(1-a)] ln a,  a = ((A-1)/(A+1))^2  (e-folds/coll)")
    print("      n_coll  = u / xi")
    print("      tau_s   = INT dE /(E xi Sigma_e v),  v = sqrt(2E/m_n)")
    print("      l_s     = sqrt(Fermi age),  age = INT [1/(3 Sigma_tr)]/(xi Sigma_e) dE/E")
    print("\n  Caveats (all quantities are elastic-only, classical slowing-down):")
    print("    * Inelastic (n,n') scattering is IGNORED.  For Pb it dominates the")
    print("      energy loss above ~1 MeV, so the Pb tau_s / l_s / n_coll here are")
    print("      large OVER-estimates; the F-19 and Li-7 top decade is affected too.")
    print("    * mu0 = 2/(3A) is the low-energy value; real scattering is forward-")
    print("      peaked at MeV energies (l_s slightly under-estimated up high).")
    print("    * Absorption is ignored (small for these nuclides over this range).")
    print("    * E_final = 1 keV is well above thermal, so the target-at-rest")
    print("      slowing-down formulas are valid across the whole interval.")

    # ---- main table --------------------------------------------------------
    print("\n" + "-" * 85)
    hdr = (f"{'Target':<7}{'Host':<7} {'N':>10} {'<Sigma_e>':>10} {'xi':>8} "
           f"{'n_coll':>8} {'tau_s [s]':>11} {'l_s [cm]':>10}")
    print(hdr)
    print(f"{'':<7}{'':<7} {'[a/b-cm]':>10} {'[1/cm]':>10} {'':>8} "
          f"{'':>8} {'':>11} {'':>10}")
    print("-" * 85)
    for label, blanket, res in results:
        if res is None:
            print(f"{label:<7}{blanket:<7} {'(no data / zero density)':>55}")
            continue
        print(f"{label:<7}{blanket:<7} {res['N']:>10.3e} {res['Sig_avg']:>10.3e} "
              f"{res['xi']:>8.4f} {res['n_coll']:>8.1f} "
              f"{res['tau_s']:>11.3e} {res['L_s']:>10.3e}")
    print("-" * 85)

    # ---- Sigma_e energy dependence + rms displacement ----------------------
    print(f"\nElastic Sigma_e [1/cm] at the interval endpoints "
          f"(shows energy dependence):")
    print(f"  {'Target':<7}{'Host':<7} {'@ '+f'{E0:.2e} eV':>16} {'@ '+f'{E_FINAL:.2e} eV':>16}"
          f" {'mfp <1/Sig_e> [cm]':>20}")
    for label, blanket, res in results:
        if res is None:
            continue
        print(f"  {label:<7}{blanket:<7} {res['Sig_E0']:>16.3e} {res['Sig_Ef']:>16.3e}"
              f" {res['lambda_s']:>20.3e}")

    print(f"\nNote: l_s above is the slowing-down LENGTH sqrt(age).  The rms "
          f"crow-flight\n      distance from birth to {E_FINAL:.0e} eV is "
          f"sqrt(6)*l_s = {math.sqrt(6):.3f} * l_s.")
    print()

    # ---- moderation power breakdown per blanket ----------------------------
    print_moderation_breakdown(host_dens, lib)


if __name__ == "__main__":
    main()