#!/usr/bin/env python3

"""Resonance integral analytical model: NR, NRM, IR, IRM approximations.

Approximations:
  NR  - Narrow Resonance (no multiplication)
  NRM - NR with subcritical Multiplication
  IR  - Intermediate Resonance (no multiplication)
  IRM - IR with subcritical Multiplication

Generalized for arbitrary fertile (U238, Th232) and fissile (U235, None) pairs.
"""

from pathlib import Path
import argparse, os, re, sys, csv
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import openmc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "Figures/Data"
OUTS = HERE / "Data"
FIGS = HERE / "Figures_All"

os.chdir(ROOT)
sys.path.insert(0, str(ROOT))
from Python.parameters import *
from Python.utilities import *

E_MIN = {'U238': 3.0, 'Th232': 10.0}
E_MAX = 15e6
E_C = {'U238': 2e4, 'Th232': 4e3}
R2_MIN = 0.97
SIGMA_P = {'U238': 11.29, 'Th232': 12.6}
NU_BAR = {
    'Th232': (np.array([1e-5, 2.53e-2, 1e6, 3e6, 4e6, 5.7e6, 7e6, 10e6, 14.7e6, 20e6]),
              np.array([2.1047, 2.1047, 2.1817, 2.3177, 2.457, 2.6807, 3.026, 3.4, 4.002, 4.82])),
    'U238':  (np.array([3.42e6, 10.15e6, 13.15e6, 20e6]),
              np.array([2.3117, 3.1395, 3.4969, 4.3394])),
    'U235':  (np.array([3.25e6, 5.29788e6, 9.29788e6, 12.2979e6, 20e6]),
              np.array([2.0523, 2.3247, 2.8567, 3.2517, 4.276])),
}

METHODS = ('NR', 'NRM', 'IR', 'IRM')

METHOD_LINESTYLES = {'NR': '-.', 'IR': '--', 'NRM': '-', 'IRM': ':'}
METHOD_MARKERS = {'NR': 'o', 'IR': 's', 'NRM': '^', 'IRM': 'D'}
BLANKET_COLORS = {'FLiBe': '#66b420', 'HCPB': '#b41f24', 'DCLL': '#0047ba'}
BLANKET_MARKERS = {'FLiBe': 'o', 'HCPB': 's', 'DCLL': '^'}


def nu_bar(isotope, energy_ev):
    e_tab, nu_tab = NU_BAR[isotope]
    return np.interp(energy_ev, e_tab, nu_tab, left=nu_tab[0], right=nu_tab[-1])


LOADINGS = [0.0, 0.01, 0.10, 0.50, 1, 10, 25, 50, 75, 100, 150, 250, 500, 750, 1000, 2000, 3000, 4000, 4499]

# =====================================================================
# OpenMC comparison data paths
# =====================================================================
CSV_FLUX = {
    ('FLiBe', 'U238'):  DATA / "FLiBe_900K_Li07.5_U238_flux.csv",
    ('FLiBe', 'Th232'): DATA / "FLiBe_900K_Li07.5_Th232_flux.csv",
    ('HCPB', 'U238'):   DATA / "HCPB_900K_Li60.0_U238_flux.csv",
    ('HCPB', 'Th232'):  DATA / "HCPB_900K_Li60.0_Th232_flux.csv",
    ('DCLL', 'U238'):   DATA / "DCLL_900K_Li90.0_U238_flux.csv",
    ('DCLL', 'Th232'):  DATA / "DCLL_900K_Li90.0_Th232_flux.csv",
}
CSV_RXNS = {
    ('FLiBe', 'U238'):  DATA / "FLiBe_900K_Li07.5_U238_summary.csv",
    ('FLiBe', 'Th232'): DATA / "FLiBe_900K_Li07.5_Th232_summary.csv",
    ('HCPB', 'U238'):   DATA / "HCPB_900K_Li60.0_U238_summary.csv",
    ('HCPB', 'Th232'):  DATA / "HCPB_900K_Li60.0_Th232_summary.csv",
    ('DCLL', 'U238'):   DATA / "DCLL_900K_Li90.0_U238_summary.csv",
    ('DCLL', 'Th232'):  DATA / "DCLL_900K_Li90.0_Th232_summary.csv",
}

COMPOUNDS = {
    ('FLiBe', 'U238'): 'uf4',   ('FLiBe', 'Th232'): 'thf4',
    ('HCPB', 'U238'):  'uo2',   ('HCPB', 'Th232'):  'tho2',
    ('DCLL', 'U238'):  'uo2',   ('DCLL', 'Th232'):  'tho2',
}
BLANKET_LOWER = {'FLiBe': 'flibe', 'HCPB': 'hcpb', 'DCLL': 'dcll'}


def out_csv_path(method, blanket, fertile):
    return OUTS / f"integral_{method.lower()}_{BLANKET_LOWER[blanket]}_{COMPOUNDS[(blanket, fertile)]}.csv"


OUT_COMPARE_LINLOG = FIGS / "integral_vs_openmc_linlog"
OUT_COMPARE_LOGLOG = FIGS / "integral_vs_openmc_loglog"
OUT_COMPARE_LINEAR = FIGS / "integral_vs_openmc_linlin"
OUT_PSI_IEFF       = FIGS / "integral_psi_vs_ieff"


# =====================================================================
# Helpers
# =====================================================================

def temp_label(labels):
    def as_k(label):
        match = re.search(r"([0-9]+(?:\.[0-9]*)?)\s*K", str(label))
        return float(match.group(1)) if match else np.inf
    return min(labels, key=lambda label: abs(as_k(label) - TEMP_K))


def hdf5_path(library, nuclide):
    for entry in library.libraries:
        materials = entry.get("materials", ())
        if isinstance(materials, str):
            materials = (materials,)
        if entry.get("type") == "neutron" and nuclide in materials:
            path = Path(entry["path"])
            xs_dir = Path(openmc.config["cross_sections"]).parent
            return path if path.is_absolute() else xs_dir / path
    raise RuntimeError(f"{nuclide} not found in OpenMC cross section library")


def eval_xs(data, label, mt, energy_ev):
    rx = data[mt]
    label = label if label in rx.xs else temp_label(rx.xs.keys())
    return np.asarray(rx.xs[label](energy_ev))


def coarse_to_density(energy_coarse, phi_u_coarse):
    """Convert bin-integrated tally values to flux density [per eV]."""
    log_e = np.log10(energy_coarse)
    log_edges_mid = 0.5 * (log_e[:-1] + log_e[1:])
    log_edge_lo = log_e[0] - (log_edges_mid[0] - log_e[0])
    log_edge_hi = log_e[-1] + (log_e[-1] - log_edges_mid[-1])
    log_edges = np.concatenate([[log_edge_lo], log_edges_mid, [log_edge_hi]])
    bin_width = np.diff(10**log_edges)
    return np.where(bin_width > 0, phi_u_coarse / bin_width, 0.0)


def safe_eval_xs(data, label, mt, energy):
    """eval_xs returning zeros on KeyError."""
    try:
        return eval_xs(data, label, mt, energy)
    except KeyError:
        return np.zeros_like(energy)


# =====================================================================
# Material builders
# =====================================================================

def make_flibe(fertile_kgm3, fertile_isotope='U238'):
    flibe = openmc.Material(name="2(LiF)-BeF2", temperature=TEMP_K)
    flibe.set_density("g/cm3", DENSITY_FLIBE)
    flibe.add_elements_from_formula("F4Li2Be", "ao", enrichment_target="Li6",
                                    enrichment_type="ao", enrichment=ENRICH_FLIBE)

    if fertile_isotope == 'U238':
        xf4 = openmc.Material(name='UF4', temperature=TEMP_K)
        xf4.add_elements_from_formula('UF4', 'ao', enrichment=ENRICH_U)
        xf4.set_density('g/cm3', DENSITY_UF4)
    elif fertile_isotope == 'Th232':
        xf4 = openmc.Material(name='ThF4', temperature=TEMP_K)
        xf4.add_elements_from_formula('ThF4', 'ao')
        xf4.set_density('g/cm3', DENSITY_ThF4)

    vf_xf4, vf_flibe = calc_xf4_vol_fracs(fertile_kgm3, fertile_isotope)
    mix = openmc.Material.mix_materials([xf4, flibe], [vf_xf4, vf_flibe], "vo")
    mix.temperature = TEMP_K
    return mix


def make_dcll(fertile_kgm3, fertile_isotope='U238'):
    pbli = openmc.Material(name='Pb-17Li', temperature=TEMP_K)
    pbli.set_density('g/cm3', DENSITY_DCLL)
    pbli.add_element('Pb', 0.83, percent_type='ao')
    pbli.add_element('Li', 0.17, percent_type='ao', enrichment_target='Li6',
                     enrichment_type='ao', enrichment=ENRICH_DCLL)

    if fertile_isotope == 'U238':
        kernel = openmc.Material(name='UO2', temperature=TEMP_K)
        kernel.add_elements_from_formula('UO2', enrichment=ENRICH_U)
        kernel.set_density('g/cm3', DENSITY_UO2)
    elif fertile_isotope == 'Th232':
        kernel = openmc.Material(name='ThO2', temperature=TEMP_K)
        kernel.add_elements_from_formula('ThO2')
        kernel.set_density('g/cm3', DENSITY_ThO2)

    sic = openmc.Material(name='SiC', temperature=TEMP_K)
    sic.add_elements_from_formula('SiC')
    sic.set_density('g/cm3', 3.2)

    biso = openmc.Material.mix_materials([kernel, sic],
                                         [BISO_KERNEL_VOL_FRAC, BISO_COAT_VOL_FRAC], 'vo')

    f82h = openmc.Material(name='F82H', temperature=TEMP_K)
    f82h.set_density('g/cm3', 7.78)
    f82h.add_element('Fe', 89.3686, percent_type='wo')
    f82h.add_element('C',  0.1000, percent_type='wo')
    f82h.add_element('Si', 0.1000, percent_type='wo')
    f82h.add_element('Mn', 0.1300, percent_type='wo')
    f82h.add_element('Cr', 8.1600, percent_type='wo')
    f82h.add_element('W',  1.9400, percent_type='wo')
    f82h.add_element('V',  0.2000, percent_type='wo')
    f82h.add_element('N',  0.0014, percent_type='wo')

    he = openmc.Material(name='He')
    he.set_density('g/cm3', 0.004)
    he.add_element('He', 1)

    vf_biso, vf_breeder, _ = calc_biso_vol_fracs(fertile_kgm3, fertile_isotope=fertile_isotope)
    breeder = openmc.Material.mix_materials([biso, pbli], [vf_biso, vf_breeder], 'vo')
    blanket = openmc.Material.mix_materials(
        [breeder, f82h, sic, he],
        [DCLL_VF_LL_NOM, DCLL_VF_FS_NOM, DCLL_VF_SI_NOM, DCLL_VF_HE_NOM], 'vo')
    blanket.temperature = TEMP_K
    return blanket


def make_hcpb(fertile_kgm3, fertile_isotope='U238'):
    li4sio4 = openmc.Material(name='Li4SiO4', temperature=TEMP_K)
    li4sio4.set_density('g/cm3', DENSITY_LI4SIO4)
    li4sio4.add_elements_from_formula('Li4SiO4', enrichment_target='Li6',
                                      enrichment_type='ao', enrichment=ENRICH_HCPB)

    be = openmc.Material(name='Be', temperature=TEMP_K)
    be.set_density('g/cm3', DENSITY_BE)
    be.add_element('Be', 1, percent_type='wo')

    if fertile_isotope == 'U238':
        kernel = openmc.Material(name='UO2', temperature=TEMP_K)
        kernel.add_elements_from_formula('UO2', enrichment=ENRICH_U)
        kernel.set_density('g/cm3', DENSITY_UO2)
    elif fertile_isotope == 'Th232':
        kernel = openmc.Material(name='ThO2', temperature=TEMP_K)
        kernel.add_elements_from_formula('ThO2')
        kernel.set_density('g/cm3', DENSITY_ThO2)

    sic = openmc.Material(name='SiC', temperature=TEMP_K)
    sic.add_elements_from_formula('SiC')
    sic.set_density('g/cm3', 3.2)

    biso = openmc.Material.mix_materials([kernel, sic],
                                         [BISO_KERNEL_VOL_FRAC, BISO_COAT_VOL_FRAC], 'vo')

    eurofer = openmc.Material(name='Eurofer', temperature=TEMP_K)
    eurofer.set_density('g/cm3', 7.8)
    eurofer.add_element('Fe', 89.36, percent_type='wo')
    eurofer.add_element('C',  0.11, percent_type='wo')
    eurofer.add_element('Cr', 9.00, percent_type='wo')
    eurofer.add_element('W',  1.10, percent_type='wo')
    eurofer.add_element('Mn', 0.40, percent_type='wo')
    eurofer.add_element('N',  0.03, percent_type='wo')

    he = openmc.Material(name='He')
    he.set_density('atom/b-cm', 0.00049800000)
    he.add_element('He', 1)

    vf_biso, vf_breeder, _ = calc_biso_vol_fracs(fertile_kgm3, fertile_isotope=fertile_isotope)
    breeder = openmc.Material.mix_materials(
        [biso, li4sio4, be],
        [vf_biso, vf_breeder * HCPB_VF_LI_IN_BREEDER_NOM,
         vf_breeder * HCPB_VF_BE_IN_BREEDER_NOM], 'vo')
    blanket = openmc.Material.mix_materials(
        [breeder, eurofer, he],
        [HCPB_VF_BREEDER_NOM, HCPB_VF_EU_NOM, HCPB_VF_HE_NOM], 'vo')
    blanket.temperature = TEMP_K
    return blanket


MAKE_MAT = {'FLiBe': make_flibe, 'HCPB': make_hcpb, 'DCLL': make_dcll}


# =====================================================================
# Core computation
# =====================================================================

def run_model(method, case_name, make_mat, loadings, out_csv,
              flux_csv=None, flux_loading=0.0, fertile='U238', fissile='U235',
              blanket='FLiBe'):
    """Compute resonance integral for a given approximation method.

    method: 'NR', 'NRM', 'IR', or 'IRM'
      - NR/IR controls self-shielding (narrow resonance vs intermediate resonance)
      - M suffix enables subcritical multiplication factor
    """
    use_ir = method in ('IR', 'IRM')
    use_mult = method in ('NRM', 'IRM')

    library = openmc.data.DataLibrary.from_xml()

    max_loading = max(l for l in loadings if l > 0)
    all_nuclides = sorted(make_mat(max_loading).get_nuclide_atom_densities().keys())
    print(f"\n{case_name} [{method}]: {len(all_nuclides)} nuclides: {', '.join(all_nuclides)}")

    data = {}
    for n in all_nuclides:
        try:
            data[n] = openmc.data.IncidentNeutron.from_hdf5(str(hdf5_path(library, n)))
        except RuntimeError:
            print(f"  Warning: no cross-section data for {n}, skipping")
    labels = {n: temp_label(data[n].temperatures) for n in data}

    # =================================================================
    # RESONANCE GROUP (E_MIN -- E_C) on fertile nuclide's fine ENDF grid
    # =================================================================
    energy_coarse, phi_u_coarse = read_flux_generic(flux_csv, flux_loading)
    order = np.argsort(energy_coarse)
    energy_coarse, phi_u_coarse = energy_coarse[order], phi_u_coarse[order]

    phi_u_density_coarse = coarse_to_density(energy_coarse, phi_u_coarse)

    energy_fine_full = data[fertile].energy[labels[fertile]]
    fine_mask = (energy_fine_full >= E_MIN[fertile]) & (energy_fine_full <= E_C[fertile])
    energy = energy_fine_full[fine_mask]

    pos = phi_u_density_coarse > 0
    log_density_fine = np.interp(np.log(energy), np.log(energy_coarse[pos]),
                                 np.log(phi_u_density_coarse[pos]))
    phi_u_density = np.exp(log_density_fine)

    # Fertile cross sections on resonance grid
    sigma_tot_fert = eval_xs(data[fertile], labels[fertile], 1, energy)
    sigma_gamma_fert = eval_xs(data[fertile], labels[fertile], 102, energy)
    sigma_elastic_fert = eval_xs(data[fertile], labels[fertile], 2, energy)
    sigma_fis_fert = safe_eval_xs(data[fertile], labels[fertile], 18, energy)
    sigma_abs_fert = sigma_gamma_fert + sigma_fis_fert
    sigma_n2n_fert = safe_eval_xs(data[fertile], labels[fertile], 16, energy)
    sigma_n3n_fert = safe_eval_xs(data[fertile], labels[fertile], 17, energy)

    # Background xs on resonance grid
    total_xs_bg = {}
    elastic_xs_bg = {}
    inelastic_xs_bg = {}
    for n in data:
        if n != fertile:
            total_xs_bg[n] = eval_xs(data[n], labels[n], 1, energy)
            try:
                elastic_xs_bg[n] = eval_xs(data[n], labels[n], 2, energy)
            except KeyError:
                elastic_xs_bg[n] = total_xs_bg[n]
            inelastic_xs_bg[n] = safe_eval_xs(data[n], labels[n], 4, energy)

    # Fissile companion fission on resonance grid
    sigma_fis_fissile = np.zeros_like(energy)
    if fissile and fissile in data:
        sigma_fis_fissile = safe_eval_xs(data[fissile], labels[fissile], 18, energy)

    nu_fert_res = nu_bar(fertile, energy)
    nu_fissile_res = nu_bar(fissile, energy) if fissile else np.zeros_like(energy)

    Phi_u = np.trapezoid(phi_u_density, energy)
    res_mask = (energy_coarse >= E_MIN[fertile]) & (energy_coarse <= E_C[fertile])
    coarse_res_sum = np.sum(phi_u_coarse[res_mask])
    ratio = Phi_u / coarse_res_sum if coarse_res_sum > 0 else float('inf')
    print(f"  Phi_u (intermediate region): fine-grid trapz = {Phi_u:.6g}, "
          f"coarse-tally sum = {coarse_res_sum:.6g} (ratio {ratio:.4f})")
    if abs(ratio - 1.0) > 0.05:
        print(f"  WARNING: fine/coarse Phi_u ratio off by >5%")

    i_infinite = np.trapezoid(sigma_gamma_fert * phi_u_density, energy) / Phi_u

    # IR-specific constants
    if use_ir:
        A_fert = int(re.search(r'\d+', fertile).group())
        lambda_fert = LAMBDA[A_fert]
        sigma_p_fert = SIGMA_P[fertile]

    # =================================================================
    # FAST GROUP (E_C -- E_MAX) on coarse tally grid
    # =================================================================
    energy_fast_c, phi_u_fast_c = read_flux_generic(
        flux_csv, flux_loading, e_min=E_C[fertile], e_max=E_MAX)
    has_fast = len(energy_fast_c) > 0

    if has_fast:
        fo = np.argsort(energy_fast_c)
        energy_fast_c, phi_u_fast_c = energy_fast_c[fo], phi_u_fast_c[fo]
        phi_u_density_fast = coarse_to_density(energy_fast_c, phi_u_fast_c)

        sigma_gamma_fert_fast = eval_xs(data[fertile], labels[fertile], 102, energy_fast_c)
        sigma_fis_fert_fast = safe_eval_xs(data[fertile], labels[fertile], 18, energy_fast_c)
        sigma_abs_fert_fast = sigma_gamma_fert_fast + sigma_fis_fert_fast

        sigma_fis_fissile_fast = np.zeros_like(energy_fast_c)
        if fissile and fissile in data:
            sigma_fis_fissile_fast = safe_eval_xs(data[fissile], labels[fissile], 18, energy_fast_c)

        nu_fert_fast = nu_bar(fertile, energy_fast_c)
        nu_fissile_fast = nu_bar(fissile, energy_fast_c) if fissile else np.zeros_like(energy_fast_c)

        sigma_n2n_fert_fast = safe_eval_xs(data[fertile], labels[fertile], 16, energy_fast_c)
        sigma_n3n_fert_fast = safe_eval_xs(data[fertile], labels[fertile], 17, energy_fast_c)

        total_xs_bg_fast = {}
        elastic_xs_bg_fast = {}
        inelastic_xs_bg_fast = {}
        for n in data:
            if n != fertile:
                total_xs_bg_fast[n] = eval_xs(data[n], labels[n], 1, energy_fast_c)
                try:
                    elastic_xs_bg_fast[n] = eval_xs(data[n], labels[n], 2, energy_fast_c)
                except KeyError:
                    elastic_xs_bg_fast[n] = total_xs_bg_fast[n]
                inelastic_xs_bg_fast[n] = safe_eval_xs(data[n], labels[n], 4, energy_fast_c)

        R_F_per_nf = np.trapezoid(sigma_gamma_fert_fast * phi_u_density_fast, energy_fast_c)
        Phi_u_fast = np.trapezoid(phi_u_density_fast, energy_fast_c)
        fast_mask = (energy_coarse >= E_C[fertile]) & (energy_coarse <= E_MAX)
        coarse_fast_sum = np.sum(phi_u_coarse[fast_mask])
        ratio_fast = Phi_u_fast / coarse_fast_sum if coarse_fast_sum > 0 else float('inf')
        print(f"  Phi_u (fast region): trapz = {Phi_u_fast:.6g}, "
              f"coarse-tally sum = {coarse_fast_sum:.6g} (ratio {ratio_fast:.4f})")
        if abs(ratio_fast - 1.0) > 0.05:
            print(f"  WARNING: fast fine/coarse Phi_u ratio off by >5%")
        print(f"  fast group ({E_C[fertile]:.0f}-{E_MAX:.0e} eV): "
              f"R_F/n_fert = {R_F_per_nf:.6g},  {len(energy_fast_c)} tally bins")
    else:
        R_F_per_nf = 0.0

    # =================================================================
    # LOADING LOOP
    # =================================================================
    rows = []
    for loading in loadings:
        material = make_mat(loading)
        atom_density = material.get_nuclide_atom_densities()
        n_fert = atom_density.get(fertile, 0.0)
        n_fiss = atom_density.get(fissile, 0.0) if fissile else 0.0

        if n_fert > 0.0:
            # --- self-shielding (resonance group) ---
            if use_ir:
                sigma_bg_lambda = np.zeros_like(energy)
                for nuclide, density in atom_density.items():
                    if nuclide != fertile and nuclide in elastic_xs_bg:
                        A = int(re.search(r'\d+', nuclide).group())
                        lam = LAMBDA[A]
                        sigma_abs_bg = np.maximum(
                            total_xs_bg[nuclide] - elastic_xs_bg[nuclide], 0.0)
                        sigma_bg_lambda += density * (sigma_abs_bg + lam * elastic_xs_bg[nuclide])
                psi_0 = sigma_bg_lambda / n_fert
                f_shield = ((lambda_fert * sigma_p_fert + psi_0)
                            / (sigma_abs_fert + lambda_fert * sigma_elastic_fert + psi_0))
            else:
                sigma_bg_tot = np.zeros_like(energy)
                for nuclide, density in atom_density.items():
                    if nuclide != fertile and nuclide in total_xs_bg:
                        sigma_bg_tot += density * total_xs_bg[nuclide]
                psi_0 = sigma_bg_tot / n_fert
                f_shield = psi_0 / (sigma_tot_fert + psi_0)

            i_eff = np.trapezoid(sigma_gamma_fert * f_shield * phi_u_density, energy) / Phi_u
            psi = np.trapezoid(psi_0 * phi_u_density, energy) / Phi_u
            R_I = n_fert * i_eff * Phi_u
            R_F = n_fert * R_F_per_nf

            # --- subcritical multiplication (both groups) ---
            if use_mult:
                phi_shielded = f_shield * phi_u_density

                prod_res = np.trapezoid(
                    (nu_fert_res * n_fert * sigma_fis_fert
                     + nu_fissile_res * n_fiss * sigma_fis_fissile
                     + 2 * n_fert * sigma_n2n_fert
                     + 3 * n_fert * sigma_n3n_fert)
                    * phi_shielded, energy)

                abs_res = np.trapezoid(
                    n_fert * (sigma_abs_fert + sigma_n2n_fert + sigma_n3n_fert)
                    * phi_shielded, energy)
                for nuclide, density in atom_density.items():
                    if nuclide != fertile and nuclide in total_xs_bg:
                        sigma_abs_bg = (total_xs_bg[nuclide] - elastic_xs_bg[nuclide]
                                        - inelastic_xs_bg.get(nuclide, 0.0))
                        abs_res += np.trapezoid(
                            density * np.maximum(sigma_abs_bg, 0.0) * phi_shielded,
                            energy)

                prod_fast = 0.0
                abs_fast = 0.0
                if has_fast:
                    prod_fast = np.trapezoid(
                        (nu_fert_fast * n_fert * sigma_fis_fert_fast
                         + nu_fissile_fast * n_fiss * sigma_fis_fissile_fast
                         + 2 * n_fert * sigma_n2n_fert_fast
                         + 3 * n_fert * sigma_n3n_fert_fast)
                        * phi_u_density_fast, energy_fast_c)

                    abs_fast = np.trapezoid(
                        n_fert * (sigma_abs_fert_fast + sigma_n2n_fert_fast
                                  + sigma_n3n_fert_fast)
                        * phi_u_density_fast, energy_fast_c)
                    for nuclide, density in atom_density.items():
                        if nuclide != fertile and nuclide in total_xs_bg_fast:
                            sigma_abs_bg_f = (total_xs_bg_fast[nuclide]
                                              - elastic_xs_bg_fast[nuclide]
                                              - inelastic_xs_bg_fast.get(nuclide, 0.0))
                            abs_fast += np.trapezoid(
                                density * np.maximum(sigma_abs_bg_f, 0.0)
                                * phi_u_density_fast, energy_fast_c)

                production = prod_res + prod_fast
                absorption = abs_res + abs_fast
                k_inf = production / absorption if absorption > 0 else 0.0
                k_eff = (1.0 - LEAKAGE[blanket]) * k_inf
                if k_eff >= 1.0:
                    print(f"  WARNING: k_eff={k_eff:.4f} >= 1 at {loading} kg/m3")
                M = 1.0 / (1.0 - k_eff) if k_eff < 1.0 else np.inf
            else:
                k_inf = 0.0
                k_eff = 0.0
                M = 1.0

        else:
            psi = np.inf
            i_eff = i_infinite
            k_inf = 0.0
            k_eff = 0.0
            M = 1.0
            R_I = 0.0
            R_F = 0.0

        rows.append({
            "kg_fertile_per_m3_breeder": loading,
            "N_fertile_atom_per_b_cm": n_fert,
            "sigma_bg": psi,
            "I_infinite_b": i_infinite,
            "I_eff_b": i_eff,
            "I_eff_over_I_infinite": i_eff / i_infinite,
            "R_I": R_I,
            "R_F": R_F,
            "R_per_source_neutron": M * (R_I + R_F),
            "k_inf": k_inf,
            "k_eff": k_eff,
            "M_subcrit": M,
        })

    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out_csv}")

    return rows


# =====================================================================
# I/O helpers
# =====================================================================

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


def _parse_stats(filename):
    """Extract total statistics (particles * cycles) from a folder name like '...4e4x25'."""
    suffix = filename.rstrip("/").split("_")[-1]
    parts = suffix.lower().split("x")
    return float(parts[0]) * int(parts[1])


def read_flux_generic(flux_csv, loading=0.0, e_min=0.0, e_max=E_MAX):
    candidates = {}
    with flux_csv.open() as f:
        for row in csv.DictReader(f):
            if float(row["fertile_kg/m3"]) == loading:
                fn = row["filename"]
                if fn not in candidates:
                    candidates[fn] = {"stats": _parse_stats(fn), "energy": [], "flux": []}
                e = float(row["energy mid [eV]"])
                if e_min <= e <= e_max:
                    candidates[fn]["energy"].append(e)
                    candidates[fn]["flux"].append(float(row["mean"]))
    if not candidates:
        return np.array([]), np.array([])
    best = max(candidates, key=lambda fn: candidates[fn]["stats"])
    return np.array(candidates[best]["energy"]), np.array(candidates[best]["flux"])


def filter_rows(rows, loadings):
    loading_set = set(loadings)
    out = [r for r in rows if r["kg_fertile_per_m3_breeder"] in loading_set]
    missing = loading_set - {r["kg_fertile_per_m3_breeder"] for r in out}
    if missing:
        raise ValueError(f"loadings {sorted(missing)} not found in cached CSV")
    return out


# =====================================================================
# Plotting
# =====================================================================

def fit_power_law(x, y):
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    c, log_a = np.polyfit(np.log(x[mask]), np.log(y[mask]), 1)
    return np.exp(log_a), c


def fit_power_law_progressive(x, y, r2_min=R2_MIN):
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

    n_used = min(2, len(x))
    a, c = fit_power_law(x[:n_used], y[:n_used])

    for k in range(3, len(x) + 1):
        a_k, c_k = fit_power_law(x[:k], y[:k])
        if r_squared(x[:k], y[:k], a_k, c_k) < r2_min:
            break
        a, c, n_used = a_k, c_k, k

    r2 = r_squared(x[:n_used], y[:n_used], a, c)
    assert n_used <= 2 or r2 >= r2_min
    return a, c, r2, n_used, x[n_used - 1]


def plot_compare_all(cases, scale='linear', suffix=''):
    """All methods vs OpenMC for multiple blankets of one fertile isotope.

    cases: list of dicts with keys:
        'label': blanket display name (e.g. 'FLiBe')
        'x_key': row column for x-axis
        'openmc_csv': path to OpenMC summary CSV
        'ng_col': column name for (n,g) rate
        'rows': dict {method: [row_dicts]}
    """
    fig, ax = plt.subplots(figsize=(4.5, 3.5))

    if scale in ('linlog', 'loglog'):
        ax.set_xscale("log")
        if suffix.endswith('_full'):
            ax.set_xlim(10**(-1 - 0.03*5), 10**(4 + 0.03*5))
        else:
            ax.set_xlim(10**(-1 - 0.03*4), 10**(3 + 0.03*4))
        if scale == 'loglog':
            ax.set_yscale("log")
        out = OUT_COMPARE_LOGLOG if scale == 'loglog' else OUT_COMPARE_LINLOG
    else:
        max_load = max(r[cases[0]['x_key']] for c in cases
                       for m in METHODS for r in c['rows'][m])
        ax.set_xlim(-0.03 * max_load, max_load + 0.03 * max_load)
        out = OUT_COMPARE_LINEAR
    out = Path(str(out) + suffix)

    def mask_for(x, y):
        return (x > 0) & (y > 0) if scale != 'linear' else np.ones_like(x, dtype=bool)

    for case in cases:
        blanket = case['label']
        color = BLANKET_COLORS[blanket]
        mkr = BLANKET_MARKERS[blanket]

        omc = read_openmc_generic(case['openmc_csv'], ng_col=case['ng_col'])
        omc_mask = mask_for(omc["x"], omc["y"])
        ax.plot(omc["x"][omc_mask], omc["y"][omc_mask],
                marker=mkr, color=color, lw=0, markersize=2.5,
                linestyle='none', label='_nolegend_')

        for method in METHODS:
            ls = METHOD_LINESTYLES[method]
            rows = case['rows'][method]
            x_m = np.array([r[case['x_key']] for r in rows])
            y_m = np.array([r["R_per_source_neutron"] for r in rows])
            m = mask_for(x_m, y_m)
            ax.plot(x_m[m], y_m[m], ls=ls, color=color, lw=0.75,
                    markersize=1, label='_nolegend_')

    ax.set_xlabel(r"Fertile isotope density [kg$/$m$^3$]")
    ax.set_ylabel(r"Fertile isotope (n,$\gamma$) rate [rxns$/$src-n]")

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    if scale != 'loglog':
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    if scale == 'linear' and not suffix.endswith('_full'):
        ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
    ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

    # Factored legend: blankets (color+marker) and methods (linestyle)
    for blanket in ['FLiBe', 'HCPB', 'DCLL']:
        ax.plot([], [], marker=BLANKET_MARKERS[blanket], color=BLANKET_COLORS[blanket],
                ls='-', lw=0.75, markersize=2.5, label=blanket)
    for method in METHODS:
        ax.plot([], [], ls=METHOD_LINESTYLES[method], color='0.4', lw=0.75, label=method)
    ax.plot([], [], 'ko', markersize=2.5, lw=0, label='OpenMC')

    if not suffix.endswith('_full'):
        if scale == 'loglog':
            ax.set_ylim(10**(-5 - 0.03*5), 10**(0 + 0.03*5))
        else:
            ax.set_ylim(-0.03 * 0.5, 0.5 + 0.03 * 0.5)

    fig.tight_layout()
    leg = ax.legend(loc='upper left', fancybox=False, edgecolor='none',
                    frameon=True, framealpha=0.0, ncol=2, fontsize=6.5)
    leg.get_frame().set_linewidth(0.5)

    fig.savefig(f'{out}.png', bbox_inches='tight', pad_inches=0.01, dpi=300)
    fig.savefig(f'{out}.pdf', bbox_inches='tight', pad_inches=0.01, dpi=300)
    print(f"wrote {out}")
    plt.close(fig)


def plot_psi_vs_ieff(cases, r2_min=R2_MIN, suffix=''):
    """psi vs I_eff for all methods and blankets (one fertile isotope).

    cases: list of dicts with keys:
        'label': blanket display name
        'rows': dict {method: [row_dicts]}
    """
    fig, ax = plt.subplots(figsize=(4.5, 3.5))

    for case in cases:
        blanket = case['label']
        color = BLANKET_COLORS[blanket]

        for method in METHODS:
            mkr = METHOD_MARKERS[method]
            ls = METHOD_LINESTYLES[method]
            rows = case['rows'][method]

            psi = np.array([r["sigma_bg"] for r in rows])
            i_eff = np.array([r["I_eff_b"] for r in rows])
            mask = np.isfinite(psi) & (psi > 0) & (i_eff > 0)

            ax.plot(psi[mask], i_eff[mask], marker=mkr, color=color, markersize=3,
                    markeredgewidth=0.5, linestyle="none", label='_nolegend_')

            a, c, r2, n_used, psi_boundary = fit_power_law_progressive(
                psi[mask], i_eff[mask], r2_min)
            n_total = int(mask.sum())
            print(f"{blanket} {method}: I_eff = {a:.4g} * psi^{c:.4g}   "
                  f"(R2={r2:.4f} using {n_used}/{n_total} points up to "
                  f"psi={psi_boundary:.3g} b)")

            x_fit = np.geomspace(psi[mask].min(), psi_boundary, 100)
            ax.plot(x_fit, a * x_fit**c, ls=ls, color=color, lw=0.75,
                    label='_nolegend_')

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\psi_o$ [b]")
    ax.set_ylabel(r"$I_\mathrm{eff}$ [b]")
    ax.set_xlim(10**(1 - 0.03*7), 10**(8 + 0.03*7))
    if not suffix.endswith('_full'):
        ax.set_ylim(10**(-1 - 0.03*2), 10**(1 + 0.03*2))

    # Factored legend
    for blanket in ['FLiBe', 'HCPB', 'DCLL']:
        ax.plot([], [], 'o', color=BLANKET_COLORS[blanket], markersize=3, label=blanket)
    for method in METHODS:
        ax.plot([], [], marker=METHOD_MARKERS[method], ls=METHOD_LINESTYLES[method],
                color='0.4', lw=0.75, markersize=3, markeredgewidth=0.5, label=method)

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
    plt.close(fig)


# =====================================================================
# Main
# =====================================================================

BLANKET_NAMES = ['FLiBe', 'HCPB', 'DCLL']
FERTILE_CONFIGS = [('U238', 'U235'), ('Th232', None)]
LEAKAGE = {'FLiBe': 0.05, 'HCPB': 0.11, 'DCLL': 0.09}

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

    os.makedirs(OUTS, exist_ok=True)

    # Run all method x blanket x fertile combos (always with full LOADINGS)
    results = {}
    for method in METHODS:
        for blanket in BLANKET_NAMES:
            for fertile, fissile in FERTILE_CONFIGS:
                key = (method, blanket, fertile)
                csv_path = out_csv_path(method, blanket, fertile)

                case_name = f"{blanket} ({COMPOUNDS[(blanket, fertile)].upper()})"
                make_fn = MAKE_MAT[blanket]
                mat_fn = (lambda l, _fn=make_fn, _fi=fertile: _fn(l, _fi))

                cached = read_cached_rows(csv_path)
                if cached:
                    results[key] = cached
                else:
                    results[key] = run_model(
                        method, case_name, mat_fn, LOADINGS, csv_path,
                        flux_csv=CSV_FLUX[(blanket, fertile)],
                        fertile=fertile, fissile=fissile,
                        blanket=blanket)

    # Filter for plotting
    plot_data = {}
    for key, rows in results.items():
        plot_data[key] = filter_rows(rows, loadings)

    # Plot per fertile isotope
    x_key = "kg_fertile_per_m3_breeder"

    for fertile, _ in FERTILE_CONFIGS:
        ng_col = f'{fertile}(n,g)'
        fert_suffix = f'_{fertile.lower()}'

        compare_cases = []
        psi_cases = []
        for blanket in BLANKET_NAMES:
            method_rows = {m: plot_data[(m, blanket, fertile)] for m in METHODS}
            compare_cases.append({
                'label': blanket,
                'x_key': x_key,
                'openmc_csv': CSV_RXNS[(blanket, fertile)],
                'ng_col': ng_col,
                'rows': method_rows,
            })
            psi_cases.append({
                'label': blanket,
                'rows': method_rows,
            })

        for scale in ('linlog', 'loglog', 'linear'):
            plot_compare_all(compare_cases, scale=scale, suffix=fert_suffix + suffix)

        plot_psi_vs_ieff(psi_cases, r2_min=R2_MIN, suffix=fert_suffix + suffix)


if __name__ == "__main__":
    main()
