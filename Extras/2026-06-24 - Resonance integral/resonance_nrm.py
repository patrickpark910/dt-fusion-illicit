#!/usr/bin/env python3

"""resonance integral-style U-238(n,gamma) test for UF4-FLiBe."""

from pathlib import Path
import os, re, sys, csv
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import openmc

# NB. the / operator for path joining only works when the left side is a Path, and everything chained after that is a Path --ppark
HERE = Path(__file__).resolve().parent
# ROOT = Path("/Users/patrickpark/git-repos/emma-openmc")
ROOT = Path("/home/patri/projects/emma-openmc")
DATA = ROOT / "Figures/Data"
OUTS = HERE / "Outputs"

os.chdir(ROOT)
sys.path.insert(0, str(ROOT))
from Python.parameters import *
from Python.utilities import *

E_MIN = {'U238': 3.0, 'Th232': 10.0}
E_MAX = 15e6
E_C = {'U238': 2e4, 'Th232': 4e3}
R2_MIN = 0.97
SIGMA_P_U238 = 11.29  # barns
NU_BAR = {
    'Th232': (np.array([1e-5, 2.53e-2, 1e6, 3e6, 4e6, 5.7e6, 7e6, 10e6, 14.7e6, 20e6]),
              np.array([2.1047, 2.1047, 2.1817, 2.3177, 2.457, 2.6807, 3.026, 3.4, 4.002, 4.82])),
    'U238':  (np.array([3.42e6, 10.15e6, 13.15e6, 20e6]),
              np.array([2.3117, 3.1395, 3.4969, 4.3394])),
    'U235':  (np.array([3.25e6, 5.29788e6, 9.29788e6, 12.2979e6, 20e6]),
              np.array([2.0523, 2.3247, 2.8567, 3.2517, 4.276])),
}

def nu_bar(isotope, energy_ev):
    e_tab, nu_tab = NU_BAR[isotope]
    return np.interp(energy_ev, e_tab, nu_tab, left=nu_tab[0], right=nu_tab[-1])

BLANKET_COLORS = ['#66b420', '#b41f24', '#0047ba']   # FLiBe, HCPB, DCLL -- shared across all plots


# kg U-238 per m3 of real, post-deduction FLiBe.
LOADINGS = [0.0, 0.10, 0.50, 1, 10, 25, 50, 75, 100, 150, 250, 500, 750, 999.99, 2000, 3000, 4000] # FERTILE_KGM3  # from Python/parameters.py

""" OpenMC comparison data to read from Figures/Data """
# Flux spectrum by energy by blanket (U-238)
CSV_FLUX_FLIBE_U = DATA / "FLiBe_900K_Li07.5_U238_flux.csv"
CSV_FLUX_HCPB_U = DATA / "HCPB_900K_Li60.0_U238_flux.csv"
CSV_FLUX_DCLL_U = DATA / "DCLL_900K_Li90.0_U238_flux.csv"
# Flux spectrum by energy by blanket (Th-232)
CSV_FLUX_FLIBE_TH = DATA / "FLiBe_900K_Li07.5_Th232_flux.csv"
CSV_FLUX_HCPB_TH = DATA / "HCPB_900K_Li60.0_Th232_flux.csv"
CSV_FLUX_DCLL_TH = DATA / "DCLL_900K_Li90.0_Th232_flux.csv"

# Reaction rates summary by blanket (U-238)
CSV_RXNS_FLIBE_U = DATA / "FLiBe_900K_Li07.5_U238_summary.csv"
CSV_RXNS_HCPB_U = DATA / "HCPB_900K_Li60.0_U238_summary.csv"
CSV_RXNS_DCLL_U = DATA / "DCLL_900K_Li90.0_U238_summary.csv"
# Reaction rates summary by blanket (Th-232)
CSV_RXNS_FLIBE_TH = DATA / "FLiBe_900K_Li07.5_Th232_summary.csv"
CSV_RXNS_HCPB_TH = DATA / "HCPB_900K_Li60.0_Th232_summary.csv"
CSV_RXNS_DCLL_TH = DATA / "DCLL_900K_Li90.0_Th232_summary.csv"


""" Names of various output files"""
OUT_COMPARE_LINLOG = HERE / "integral_vs_openmc_linlog"
OUT_COMPARE_LOGLOG = HERE / "integral_vs_openmc_loglog"
OUT_COMPARE_LINEAR = HERE / "integral_vs_openmc_linlin"
OUT_PSI_IEFF = HERE / "integral_psi_vs_ieff"

OUT_CSV_FLIBE_U_NRM = OUTS / "integral_nrm_flibe_uf4.csv"
OUT_CSV_DCLL_U_NRM  = OUTS / "integral_nrm_dcll_uo2.csv"
OUT_CSV_HCPB_U_NRM  = OUTS / "integral_nrm_hcpb_uo2.csv"
OUT_CSV_FLIBE_TH_NRM = OUTS / "integral_nrm_flibe_thf4.csv"
OUT_CSV_DCLL_TH_NRM  = OUTS / "integral_nrm_dcll_tho2.csv"
OUT_CSV_HCPB_TH_NRM  = OUTS / "integral_nrm_hcpb_tho2.csv"


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


def uranium_mix():
    w25 = ENRICH_U / 100.0
    n25 = w25 / AMU_U235
    n28 = (1.0 - w25) / AMU_U238
    x25 = n25 / (n25 + n28)
    x28 = 1.0 - x25
    amu_u = x25 * AMU_U235 + x28 * AMU_U238
    uf4_per_u238 = (amu_u + 4.0 * AMU_F19) / (x28 * AMU_U238)
    return x25, x28, uf4_per_u238


def make_flibe(kg_u238_per_m3_breeder):
    x25, x28, uf4_per_u238 = uranium_mix()

    flibe = openmc.Material(name="FLiBe", temperature=TEMP_K)
    flibe.set_density("g/cm3", DENSITY_FLIBE)
    flibe.add_elements_from_formula("F4Li2Be", "ao", enrichment_target="Li6", enrichment_type="ao", enrichment=ENRICH_FLIBE )

    uf4 = openmc.Material(name="UF4", temperature=TEMP_K)
    uf4.set_density("g/cm3", DENSITY_UF4)
    uf4.add_nuclide("U235", x25, "ao")
    uf4.add_nuclide("U238", x28, "ao")
    uf4.add_element("F", 4.0, "ao")

    v_uf4_per_v_flibe = kg_u238_per_m3_breeder * uf4_per_u238 / (DENSITY_UF4 * 1000.0)
    v_mix_per_v_flibe = 1.0 + v_uf4_per_v_flibe
    vf_flibe = 1.0 / v_mix_per_v_flibe
    vf_uf4 = v_uf4_per_v_flibe / v_mix_per_v_flibe

    if vf_uf4 == 0.0:
        return flibe

    mix = openmc.Material.mix_materials([flibe, uf4], [vf_flibe, vf_uf4], "vo")
    mix.temperature = TEMP_K
    return mix


def make_dcll(kg_u238_per_m3_breeder):
    """Build DCLL blanket material matching Python/dcll.py."""
    pbli = openmc.Material(name='Pb-17Li', temperature=TEMP_K)
    pbli.set_density('g/cm3', DENSITY_DCLL)
    pbli.add_element('Pb', 0.83, percent_type='ao')
    pbli.add_element('Li', 0.17, percent_type='ao',
                     enrichment_target='Li6', enrichment_type='ao',
                     enrichment=ENRICH_DCLL)

    kernel = openmc.Material(name='UO2', temperature=TEMP_K)
    kernel.add_elements_from_formula('UO2', enrichment=ENRICH_U)
    kernel.set_density('g/cm3', DENSITY_UO2)

    sic = openmc.Material(name='SiC', temperature=TEMP_K)
    sic.add_elements_from_formula('SiC')
    sic.set_density('g/cm3', 3.2)

    biso = openmc.Material.mix_materials(
        [kernel, sic], [BISO_KERNEL_VOL_FRAC, BISO_COAT_VOL_FRAC], 'vo')

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

    vf_biso_br, vf_pbli_br, _ = calc_biso_breeder_vol_fracs(kg_u238_per_m3_breeder, fertile_isotope='U238')
    vf_biso_bl = vf_biso_br * DCLL_VF_LL_NOM
    vf_pbli_bl = vf_pbli_br * DCLL_VF_LL_NOM

    if kg_u238_per_m3_breeder == 0.0:
        mix = openmc.Material.mix_materials([pbli, f82h, sic, he], [DCLL_VF_LL_NOM, DCLL_VF_FS_NOM, DCLL_VF_SI_NOM, DCLL_VF_HE_NOM], 'vo')
    else:
        mix = openmc.Material.mix_materials([biso, pbli, f82h, sic, he], [vf_biso_bl, vf_pbli_bl, DCLL_VF_FS_NOM, DCLL_VF_SI_NOM, DCLL_VF_HE_NOM], 'vo')
    mix.temperature = TEMP_K
    return mix


def make_hcpb(kg_u238_per_m3_breeder):
    """Build HCPB blanket material matching Python/hcpb.py."""
    li4sio4 = openmc.Material(name='Li4SiO4', temperature=TEMP_K)
    li4sio4.set_density('g/cm3', DENSITY_LI4SIO4)
    li4sio4.add_elements_from_formula(
        'Li4SiO4', enrichment_target='Li6', enrichment_type='ao',
        enrichment=ENRICH_HCPB)

    be_mat = openmc.Material(name='Be', temperature=TEMP_K)
    be_mat.set_density('g/cm3', DENSITY_BE)
    be_mat.add_element('Be', 1, percent_type='wo')

    kernel = openmc.Material(name='UO2', temperature=TEMP_K)
    kernel.add_elements_from_formula('UO2', enrichment=ENRICH_U)
    kernel.set_density('g/cm3', DENSITY_UO2)

    sic = openmc.Material(name='SiC', temperature=TEMP_K)
    sic.add_elements_from_formula('SiC')
    sic.set_density('g/cm3', 3.2)

    biso = openmc.Material.mix_materials(
        [kernel, sic], [BISO_KERNEL_VOL_FRAC, BISO_COAT_VOL_FRAC], 'vo')

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

    vf_biso_br, vf_libe_br, _ = calc_biso_breeder_vol_fracs(
        kg_u238_per_m3_breeder, fertile_isotope='U238')
    vf_biso_bl = vf_biso_br * (HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)
    vf_libe_bl = vf_libe_br * (HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)
    vf_li_bl = vf_libe_bl * HCPB_VF_LI_NOM / (HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)
    vf_be_bl = vf_libe_bl * HCPB_VF_BE_NOM / (HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)

    if kg_u238_per_m3_breeder == 0.0:
        mix = openmc.Material.mix_materials(
            [li4sio4, be_mat, eurofer, he],
            [HCPB_VF_LI_NOM, HCPB_VF_BE_NOM, HCPB_VF_EU_NOM, HCPB_VF_HE_NOM],
            'vo')
    else:
        mix = openmc.Material.mix_materials(
            [biso, li4sio4, be_mat, eurofer, he],
            [vf_biso_bl, vf_li_bl, vf_be_bl, HCPB_VF_EU_NOM, HCPB_VF_HE_NOM],
            'vo')
    mix.temperature = TEMP_K
    return mix


def make_flibe_th(kg_th232_per_m3_breeder):
    flibe = openmc.Material(name="FLiBe", temperature=TEMP_K)
    flibe.set_density("g/cm3", DENSITY_FLIBE)
    flibe.add_elements_from_formula("F4Li2Be", "ao", enrichment_target="Li6", enrichment_type="ao", enrichment=ENRICH_FLIBE)

    thf4 = openmc.Material(name="ThF4", temperature=TEMP_K)
    thf4.set_density("g/cm3", DENSITY_ThF4)
    thf4.add_nuclide("Th232", 1.0, "ao")
    thf4.add_element("F", 4.0, "ao")

    xf4_x_ratio = AMU_ThF4 / AMU_Th232
    v_thf4_per_v_flibe = kg_th232_per_m3_breeder * xf4_x_ratio / (DENSITY_ThF4 * 1000.0)
    v_mix_per_v_flibe = 1.0 + v_thf4_per_v_flibe
    vf_flibe = 1.0 / v_mix_per_v_flibe
    vf_thf4 = v_thf4_per_v_flibe / v_mix_per_v_flibe

    if vf_thf4 == 0.0:
        return flibe

    mix = openmc.Material.mix_materials([flibe, thf4], [vf_flibe, vf_thf4], "vo")
    mix.temperature = TEMP_K
    return mix


def make_dcll_th(kg_th232_per_m3_breeder):
    pbli = openmc.Material(name='Pb-17Li', temperature=TEMP_K)
    pbli.set_density('g/cm3', DENSITY_DCLL)
    pbli.add_element('Pb', 0.83, percent_type='ao')
    pbli.add_element('Li', 0.17, percent_type='ao',
                     enrichment_target='Li6', enrichment_type='ao',
                     enrichment=ENRICH_DCLL)

    kernel = openmc.Material(name='ThO2', temperature=TEMP_K)
    kernel.add_elements_from_formula('ThO2')
    kernel.set_density('g/cm3', DENSITY_ThO2)

    sic = openmc.Material(name='SiC', temperature=TEMP_K)
    sic.add_elements_from_formula('SiC')
    sic.set_density('g/cm3', 3.2)

    biso = openmc.Material.mix_materials(
        [kernel, sic], [BISO_KERNEL_VOL_FRAC, BISO_COAT_VOL_FRAC], 'vo')

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

    vf_biso_br, vf_pbli_br, _ = calc_biso_breeder_vol_fracs(kg_th232_per_m3_breeder, fertile_isotope='Th232')
    vf_biso_bl = vf_biso_br * DCLL_VF_LL_NOM
    vf_pbli_bl = vf_pbli_br * DCLL_VF_LL_NOM

    if kg_th232_per_m3_breeder == 0.0:
        mix = openmc.Material.mix_materials([pbli, f82h, sic, he], [DCLL_VF_LL_NOM, DCLL_VF_FS_NOM, DCLL_VF_SI_NOM, DCLL_VF_HE_NOM], 'vo')
    else:
        mix = openmc.Material.mix_materials([biso, pbli, f82h, sic, he], [vf_biso_bl, vf_pbli_bl, DCLL_VF_FS_NOM, DCLL_VF_SI_NOM, DCLL_VF_HE_NOM], 'vo')
    mix.temperature = TEMP_K
    return mix


def make_hcpb_th(kg_th232_per_m3_breeder):
    li4sio4 = openmc.Material(name='Li4SiO4', temperature=TEMP_K)
    li4sio4.set_density('g/cm3', DENSITY_LI4SIO4)
    li4sio4.add_elements_from_formula(
        'Li4SiO4', enrichment_target='Li6', enrichment_type='ao',
        enrichment=ENRICH_HCPB)

    be_mat = openmc.Material(name='Be', temperature=TEMP_K)
    be_mat.set_density('g/cm3', DENSITY_BE)
    be_mat.add_element('Be', 1, percent_type='wo')

    kernel = openmc.Material(name='ThO2', temperature=TEMP_K)
    kernel.add_elements_from_formula('ThO2')
    kernel.set_density('g/cm3', DENSITY_ThO2)

    sic = openmc.Material(name='SiC', temperature=TEMP_K)
    sic.add_elements_from_formula('SiC')
    sic.set_density('g/cm3', 3.2)

    biso = openmc.Material.mix_materials(
        [kernel, sic], [BISO_KERNEL_VOL_FRAC, BISO_COAT_VOL_FRAC], 'vo')

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

    vf_biso_br, vf_libe_br, _ = calc_biso_breeder_vol_fracs(
        kg_th232_per_m3_breeder, fertile_isotope='Th232')
    vf_biso_bl = vf_biso_br * (HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)
    vf_libe_bl = vf_libe_br * (HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)
    vf_li_bl = vf_libe_bl * HCPB_VF_LI_NOM / (HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)
    vf_be_bl = vf_libe_bl * HCPB_VF_BE_NOM / (HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)

    if kg_th232_per_m3_breeder == 0.0:
        mix = openmc.Material.mix_materials(
            [li4sio4, be_mat, eurofer, he],
            [HCPB_VF_LI_NOM, HCPB_VF_BE_NOM, HCPB_VF_EU_NOM, HCPB_VF_HE_NOM],
            'vo')
    else:
        mix = openmc.Material.mix_materials(
            [biso, li4sio4, be_mat, eurofer, he],
            [vf_biso_bl, vf_li_bl, vf_be_bl, HCPB_VF_EU_NOM, HCPB_VF_HE_NOM],
            'vo')
    mix.temperature = TEMP_K
    return mix



def run_nrm(case_name, make_mat, loadings, out_csv,
            flux_csv=None, flux_loading=0.0,
            fertile='U238', fissile='U235'):
    """NR with subcritical multiplication: R_NRM = M * R_NR, M = 1/(1-k_inf).

    Generalized for arbitrary fertile/fissile pair.  Pass fissile=None for
    pure fertile dopants (e.g. ThF4, ThO2) that carry no companion fissile.
    """
    library = openmc.data.DataLibrary.from_xml()

    max_loading = max(l for l in loadings if l > 0)
    all_nuclides = sorted(make_mat(max_loading).get_nuclide_atom_densities().keys())
    print(f"\n{case_name} [NRM]: {len(all_nuclides)} nuclides: {', '.join(all_nuclides)}")

    data = {}
    for n in all_nuclides:
        try:
            data[n] = openmc.data.IncidentNeutron.from_hdf5(
                str(hdf5_path(library, n)))
        except RuntimeError:
            print(f"  Warning: no cross-section data for {n}, skipping")
    labels = {n: temp_label(data[n].temperatures) for n in data}

    # =====================================================================
    # RESONANCE GROUP (E_MIN -- E_C) on fertile nuclide's fine ENDF grid
    # =====================================================================
    energy_coarse, phi_u_coarse = read_flux_generic(flux_csv, flux_loading)
    order = np.argsort(energy_coarse)
    energy_coarse, phi_u_coarse = energy_coarse[order], phi_u_coarse[order]

    log_e = np.log10(energy_coarse)
    log_edges_mid = 0.5 * (log_e[:-1] + log_e[1:])
    log_edge_lo = log_e[0] - (log_edges_mid[0] - log_e[0])
    log_edge_hi = log_e[-1] + (log_e[-1] - log_edges_mid[-1])
    log_edges = np.concatenate([[log_edge_lo], log_edges_mid, [log_edge_hi]])
    bin_width = np.diff(10**log_edges)
    phi_u_density_coarse = phi_u_coarse / bin_width

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
    try:
        sigma_fis_fert = eval_xs(data[fertile], labels[fertile], 18, energy)
    except KeyError:
        sigma_fis_fert = np.zeros_like(energy)
    sigma_a_fert = sigma_gamma_fert + sigma_fis_fert

    try:
        sigma_n2n_fert = eval_xs(data[fertile], labels[fertile], 16, energy)
    except KeyError:
        sigma_n2n_fert = np.zeros_like(energy)
    try:
        sigma_n3n_fert = eval_xs(data[fertile], labels[fertile], 17, energy)
    except KeyError:
        sigma_n3n_fert = np.zeros_like(energy)

    # Background xs on resonance grid (total and elastic, for NR and k_inf)
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
            try:
                inelastic_xs_bg[n] = eval_xs(data[n], labels[n], 4, energy)
            except KeyError:
                inelastic_xs_bg[n] = np.zeros_like(energy)

    # Fissile companion fission on resonance grid
    sigma_fis_fissile = np.zeros_like(energy)
    if fissile and fissile in data:
        try:
            sigma_fis_fissile = eval_xs(data[fissile], labels[fissile], 18, energy)
        except KeyError:
            pass

    nu_fert_res = nu_bar(fertile, energy)
    nu_fissile_res = nu_bar(fissile, energy) if fissile else np.zeros_like(energy)

    Phi_u = np.trapezoid(phi_u_density, energy)
    res_mask = (energy_coarse >= E_MIN[fertile]) & (energy_coarse <= E_C[fertile])
    coarse_res_sum = np.sum(phi_u_coarse[res_mask])
    ratio = Phi_u / coarse_res_sum if coarse_res_sum > 0 else float('inf')
    print(f"  Phi_u (intermediate region): fine-grid trapz = {Phi_u:.6g}, coarse-tally sum = {coarse_res_sum:.6g} (ratio {ratio:.4f})")
    if abs(ratio - 1.0) > 0.05:
        print(f"  WARNING: fine/coarse Phi_u ratio off by >5%")

    i_infinite = np.trapezoid(sigma_gamma_fert * phi_u_density, energy) / Phi_u

    # =====================================================================
    # FAST GROUP (E_C -- E_MAX) on coarse tally grid
    # =====================================================================
    energy_fast_c, phi_u_fast_c = read_flux_generic(
        flux_csv, flux_loading, e_min=E_C[fertile], e_max=E_MAX)
    has_fast = len(energy_fast_c) > 0

    if has_fast:
        fo = np.argsort(energy_fast_c)
        energy_fast_c, phi_u_fast_c = energy_fast_c[fo], phi_u_fast_c[fo]

        log_ef = np.log10(energy_fast_c)
        log_efm = 0.5 * (log_ef[:-1] + log_ef[1:])
        log_efl = log_ef[0] - (log_efm[0] - log_ef[0])
        log_efh = log_ef[-1] + (log_ef[-1] - log_efm[-1])
        log_ef_edges = np.concatenate([[log_efl], log_efm, [log_efh]])
        bw_fast = np.diff(10**log_ef_edges)
        phi_u_density_fast = phi_u_fast_c / bw_fast

        # Fertile xs on fast grid
        sigma_gamma_fert_fast = eval_xs(data[fertile], labels[fertile], 102, energy_fast_c)
        try:
            sigma_fis_fert_fast = eval_xs(data[fertile], labels[fertile], 18, energy_fast_c)
        except KeyError:
            sigma_fis_fert_fast = np.zeros_like(energy_fast_c)
        sigma_a_fert_fast = sigma_gamma_fert_fast + sigma_fis_fert_fast

        # Fissile companion fission on fast grid
        sigma_fis_fissile_fast = np.zeros_like(energy_fast_c)
        if fissile and fissile in data:
            try:
                sigma_fis_fissile_fast = eval_xs(data[fissile], labels[fissile], 18, energy_fast_c)
            except KeyError:
                pass

        nu_fert_fast = nu_bar(fertile, energy_fast_c)
        nu_fissile_fast = nu_bar(fissile, energy_fast_c) if fissile else np.zeros_like(energy_fast_c)

        try:
            sigma_n2n_fert_fast = eval_xs(data[fertile], labels[fertile], 16, energy_fast_c)
        except KeyError:
            sigma_n2n_fert_fast = np.zeros_like(energy_fast_c)
        try:
            sigma_n3n_fert_fast = eval_xs(data[fertile], labels[fertile], 17, energy_fast_c)
        except KeyError:
            sigma_n3n_fert_fast = np.zeros_like(energy_fast_c)

        # Background xs on fast grid (for k_inf removal)
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
                try:
                    inelastic_xs_bg_fast[n] = eval_xs(data[n], labels[n], 4, energy_fast_c)
                except KeyError:
                    inelastic_xs_bg_fast[n] = np.zeros_like(energy_fast_c)

        R_F_per_nf = np.trapezoid(sigma_gamma_fert_fast * phi_u_density_fast, energy_fast_c)
        Phi_u_fast = np.trapezoid(phi_u_density_fast, energy_fast_c)
        fast_mask = (energy_coarse >= E_C[fertile]) & (energy_coarse <= E_MAX)
        coarse_fast_sum = np.sum(phi_u_coarse[fast_mask])
        ratio_fast = Phi_u_fast / coarse_fast_sum if coarse_fast_sum > 0 else float('inf')
        print(f"  Phi_u (fast region): trapz = {Phi_u_fast:.6g}, coarse-tally sum = {coarse_fast_sum:.6g} (ratio {ratio_fast:.4f})")
        if abs(ratio_fast - 1.0) > 0.05:
            print(f"  WARNING: fast fine/coarse Phi_u ratio off by >5%")
        print(f"  fast group ({E_C[fertile]:.0f}–{E_MAX:.0e} eV): "
              f"R_F/n_fert = {R_F_per_nf:.6g},  {len(energy_fast_c)} tally bins")
    else:
        R_F_per_nf = 0.0

    # =====================================================================
    # LOADING LOOP
    # =====================================================================
    rows = []
    for loading in loadings:
        # print(loading)
        material = make_mat(loading)
        atom_density = material.get_nuclide_atom_densities()
        n_fert = atom_density.get(fertile, 0.0)
        n_fiss = atom_density.get(fissile, 0.0) if fissile else 0.0

        if n_fert > 0.0:
            # --- NR self-shielding (resonance group) ---
            sigma_bg_tot = np.zeros_like(energy)
            for nuclide, density in atom_density.items():
                if nuclide != fertile and nuclide in total_xs_bg:
                    sigma_bg_tot += density * total_xs_bg[nuclide]
            psi_0 = sigma_bg_tot / n_fert
            f_shield = psi_0 / (sigma_tot_fert + psi_0)

            i_eff = np.trapezoid(sigma_gamma_fert * f_shield * phi_u_density, energy) / Phi_u
            psi = np.trapezoid(psi_0 * phi_u_density, energy) / Phi_u
            R_I_base = n_fert * i_eff * Phi_u
            R_F_base = n_fert * R_F_per_nf

            # --- subcritical multiplication (both groups) ---
            phi_shielded = f_shield * phi_u_density

            # Resonance-group production & absorption
            prod_res = np.trapezoid(
                (nu_fert_res * n_fert * sigma_fis_fert
                 + nu_fissile_res * n_fiss * sigma_fis_fissile
                 + 2 * n_fert * sigma_n2n_fert
                 + 3 * n_fert * sigma_n3n_fert)
                * phi_shielded, energy)

            abs_res = np.trapezoid(
                n_fert * (sigma_a_fert + sigma_n2n_fert + sigma_n3n_fert) * phi_shielded, energy)
            for nuclide, density in atom_density.items():
                if nuclide != fertile and nuclide in total_xs_bg:
                    sigma_abs_bg = (total_xs_bg[nuclide] - elastic_xs_bg[nuclide]
                                   - inelastic_xs_bg.get(nuclide, 0.0))
                    abs_res += np.trapezoid(
                        density * np.maximum(sigma_abs_bg, 0.0) * phi_shielded,
                        energy)

            # Fast-group production & absorption (unshielded flux)
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
                    n_fert * (sigma_a_fert_fast + sigma_n2n_fert_fast + sigma_n3n_fert_fast) * phi_u_density_fast, energy_fast_c)
                for nuclide, density in atom_density.items():
                    if nuclide != fertile and nuclide in total_xs_bg_fast:
                        sigma_abs_bg_f = (total_xs_bg_fast[nuclide] - elastic_xs_bg_fast[nuclide]
                                         - inelastic_xs_bg_fast.get(nuclide, 0.0))
                        abs_fast += np.trapezoid(
                            density * np.maximum(sigma_abs_bg_f, 0.0) * phi_u_density_fast,
                            energy_fast_c)

            production = prod_res + prod_fast
            absorption = abs_res + abs_fast
            k_inf = production / absorption if absorption > 0 else 0.0
            if k_inf >= 1.0:
                print(f"  WARNING: k_inf={k_inf:.4f} >= 1 at {loading} kg/m3")
            M = 1.0 / (1.0 - k_inf) if k_inf < 1.0 else np.inf

            R_I = M * R_I_base
            R_F = M * R_F_base
        else:
            psi = np.inf
            i_eff = i_infinite
            k_inf = 0.0
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
            "R_per_source_neutron": R_I + R_F,
            "k_inf": k_inf,
            "M_subcrit": M,
        })

    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out_csv}")

    return rows



def read_cached_rows(csv_path):
    if not csv_path.exists():
        return None
    print(f"  cached: {csv_path}")
    with csv_path.open() as f:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]
    # normalize legacy column name
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


def read_flux_generic(flux_csv, loading=0.0, e_min=0.0, e_max=E_MAX):
    energy, flux = [], []
    with flux_csv.open() as f:
        for row in csv.DictReader(f):
            if float(row["fertile_kg/m3"]) == loading:
                e = float(row["energy mid [eV]"])
                if e_min <= e <= e_max:
                    energy.append(e)
                    flux.append(float(row["mean"]))
    return np.array(energy), np.array(flux)


def plot_compare_all(cases, scale='linear'):
    """Combined NRM and OpenMC for all blankets.

    scale: 'linear', 'linlog', or 'loglog'
    cases: list of (label, nrm_rows, x_key, openmc_csv_path, ng_col, color, linestyle, marker, fillstyle)
    """
    fig, ax = plt.subplots(figsize=(3.5, 3.0))

    if scale == 'linlog':
        ax.set_xscale("log")
        ax.set_yscale("linear")
        ax.set_xlim(10**(-1 - 0.03*5), 10**(4 + 0.03*5))
        out = OUT_COMPARE_LINLOG
    elif scale == 'loglog':
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(10**(-1 - 0.03*5), 10**(4 + 0.03*5))
        out = OUT_COMPARE_LOGLOG
    else:
        max_loading = max(LOADINGS)
        ax.set_xlim(-0.03 * max_loading, max_loading + 0.03 * max_loading)
        out = OUT_COMPARE_LINEAR

    def mask_for(x, y):
        return (x > 0) & (y > 0) if scale != 'linear' else np.ones_like(x, dtype=bool)

    dummy_cases = []
    for label, nrm_rows, x_key, openmc_csv, ng_col, color, linestyle, mkr, fillstyle in cases:
        omc = read_openmc_generic(openmc_csv, ng_col=ng_col)
        omc_mask = mask_for(omc["x"], omc["y"])
        ax.plot(omc["x"][omc_mask], omc["y"][omc_mask],
                marker=mkr, color=color, lw=0.75, markersize=2,
                fillstyle=fillstyle, linestyle='none',
                label='_nolegend_')

        x_m = np.array([r[x_key] for r in nrm_rows])
        y_m = np.array([r["R_per_source_neutron"] for r in nrm_rows])
        mask = mask_for(x_m, y_m)
        ax.plot(x_m[mask], y_m[mask], ls=linestyle, color=color, lw=0.75,
                markersize=1, label='_nolegend_')

        dummy_cases.append((label, color, linestyle, mkr, fillstyle))

    ax.set_xlabel(r"Fertile isotope density [kg$/$m³]")
    ax.set_ylabel(r"Fertile isotope (n,$\gamma$) rate [rxns$/$src-n]")
    ax.set_ylim(-0.03, 1 + 0.03)

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
    ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

    # Dummy plots for legend — Th (left col) then U (right col)
    th = [d for d in dummy_cases if d[4] == 'none']
    u  = [d for d in dummy_cases if d[4] == 'full']
    for label, color, linestyle, mkr, fillstyle in th + u:
        ax.plot([9e8, 9e9], [9e8, 9e9], marker=mkr, color=color, markersize=2,
                fillstyle=fillstyle, linestyle=linestyle, linewidth=0.75,
                label=label)

    fig.tight_layout()
    leg = ax.legend(loc='lower right', fancybox=False, edgecolor='none',
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


def plot_psi_vs_ieff(cases, r2_min=R2_MIN):
    """psi vs I_eff for all blankets, log-log.

    cases: list of (label, nrm_rows, color, marker, fillstyle, fit_ls)
    """
    fig, ax = plt.subplots(figsize=(3.5, 3.0))

    dummy_fits = []
    for label, nrm_rows, color, marker, fillstyle, fit_ls in cases:
        psi = np.array([r["sigma_bg"] for r in nrm_rows])
        i_eff = np.array([r["I_eff_b"] for r in nrm_rows])
        mask = np.isfinite(psi) & (psi > 0) & (i_eff > 0)

        ax.plot(psi[mask], i_eff[mask], marker, color=color, markersize=3,
                fillstyle=fillstyle, linestyle="none", label='_nolegend_')

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
    ax.set_xlim(10**(1 - 0.03*6), 10**(7 + 0.03*6))
    ax.set_ylim(10**(-1 - 0.03*2), 10**(1 + 0.03*2))

    for label, a, c, color, mkr, fillstyle, fit_ls in th + u:
        ax.plot([9e8, 9e9], [9e8, 9e9], marker=mkr, color=color, markersize=3,
                fillstyle=fillstyle, linestyle=fit_ls, linewidth=0.75,
                label=rf"{label}: ${a:.3g}\,\psi_o^{{{c:.3g}}}$")

    from matplotlib.ticker import ScalarFormatter
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_major_formatter(lambda x, _: f'{x:g}')
    ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
    ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

    fig.tight_layout()
    leg = ax.legend(loc='lower right', fancybox=False, edgecolor='none',
                    frameon=True, framealpha=0.0, ncol=2, fontsize=6.5)
    leg.get_frame().set_linewidth(0.5)

    fig.savefig(f'{OUT_PSI_IEFF}.png', bbox_inches='tight', pad_inches=0.01, dpi=300)
    fig.savefig(f'{OUT_PSI_IEFF}.pdf', bbox_inches='tight', pad_inches=0.01, dpi=300)
    print(f"wrote {OUT_PSI_IEFF}")


def main():

    # U-238 cases
    flibe_u_nrm = read_cached_rows(OUT_CSV_FLIBE_U_NRM) or run_nrm("FLiBe (UF4)", make_flibe, LOADINGS, OUT_CSV_FLIBE_U_NRM, flux_csv=CSV_FLUX_FLIBE_U, fertile='U238', fissile='U235')
    hcpb_u_nrm  = read_cached_rows(OUT_CSV_HCPB_U_NRM)  or run_nrm("HCPB (UO2-SiC)", make_hcpb, LOADINGS, OUT_CSV_HCPB_U_NRM, flux_csv=CSV_FLUX_HCPB_U, fertile='U238', fissile='U235')
    dcll_u_nrm  = read_cached_rows(OUT_CSV_DCLL_U_NRM)  or run_nrm("DCLL (UO2-SiC)", make_dcll, LOADINGS, OUT_CSV_DCLL_U_NRM, flux_csv=CSV_FLUX_DCLL_U, fertile='U238', fissile='U235')

    # Th-232 cases
    flibe_th_nrm = read_cached_rows(OUT_CSV_FLIBE_TH_NRM) or run_nrm("FLiBe (ThF4)", make_flibe_th, LOADINGS, OUT_CSV_FLIBE_TH_NRM, flux_csv=CSV_FLUX_FLIBE_TH, fertile='Th232', fissile=None)
    hcpb_th_nrm  = read_cached_rows(OUT_CSV_HCPB_TH_NRM)  or run_nrm("HCPB (ThO2-SiC)", make_hcpb_th, LOADINGS, OUT_CSV_HCPB_TH_NRM, flux_csv=CSV_FLUX_HCPB_TH, fertile='Th232', fissile=None)
    dcll_th_nrm  = read_cached_rows(OUT_CSV_DCLL_TH_NRM)  or run_nrm("DCLL (ThO2-SiC)", make_dcll_th, LOADINGS, OUT_CSV_DCLL_TH_NRM, flux_csv=CSV_FLUX_DCLL_TH, fertile='Th232', fissile=None)

    x_key = "kg_fertile_per_m3_breeder"
    compare_cases = [
        (r"FLiBe-UF$_4$",  flibe_u_nrm,  x_key, CSV_RXNS_FLIBE_U,  'U238(n,g)',   '#66b420', '-',        'o', 'full'),
        (r"FLiBe-ThF$_4$", flibe_th_nrm, x_key, CSV_RXNS_FLIBE_TH, 'Th232(n,g)',  '#66b420', '--',  'o', 'none'),
        (r"HCPB-UO$_2$",   hcpb_u_nrm,   x_key, CSV_RXNS_HCPB_U,   'U238(n,g)',   '#b41f24', '-',        's', 'full'),
        (r"HCPB-ThO$_2$",  hcpb_th_nrm,  x_key, CSV_RXNS_HCPB_TH,  'Th232(n,g)',  '#b41f24', '--',  's', 'none'),
        (r"DCLL-UO$_2$",   dcll_u_nrm,   x_key, CSV_RXNS_DCLL_U,   'U238(n,g)',   '#0047ba', '-',        '^', 'full'),
        (r"DCLL-ThO$_2$",  dcll_th_nrm,  x_key, CSV_RXNS_DCLL_TH,  'Th232(n,g)',  '#0047ba', '--',  '^', 'none'),
    ]
    plot_compare_all(compare_cases, scale='linlog')
    plot_compare_all(compare_cases, scale='loglog')
    plot_compare_all(compare_cases, scale='linear')

    psi_cases = [
        (r"FLiBe-UF$_4$",  flibe_u_nrm,  '#66b420', 'o', 'full',  '-'),
        (r"FLiBe-ThF$_4$", flibe_th_nrm, '#66b420', 'o', 'none',  '--'),
        (r"HCPB-UO$_2$",   hcpb_u_nrm,   '#b41f24', 's', 'full',  '-'),
        (r"HCPB-ThO$_2$",  hcpb_th_nrm,  '#b41f24', 's', 'none',  '--'),
        (r"DCLL-UO$_2$",   dcll_u_nrm,   '#0047ba', '^', 'full',  '-'),
        (r"DCLL-ThO$_2$",  dcll_th_nrm,  '#0047ba', '^', 'none',  '--'),
    ]
    plot_psi_vs_ieff(psi_cases, r2_min=R2_MIN)


if __name__ == "__main__":
    main()