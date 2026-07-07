#!/usr/bin/env python3

"""resonance integral-style U-238(n,gamma) test for UF4-FLiBe."""

from pathlib import Path
import os, re, sys, csv
import matplotlib.pyplot as plt
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

E_MIN_EV, E_MAX_EV = 0.01, 2e4
E_FAST_MAX_EV = 15e6
R2_MIN = 0.97
SIGMA_P_U238 = 11.29  # barns
NU_U235 = 2.43
NU_U238 = 2.80

BLANKET_COLORS = ['C0', 'C1', 'C2']   # FLiBe, HCPB, DCLL -- shared across all plots

# kg U-238 per m3 of real, post-deduction FLiBe.
LOADINGS = [0.0, 0.10, 0.50, 1, 10, 25, 50, 75, 100, 150, 250, 500, 750, 999.99, 2000, 3000, 4000] # FERTILE_KGM3  # from Python/parameters.py

""" OpenMC comparison data to read from Figures/Data """
# Flux spectrum by energy by blanket
CSV_FLUX_FLIBE = DATA / "FLiBe_900K_Li07.5_U238_flux.csv"
CSV_FLUX_HCPB = DATA / "HCPB_900K_Li60.0_U238_flux.csv"
CSV_FLUX_DCLL = DATA / "DCLL_900K_Li90.0_U238_flux.csv"

# Reaction rates summary by blanket
CSV_RXNS_FLIBE = DATA / "FLiBe_900K_Li07.5_U238_summary.csv"
CSV_RXNS_HCPB = DATA / "HCPB_900K_Li60.0_U238_summary.csv"
CSV_RXNS_DCLL = DATA / "DCLL_900K_Li90.0_U238_summary.csv"


""" Names of various output files"""
OUT_COMPARE_LOG = HERE / "integral_vs_openmc_log.png"
OUT_COMPARE_LINEAR = HERE / "integral_vs_openmc_linear.png"
OUT_PSI_IEFF = HERE / "integral_psi_vs_ieff.png"

OUT_CSV_FLIBE_NR = OUTS / "integral_nr_flibe_uf4.csv"
OUT_CSV_DCLL_NR  = OUTS / "integral_nr_dcll_uo2.csv"
OUT_CSV_HCPB_NR  = OUTS / "integral_nr_hcpb_uo2.csv"

OUT_CSV_FLIBE_IR = OUTS / "integral_ir_flibe_uf4.csv"
OUT_CSV_DCLL_IR  = OUTS / "integral_ir_dcll_uo2.csv"
OUT_CSV_HCPB_IR  = OUTS / "integral_ir_hcpb_uo2.csv"

OUT_CSV_FLIBE_IRM = OUTS / "integral_irm_flibe_uf4.csv"
OUT_CSV_DCLL_IRM  = OUTS / "integral_irm_dcll_uo2.csv"
OUT_CSV_HCPB_IRM  = OUTS / "integral_irm_hcpb_uo2.csv"

OUT_CSV_FLIBE_NRM = OUTS / "integral_nrm_flibe_uf4.csv"
OUT_CSV_DCLL_NRM  = OUTS / "integral_nrm_dcll_uo2.csv"
OUT_CSV_HCPB_NRM  = OUTS / "integral_nrm_hcpb_uo2.csv"


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


def compute_fast_group(data, labels, flux_csv, flux_loading=0.0):
    """Naive fast-group integral (no self-shielding):

        R_F_per_n28 = int_{E_MAX_EV}^{E_FAST_MAX_EV} sigma_gamma^{28}(E) phi_u(E) dE

    Uses the coarse tally grid directly (cross sections are smooth above
    20 keV, so the ~50 pts/decade tally grid is adequate).  The result is
    loading-independent: R_F = n28 * R_F_per_n28.
    """
    energy_coarse, phi_u_coarse = read_flux_generic(
        flux_csv, flux_loading, e_min=E_MAX_EV, e_max=E_FAST_MAX_EV)
    if len(energy_coarse) == 0:
        return 0.0
    order = np.argsort(energy_coarse)
    energy_coarse = energy_coarse[order]
    phi_u_coarse = phi_u_coarse[order]

    # Reconstruct bin edges and convert bin-integrated tally to density
    log_e = np.log10(energy_coarse)
    log_edges_mid = 0.5 * (log_e[:-1] + log_e[1:])
    log_edge_lo = log_e[0] - (log_edges_mid[0] - log_e[0])
    log_edge_hi = log_e[-1] + (log_e[-1] - log_edges_mid[-1])
    log_edges = np.concatenate([[log_edge_lo], log_edges_mid, [log_edge_hi]])
    bin_width = np.diff(10**log_edges)
    phi_u_density = phi_u_coarse / bin_width

    sigma_gamma = eval_xs(data["U238"], labels["U238"], 102, energy_coarse)

    R_F_per_n28 = np.trapezoid(sigma_gamma * phi_u_density, energy_coarse)
    print(f"  fast group ({E_MAX_EV:.0f}–{E_FAST_MAX_EV:.0e} eV): "
          f"R_F/n28 = {R_F_per_n28:.6g},  {len(energy_coarse)} tally bins")
    return R_F_per_n28


def run_nr(case_name, make_mat, loadings, out_csv,
                   flux_csv=None, flux_loading=0.0):
    """I_eff from the nominal (undoped) OpenMC flux phi_u(E), per the memo:

        I_eff = (1/Phi_u) * int_E sigma_gamma^28(E) * f(E,n28) * phi_u(E) dE
        f(E,n28) = psi_0(E,n28) / (sigma_tot^28(E) + psi_0(E,n28))
        psi_0(E,n28) = (1/n28) * sum_{i != 28} N_i * sigma_tot^i(E)

    phi_u(E) is read from flux_csv as bin-integrated OpenMC tally values on
    the coarse tally grid (energy_filter = logspace_per_decade(..., 50) =>
    ~500 points over 0.01 eV-14 MeV). That grid is far too coarse to resolve
    U-238's resolved resonances (~500 resonances between 1 eV-10 keV alone,
    so well under one tally point per resonance on average) -- but it's also
    the wrong thing to fix by re-tallying finer, since there's no U-238 in
    the *undoped* run this flux comes from: phi_u(E) itself has no resonance
    structure to lose, it's smooth.

    So instead: convert phi_u to a flux density, interpolate that (smooth)
    density onto U-238's own native ENDF grid (already resonance-resolved by
    construction -- every reaction for a nuclide/temperature shares one grid
    in OpenMC's HDF5 format), and do the actual self-shielding integral
    there. Sums-over-tally-bins become proper trapezoidal integrals now that
    this is a real pointwise grid rather than raw tally bins.
    """
    library = openmc.data.DataLibrary.from_xml()

    max_loading = max(l for l in loadings if l > 0)
    all_nuclides = sorted(make_mat(max_loading).get_nuclide_atom_densities().keys())
    print(f"\n{case_name}: {len(all_nuclides)} nuclides: {', '.join(all_nuclides)}")

    data = {}
    for n in all_nuclides:
        try:
            data[n] = openmc.data.IncidentNeutron.from_hdf5(
                str(hdf5_path(library, n)))
        except RuntimeError:
            print(f"  Warning: no cross-section data for {n}, skipping")
    labels = {n: temp_label(data[n].temperatures) for n in data}

    # --- coarse tally grid: phi_u as bin-integrated values, as tallied ---
    energy_coarse, phi_u_coarse = read_flux_generic(flux_csv, flux_loading)
    order = np.argsort(energy_coarse)
    energy_coarse, phi_u_coarse = energy_coarse[order], phi_u_coarse[order]

    # Recover bin edges from bin *midpoints* (log-space midpoint rule) so we
    # can convert bin-integrated tally values into a flux density [per eV].
    # Reconstructed purely from energy_coarse itself, so this is correct
    # regardless of exactly how "energy mid" was computed upstream, and it
    # automatically respects whatever E_MIN_EV/E_MAX_EV window was already
    # applied in read_flux_generic.
    log_e = np.log10(energy_coarse)
    log_edges_mid = 0.5 * (log_e[:-1] + log_e[1:])
    log_edge_lo = log_e[0] - (log_edges_mid[0] - log_e[0])
    log_edge_hi = log_e[-1] + (log_e[-1] - log_edges_mid[-1])
    log_edges = np.concatenate([[log_edge_lo], log_edges_mid, [log_edge_hi]])
    bin_width = np.diff(10**log_edges)
    phi_u_density_coarse = phi_u_coarse / bin_width   # flux per eV, coarse grid

    # --- fine grid: U-238's own native (resonance-resolved) energy grid ---
    energy_fine_full = data["U238"].energy[labels["U238"]]
    fine_mask = (energy_fine_full >= E_MIN_EV) & (energy_fine_full <= E_MAX_EV)
    energy = energy_fine_full[fine_mask]

    # Interpolate the *density* (log-log, since it spans ~10 decades) onto
    # the fine grid -- phi_u never had resonance structure to lose here.
    pos = phi_u_density_coarse > 0
    log_density_fine = np.interp(np.log(energy), np.log(energy_coarse[pos]),
                                  np.log(phi_u_density_coarse[pos]))
    phi_u_density = np.exp(log_density_fine)

    sigma_tot_u238 = eval_xs(data["U238"], labels["U238"], 1, energy)
    sigma_gamma_u238 = eval_xs(data["U238"], labels["U238"], 102, energy)
    total_xs_bg = {n: eval_xs(data[n], labels[n], 1, energy)
                   for n in data if n != "U238"}

    Phi_u = np.trapezoid(phi_u_density, energy)
    ratio = Phi_u / np.sum(phi_u_coarse)
    print(f"  Phi_u: fine-grid trapz = {Phi_u:.6g}, coarse-tally sum = "
          f"{np.sum(phi_u_coarse):.6g} (ratio {ratio:.4f})")
    if abs(ratio - 1.0) > 0.05:
        print(f"  WARNING: fine/coarse Phi_u ratio off by >5% -- check the "
              f"density conversion/interpolation before trusting these results")

    # infinite-dilution limit: n28 -> 0 => psi_0 -> infinity => f -> 1
    i_infinite = np.trapezoid(sigma_gamma_u238 * phi_u_density, energy) / Phi_u

    # --- fast group: naive sigma*phi above E_MAX_EV (no self-shielding) ---
    R_F_per_n28 = compute_fast_group(data, labels, flux_csv, flux_loading)

    rows = []
    for loading in loadings:
        print(loading)
        material = make_mat(loading)
        atom_density = material.get_nuclide_atom_densities()
        n28 = atom_density.get("U238", 0.0)

        if n28 > 0.0:
            sigma_bg_tot = np.zeros_like(energy)
            for nuclide, density in atom_density.items():
                if nuclide != "U238" and nuclide in total_xs_bg:
                    sigma_bg_tot += density * total_xs_bg[nuclide]
            psi_0 = sigma_bg_tot / n28
            f_shield = psi_0 / (sigma_tot_u238 + psi_0)

            i_eff = np.trapezoid(sigma_gamma_u238 * f_shield * phi_u_density, energy) / Phi_u
            psi = np.trapezoid(psi_0 * phi_u_density, energy) / Phi_u   # flux-weighted avg background xs
            R_I = n28 * i_eff * Phi_u
            R_F = n28 * R_F_per_n28
        else:
            psi = np.inf
            i_eff = i_infinite
            R_I = 0.0
            R_F = 0.0

        rows.append({
            "kgU238_per_m3_breeder": loading,
            "N_U238_atom_per_b_cm": n28,
            "sigma_bg": psi,
            "I_infinite_b": i_infinite,
            "I_eff_b": i_eff,
            "I_eff_over_I_infinite": i_eff / i_infinite,
            "R_I": R_I,
            "R_F": R_F,
            "R_per_source_neutron": R_I + R_F,
        })

    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out_csv}")

    return rows


def run_nrm(case_name, make_mat, loadings, out_csv,
            flux_csv=None, flux_loading=0.0):
    """NR with subcritical multiplication: R_NRM = M * R_NR, M = 1/(1-k_inf).

    Like run_irm but with NR self-shielding, and with the k_inf calculation
    extended to cover the fast group (20 keV-15 MeV) where U-238 fast
    fission and (n,Xn) live.  The resonance-only k_inf in run_irm misses
    these channels entirely because sigma_fis^{28} is zero below ~1 MeV.

    nu-bar is evaluated from ENDF data (MT=452) where available, falling
    back to the constants NU_U235 and NU_U238.
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
    # RESONANCE GROUP (E_MIN_EV -- E_MAX_EV) on U-238's fine ENDF grid
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

    energy_fine_full = data["U238"].energy[labels["U238"]]
    fine_mask = (energy_fine_full >= E_MIN_EV) & (energy_fine_full <= E_MAX_EV)
    energy = energy_fine_full[fine_mask]

    pos = phi_u_density_coarse > 0
    log_density_fine = np.interp(np.log(energy), np.log(energy_coarse[pos]),
                                  np.log(phi_u_density_coarse[pos]))
    phi_u_density = np.exp(log_density_fine)

    # U-238 cross sections on resonance grid
    sigma_tot_u238 = eval_xs(data["U238"], labels["U238"], 1, energy)
    sigma_gamma_u238 = eval_xs(data["U238"], labels["U238"], 102, energy)
    sigma_elastic_u238 = eval_xs(data["U238"], labels["U238"], 2, energy)
    try:
        sigma_fis_u238 = eval_xs(data["U238"], labels["U238"], 18, energy)
    except KeyError:
        sigma_fis_u238 = np.zeros_like(energy)
    sigma_a_u238 = sigma_gamma_u238 + sigma_fis_u238

    # U-238 (n,2n) and (n,3n) on resonance grid (zero below thresholds,
    # but included for completeness and consistency with the fast group)
    try:
        sigma_n2n_u238 = eval_xs(data["U238"], labels["U238"], 16, energy)
    except KeyError:
        sigma_n2n_u238 = np.zeros_like(energy)
    try:
        sigma_n3n_u238 = eval_xs(data["U238"], labels["U238"], 17, energy)
    except KeyError:
        sigma_n3n_u238 = np.zeros_like(energy)

    # Background xs on resonance grid (total and elastic, for NR and k_inf)
    total_xs_bg = {}
    elastic_xs_bg = {}
    inelastic_xs_bg = {}
    for n in data:
        if n != "U238":
            total_xs_bg[n] = eval_xs(data[n], labels[n], 1, energy)
            try:
                elastic_xs_bg[n] = eval_xs(data[n], labels[n], 2, energy)
            except KeyError:
                elastic_xs_bg[n] = total_xs_bg[n]
            # inelastic (MT=4) for corrected k_inf removal
            try:
                inelastic_xs_bg[n] = eval_xs(data[n], labels[n], 4, energy)
            except KeyError:
                inelastic_xs_bg[n] = np.zeros_like(energy)

    # U-235 fission on resonance grid
    sigma_fis_u235 = np.zeros_like(energy)
    if "U235" in data:
        try:
            sigma_fis_u235 = eval_xs(data["U235"], labels["U235"], 18, energy)
        except KeyError:
            pass

    # nu-bar on resonance grid (constant is fine here: U-238 fission is
    # zero below ~1 MeV, and nu_25 ~ 2.43 is nearly flat in the resonance range)
    nu_u238_res = NU_U238 * np.ones_like(energy)
    nu_u235_res = NU_U235 * np.ones_like(energy)

    Phi_u = np.trapezoid(phi_u_density, energy)
    ratio = Phi_u / np.sum(phi_u_coarse)
    print(f"  Phi_u (res): fine-grid trapz = {Phi_u:.6g}, coarse-tally sum = "
          f"{np.sum(phi_u_coarse):.6g} (ratio {ratio:.4f})")
    if abs(ratio - 1.0) > 0.05:
        print(f"  WARNING: fine/coarse Phi_u ratio off by >5%")

    i_infinite = np.trapezoid(sigma_gamma_u238 * phi_u_density, energy) / Phi_u

    # =====================================================================
    # FAST GROUP (E_MAX_EV -- E_FAST_MAX_EV) on coarse tally grid
    # =====================================================================
    energy_fast_c, phi_u_fast_c = read_flux_generic(
        flux_csv, flux_loading, e_min=E_MAX_EV, e_max=E_FAST_MAX_EV)
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

        # U-238 xs on fast grid
        sigma_gamma_u238_fast = eval_xs(data["U238"], labels["U238"], 102, energy_fast_c)
        try:
            sigma_fis_u238_fast = eval_xs(data["U238"], labels["U238"], 18, energy_fast_c)
        except KeyError:
            sigma_fis_u238_fast = np.zeros_like(energy_fast_c)
        sigma_a_u238_fast = sigma_gamma_u238_fast + sigma_fis_u238_fast

        # U-235 fission on fast grid
        sigma_fis_u235_fast = np.zeros_like(energy_fast_c)
        if "U235" in data:
            try:
                sigma_fis_u235_fast = eval_xs(data["U235"], labels["U235"], 18, energy_fast_c)
            except KeyError:
                pass

        # nu-bar on fast grid -- try ENDF (MT=452), fall back to constants
        # (nu_28 rises from ~2.5 near threshold to ~4.5 at 14 MeV;
        # the constant 2.80 materially underweights 14-MeV-induced fission)
        try:
            nu_u238_fast = eval_xs(data["U238"], labels["U238"], 452, energy_fast_c)
            print(f"  nu_U238: using ENDF MT=452 (range {nu_u238_fast.min():.2f}–{nu_u238_fast.max():.2f})")
        except (KeyError, Exception):
            nu_u238_fast = NU_U238 * np.ones_like(energy_fast_c)
            print(f"  nu_U238: MT=452 not available, using constant {NU_U238}")
        try:
            nu_u235_fast = eval_xs(data["U235"], labels["U235"], 452, energy_fast_c)
        except (KeyError, Exception):
            nu_u235_fast = NU_U235 * np.ones_like(energy_fast_c)

        # U-238 (n,2n) and (n,3n) on fast grid
        try:
            sigma_n2n_u238_fast = eval_xs(data["U238"], labels["U238"], 16, energy_fast_c)
        except KeyError:
            sigma_n2n_u238_fast = np.zeros_like(energy_fast_c)
        try:
            sigma_n3n_u238_fast = eval_xs(data["U238"], labels["U238"], 17, energy_fast_c)
        except KeyError:
            sigma_n3n_u238_fast = np.zeros_like(energy_fast_c)

        # Background xs on fast grid (for k_inf removal)
        total_xs_bg_fast = {}
        elastic_xs_bg_fast = {}
        inelastic_xs_bg_fast = {}
        for n in data:
            if n != "U238":
                total_xs_bg_fast[n] = eval_xs(data[n], labels[n], 1, energy_fast_c)
                try:
                    elastic_xs_bg_fast[n] = eval_xs(data[n], labels[n], 2, energy_fast_c)
                except KeyError:
                    elastic_xs_bg_fast[n] = total_xs_bg_fast[n]
                try:
                    inelastic_xs_bg_fast[n] = eval_xs(data[n], labels[n], 4, energy_fast_c)
                except KeyError:
                    inelastic_xs_bg_fast[n] = np.zeros_like(energy_fast_c)

        R_F_per_n28 = np.trapezoid(sigma_gamma_u238_fast * phi_u_density_fast, energy_fast_c)
        print(f"  fast group ({E_MAX_EV:.0f}–{E_FAST_MAX_EV:.0e} eV): "
              f"R_F/n28 = {R_F_per_n28:.6g},  {len(energy_fast_c)} tally bins")
    else:
        R_F_per_n28 = 0.0

    # =====================================================================
    # LOADING LOOP
    # =====================================================================
    rows = []
    for loading in loadings:
        print(loading)
        material = make_mat(loading)
        atom_density = material.get_nuclide_atom_densities()
        n28 = atom_density.get("U238", 0.0)
        n25 = atom_density.get("U235", 0.0)

        if n28 > 0.0:
            # --- NR self-shielding (resonance group) ---
            sigma_bg_tot = np.zeros_like(energy)
            for nuclide, density in atom_density.items():
                if nuclide != "U238" and nuclide in total_xs_bg:
                    sigma_bg_tot += density * total_xs_bg[nuclide]
            psi_0 = sigma_bg_tot / n28
            f_shield = psi_0 / (sigma_tot_u238 + psi_0)

            i_eff = np.trapezoid(sigma_gamma_u238 * f_shield * phi_u_density, energy) / Phi_u
            psi = np.trapezoid(psi_0 * phi_u_density, energy) / Phi_u
            R_I_base = n28 * i_eff * Phi_u
            R_F_base = n28 * R_F_per_n28

            # --- subcritical multiplication (both groups) ---
            phi_shielded = f_shield * phi_u_density

            # Resonance-group production & absorption
            prod_res = np.trapezoid(
                (nu_u238_res * n28 * sigma_fis_u238
                 + nu_u235_res * n25 * sigma_fis_u235
                 + 2 * n28 * sigma_n2n_u238
                 + 3 * n28 * sigma_n3n_u238)
                * phi_shielded, energy)

            abs_res = np.trapezoid(
                n28 * (sigma_a_u238 + sigma_n2n_u238 + sigma_n3n_u238) * phi_shielded, energy)
            for nuclide, density in atom_density.items():
                if nuclide != "U238" and nuclide in total_xs_bg:
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
                    (nu_u238_fast * n28 * sigma_fis_u238_fast
                     + nu_u235_fast * n25 * sigma_fis_u235_fast
                     + 2 * n28 * sigma_n2n_u238_fast
                     + 3 * n28 * sigma_n3n_u238_fast)
                    * phi_u_density_fast, energy_fast_c)

                abs_fast = np.trapezoid(
                    n28 * (sigma_a_u238_fast + sigma_n2n_u238_fast + sigma_n3n_u238_fast) * phi_u_density_fast, energy_fast_c)
                for nuclide, density in atom_density.items():
                    if nuclide != "U238" and nuclide in total_xs_bg_fast:
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
            "kgU238_per_m3_breeder": loading,
            "N_U238_atom_per_b_cm": n28,
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


def run_ir(case_name, make_mat, loadings, out_csv,
           flux_csv=None, flux_loading=0.0):
    """I_eff from the intermediate resonance (IR) approximation:

        phi(E, N28) = phi_u(E) * [lambda_28*sigma_p^28 + psi_lambda(E)]
                                / [sigma_a^28(E) + lambda_28*sigma_s^28(E) + psi_lambda(E)]

        psi_lambda(E) = (1/N28) * sum_{i != 28} N_i * [sigma_a^i(E) + lambda_i * sigma_s^i(E)]

    Background absorption (e.g. Li-6(n,t)) always removes neutrons regardless
    of scattering kinematics, so it enters psi_lambda at full weight; only the
    scattering channel is interpolated by lambda_i.
    """
    library = openmc.data.DataLibrary.from_xml()

    max_loading = max(l for l in loadings if l > 0)
    all_nuclides = sorted(make_mat(max_loading).get_nuclide_atom_densities().keys())
    print(f"\n{case_name} [IR]: {len(all_nuclides)} nuclides: {', '.join(all_nuclides)}")

    data = {}
    for n in all_nuclides:
        try:
            data[n] = openmc.data.IncidentNeutron.from_hdf5(
                str(hdf5_path(library, n)))
        except RuntimeError:
            print(f"  Warning: no cross-section data for {n}, skipping")
    labels = {n: temp_label(data[n].temperatures) for n in data}

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

    energy_fine_full = data["U238"].energy[labels["U238"]]
    fine_mask = (energy_fine_full >= E_MIN_EV) & (energy_fine_full <= E_MAX_EV)
    energy = energy_fine_full[fine_mask]

    pos = phi_u_density_coarse > 0
    log_density_fine = np.interp(np.log(energy), np.log(energy_coarse[pos]),
                                  np.log(phi_u_density_coarse[pos]))
    phi_u_density = np.exp(log_density_fine)

    sigma_gamma_u238 = eval_xs(data["U238"], labels["U238"], 102, energy)
    sigma_elastic_u238 = eval_xs(data["U238"], labels["U238"], 2, energy)
    try:
        sigma_fis_u238 = eval_xs(data["U238"], labels["U238"], 18, energy)
    except KeyError:
        sigma_fis_u238 = np.zeros_like(energy)
    sigma_a_u238 = sigma_gamma_u238 + sigma_fis_u238

    lambda_28 = LAMBDA[238]

    # Background cross sections: elastic (for lambda-weighted scattering) and
    # total (so we can extract absorption = total - elastic)
    elastic_xs_bg = {}
    total_xs_bg = {}
    for n in data:
        if n != "U238":
            try:
                elastic_xs_bg[n] = eval_xs(data[n], labels[n], 2, energy)
            except KeyError:
                elastic_xs_bg[n] = eval_xs(data[n], labels[n], 1, energy)
            total_xs_bg[n] = eval_xs(data[n], labels[n], 1, energy)

    Phi_u = np.trapezoid(phi_u_density, energy)
    ratio = Phi_u / np.sum(phi_u_coarse)
    print(f"  Phi_u: fine-grid trapz = {Phi_u:.6g}, coarse-tally sum = "
          f"{np.sum(phi_u_coarse):.6g} (ratio {ratio:.4f})")
    if abs(ratio - 1.0) > 0.05:
        print(f"  WARNING: fine/coarse Phi_u ratio off by >5%")

    i_infinite = np.trapezoid(sigma_gamma_u238 * phi_u_density, energy) / Phi_u

    # --- fast group: naive sigma*phi above E_MAX_EV (no self-shielding) ---
    R_F_per_n28 = compute_fast_group(data, labels, flux_csv, flux_loading)

    rows = []
    for loading in loadings:
        print(loading)
        material = make_mat(loading)
        atom_density = material.get_nuclide_atom_densities()
        n28 = atom_density.get("U238", 0.0)

        if n28 > 0.0:
            sigma_bg_lambda = np.zeros_like(energy)
            for nuclide, density in atom_density.items():
                if nuclide != "U238" and nuclide in elastic_xs_bg:
                    A = int(re.search(r'\d+', nuclide).group())
                    lam = LAMBDA[A]
                    print(f"  {nuclide}: A={A}, lambda={lam}")
                    sigma_a_bg = np.maximum(total_xs_bg[nuclide] - elastic_xs_bg[nuclide], 0.0)
                    sigma_bg_lambda += density * (sigma_a_bg + lam * elastic_xs_bg[nuclide])
            psi_lambda = sigma_bg_lambda / n28

            f_num = lambda_28 * SIGMA_P_U238 + psi_lambda
            f_den = sigma_a_u238 + lambda_28 * sigma_elastic_u238 + psi_lambda
            f_shield = f_num / f_den

            i_eff = np.trapezoid(sigma_gamma_u238 * f_shield * phi_u_density, energy) / Phi_u
            psi = np.trapezoid(psi_lambda * phi_u_density, energy) / Phi_u
            R_I = n28 * i_eff * Phi_u
            R_F = n28 * R_F_per_n28
        else:
            psi = np.inf
            i_eff = i_infinite
            R_I = 0.0
            R_F = 0.0

        rows.append({
            "kgU238_per_m3_breeder": loading,
            "N_U238_atom_per_b_cm": n28,
            "sigma_bg": psi,
            "I_infinite_b": i_infinite,
            "I_eff_b": i_eff,
            "I_eff_over_I_infinite": i_eff / i_infinite,
            "R_I": R_I,
            "R_F": R_F,
            "R_per_source_neutron": R_I + R_F,
        })

    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out_csv}")

    return rows


def run_irm(case_name, make_mat, loadings, out_csv,
            flux_csv=None, flux_loading=0.0):
    """IR with subcritical multiplication: R_IRM = M * R_IR, M = 1/(1-k_inf).

    k_inf is computed over both the resonance group (0.01-20 keV, using IR-
    shielded flux) and the fast group (20 keV-15 MeV, using unshielded phi_u),
    so that U-238 fast fission and (n,Xn) are captured.

    nu-bar is evaluated from ENDF data (MT=452) where available, falling
    back to the constants NU_U235 and NU_U238.
    """
    library = openmc.data.DataLibrary.from_xml()

    max_loading = max(l for l in loadings if l > 0)
    all_nuclides = sorted(make_mat(max_loading).get_nuclide_atom_densities().keys())
    print(f"\n{case_name} [IRM]: {len(all_nuclides)} nuclides: {', '.join(all_nuclides)}")

    data = {}
    for n in all_nuclides:
        try:
            data[n] = openmc.data.IncidentNeutron.from_hdf5(
                str(hdf5_path(library, n)))
        except RuntimeError:
            print(f"  Warning: no cross-section data for {n}, skipping")
    labels = {n: temp_label(data[n].temperatures) for n in data}

    # =====================================================================
    # RESONANCE GROUP (E_MIN_EV -- E_MAX_EV) on U-238's fine ENDF grid
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

    energy_fine_full = data["U238"].energy[labels["U238"]]
    fine_mask = (energy_fine_full >= E_MIN_EV) & (energy_fine_full <= E_MAX_EV)
    energy = energy_fine_full[fine_mask]

    pos = phi_u_density_coarse > 0
    log_density_fine = np.interp(np.log(energy), np.log(energy_coarse[pos]),
                                  np.log(phi_u_density_coarse[pos]))
    phi_u_density = np.exp(log_density_fine)

    # U-238 cross sections
    sigma_gamma_u238 = eval_xs(data["U238"], labels["U238"], 102, energy)
    sigma_elastic_u238 = eval_xs(data["U238"], labels["U238"], 2, energy)
    sigma_tot_u238 = eval_xs(data["U238"], labels["U238"], 1, energy)
    try:
        sigma_fis_u238 = eval_xs(data["U238"], labels["U238"], 18, energy)
    except KeyError:
        sigma_fis_u238 = np.zeros_like(energy)
    sigma_a_u238 = sigma_gamma_u238 + sigma_fis_u238

    lambda_28 = LAMBDA[238]

    # Background cross sections: elastic (for psi_lambda), total, and
    # inelastic (for corrected k_inf removal = nonelastic - inelastic)
    elastic_xs_bg = {}
    total_xs_bg = {}
    inelastic_xs_bg = {}
    for n in data:
        if n != "U238":
            try:
                elastic_xs_bg[n] = eval_xs(data[n], labels[n], 2, energy)
            except KeyError:
                elastic_xs_bg[n] = eval_xs(data[n], labels[n], 1, energy)
            total_xs_bg[n] = eval_xs(data[n], labels[n], 1, energy)
            try:
                inelastic_xs_bg[n] = eval_xs(data[n], labels[n], 4, energy)
            except KeyError:
                inelastic_xs_bg[n] = np.zeros_like(energy)

    # U-238 (n,2n) and (n,3n) on resonance grid
    try:
        sigma_n2n_u238 = eval_xs(data["U238"], labels["U238"], 16, energy)
    except KeyError:
        sigma_n2n_u238 = np.zeros_like(energy)
    try:
        sigma_n3n_u238 = eval_xs(data["U238"], labels["U238"], 17, energy)
    except KeyError:
        sigma_n3n_u238 = np.zeros_like(energy)

    # U-235 fission (for k_inf production term)
    sigma_fis_u235 = np.zeros_like(energy)
    if "U235" in data:
        try:
            sigma_fis_u235 = eval_xs(data["U235"], labels["U235"], 18, energy)
        except KeyError:
            pass

    # nu-bar on resonance grid (constants are fine: U-238 fission is zero
    # below ~1 MeV, and nu_25 ~ 2.43 is nearly flat in the resonance range)
    nu_u238_res = NU_U238 * np.ones_like(energy)
    nu_u235_res = NU_U235 * np.ones_like(energy)

    Phi_u = np.trapezoid(phi_u_density, energy)
    ratio = Phi_u / np.sum(phi_u_coarse)
    print(f"  Phi_u (res): fine-grid trapz = {Phi_u:.6g}, coarse-tally sum = "
          f"{np.sum(phi_u_coarse):.6g} (ratio {ratio:.4f})")
    if abs(ratio - 1.0) > 0.05:
        print(f"  WARNING: fine/coarse Phi_u ratio off by >5%")

    i_infinite = np.trapezoid(sigma_gamma_u238 * phi_u_density, energy) / Phi_u

    # =====================================================================
    # FAST GROUP (E_MAX_EV -- E_FAST_MAX_EV) on coarse tally grid
    # =====================================================================
    energy_fast_c, phi_u_fast_c = read_flux_generic(
        flux_csv, flux_loading, e_min=E_MAX_EV, e_max=E_FAST_MAX_EV)
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

        # U-238 xs on fast grid
        sigma_gamma_u238_fast = eval_xs(data["U238"], labels["U238"], 102, energy_fast_c)
        try:
            sigma_fis_u238_fast = eval_xs(data["U238"], labels["U238"], 18, energy_fast_c)
        except KeyError:
            sigma_fis_u238_fast = np.zeros_like(energy_fast_c)
        sigma_a_u238_fast = sigma_gamma_u238_fast + sigma_fis_u238_fast

        # U-235 fission on fast grid
        sigma_fis_u235_fast = np.zeros_like(energy_fast_c)
        if "U235" in data:
            try:
                sigma_fis_u235_fast = eval_xs(data["U235"], labels["U235"], 18, energy_fast_c)
            except KeyError:
                pass

        # nu-bar on fast grid -- try ENDF (MT=452), fall back to constants
        try:
            nu_u238_fast = eval_xs(data["U238"], labels["U238"], 452, energy_fast_c)
            print(f"  nu_U238: using ENDF MT=452 (range {nu_u238_fast.min():.2f}–{nu_u238_fast.max():.2f})")
        except (KeyError, Exception):
            nu_u238_fast = NU_U238 * np.ones_like(energy_fast_c)
            print(f"  nu_U238: MT=452 not available, using constant {NU_U238}")
        try:
            nu_u235_fast = eval_xs(data["U235"], labels["U235"], 452, energy_fast_c)
        except (KeyError, Exception):
            nu_u235_fast = NU_U235 * np.ones_like(energy_fast_c)

        # U-238 (n,2n) and (n,3n) on fast grid
        try:
            sigma_n2n_u238_fast = eval_xs(data["U238"], labels["U238"], 16, energy_fast_c)
        except KeyError:
            sigma_n2n_u238_fast = np.zeros_like(energy_fast_c)
        try:
            sigma_n3n_u238_fast = eval_xs(data["U238"], labels["U238"], 17, energy_fast_c)
        except KeyError:
            sigma_n3n_u238_fast = np.zeros_like(energy_fast_c)

        # Background xs on fast grid (for k_inf removal)
        total_xs_bg_fast = {}
        elastic_xs_bg_fast = {}
        inelastic_xs_bg_fast = {}
        for n in data:
            if n != "U238":
                total_xs_bg_fast[n] = eval_xs(data[n], labels[n], 1, energy_fast_c)
                try:
                    elastic_xs_bg_fast[n] = eval_xs(data[n], labels[n], 2, energy_fast_c)
                except KeyError:
                    elastic_xs_bg_fast[n] = total_xs_bg_fast[n]
                try:
                    inelastic_xs_bg_fast[n] = eval_xs(data[n], labels[n], 4, energy_fast_c)
                except KeyError:
                    inelastic_xs_bg_fast[n] = np.zeros_like(energy_fast_c)

        R_F_per_n28 = np.trapezoid(sigma_gamma_u238_fast * phi_u_density_fast, energy_fast_c)
        print(f"  fast group ({E_MAX_EV:.0f}–{E_FAST_MAX_EV:.0e} eV): "
              f"R_F/n28 = {R_F_per_n28:.6g},  {len(energy_fast_c)} tally bins")
    else:
        R_F_per_n28 = 0.0

    # =====================================================================
    # LOADING LOOP
    # =====================================================================
    rows = []
    for loading in loadings:
        print(loading)
        material = make_mat(loading)
        atom_density = material.get_nuclide_atom_densities()
        n28 = atom_density.get("U238", 0.0)
        n25 = atom_density.get("U235", 0.0)

        if n28 > 0.0:
            # --- IR self-shielding (same as run_ir) ---
            sigma_bg_lambda = np.zeros_like(energy)
            for nuclide, density in atom_density.items():
                if nuclide != "U238" and nuclide in elastic_xs_bg:
                    A = int(re.search(r'\d+', nuclide).group())
                    lam = LAMBDA[A]
                    print(f"  {nuclide}: A={A}, lambda={lam}")
                    sigma_a_bg = np.maximum(total_xs_bg[nuclide] - elastic_xs_bg[nuclide], 0.0)
                    sigma_bg_lambda += density * (sigma_a_bg + lam * elastic_xs_bg[nuclide])
            psi_lambda = sigma_bg_lambda / n28

            f_num = lambda_28 * SIGMA_P_U238 + psi_lambda
            f_den = sigma_a_u238 + lambda_28 * sigma_elastic_u238 + psi_lambda
            f_shield = f_num / f_den

            i_eff_ir = np.trapezoid(sigma_gamma_u238 * f_shield * phi_u_density, energy) / Phi_u
            psi = np.trapezoid(psi_lambda * phi_u_density, energy) / Phi_u
            R_I_base = n28 * i_eff_ir * Phi_u
            R_F_base = n28 * R_F_per_n28

            # --- subcritical multiplication (both groups) ---
            phi_shielded = f_shield * phi_u_density

            # Resonance-group production & absorption
            prod_res = np.trapezoid(
                (nu_u238_res * n28 * sigma_fis_u238
                 + nu_u235_res * n25 * sigma_fis_u235
                 + 2 * n28 * sigma_n2n_u238
                 + 3 * n28 * sigma_n3n_u238)
                * phi_shielded, energy)

            abs_res = np.trapezoid(
                n28 * (sigma_a_u238 + sigma_n2n_u238 + sigma_n3n_u238) * phi_shielded, energy)
            for nuclide, density in atom_density.items():
                if nuclide != "U238" and nuclide in total_xs_bg:
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
                    (nu_u238_fast * n28 * sigma_fis_u238_fast
                     + nu_u235_fast * n25 * sigma_fis_u235_fast
                     + 2 * n28 * sigma_n2n_u238_fast
                     + 3 * n28 * sigma_n3n_u238_fast)
                    * phi_u_density_fast, energy_fast_c)

                abs_fast = np.trapezoid(
                    n28 * (sigma_a_u238_fast + sigma_n2n_u238_fast + sigma_n3n_u238_fast) * phi_u_density_fast, energy_fast_c)
                for nuclide, density in atom_density.items():
                    if nuclide != "U238" and nuclide in total_xs_bg_fast:
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
            i_eff_ir = i_infinite
            k_inf = 0.0
            M = 1.0
            R_I = 0.0
            R_F = 0.0

        rows.append({
            "kgU238_per_m3_breeder": loading,
            "N_U238_atom_per_b_cm": n28,
            "sigma_bg": psi,
            "I_infinite_b": i_infinite,
            "I_eff_b": i_eff_ir,
            "I_eff_over_I_infinite": i_eff_ir / i_infinite,
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
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


def read_openmc_generic(csv_path):
    with csv_path.open() as f:
        rows = list(csv.DictReader(f))
    return {
        "x": np.array([float(r["fertile_kg/m3"]) for r in rows]),
        "y": np.array([float(r["U238(n,g)"]) for r in rows]),
        "sd": np.array([float(r["U238(n,g)_sd"]) for r in rows]),
    }


def read_flux_generic(flux_csv, loading=0.0, e_min=None, e_max=None):
    if e_min is None:
        e_min = E_MIN_EV
    if e_max is None:
        e_max = E_MAX_EV
    energy, flux = [], []
    with flux_csv.open() as f:
        for row in csv.DictReader(f):
            if float(row["fertile_kg/m3"]) == loading:
                e = float(row["energy mid [eV]"])
                if e_min <= e <= e_max:
                    energy.append(e)
                    flux.append(float(row["mean"]))
    return np.array(energy), np.array(flux)


def plot_compare_all(cases, log_scale):
    """Combined NR, NRM, IR, IRM, and OpenMC for all blankets.

    cases: list of (label, nr_rows, nrm_rows, ir_rows, irm_rows, x_key, openmc_csv_path)
    """
    fig, ax = plt.subplots(figsize=(5.8, 3.8))

    if log_scale:
        ax.set_xscale("log")
        ax.set_yscale("log")
        out = OUT_COMPARE_LOG
    else:
        out = OUT_COMPARE_LINEAR

    def mask_for(x, y):
        return (x > 0) & (y > 0) if log_scale else np.ones_like(x, dtype=bool)

    for idx, (label, nr_rows, nrm_rows, ir_rows, irm_rows, x_key, openmc_csv) in enumerate(cases):
        c = BLANKET_COLORS[idx]

        omc = read_openmc_generic(openmc_csv)
        omc_mask = mask_for(omc["x"], omc["y"])
        ax.errorbar(omc["x"][omc_mask], omc["y"][omc_mask],
                     yerr=omc["sd"][omc_mask],
                     fmt="o-", color=c, lw=0.75, markersize=2,
                     capsize=2, label=f"{label} OpenMC")

        for tag, rows, ls in [("IRM", irm_rows, ":"),
                               ("NRM", nrm_rows, (0, (3, 1, 1, 1))),
                               ("IR",  ir_rows,  "--"),
                               ("NR",  nr_rows,  "-.")]:
            x_m = np.array([r[x_key] for r in rows])
            y_m = np.array([r["R_per_source_neutron"] for r in rows])
            mask = mask_for(x_m, y_m)
            ax.plot(x_m[mask], y_m[mask], ls=ls, color=c, lw=0.75,
                    markersize=1, label=f"{label} {tag}")

    ax.set_xlabel("U-238 loading [kg/m³ breeder]")
    ax.set_ylabel("U-238(n,γ) [reactions/source neutron]")
    ax.grid(True, which="both", alpha=0.35)
    ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
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
    """psi vs I_eff for all blankets, log-log, showing NR, NRM, IR, and IRM.

    cases: list of (label, nr_rows, nrm_rows, ir_rows, irm_rows)
    """
    fig, ax = plt.subplots(figsize=(5.8, 3.8))

    for idx, (label, nr_rows, nrm_rows, ir_rows, irm_rows) in enumerate(cases):
        color = BLANKET_COLORS[idx]

        for tag, rows, marker, ls in [("IRM", irm_rows, "D", ":"),
                                       ("NRM", nrm_rows, "^", (0, (3, 1, 1, 1))),
                                       ("IR",  ir_rows,  "s", "--"),
                                       ("NR",  nr_rows,  "o", "-.")]:
            psi = np.array([r["sigma_bg"] for r in rows])
            i_eff = np.array([r["I_eff_b"] for r in rows])
            mask = np.isfinite(psi) & (psi > 0) & (i_eff > 0)

            ax.plot(psi[mask], i_eff[mask], marker, color=color, markersize=3,
                    linestyle="none", label=f"{label} {tag}")

            a, c, r2, n_used, psi_boundary = fit_power_law_progressive(
                psi[mask], i_eff[mask], r2_min)
            n_total = int(mask.sum())
            print(f"{label} {tag}: I_eff = {a:.4g} * psi^{c:.4g}   "
                  f"(R2={r2:.4f} using {n_used}/{n_total} points up to "
                  f"psi={psi_boundary:.3g} b)")

            x_fit = np.geomspace(psi[mask].min(), psi_boundary, 100)
            ax.plot(x_fit, a * x_fit**c, ls=ls, color=color, lw=1.0,
                    label=rf"{label} {tag} fit: $I_\mathrm{{eff}} = {a:.3g}\,\psi^{{{c:.3g}}}$")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\psi$ [b]")
    ax.set_ylabel(r"$I_\mathrm{eff}$ [b]")
    ax.grid(True, which="both", alpha=0.35)
    ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(OUT_PSI_IEFF, dpi=300)
    print(f"wrote {OUT_PSI_IEFF}")


def main():
    flibe_nr = read_cached_rows(OUT_CSV_FLIBE_NR) or run_nr("FLiBe (UF4)", make_flibe, LOADINGS, OUT_CSV_FLIBE_NR, flux_csv=CSV_FLUX_FLIBE)
    hcpb_nr  = read_cached_rows(OUT_CSV_HCPB_NR)  or run_nr("HCPB (Li4SiO4-Be + UO2-SiC)", make_hcpb, LOADINGS, OUT_CSV_HCPB_NR, flux_csv=CSV_FLUX_HCPB)
    dcll_nr  = read_cached_rows(OUT_CSV_DCLL_NR)  or run_nr("DCLL (Pb-17Li + UO2-SiC)", make_dcll, LOADINGS, OUT_CSV_DCLL_NR, flux_csv=CSV_FLUX_DCLL)

    flibe_nrm = read_cached_rows(OUT_CSV_FLIBE_NRM) or run_nrm("FLiBe (UF4)", make_flibe, LOADINGS, OUT_CSV_FLIBE_NRM, flux_csv=CSV_FLUX_FLIBE)
    hcpb_nrm  = read_cached_rows(OUT_CSV_HCPB_NRM)  or run_nrm("HCPB (Li4SiO4-Be + UO2-SiC)", make_hcpb, LOADINGS, OUT_CSV_HCPB_NRM, flux_csv=CSV_FLUX_HCPB)
    dcll_nrm  = read_cached_rows(OUT_CSV_DCLL_NRM)  or run_nrm("DCLL (Pb-17Li + UO2-SiC)", make_dcll, LOADINGS, OUT_CSV_DCLL_NRM, flux_csv=CSV_FLUX_DCLL)

    flibe_ir  = read_cached_rows(OUT_CSV_FLIBE_IR)  or run_ir("FLiBe (UF4)", make_flibe, LOADINGS, OUT_CSV_FLIBE_IR, flux_csv=CSV_FLUX_FLIBE)
    hcpb_ir   = read_cached_rows(OUT_CSV_HCPB_IR)  or run_ir("HCPB (Li4SiO4-Be + UO2-SiC)", make_hcpb, LOADINGS, OUT_CSV_HCPB_IR, flux_csv=CSV_FLUX_HCPB)
    dcll_ir   = read_cached_rows(OUT_CSV_DCLL_IR)  or run_ir("DCLL (Pb-17Li + UO2-SiC)", make_dcll, LOADINGS, OUT_CSV_DCLL_IR, flux_csv=CSV_FLUX_DCLL)

    flibe_irm = read_cached_rows(OUT_CSV_FLIBE_IRM) or run_irm("FLiBe (UF4)", make_flibe, LOADINGS, OUT_CSV_FLIBE_IRM, flux_csv=CSV_FLUX_FLIBE)
    hcpb_irm  = read_cached_rows(OUT_CSV_HCPB_IRM)  or run_irm("HCPB (Li4SiO4-Be + UO2-SiC)", make_hcpb, LOADINGS, OUT_CSV_HCPB_IRM, flux_csv=CSV_FLUX_HCPB)
    dcll_irm  = read_cached_rows(OUT_CSV_DCLL_IRM)  or run_irm("DCLL (Pb-17Li + UO2-SiC)", make_dcll, LOADINGS, OUT_CSV_DCLL_IRM, flux_csv=CSV_FLUX_DCLL)

    cases = [
        ("FLiBe", flibe_nr, flibe_nrm, flibe_ir, flibe_irm, "kgU238_per_m3_breeder", CSV_RXNS_FLIBE),
        ("HCPB",  hcpb_nr,  hcpb_nrm,  hcpb_ir,  hcpb_irm,  "kgU238_per_m3_breeder", CSV_RXNS_HCPB),
        ("DCLL",  dcll_nr,  dcll_nrm,  dcll_ir,  dcll_irm,  "kgU238_per_m3_breeder", CSV_RXNS_DCLL),
    ]
    plot_compare_all(cases, log_scale=True)
    plot_compare_all(cases, log_scale=False)

    plot_psi_vs_ieff([
        ("FLiBe", flibe_nr, flibe_nrm, flibe_ir, flibe_irm),
        ("HCPB",  hcpb_nr,  hcpb_nrm,  hcpb_ir,  hcpb_irm),
        ("DCLL",  dcll_nr,  dcll_nrm,  dcll_ir,  dcll_irm),
    ], r2_min=R2_MIN)


if __name__ == "__main__":
    main()