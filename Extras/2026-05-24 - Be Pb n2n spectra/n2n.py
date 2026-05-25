"""
Thin-target secondary neutron spectra for Be-9, Pb, Li-7, natural U, and Th-232.

Runs separate OpenMC fixed-source calculations (one per target) in
dedicated subdirectories, then overlays the results on the same plots.

Each target tallies reaction-specific outgoing spectra:
  - Be-9, Pb, Li-7 : (n,2n) and inelastic (MT 51-91)
  - U, Th          : (n,2n), inelastic (MT 51-91), and nu-fission
  Discrete inelastic levels are summed in post-processing.
"""

import openmc
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import os
from pathlib import Path
import pandas as pd
from params import *

# ============================================================
# USER SETTINGS
# ============================================================
TARGETS = ['Be', 'Pb', 'Li7', 'U', 'Th']
THICKNESS = 0.1        # cm — thin shell (keep << MFP)
INNER_R   = 1.0        # cm
PARTICLES = int(1e8)
BATCHES = int(1e2)

# Mask plot bins whose relative statistical error exceeds this value.
# Set to None or float('inf') to disable masking.
REL_ERR_THRESHOLD = 0.15

# Fine energy bins: 1 keV to 16 MeV, 1000 log-spaced bins (eV)
E_BINS = np.logspace(3, 7.21, 1000)
E_CENTERS = 0.5 * (E_BINS[:-1] + E_BINS[1:])
E_WIDTHS = np.diff(E_BINS)

ORIG_DIR = Path.cwd()
RUN_DIR = Path('./')
RUN_DIR.mkdir(exist_ok=True)

# ============================================================
# Per-target scores for the outgoing-energy spectrum tally
# ============================================================
# Individual discrete inelastic levels (MT=51 through MT=91)
INELASTIC_MTS = [str(mt) for mt in range(51, 92)]

SPECIFIC_SCORES = {
    'Be':  ['(n,2n)'] + INELASTIC_MTS,
    'Pb':  ['(n,2n)'] + INELASTIC_MTS,
    'Li7': ['(n,2n)'] + INELASTIC_MTS,
    'U':   ['(n,2n)', 'nu-fission'] + INELASTIC_MTS,
    'Th':  ['(n,2n)', 'nu-fission'] + INELASTIC_MTS,
}

# Which curves to show on the plot
PLOT_CURVES = {
    'Be':  ['(n,2n)'],
    'Pb':  ['(n,2n)', "(n,n')"],
    'Li7': ['(n,2n)', "(n,n')"],
    'U':   ['nu-fission', '(n,2n)', "(n,n')"],
    'Th':  ['nu-fission', '(n,2n)', "(n,n')"],
}

def log_buffer(lo, hi, buf=0.03):
    f = (hi / lo) ** buf
    return (lo / f, hi * f)

# ============================================================
# Extract results from a statepoint file
# ============================================================
def _extract_results(target_name, sp_file):
    """Read tallies from an existing statepoint and return results dict."""
    sp = openmc.StatePoint(str(sp_file))

    # ---- Reaction rates ----
    t = sp.get_tally(name='reaction_rates')
    rxn_rates = {}
    # Known non-inelastic scores (OpenMC preserves these names)
    non_inelastic = {'(n,2n)', '(n,3n)', '(n,elastic)',
                     'absorption', '(n,gamma)', '(n,Xt)', 'fission'}

    # Extract non-inelastic scores individually
    for score in t.scores:
        if score in non_inelastic:
            val = t.get_values(scores=[score]).flatten()[0]
            err = t.get_values(scores=[score], value='std_dev').flatten()[0]
            rxn_rates[score] = (val, err)

    # Sum discrete inelastic levels (everything else) into one entry
    inelastic_val = 0.0
    inelastic_err_sq = 0.0
    for score in t.scores:
        if score not in non_inelastic:
            val = t.get_values(scores=[score]).flatten()[0]
            err = t.get_values(scores=[score], value='std_dev').flatten()[0]
            inelastic_val += val
            inelastic_err_sq += err**2
    rxn_rates["(n,n')"] = (inelastic_val, np.sqrt(inelastic_err_sq))

    print(f"\nReaction rates for {target_name} (per source neutron):")
    for score, (val, err) in rxn_rates.items():
        print(f"  {score:15s}: {val:.6e} +/- {err:.6e}")

    # ---- Reaction-specific outgoing spectra ----
    t_ns = sp.get_tally(name='specific_scatter_eout')
    spectra = {}

    # Scores to extract individually (not summed into inelastic)
    individual_scores = {'(n,2n)', 'nu-fission'}
    for score in t_ns.scores:
        if score in individual_scores:
            vals = t_ns.get_values(scores=[score]).flatten()
            errs = t_ns.get_values(scores=[score], value='std_dev').flatten()
            spectra[score] = {
                'vals': vals / E_WIDTHS,
                'errs': errs / E_WIDTHS,
            }

    # Sum all remaining scores (the discrete inelastic levels MT 51-91,
    # which OpenMC renames to (n,n1), (n,n2), … , (n,nc))
    inelastic_scores = [s for s in t_ns.scores if s not in individual_scores]
    vals_sum = np.zeros_like(E_CENTERS)
    errs_sq  = np.zeros_like(E_CENTERS)
    for s in inelastic_scores:
        v = t_ns.get_values(scores=[s]).flatten()
        e = t_ns.get_values(scores=[s], value='std_dev').flatten()
        vals_sum += v
        errs_sq  += e**2
    spectra["(n,n')"] = {
        'vals': vals_sum / E_WIDTHS,
        'errs': np.sqrt(errs_sq) / E_WIDTHS,
    }

    sp.close()
    os.chdir(ORIG_DIR)

    return {
        'rxn_rates': rxn_rates,
        'spectra':   spectra,
    }


# ============================================================
# Build, run, and extract results for one target
# ============================================================
def run_target(target_name):
    """Build model, run OpenMC, return dict of results."""

    openmc.reset_auto_ids()

    sub_dir = RUN_DIR / target_name
    sub_dir.mkdir(exist_ok=True)
    os.chdir(sub_dir)
    print(f"\n{'='*60}")
    print(f"  Target: {target_name}  |  Directory: {sub_dir}")
    print(f"{'='*60}")

    sp_file = Path(f'statepoint.{BATCHES}.h5')
    if sp_file.exists():
        print(f"  Statepoint already exists — skipping build & run.")
        return _extract_results(target_name, sp_file)

    # --- Material ---
    mat = openmc.Material(name=target_name)
    if target_name == 'Be':
        mat.add_element('Be', 1.0)
        mat.set_density('g/cm3', 1.85)
    elif target_name == 'Pb':
        mat.add_element('Pb', 1.0)
        mat.set_density('g/cm3', 11.35)
    elif target_name == 'Li7':
        mat.add_nuclide('Li7', 1.0)
        mat.set_density('g/cm3', 0.534)
    elif target_name == 'U':
        mat.add_element('U', 1.0)
        mat.set_density('g/cm3', 19.1)
    elif target_name == 'Th':
        mat.add_element('Th', 1.0)
        mat.set_density('g/cm3', 11.7)
    else:
        raise ValueError(f"Unknown target: {target_name}")

    materials = openmc.Materials([mat])

    # --- Geometry ---
    outer_r = INNER_R + THICKNESS
    s_inner = openmc.Sphere(r=INNER_R)
    s_outer = openmc.Sphere(r=outer_r)
    s_bound = openmc.Sphere(r=outer_r + 5.0, boundary_type='vacuum')

    cell_void_in  = openmc.Cell(fill=None, region=-s_inner)
    cell_target   = openmc.Cell(fill=mat,  region=+s_inner & -s_outer)
    cell_void_out = openmc.Cell(fill=None, region=+s_outer & -s_bound)

    geometry = openmc.Geometry([cell_void_in, cell_target, cell_void_out])

    # --- Source ---
    source = openmc.IndependentSource()
    source.space = openmc.stats.Point((0, 0, 0))
    source.energy = openmc.stats.Discrete([14.1e6], [1.0])
    source.particle = 'neutron'

    settings = openmc.Settings()
    settings.source = [source]
    settings.run_mode = 'fixed source'
    settings.batches = BATCHES
    settings.particles = PARTICLES

    # --- Tallies ---
    cell_filter  = openmc.CellFilter(cell_target)
    e_in_filter  = openmc.EnergyFilter([14.0e6, 14.2e6])
    e_out_filter = openmc.EnergyoutFilter(E_BINS)

    # 1) Integral reaction rates in the target cell
    t_rates = openmc.Tally(name='reaction_rates')
    t_rates.filters = [cell_filter]
    t_rates.scores = [
        '(n,2n)', '(n,3n)', '(n,elastic)',
        'absorption', '(n,gamma)', '(n,Xt)', 'fission',
    ] + INELASTIC_MTS

    # 2) Reaction-specific outgoing neutron spectra
    #    Each target gets only the scores relevant to it.
    t_specific = openmc.Tally(name='specific_scatter_eout')
    t_specific.filters = [cell_filter, e_in_filter, e_out_filter]
    t_specific.scores = SPECIFIC_SCORES[target_name]

    tallies = openmc.Tallies([t_rates, t_specific])

    # --- Export & Run ---
    materials.export_to_xml()
    geometry.export_to_xml()
    settings.export_to_xml()
    tallies.export_to_xml()

    openmc.run(cwd='.', output=True)

    return _extract_results(target_name, Path(f'statepoint.{BATCHES}.h5'))


# ============================================================
# Run all targets
# ============================================================
results = {}
for tgt in TARGETS:
    results[tgt] = run_target(tgt)


# ============================================================
# Combined plots
# ============================================================
os.chdir(RUN_DIR)

# --- Single plot: reaction-specific outgoing spectra ---
fig1, ax1 = plt.subplots(figsize=(7.0, 3.0))

color_idx = 0
for tgt in TARGETS:
    for score_name in PLOT_CURVES[tgt]:
        sdata = results[tgt]['spectra'][score_name]
        spec = sdata['vals']
        errs = sdata['errs']

        # Mask: keep bins with positive value AND relative error below threshold
        mask = spec > 0
        if REL_ERR_THRESHOLD is not None:
            with np.errstate(divide='ignore', invalid='ignore'):
                rel_err = np.where(spec > 0, errs / spec, np.inf)
            mask &= rel_err < REL_ERR_THRESHOLD

        if mask.any():
            label = f"{tgt} {score_name}"
            ax1.step(E_CENTERS[mask], spec[mask], where='mid', lw=0.75,
                     color=f'C{color_idx}', label=label)
        color_idx += 1

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('Outgoing neutron energy [eV]')
ax1.set_ylabel(r'Produced neutrons$/$src-n$/$eV')
leg = ax1.legend(loc='lower left', fontsize=7, ncol=2,
                 fancybox=False, edgecolor='black',
                 frameon=True, framealpha=0.75)
leg.get_frame().set_linewidth(0.5)
ax1.set_xlim(log_buffer(1e4, 14e6, buf=0.03))
ax1.set_ylim(log_buffer(1e-11, 1e-8, buf=0.03))


# Format tick labels as 1eX scientific notation on both axes
sci_fmt = FuncFormatter(lambda x, _: f'1e{int(round(np.log10(x)))}' if x > 0 else '0')
ax1.xaxis.set_major_formatter(sci_fmt)
ax1.yaxis.set_major_formatter(sci_fmt)
fig1.tight_layout()
fig1.savefig('compare.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
fig1.savefig('compare.png', bbox_inches='tight', pad_inches=0.01, format='png', dpi=300)

print(f"\nSaved {RUN_DIR / 'compare.png'}")

# --- Export all spectra to CSV ---
csv_data = {}
for tgt in TARGETS:
    for score_name, sdata in results[tgt]['spectra'].items():
        csv_data[f"{tgt} {score_name}"] = sdata['vals']
        csv_data[f"{tgt} {score_name} err"] = sdata['errs']

df = pd.DataFrame(csv_data, index=E_CENTERS)
df.index.name = 'energy_mid_eV'
df.to_csv('spectra.csv')
print(f"Saved {RUN_DIR / 'spectra.csv'}")

# --- Print reaction rate comparison ---
print(f"\n{'='*120}")
print(f"{'Reaction':<16s}", end='')
for tgt in TARGETS:
    print(f"  {tgt:>16s}", end='')
print()
print('-' * 120)

all_scores = list(results[TARGETS[0]]['rxn_rates'].keys())
for score in all_scores:
    print(f"{score:<16s}", end='')
    for tgt in TARGETS:
        try:
            val, err = results[tgt]['rxn_rates'][score]
            print(f"  {val:>12.4e} ± {err:.0e}", end='')
        except KeyError:
            print(f"  {'N/A':>16s}", end='')
    print()

os.chdir(ORIG_DIR)
print(f"\nDone. All files in {RUN_DIR}/")