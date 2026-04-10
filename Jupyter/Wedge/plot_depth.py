"""
plot_depth_profile.py

4x2 reaction-rate vs depth profiles through the outboard breeding region.

Layout:
    Columns : 0.1 kg/m³ (left)  |  1000 kg/m³ (right)
    Rows    : TBR  |  FPR  |  fission  |  (n,2n)
    Lines   : HCPB homog., HCPB lattice, DCLL homog., DCLL lattice

Usage:
    python plot_depth_profile.py
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from prism_utilities import *


# ──────────────────────────────────────────────────────────────────────────────
# USER CONFIGURATION
# ──────────────────────────────────────────────────────────────────────────────

ISOTOPE = 'U238'
FERTILE_LOADINGS = [0.1, 999.99]
CASES = ['A', 'C']
OPENMC_BASE = './OpenMC'

BLANKETS = {
    'DCLL': {'prefix': 'dcll', 'enrich': '90.0'},
    'HCPB': {'prefix': 'hcpb', 'enrich': '60.0'},
}

# Depth tally names (must match tallies() in model scripts)
TALLY_NUCLIDE_DEPTH = 'Reactions per nuclide-depth'    # has nuclide filter
TALLY_TOTAL_DEPTH   = 'Reactions per depth'            # no nuclide filter, has 'fission'

# Z-bounds of the outboard breeding region (cm) — used to compute depth from FW
Z_BOUNDS = {
    'DCLL': (480.7, 550.1),
    'HCPB': (450.4, 532.4),
}

# Per-blanket (n,2n) nuclides
N2N_NUCLIDES = {
    'DCLL': ['Pb204', 'Pb206', 'Pb207', 'Pb208'],
    'HCPB': ['Be9'],
}


# ──────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

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
        raise FileNotFoundError(f"No statepoint found for pattern:\n  {pattern}")
    return matches[-1]


def extract_depth_profile(sp_path, tally_name, score, nuclides=None):
    """
    Extract a 1-D depth profile (reaction rate vs mesh z-bin) from a
    depth-mesh tally.

    If nuclides is not None, slices by each nuclide and sums.
    If nuclides is None, just slices by score (for total tallies).

    Returns:
        depth_mid  – 1-D array of mesh-bin z-midpoints (cm, absolute)
        mean       – 1-D array of reaction rates per depth bin
    """
    import openmc

    sp = openmc.StatePoint(sp_path)
    tally = sp.get_tally(name=tally_name)

    # Find the mesh filter to get the z-bin structure
    mesh_filter = None
    for f in tally.filters:
        if isinstance(f, openmc.MeshFilter):
            mesh_filter = f
            break

    mesh = mesh_filter.mesh
    n_z = mesh.dimension[2]
    z_lo = mesh.lower_left[2]
    z_hi = mesh.upper_right[2]
    z_edges = np.linspace(z_lo, z_hi, n_z + 1)
    depth_mid = 0.5 * (z_edges[:-1] + z_edges[1:])

    if nuclides is not None:
        total = np.zeros(n_z)
        for nuc in nuclides:
            sl = tally.get_slice(scores=[score], nuclides=[nuc])
            vals = sl.mean.flatten()
            # With mesh dim [1,1,N], vals length is N (one per z-bin)
            # but if there are extra dimensions, reshape
            n_other = len(vals) // n_z
            total += vals.reshape(n_other, n_z).sum(axis=0)
    else:
        sl = tally.get_slice(scores=[score])
        vals = sl.mean.flatten()
        n_other = len(vals) // n_z
        total = vals.reshape(n_other, n_z).sum(axis=0)

    sp.close()
    return depth_mid, total


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():

    # ── Nuclide / label config ────────────────────────────────────────────
    if ISOTOPE == 'U238':
        fpr_nuclides = ['U238']
        fis_nuclides = ['U235', 'U238']
        XO2 = r'UO$_2$'
        X   = r'U'
    else:
        fpr_nuclides = ['Th232']
        fis_nuclides = ['Th232']
        XO2 = r'ThO$_2$'
        X   = r'Th'

    tbr_nuclides = ['Li6', 'Li7']

    # ── Line style map ────────────────────────────────────────────────────
    line_styles = {
        ('HCPB', 'A'): {'color': '#ff1f5b', 'label': rf'HCPB-{XO2} - homog.'},
        ('HCPB', 'C'): {'color': '#f48628', 'label': rf'HCPB-{XO2} - lattice'},
        ('DCLL', 'A'): {'color': '#04cc6c', 'label': rf'DCLL-{XO2} - homog.'},
        ('DCLL', 'C'): {'color': '#0c9edd', 'label': rf'DCLL-{XO2} - lattice'},
    }

    # ── Row config ────────────────────────────────────────────────────────
    # (row, qty_key, tally_name, score, nuclides, title_suffix)
    # For FIS: use total tally (has 'fission' score without nuclide filter)
    # For N2N: nuclides resolved per-blanket below
    row_config = [
        (0, 'TBR', TALLY_NUCLIDE_DEPTH, '(n,Xt)',    tbr_nuclides, rf'Li(n,Xt)'),
        (1, 'FPR', TALLY_NUCLIDE_DEPTH, '(n,gamma)', fpr_nuclides, rf'{X}(n,$\gamma$)'),
        (2, 'FIS', TALLY_TOTAL_DEPTH,   'fission',   None,         rf'{X}(n,fis)'),
        (3, 'N2N', TALLY_NUCLIDE_DEPTH, '(n,2n)',    None,         r'Total (n,2n)'),
    ]

    # ── Create 4×2 figure ─────────────────────────────────────────────────
    col_loadings = {0: 0.1, 1: 999.99}
    col_labels   = {0: r'0.1 kg$/$m³', 1: r'1000 kg$/$m³'}

    fig, axes = plt.subplots(4, 2, figsize=(16, 24), sharex=False, sharey=False)

    for row, qty, tally_name, score, nuclides, title_suffix in row_config:
        for col, loading in col_loadings.items():
            ax = axes[row, col]
            title = rf'{title_suffix} — {col_labels[col]}'

            for blanket_key in ['HCPB', 'DCLL']:

                # Resolve per-blanket nuclides for N2N
                if qty == 'N2N':
                    nucs = N2N_NUCLIDES[blanket_key]
                else:
                    nucs = nuclides

                z_start = Z_BOUNDS[blanket_key][0]

                for case in CASES:
                    sp_path = find_statepoint(blanket_key, case, loading)
                    print(f"  Loading {blanket_key} {qty} | case {case} | "
                          f"{loading} kg/m³ → {sp_path}")

                    depth_mid, mean = extract_depth_profile(
                        sp_path, tally_name, score, nucs
                    )

                    # Convert absolute z to depth from first wall
                    depth = depth_mid - z_start

                    style = line_styles[(blanket_key, case)]
                    ax.step(depth, mean, where='mid', linewidth=1.2,
                            color=style['color'], label=style['label'])

            # ── Axis formatting ───────────────────────────────────────────
            ax.set_xlim(0, None)
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
            ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

            leg = ax.legend(title=title, title_fontsize=13, fontsize=9,
                            fancybox=False, edgecolor='black', frameon=True,
                            framealpha=0.75, ncol=1, loc='best')
            leg.get_frame().set_linewidth(0.5)

    # ── Shared axis labels ────────────────────────────────────────────────
    for ax in axes[3, :]:
        ax.set_xlabel('Depth from first wall [cm]', fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel(r'Reaction rate [rxn$/$src-n]', fontsize=14)

    fig.tight_layout()

    # ── Save ──────────────────────────────────────────────────────────────
    os.makedirs('./Figures/PDF', exist_ok=True)
    os.makedirs('./Figures/PNG', exist_ok=True)
    plt.savefig(f'./Figures/PDF/fig_depth_profile_{ISOTOPE}.pdf', bbox_inches='tight')
    plt.savefig(f'./Figures/PNG/fig_depth_profile_{ISOTOPE}.png', bbox_inches='tight')
    print(f"\nSaved to ./Figures/PDF/ and ./Figures/PNG/")
    plt.show()


if __name__ == '__main__':
    main()