import os, sys
import openmc
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
from matplotlib.ticker import MultipleLocator, ScalarFormatter, FuncFormatter
import numpy as np


current_dir = os.getcwd()
print(current_dir)

def extract_tallies(case='C', fertile_isotope='U238', fertile_loading=30): 
    """ Load tallies """
    start_dir = os.getcwd()
    path = os.path.join(current_dir, "Jupyter/Wedge/OpenMC", f"prism_{case}_{fertile_isotope}_{fertile_loading:06.2f}kgm3")
    os.chdir(path)
    try:
        sp_path = os.path.join(path, "statepoint.100.h5")
        if not os.path.isfile(sp_path):
            sp_path = os.path.join(path, "statepoint.10.h5")
        print(f"Loading statepoint: {sp_path}")
        
        try:
            sp = openmc.StatePoint(sp_path) 
            print([t.name for t in sp.tallies.values()])
        except Exception as e:
            print(f"\n{e}\n")
            sys.exit(f"can't read the sp")

        print(f"Reading tallies...")

        # fertile      = sp.get_tally(name=f'Fertile rxn rates').get_pandas_dataframe()
        fertile_spec = sp.get_tally(name=f'Fertile rxn rates spectrum').get_pandas_dataframe() # or U
        flux_spec    = sp.get_tally(name=f'flux spectrum').get_pandas_dataframe()

        # Add new column for energy bin midpoint (for plotting)
        for df in [fertile_spec, flux_spec]:
            df['energy mid [eV]'] = (df['energy low [eV]'] + df['energy high [eV]'])/ 2

        U238_ng_Ebin_df  = fertile_spec[(fertile_spec['nuclide'] == fertile_isotope) & (fertile_spec['score'] == '(n,gamma)')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        U238_ng_Ebin_df.to_csv(f'{path}/{fertile_isotope}_n-gamma_Ebins.csv', index=False)

        U238_fis_Ebin_df = fertile_spec[(fertile_spec['nuclide'] == fertile_isotope) & (fertile_spec['score'] == 'fission')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        U238_fis_Ebin_df.to_csv(f'{path}/{fertile_isotope}_fission_Ebins.csv', index=False)

        # Flux in blanket cells only (c13, c23 → ids 13, 23)
        flux_Ebin_df = flux_spec[flux_spec['cell'].isin([13, 23])][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']].copy()
        flux_Ebin_df.to_csv(f'{path}/{fertile_isotope}_flux_Ebins.csv', index=False)

        print(f"Extracted tallies for {fertile_loading} kg/m³")

    finally:
        os.chdir(start_dir)
#############################################################
# 30, 60, 90, 120, 150
loadings = [0.1, 0.5, 1.5, 15.0, 30.0, 60.0, 90.0, 120.0, 150.0, 250.0, 500.0, 750.0, 999.99]

#for L in loadings:
#    extract_tallies(case='A', fertile_isotope='U238', fertile_loading=L)
#    print(f"Extracted tallies for {L} kg/m³")   
#for L in loadings:
#    extract_tallies(case='B', fertile_isotope='U238', fertile_loading=L)
#    print(f"Extracted tallies for {L} kg/m³")  
#for fertile_isotope in ['U238', 'Th232']:
#    for L in loadings:
#        extract_tallies(case='C', fertile_isotope=fertile_isotope, fertile_loading=L)
#        print(f"Extracted tallies for {L} kg/m³")  
# COLLATE TALLIES ###########################################
# ... existing code ...

# COLLATE TALLIES ###########################################
def collate_ng_tallies(case='C'):
    # Fix: Search for both U238 and Th232 folders
    tally_folders_u238 = [x for x in os.listdir("./Jupyter/Wedge/OpenMC/") if (x.startswith(f"prism_{case}_U238_"))]
    tally_folders_th232 = [x for x in os.listdir("./Jupyter/Wedge/OpenMC/") if (x.startswith(f"prism_{case}_Th232_"))]
    tally_folders = sorted(set(tally_folders_u238 + tally_folders_th232))
    print(tally_folders)
    ng_list = []
    fertile_isotopes= ['U238', 'Th232']
    for fertile_isotope in fertile_isotopes:
        # Filter folders for this isotope
        isotope_folders = [f for f in tally_folders if f.startswith(f"prism_{case}_{fertile_isotope}_")]
        for folder in isotope_folders:
            part = folder.split("_")[-1]
            fertile = float(part.replace("kgm3", ""))
            
            tally_ng      = f"./Jupyter/Wedge/OpenMC/{folder}/{fertile_isotope}_n-gamma_Ebins.csv"

            try:
                df_ng = pd.read_csv(tally_ng)
            except FileNotFoundError:
                print(f"File '{fertile_isotope}_n-gamma_Ebins.csv' not found in {folder}, skipping...")
                continue

            cols = ["energy low [eV]", "energy high [eV]", "energy mid [eV]", "mean"]  
            ng =  (df_ng[cols].groupby("energy mid [eV]", as_index=False)
                                .agg(**{"energy low [eV]" : ("energy low [eV]", "first"),
                                        "energy high [eV]": ("energy high [eV]", "first"),
                                        "mean": ("mean", "sum"),} ) )
            
            ng['filename'],  ng['fertile_kg/m3'] = folder, fertile
            ng_list.append(ng)
            

        df_ng_collated = pd.concat(ng_list, ignore_index=True) if ng_list else pd.DataFrame()
        df_ng_collated = df_ng_collated[['filename', 'fertile_kg/m3',
                                        'energy mid [eV]', 'mean', 'energy low [eV]', 'energy high [eV]']]

        # Fix: Remove braces around 900
        df_ng_collated.to_csv(f"./Figures/Data/prism_{case}_{fertile_isotope}_n-gamma_900K.csv",index=False)

# COLLATE FLUX TALLIES #######################################
def collate_flux_tallies(case='C'):
    # Fix: Search for both U238 and Th232 folders
    tally_folders_u238 = [x for x in os.listdir("./Jupyter/Wedge/OpenMC/") if (x.startswith(f"prism_{case}_U238_"))]
    tally_folders_th232 = [x for x in os.listdir("./Jupyter/Wedge/OpenMC/") if (x.startswith(f"prism_{case}_Th232_"))]
    tally_folders = sorted(set(tally_folders_u238 + tally_folders_th232))
    print(tally_folders)
    flux_list = []
    fertile_isotopes = ['U238', 'Th232']
    for fertile_isotope in fertile_isotopes:
        # Filter folders for this isotope
        isotope_folders = [f for f in tally_folders if f.startswith(f"prism_{case}_{fertile_isotope}_")]
        for folder in isotope_folders:
            part = folder.split("_")[-1]
            fertile = float(part.replace("kgm3", ""))
            tally_flux = f"./Jupyter/Wedge/OpenMC/{folder}/{fertile_isotope}_flux_Ebins.csv"
            try:
                df_flux = pd.read_csv(tally_flux)
            except FileNotFoundError:
                print(f"File '{fertile_isotope}_flux_Ebins.csv' not found in {folder}, skipping...")
                continue
            cols = ["energy low [eV]", "energy high [eV]", "energy mid [eV]", "mean"]
            flux = (df_flux[cols].groupby("energy mid [eV]", as_index=False)
                        .agg(**{"energy low [eV]": ("energy low [eV]", "first"),
                                "energy high [eV]": ("energy high [eV]", "first"),
                                "mean": ("mean", "sum")}))
            flux['filename'] = folder
            flux['fertile_kg/m3'] = fertile
            flux_list.append(flux)
        df_flux_collated = pd.concat(flux_list, ignore_index=True) if flux_list else pd.DataFrame()
        if not df_flux_collated.empty:
            df_flux_collated = df_flux_collated[['filename', 'fertile_kg/m3', 'energy mid [eV]', 'mean', 'energy low [eV]', 'energy high [eV]']]
            dst = f"./Figures/Data/prism_{case}_{fertile_isotope}_flux_900K.csv"
            os.makedirs("./Figures/Data", exist_ok=True)
            df_flux_collated.to_csv(dst, index=False)
            print(f"Collated wedge {case} flux tallies to: {dst}")
        else:
            print("No flux Ebins files found; nothing written.")

# COLLATE TALLIES ###########################################
def collate_fission_tallies(case='C'):
    # Fix: Search for both U238 and Th232 folders
    tally_folders_u238 = [x for x in os.listdir("./Jupyter/Wedge/OpenMC/") if (x.startswith(f"prism_{case}_U238_"))]
    tally_folders_th232 = [x for x in os.listdir("./Jupyter/Wedge/OpenMC/") if (x.startswith(f"prism_{case}_Th232_"))]
    tally_folders = sorted(set(tally_folders_u238 + tally_folders_th232))
    print(tally_folders)
    fission_list = []
    fertile_isotopes= ['U238', 'Th232']
    for fertile_isotope in fertile_isotopes:
        # Filter folders for this isotope
        isotope_folders = [f for f in tally_folders if f.startswith(f"prism_{case}_{fertile_isotope}_")]
        for folder in isotope_folders:
            part = folder.split("_")[-1]
            fertile = float(part.replace("kgm3", ""))
            
            tally_fission      = f"./Jupyter/Wedge/OpenMC/{folder}/{fertile_isotope}_fission_Ebins.csv"

            try:
                df_fission = pd.read_csv(tally_fission)
            except FileNotFoundError:
                print(f"File '{fertile_isotope}_fission_Ebins.csv' not found in {folder}, skipping...")
                continue

            cols = ["energy low [eV]", "energy high [eV]", "energy mid [eV]", "mean"]  
            fission =  (df_fission[cols].groupby("energy mid [eV]", as_index=False)
                                .agg(**{"energy low [eV]" : ("energy low [eV]", "first"),
                                        "energy high [eV]": ("energy high [eV]", "first"),
                                        "mean": ("mean", "sum"),} ) )
            
            fission['filename'],  fission['fertile_kg/m3'] = folder, fertile
            fission_list.append(fission)
            
        df_fission_collated = pd.concat(fission_list, ignore_index=True) if fission_list else pd.DataFrame()
        df_fission_collated = df_fission_collated[['filename', 'fertile_kg/m3',
                                        'energy mid [eV]', 'mean', 'energy low [eV]', 'energy high [eV]']]

        df_fission_collated.to_csv(f"./Figures/Data/prism_{case}_{fertile_isotope}_fission_900K.csv",index=False)


# PLOT ######################################################
def readtxtFile(path): 
    energy, microxs = [], []

    with open(path, 'r') as file:
        file.readline()
        file.readline()

        for line in file:
            values = line.split()
            energy.append(float(values[0]))
            microxs.append(float(values[1]))

    return np.array(energy), np.log(np.array(microxs))


def plot_cum_norm_histogram(fertile_isotope='U238'):
    """
    Left: Prism C cum-norm (n,gamma) vs energy
    Right: HCPB cum-norm (n,gamma) vs energy
    """
    print(f"\nPlotting cumulative, normalized fissile production vs. energy (Prism C vs HCPB) for {fertile_isotope}...")

    # Optional U-238 reference curve
    u_path = "./Figures/XSPlot/U238gamma.txt"
    if os.path.isfile(u_path):
        u238_energy, u238_mxs = readtxtFile(u_path)
        u238_mxs_shifted = (u238_mxs - np.min(u238_mxs)) * 0.8 / (np.max(u238_mxs) - np.min(u238_mxs)) + 0.1
    else:
        u238_energy = u238_mxs_shifted = None  # not used

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    loadings_wanted = [15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]
    colors = ['black', 'red', 'orange', 'yellow', 'green', 'blue',
              'purple', 'deeppink', 'lightpink', 'aqua']
    lw = 0.5

    def _plot_one(ax, csv_path, title):
        if not os.path.isfile(csv_path):
            ax.set_title(f"{title} (no data)")
            return

        df = pd.read_csv(csv_path)
        ax.set_xlabel("Incident neutron energy [eV]")
        ax.set_ylabel(r"Cumulative fraction of fertile $($n,$\gamma)$ rxns")
        ax.set_title(title)

        # Optional U-238 reference (only if file exists)
        if u238_energy is not None and u238_mxs_shifted is not None:
            uranium, = ax.plot(
                u238_energy,
                u238_mxs_shifted,
                color='gray',
                linewidth=0.7,
                alpha=0.4,
                label=r'U238 (n, $\gamma$)'
            )
            ax.add_artist(
                ax.legend(
                    handles=[uranium],
                    loc='upper left',
                    edgecolor='gray',
                    frameon=True,
                    framealpha=.75
                )
            )

        df = df[df["fertile_kg/m3"].isin(loadings_wanted)]
        if df.empty:
            return

        # Normalize by total (n,gamma) per loading
        df_mean = df.groupby("fertile_kg/m3")["mean"].sum().to_frame()
        df_mean.columns = ["sum"]
        df = df.merge(df_mean, on="fertile_kg/m3", how="left")
        df['norm_mean'] = df['mean'] / df['sum']

        bins = np.sort(df['energy mid [eV]'].unique())
        legend_handles = []

        for load, color in zip(loadings_wanted, colors):
            sub = df[df['fertile_kg/m3'] == load]
            if sub.empty:
                continue
            label = 1000 if load == 999.99 else load
            _, _, patches = ax.hist(
                sub['energy mid [eV]'],
                bins=bins,
                weights=sub['norm_mean'],
                cumulative=True,
                histtype='step',
                color=color,
                label=fr'{label} kg$/$m$^3$',
                lw=lw
            )
            legend_handles.append(patches[0])

        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)
        ax.set_xscale('log')
        ax.set_xlim(0.67 * 1e0, 1.5 * 1e7)
        ax.set_ylim(-0.03, 1.03)

        if legend_handles:
            leg = ax.legend(
                handles=legend_handles,
                title=title,
                fancybox=False,
                edgecolor='black',
                frameon=True,
                framealpha=.75,
                ncol=1,
                loc="lower right"
            )
            leg.get_frame().set_linewidth(0.5)

    # Left: Prism C
    prism_csv = f"./Figures/Data/prism_C_{fertile_isotope}_n-gamma_900K.csv"
    _plot_one(axes[0], prism_csv, "Prism C")

    # Right: HCPB – choose U vs Th filename based on isotope
    hcpb_iso = "U" if fertile_isotope.startswith("U") else "Th"
    hcpb_csv = f"./Figures/Data/HCPB_{hcpb_iso}_n-gamma_900K.csv"
    _plot_one(axes[1], hcpb_csv, "HCPB")

    fig.tight_layout()
    os.makedirs("./Figures/pdf", exist_ok=True)
    out = f'./Figures/pdf/fig_cum_norm_histogram_{fertile_isotope}.pdf'
    plt.savefig(out, bbox_inches='tight', format='pdf')
    print(f"Saved: {out}")
    plt.close()


def plot_flux_prism(fertile_isotope='U238'):
    """Log-log and lin-lin neutron flux spectrum for Prism C and HCPB. Saves to Figures/pdf and Figures/png."""
    print("\nPlotting flux spectrum (log-log and lin-lin) for Prism C and HCPB...")

    def sci_notation(x, pos):
        return f"{x:.0e}"

    # 2x2 grid: row 0 = Prism C (log, lin), row 1 = HCPB (log, lin)
    fig, axes = plt.subplots(2, 2, figsize=(20, 12))
    
    # Determine HCPB isotope name for CSV
    hcpb_iso = "U238" if fertile_isotope == "U238" else "Th232"
    
    cases = [
        ("C", f"./Figures/Data/prism_C_{fertile_isotope}_flux_900K.csv", "Prism C"),
        ("HCPB", f"./Figures/Data/HCPB_{hcpb_iso}_flux_900K.csv", "HCPB")
    ]
    
    for row_idx, (case_name, csv_path, title_prefix) in enumerate(cases):
        if not os.path.isfile(csv_path):
            print(f"Flux CSV not found: {csv_path}. Skipping {case_name}.")
            axes[row_idx, 0].text(0.5, 0.5, f"No data: {case_name}", ha="center", va="center", transform=axes[row_idx, 0].transAxes)
            axes[row_idx, 1].text(0.5, 0.5, f"No data: {case_name}", ha="center", va="center", transform=axes[row_idx, 1].transAxes)
            continue
            
        df = pd.read_csv(csv_path)
        loadings = sorted(df["fertile_kg/m3"].unique())
        if len(loadings) == 0:
            print(f"No loadings in flux CSV for {case_name}; skipping.")
            continue
            
        ax_log = axes[row_idx, 0]
        ax_lin = axes[row_idx, 1]
        
        for loading in loadings:
            sub = df[df["fertile_kg/m3"] == loading].sort_values("energy mid [eV]")
            sub = sub[sub["mean"] > 0]
            if sub.empty:
                continue
            lbl = f"{int(loading) if loading == int(loading) else loading} kg/m³"
            ax_log.plot(sub["energy mid [eV]"], sub["mean"], label=lbl, linewidth=0.5, alpha=0.8)
            ax_lin.plot(sub["energy mid [eV]"], sub["mean"], label=lbl, linewidth=0.5, alpha=0.8)
        
        # Log-log subplot
        ax_log.set_xscale("log")
        ax_log.set_yscale("log")
        ax_log.set_xlim(1e0, 2e7)
        ax_log.set_ylim(1e-4, 1e1)
        ax_log.set_xlabel("Energy [eV]", fontsize=12)
        ax_log.set_ylabel("Flux [n/cm²-s]", fontsize=12)
        ax_log.set_title(f"{title_prefix} flux (log-log)", fontsize=14)
        ax_log.legend(title=f"{fertile_isotope} [kg/m³]", fontsize=8, loc="best")
        ax_log.grid(True, which="both", alpha=0.3)
        
        # Lin-lin subplot
        ax_lin.set_xlim(1e5, 5e6)
        ax_lin.xaxis.set_major_formatter(FuncFormatter(sci_notation))
        ymax = df["mean"].max()
        ax_lin.set_ylim(0, 1)
        ax_lin.set_xlabel("Energy [eV]", fontsize=12)
        ax_lin.set_ylabel("Flux [n/cm²-s]", fontsize=12)
        ax_lin.set_title(rf"{title_prefix} flux, lin-lin 0.1$-$5 MeV", fontsize=14)
        ax_lin.legend(title=f"{fertile_isotope} [kg/m³]", fontsize=8, loc="best")
        ax_lin.grid(True, alpha=0.3)

    plt.tight_layout()
    for d in ["./Figures/pdf", "./Figures/png"]:
        os.makedirs(d, exist_ok=True)
    plt.savefig(f"./Figures/pdf/fig_flux_prismC_vs_HCPB_{fertile_isotope}.pdf", bbox_inches="tight", format="pdf")
    plt.savefig(f"./Figures/png/fig_flux_prismC_vs_HCPB_{fertile_isotope}.png", bbox_inches="tight", format="png")
    print(f"Saved flux plots: ./Figures/pdf/fig_flux_prismC_vs_HCPB_{fertile_isotope}.pdf, ./Figures/png/fig_flux_prismC_vs_HCPB_{fertile_isotope}.png")
    plt.close()

def plot_wedge_c_fission_vs_loading(case='C', fertile_isotope='U238'):
    """
    Plot total fertile fission reaction rate vs loading for wedge Prism {case}.

    Requires per-run files like:
      ./Jupyter/Wedge/OpenMC/prism_{case}_{fertile_isotope}_030.00kgm3/{fertile_isotope}_fission_Ebins.csv

    Uses total fission RR = sum(mean) over all cells and energy bins.
    """
    base = "./Jupyter/Wedge/OpenMC"
    if not os.path.isdir(base):
        print(f"Missing directory: {base} (run from repo root)")
        return

    # Find all prism {case} folders for this isotope
    folders = sorted(
        d for d in os.listdir(base)
        if d.startswith(f"prism_{case}_{fertile_isotope}_") and d.endswith("kgm3")
    )

    xs, ys = [], []
    for folder in folders:
        try:
            loading = float(folder.split("_")[-1].replace("kgm3", ""))
        except Exception:
            continue

        fpath = os.path.join(base, folder, f"{fertile_isotope}_fission_Ebins.csv")
        if not os.path.isfile(fpath):
            print(f"Missing fission Ebins: {fpath} (skip)")
            continue

        df = pd.read_csv(fpath)
        if "mean" not in df.columns:
            print(f"No 'mean' column in: {fpath} (skip)")
            continue

        total_fis = df["mean"].sum()
        xs.append(loading)
        ys.append(total_fis)

    if len(xs) < 2:
        print("Not enough points to plot/interpolate (need >= 2 loadings with fission Ebins).")
        return

    # Sort by loading
    xs = np.array(xs, dtype=float)
    ys = np.array(ys, dtype=float)
    order = np.argsort(xs)
    xs, ys = xs[order], ys[order]

    # Simple interpolation (PCHIP is shape-preserving; Akima is fine too)
    x_fine = np.linspace(xs.min(), xs.max(), 300)
    y_smooth = PchipInterpolator(xs, ys)(x_fine)

    plt.figure(figsize=(7.5, 5))
    ax = plt.gca()

    ax.scatter(xs, ys, s=35, color="black", label="OpenMC points")
    ax.plot(x_fine, y_smooth, linewidth=1.5, color="red", label="PCHIP interp")

    ax.set_xlabel(r"Fertile isotope density in blanket [kg/m$^3$]")
    ax.set_ylabel("Total fission reaction rate [per source particle]")
    ax.set_title(f"Wedge Prism {case}: {fertile_isotope} fission RR vs loading")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=True)

    plt.tight_layout()
    os.makedirs("./Figures/pdf", exist_ok=True)
    os.makedirs("./Figures/png", exist_ok=True)
    plt.savefig(f"./Figures/pdf/fig_wedge_prism{case}_{fertile_isotope}_fission_vs_loading.pdf",
                bbox_inches="tight", format="pdf")
    plt.savefig(f"./Figures/png/fig_wedge_prism{case}_{fertile_isotope}_fission_vs_loading.png",
                bbox_inches="tight", format="png")
    plt.close("all")

    print("Saved fission-vs-loading plots to ./Figures/pdf and ./Figures/png")

def wedge_c_fissile_per_yr(fertile_isotope='U238'):
    """
    Compute fissile [kg/yr] per loading for wedge case C from collated (n,gamma) CSV.
    Same conversion as reactor.py: total (n,gamma) per source neutron -> kg/yr.
    """
    AVO = 6.022141076e23
    SEC_PER_YR = 3600 * 24 * 365
    AMU_PU239 = 239.0521634
    AMU_U233 = 233.039635207
    NPS_FUS = 1000 * 3.546e17  # n/s
    case = 'C'
    csv_path = f"./Figures/Data/prism_{case}_{fertile_isotope}_n-gamma_900K.csv"
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"Run collate_ng_tallies(case='{case}') first; need {csv_path}")
    df = pd.read_csv(csv_path)
    # Total (n,gamma) per source neutron per loading [sum over cells and energy]
    total_ng = df.groupby("fertile_kg/m3", as_index=False)["mean"].sum()
    total_ng.columns = ["fertile_kg/m3", "total_ng_per_srcn"]
    if fertile_isotope == "U238":
        total_ng["fissile_kg/yr"] = (
            total_ng["total_ng_per_srcn"] * NPS_FUS * SEC_PER_YR * AMU_PU239 / AVO / 1e3
        )
    elif fertile_isotope == "Th232":
        total_ng["fissile_kg/yr"] = (
            total_ng["total_ng_per_srcn"] * NPS_FUS * SEC_PER_YR * AMU_U233 / AVO / 1e3
        )
    else:
        raise ValueError("fertile_isotope must be U238 or Th232")
    return total_ng[["fertile_kg/m3", "fissile_kg/yr"]]

def plot_fudge_factor(fertile_isotope='U238'):
    """
    Plot 3 subplots:
    1. Wedge C fissile production [kg/yr] vs loading
    2. HCPB fissile production [kg/yr] vs loading
    3. Factor (HCPB / Wedge C) vs loading
    """
    print(f"\nPlotting fissile production and fudge factor for {fertile_isotope}...")
    
    # Wedge C: from collated n-gamma
    wedge_df = wedge_c_fissile_per_yr(fertile_isotope=fertile_isotope)
    x_w = wedge_df["fertile_kg/m3"].values
    y_w = wedge_df["fissile_kg/yr"].values
    
    # HCPB: premade rxns CSV with fissile kg/yr column
    if fertile_isotope == "U238":
        hcpb_path = "./Figures/Data/HCPB_U238_rxns_900K.csv"
        col_kg_yr = "Pu239_kg/yr"
    else:
        hcpb_path = "./Figures/Data/HCPB_Th232_rxns_900K.csv"
        col_kg_yr = "U233_kg/yr"
    
    if not os.path.isfile(hcpb_path):
        print(f"HCPB rxns not found: {hcpb_path}")
        return
    
    hcpb_df = pd.read_csv(hcpb_path)
    x_h = hcpb_df["fertile_kg/m3"].values
    y_h = hcpb_df[col_kg_yr].values
    
    # Common x: union of loadings, interpolate both at common points for factor
    x_common = np.unique(np.concatenate([x_w, x_h]))
    x_common = x_common[(x_common >= min(x_w.min(), x_h.min())) & 
                         (x_common <= max(x_w.max(), x_h.max()))]
    x_common = np.sort(x_common)

    def _prep_xy(x, y):
        d = pd.DataFrame({"x": x, "y": y}).dropna()
        d["x"] = d["x"].astype(float)
        d["y"] = d["y"].astype(float)
        # if duplicates exist, collapse them (mean is fine; can also use first)
        d = d.groupby("x", as_index=False)["y"].mean()
        d = d.sort_values("x")
        return d["x"].to_numpy(), d["y"].to_numpy()

    # after you define x_w, y_w, x_h, y_h:
    x_w, y_w = _prep_xy(x_w, y_w)
    x_h, y_h = _prep_xy(x_h, y_h)

    # align by loading (x) so we can drop the same x from both
    w = pd.DataFrame({"x": x_w, "y_w": y_w})
    h = pd.DataFrame({"x": x_h, "y_h": y_h})
    m = w.merge(h, on="x", how="outer")  # keep all x's from either set

    def _keep_not_outlier(series, z_thresh=6.0):
        s = series.dropna().astype(float).to_numpy()
        if s.size < 3:
            return pd.Series(True, index=series.index)  # not enough points to judge
        v = np.log10(np.clip(series.astype(float).to_numpy(), 1e-30, None))
        med = np.nanmedian(v)
        mad = np.nanmedian(np.abs(v - med))
        if mad == 0 or np.isnan(mad):
            return pd.Series(True, index=series.index)
        z = 0.6745 * (v - med) / mad
        return pd.Series(np.abs(z) <= z_thresh, index=series.index)

    # if wedge OR hcpb is an outlier at that x, drop that x from BOTH
    keep = _keep_not_outlier(m["y_w"]) & _keep_not_outlier(m["y_h"])
    m = m[keep].sort_values("x")

    # back to arrays (and drop any x where a series is missing)
    w_ok = m.dropna(subset=["y_w"])
    h_ok = m.dropna(subset=["y_h"])
    x_w, y_w = w_ok["x"].to_numpy(), w_ok["y_w"].to_numpy()
    x_h, y_h = h_ok["x"].to_numpy(), h_ok["y_h"].to_numpy()

    if len(x_w) < 2 or len(x_h) < 2:
        print("Not enough points left after outlier removal.")
        return

    interp_w = PchipInterpolator(x_w, y_w)
    interp_h = PchipInterpolator(x_h, y_h)

    # Interpolate both at common x points

    y_w_common = interp_w(x_common)
    y_h_common = interp_h(x_common)
    
    # Factor at each x (avoid div by zero)
    with np.errstate(divide="ignore", invalid="ignore"):
        factor = np.where(y_w_common > 0, y_h_common / y_w_common, np.nan)
    
    # Create 3 subplots
    fig, axes = plt.subplots(3, 1, figsize=(8, 12), sharex=True)
    ax1, ax2, ax3 = axes[0], axes[1], axes[2]
    
    # Subplot 1: Wedge C fissile per year
    ax1.scatter(x_w, y_w, s=40, color="C0", label="Wedge Prism C", zorder=3)
    x_fine_w = np.linspace(x_w.min(), x_w.max(), 200)
    ax1.plot(x_fine_w, interp_w(x_fine_w), "-", color="C0", linewidth=1.5, alpha=0.7, label="PCHIP interp")
    ax1.set_ylabel("Fissile production [kg/yr]", fontsize=12)
    ax1.set_title(f"Wedge Prism C fissile production ({fertile_isotope})", fontsize=14)
    ax1.legend(frameon=True, loc="best")
    ax1.grid(True, alpha=0.3)
    
    # Subplot 2: HCPB fissile per year
    ax2.scatter(x_h, y_h, s=40, color="C1", label="HCPB", zorder=3)
    x_fine_h = np.linspace(x_h.min(), x_h.max(), 200)
    ax2.plot(x_fine_h, interp_h(x_fine_h), "-", color="C1", linewidth=1.5, alpha=0.7, label="PCHIP interp")
    ax2.set_ylabel("Fissile production [kg/yr]", fontsize=12)
    ax2.set_title(f"HCPB fissile production ({fertile_isotope})", fontsize=14)
    ax2.legend(frameon=True, loc="best")
    ax2.grid(True, alpha=0.3)
    
    # Subplot 3: Factor (HCPB / Wedge C)
    valid_mask = ~np.isnan(factor)
    if np.any(valid_mask):
        ax3.scatter(x_common[valid_mask], factor[valid_mask], s=40, color="black", 
                   label="HCPB / Wedge C", zorder=3)
        # Interpolate factor for smooth line
        factor_interp = PchipInterpolator(x_common[valid_mask], factor[valid_mask])
        x_fine_factor = np.linspace(x_common[valid_mask].min(), x_common[valid_mask].max(), 200)
        ax3.plot(x_fine_factor, factor_interp(x_fine_factor), "-", color="red", 
                linewidth=1.5, alpha=0.7, label="PCHIP interp")
    ax3.axhline(1.0, color="gray", linestyle="--", linewidth=0.8, label="Factor = 1")
    ax3.set_xlabel(r"Fertile isotope density [kg/m$^3$]", fontsize=12)
    ax3.set_ylabel("Factor (HCPB / Wedge C)", fontsize=12)
    ax3.set_title(f"Fudge factor: HCPB / Wedge C ({fertile_isotope})", fontsize=14)
    ax3.legend(frameon=True, loc="best")
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    os.makedirs("./Figures/pdf", exist_ok=True)
    os.makedirs("./Figures/png", exist_ok=True)
    base = f"fig_fudge_factor_prismC_vs_HCPB_{fertile_isotope}"
    plt.savefig(f"./Figures/pdf/{base}.pdf", bbox_inches="tight", format="pdf")
    plt.savefig(f"./Figures/png/{base}.png", bbox_inches="tight", format="png")
    plt.close()
    print(f"Saved {base}.pdf/.png")

for case in ['C']:
    collate_ng_tallies(case=case)
    print(f"Collated ng tallies for {case}")
    collate_flux_tallies(case=case)
    print(f"Collated flux tallies for {case}")
    collate_fission_tallies(case=case)
    print(f"Collated fission tallies for {case}")
    for fertile_isotope in ['U238', 'Th232']:
        plot_wedge_c_fission_vs_loading(case=case, fertile_isotope=fertile_isotope)
        print("Plotted fission vs. loading")
        wedge_c_fissile_per_yr(fertile_isotope=fertile_isotope)
        print("Computed fissile production per year for Wedge C")

for fertile_isotope in ['U238', 'Th232']:
    plot_cum_norm_histogram(fertile_isotope=fertile_isotope)
    print("Plotted cumulative, normalized fissile production vs. energy")
    plot_flux_prism(fertile_isotope=fertile_isotope)
    print("Plotted flux spectrum")
    plot_fudge_factor(fertile_isotope=fertile_isotope)
    print("Plotted fudge factor")



print("All plots completed and saved to ./Figures/pdf and ./Figures/png.")
