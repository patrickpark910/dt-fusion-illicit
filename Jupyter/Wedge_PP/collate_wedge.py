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
            
            tally_ng  = f"./Jupyter/Wedge/OpenMC/{folder}/{fertile_isotope}_n-gamma_Ebins.csv"

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


# COLLATE TBR TALLIES #######################################
def collate_tbr_tallies(case='C'):
    """
    Extract TBR (Li6(n,Xt) + Li7(n,Xt) per source neutron) from each wedge run;
    one row per loading (x value). Collate to Figures/Data/prism_{case}_{isotope}_tbr_900K.csv.
    """
    base = "./Jupyter/Wedge/OpenMC"
    if not os.path.isdir(base):
        print(f"Missing directory: {base}")
        return
    tally_folders_u238 = [x for x in os.listdir(base) if x.startswith(f"prism_{case}_U238_")]
    tally_folders_th232 = [x for x in os.listdir(base) if x.startswith(f"prism_{case}_Th232_")]
    os.makedirs("./Figures/Data", exist_ok=True)
    for fertile_isotope, folders in [("U238", tally_folders_u238), ("Th232", tally_folders_th232)]:
        tbr_list = []
        for folder in sorted(folders):
            path = os.path.join(base, folder)
            try:
                loading = float(folder.split("_")[-1].replace("kgm3", ""))
            except Exception:
                continue
            sp_path = os.path.join(path, "statepoint.100.h5")
            if not os.path.isfile(sp_path):
                sp_path = os.path.join(path, "statepoint.10.h5")
            if not os.path.isfile(sp_path):
                print(f"No statepoint in {folder}, skipping TBR")
                continue
            try:
                sp = openmc.StatePoint(sp_path)
                li_tally = sp.get_tally(name="Total Li rxn rate")
                df = li_tally.get_pandas_dataframe()
            except Exception as e:
                print(f"Could not load TBR from {folder}: {e}")
                continue
            if not all(c in df.columns for c in ["nuclide", "score", "mean"]):
                print(f"Tally in {folder} missing required columns, skip")
                continue
            # TBR = sum of Li6 (n,Xt) + Li7 (n,Xt) per source neutron
            mask = (df["nuclide"].isin(["Li6", "Li7"])) & (df["score"] == "(n,Xt)")
            tbr = float(df.loc[mask, "mean"].sum())
            tbr_list.append({"fertile_kg/m3": loading, "tbr": tbr})
        if not tbr_list:
            print(f"No TBR data for prism {case} {fertile_isotope}")
            continue
        df_tbr = pd.DataFrame(tbr_list).sort_values("fertile_kg/m3")
        dst = f"./Figures/Data/prism_{case}_{fertile_isotope}_tbr_900K.csv"
        df_tbr.to_csv(dst, index=False)
        print(f"Collated wedge {case} TBR for {fertile_isotope} to {dst} ({len(df_tbr)} points)")


# PLOTINGGGG ######################################################
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

    # Simple interpolation 
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


def _prep_xy(x, y):
    d = pd.DataFrame({"x": x, "y": y}).dropna()
    d["x"] = d["x"].astype(float)
    d["y"] = d["y"].astype(float)
    d = d.groupby("x", as_index=False)["y"].mean()
    d = d.sort_values("x")
    return d["x"].to_numpy(), d["y"].to_numpy()


def _keep_not_outlier(series, z_thresh=6.0):
    s = series.dropna().astype(float).to_numpy()
    if s.size < 3:
        return pd.Series(True, index=series.index)
    v = np.log10(np.clip(series.astype(float).to_numpy(), 1e-30, None))
    med = np.nanmedian(v)
    mad = np.nanmedian(np.abs(v - med))
    if mad == 0 or np.isnan(mad):
        return pd.Series(True, index=series.index)
    z = 0.6745 * (v - med) / mad
    return pd.Series(np.abs(z) <= z_thresh, index=series.index)


def plot_fissile_factor():
    """
    Single plot: Fissile production [kg/yr] (left axis) for HCPB and Wedge, U238 and Th232;
    Factor Wedge/HCPB (right axis). HCPB red, Wedge orange; solid U238, dashed Th232. PCHIP.
    """
    print("\nPlotting fissile production and Wedge/HCPB factor (U238 + Th232)...")
    HCPB_RED = '#b41f24'
    WEDGE_ORANGE = '#ff6600'
    data = {}
    for iso in ['U238', 'Th232']:
        col_kg_yr = "Pu239_kg/yr" if iso == "U238" else "U233_kg/yr"
        hcpb_path = f"./Figures/Data/HCPB_{iso}_rxns_900K.csv"
        if not os.path.isfile(hcpb_path):
            print(f"HCPB rxns not found: {hcpb_path}, skip {iso}")
            continue
        try:
            wedge_df = wedge_c_fissile_per_yr(fertile_isotope=iso)
        except FileNotFoundError:
            print(f"Wedge n-gamma not found for {iso}, skip")
            continue
        x_w = wedge_df["fertile_kg/m3"].values
        y_w = wedge_df["fissile_kg/yr"].values
        hcpb_df = pd.read_csv(hcpb_path)
        x_h = hcpb_df["fertile_kg/m3"].values
        y_h = hcpb_df[col_kg_yr].values
        x_w, y_w = _prep_xy(x_w, y_w)
        x_h, y_h = _prep_xy(x_h, y_h)
        w = pd.DataFrame({"x": x_w, "y_w": y_w})
        h = pd.DataFrame({"x": x_h, "y_h": y_h})
        m = w.merge(h, on="x", how="outer")
        keep = _keep_not_outlier(m["y_w"]) & _keep_not_outlier(m["y_h"])
        m = m[keep].sort_values("x")
        w_ok = m.dropna(subset=["y_w"])
        h_ok = m.dropna(subset=["y_h"])
        x_w = w_ok["x"].to_numpy()
        y_w = w_ok["y_w"].to_numpy()
        x_h = h_ok["x"].to_numpy()
        y_h = h_ok["y_h"].to_numpy()
        if len(x_w) < 2 or len(x_h) < 2:
            continue
        for _x, _y in [(x_w, y_w), (x_h, y_h)]:
            o = np.argsort(_x)
            _x[:], _y[:] = _x[o], _y[o]
        x_common = np.unique(np.concatenate([x_w, x_h]))
        x_common = x_common[(x_common >= min(x_w.min(), x_h.min())) & (x_common <= max(x_w.max(), x_h.max()))]
        x_common = np.sort(x_common)
        interp_w = PchipInterpolator(x_w, y_w)
        interp_h = PchipInterpolator(x_h, y_h)
        y_w_c = interp_w(x_common)
        y_h_c = interp_h(x_common)
        with np.errstate(divide="ignore", invalid="ignore"):
            factor = np.where(y_h_c > 0, y_w_c / y_h_c, np.nan)
        data[iso] = dict(x_w=x_w, y_w=y_w, x_h=x_h, y_h=y_h, interp_w=interp_w, interp_h=interp_h,
                         x_common=x_common, factor=factor)
                
    if not data:
        print("No fissile data to plot.")
        return

    # Save factor to CSV (computed points)
    rows = []
    for iso in data:
        d = data[iso]
        v = ~np.isnan(d["factor"])
        if np.any(v):
            for x, f in zip(d["x_common"][v], d["factor"][v]):
                rows.append({"fertile_kg/m3": float(x), "isotope": iso, "factor": float(f)})
    if rows:
        os.makedirs("./Figures/Data", exist_ok=True)
        pd.DataFrame(rows).sort_values(["isotope", "fertile_kg/m3"]).to_csv(
            "./Figures/Data/wedge_HCPB_fissile_factor_900K.csv", index=False
        )
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax2 = ax.twinx()
    x_min, x_max = np.inf, -np.inf
    y_max = -np.inf
    for iso, ls, mk_h, mk_w in [
        ('U238', '-', 's', 'o'),
        ('Th232', '--', '1', '+'),
    ]:
        if iso not in data:
            continue
        d = data[iso]
        x_w, y_w, x_h, y_h = d['x_w'], d['y_w'], d['x_h'], d['y_h']
        x_min = min(x_min, x_w.min(), x_h.min())
        x_max = max(x_max, x_w.max(), x_h.max())
        y_max = max(y_max, y_w.max(), y_h.max())
        ax.scatter(x_h, y_h, marker=mk_h, s=30, color=HCPB_RED, zorder=3)
        ax.scatter(x_w, y_w, marker=mk_w, s=30 if mk_w == 'o' else 60, color=WEDGE_ORANGE, zorder=3)
        x_f = np.linspace(d['x_common'].min(), d['x_common'].max(), 200)
        ax.plot(x_f, d['interp_h'](x_f), ls, linewidth=1, color=HCPB_RED)
        ax.plot(x_f, d['interp_w'](x_f), ls, linewidth=1, color=WEDGE_ORANGE)
    for iso, ls in [('U238', '-'), ('Th232', '--')]:
        if iso not in data:
            continue
        d = data[iso]
        v = ~np.isnan(d['factor'])
        if not np.any(v):
            continue
        x_c = d['x_common'][v]
        f = d['factor'][v]
        ax2.scatter(x_c, f, s=20, color='black', zorder=3, alpha=0.7)
        fi = PchipInterpolator(x_c, f)
        x_f = np.linspace(x_c.min(), x_c.max(), 200)
        ax2.plot(x_f, fi(x_f), ls, linewidth=0.5, color='black')
    ax2.axhline(1.0, color='gray', linestyle=':', linewidth=0.8)
    ax.set_xlabel(r'Fertile isotope density in blanket [kg/m$^3$]')
    ax.set_ylabel('Initial fissile production rate [kg/yr]')
    ax2.set_ylabel('Factor (Wedge / HCPB)')
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_label_position('right')
    ax2.set_ylim(0, 1.9)
    ax.set_xlim(max(0, x_min - 25), x_max + 25)
    ax.set_ylim(0, y_max * 1.05)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(50))
    ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
    ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

    # Legends: use ONLY dummy handles (do not label real curves; avoid axis autoscale to 1e9)
    from matplotlib.lines import Line2D

    left_handles = [
        Line2D([], [], color=HCPB_RED, linestyle='-', marker='s', linewidth=1, markersize=5, label=r'HCPB-U$^{238}$'),
        Line2D([], [], color=HCPB_RED, linestyle='--', marker='1', linewidth=1, markersize=9, label=r'HCPB-Th$^{232}$'),
        Line2D([], [], color=WEDGE_ORANGE, linestyle='-', marker='o', linewidth=1, markersize=5, label=r'Wedge-U$^{238}$'),
        Line2D([], [], color=WEDGE_ORANGE, linestyle='--', marker='+', linewidth=1, markersize=8, label=r'Wedge-Th$^{232}$'),
    ]
    right_handles = [
        Line2D([], [], color='black', linestyle='-', linewidth=1, label=r'W/HCPB-U$^{238}$'),
        Line2D([], [], color='black', linestyle='--', linewidth=1, label=r'W/HCPB-Th$^{232}$'),
    ]

    ax.legend(handles=left_handles, fancybox=False, edgecolor='black', frameon=True, framealpha=0.75, loc='lower right')
    ax2.legend(handles=right_handles, fancybox=False, edgecolor='black', frameon=True, framealpha=0.75, loc='upper left')
    for leg in [ax.get_legend(), ax2.get_legend()]:
        if leg is not None:
            leg.get_frame().set_linewidth(0.5)
    plt.tight_layout()
    os.makedirs("./Figures/pdf", exist_ok=True)
    os.makedirs("./Figures/png", exist_ok=True)
    base = "fig_fissile_factor_wedge_vs_HCPB"
    plt.savefig(f"./Figures/pdf/{base}.pdf", bbox_inches="tight", format="pdf")
    plt.savefig(f"./Figures/png/{base}.png", bbox_inches="tight", format="png")
    plt.close()
    print(f"Saved {base}.pdf/.png")


def plot_tbr_factor():
    """
    Single plot: TBR (left axis) for HCPB and Wedge, U238 and Th232; Factor Wedge/HCPB (right axis).
    HCPB red, Wedge orange; solid U238, dashed Th232. PCHIP interpolation.
    Uses one point per x (via _prep_xy); drops at most one factor-outlier x per isotope, then plots the rest.
    Legends use Line2D handles only (no labels on plot artists).
    """
    print("\nPlotting TBR and Wedge/HCPB factor (U238 + Th232)...")
    case = 'C'
    HCPB_RED = '#b41f24'
    WEDGE_ORANGE = '#ff6600'
    data = {}
    for iso in ['U238', 'Th232']:
        wp = f"./Figures/Data/prism_{case}_{iso}_tbr_900K.csv"
        hp = f"./Figures/Data/HCPB_{iso}_rxns_900K.csv"
        if not os.path.isfile(wp) or not os.path.isfile(hp):
            print(f"Skip {iso}: missing {wp} or {hp}")
            continue
        df_w = pd.read_csv(wp)
        df_h = pd.read_csv(hp)
        if "tbr" not in df_w.columns or "tbr" not in df_h.columns:
            print(f"Skip {iso}: missing 'tbr' column in wedge or HCPB CSV")
            continue
        # One point per x, sorted (proper tallies per loading)
        x_w, y_w = _prep_xy(df_w["fertile_kg/m3"].values, df_w["tbr"].values)
        x_h, y_h = _prep_xy(df_h["fertile_kg/m3"].values, df_h["tbr"].values)
        w = pd.DataFrame({"x": x_w, "y_w": y_w})
        h = pd.DataFrame({"x": x_h, "y_h": y_h})
        m = w.merge(h, on="x", how="outer").sort_values("x")

        # Drop AT MOST ONE big factor-outlier x (plot the rest; factor curve excludes this x)
        m_both = m.dropna(subset=["y_w", "y_h"]).copy()
        if len(m_both) >= 3:
            f_raw = (m_both["y_w"].astype(float) / m_both["y_h"].astype(float)).to_numpy()
            v = np.log10(np.clip(f_raw, 1e-30, None))
            med = np.nanmedian(v)
            mad = np.nanmedian(np.abs(v - med))
            if mad != 0 and not np.isnan(mad):
                z = 0.6745 * (v - med) / mad
                j = int(np.nanargmax(np.abs(z)))
                if np.abs(z[j]) > 8.0:
                    x_drop = float(m_both["x"].to_numpy()[j])
                    print(f"Drop single TBR factor outlier at x={x_drop:g} kg/m3 for {iso}; plot rest.")
                    m = m[m["x"] != x_drop]
        w_ok = m.dropna(subset=["y_w"])
        h_ok = m.dropna(subset=["y_h"])
        x_w = w_ok["x"].to_numpy().astype(float)
        y_w = w_ok["y_w"].to_numpy().astype(float)
        x_h = h_ok["x"].to_numpy().astype(float)
        y_h = h_ok["y_h"].to_numpy().astype(float)
        if len(x_w) < 2 or len(x_h) < 2:
            continue
        for _x, _y in [(x_w, y_w), (x_h, y_h)]:
            o = np.argsort(_x)
            _x[:], _y[:] = _x[o], _y[o]
        x_common = np.unique(np.concatenate([x_w, x_h]))
        x_common = x_common[(x_common >= min(x_w.min(), x_h.min())) & (x_common <= max(x_w.max(), x_h.max()))]
        x_common = np.sort(x_common)
        interp_w = PchipInterpolator(x_w, y_w)
        interp_h = PchipInterpolator(x_h, y_h)
        y_w_c = interp_w(x_common)
        y_h_c = interp_h(x_common)
        with np.errstate(divide="ignore", invalid="ignore"):
            factor = np.where(y_h_c > 0, y_w_c / y_h_c, np.nan)
        data[iso] = dict(x_w=x_w, y_w=y_w, x_h=x_h, y_h=y_h, interp_w=interp_w, interp_h=interp_h,
                         x_common=x_common, factor=factor)
                
    if not data:
        print("No TBR data to plot.")
        return

    # Save factor to CSV (computed points)
    rows = []
    for iso in data:
        d = data[iso]
        v = ~np.isnan(d["factor"])
        if np.any(v):
            for x, f in zip(d["x_common"][v], d["factor"][v]):
                rows.append({"fertile_kg/m3": float(x), "isotope": iso, "factor": float(f)})
    if rows:
        os.makedirs("./Figures/Data", exist_ok=True)
        pd.DataFrame(rows).sort_values(["isotope", "fertile_kg/m3"]).to_csv(
            "./Figures/Data/wedge_HCPB_TBR_factor_900K.csv", index=False
        )
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax2 = ax.twinx()
    x_min, x_max = np.inf, -np.inf
    tbr_min, tbr_max = np.inf, -np.inf
    # Left axis: HCPB (red), Wedge (orange); solid U238, dashed Th232
    for iso, ls, mk_h, mk_w in [
        ('U238', '-', 's', 'o'),
        ('Th232', '--', '1', '+'),
    ]:
        if iso not in data:
            continue
        d = data[iso]
        x_w, y_w, x_h, y_h = d['x_w'], d['y_w'], d['x_h'], d['y_h']
        x_min, x_max = min(x_min, x_w.min(), x_h.min()), max(x_max, x_w.max(), x_h.max())
        tbr_min = min(tbr_min, y_w.min(), y_h.min())
        tbr_max = max(tbr_max, y_w.max(), y_h.max())
        ax.scatter(x_h, y_h, marker=mk_h, s=30, color=HCPB_RED, zorder=3)
        ax.scatter(x_w, y_w, marker=mk_w, s=30 if mk_w == 'o' else 60, color=WEDGE_ORANGE, zorder=3)
        x_f = np.linspace(d['x_common'].min(), d['x_common'].max(), 200)
        ax.plot(x_f, d['interp_h'](x_f), ls, linewidth=1, color=HCPB_RED)
        ax.plot(x_f, d['interp_w'](x_f), ls, linewidth=1, color=WEDGE_ORANGE)
    # Right axis: Factor (Wedge/HCPB) for U238 and Th232
    for iso, ls in [('U238', '-'), ('Th232', '--')]:
        if iso not in data:
            continue
        d = data[iso]
        v = ~np.isnan(d['factor'])
        if not np.any(v):
            continue
        x_c = d['x_common'][v]
        f = d['factor'][v]
        ax2.scatter(x_c, f, s=20, color='black', zorder=3, alpha=0.7) 
        fi = PchipInterpolator(x_c, f)
        x_f = np.linspace(x_c.min(), x_c.max(), 200)
        ax2.plot(x_f, fi(x_f), ls, linewidth=0.6, color='black')
    ax2.axhline(1.0, color='gray', linestyle=':', linewidth=0.8)
    ax.set_xlabel(r'Fertile isotope density in blanket [kg/m$^3$]')
    ax.set_ylabel('Tritium breeding ratio')
    ax2.set_ylabel('Factor (Wedge / HCPB)')
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_label_position('right')
    ax.set_xlim(max(0, x_min - 25), x_max + 25)
    ax.set_ylim(max(0, tbr_min - 0.02), tbr_max + 0.02)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(50))
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
    ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

    # Legends: only Line2D handles (no label= on scatter/plot; same as fissile factor plot)
    from matplotlib.lines import Line2D
    left_handles = [
        Line2D([], [], color=HCPB_RED, linestyle='-', marker='s', linewidth=1, markersize=5, label=r'HCPB-U$^{238}$'),
        Line2D([], [], color=HCPB_RED, linestyle='--', marker='1', linewidth=1, markersize=9, label=r'HCPB-Th$^{232}$'),
        Line2D([], [], color=WEDGE_ORANGE, linestyle='-', marker='o', linewidth=1, markersize=5, label=r'Wedge-U$^{238}$'),
        Line2D([], [], color=WEDGE_ORANGE, linestyle='--', marker='+', linewidth=1, markersize=8, label=r'Wedge-Th$^{232}$'),
    ]
    right_handles = [
        Line2D([], [], color='black', linestyle='-', linewidth=1, label=r'W/HCPB-U$^{238}$'),
        Line2D([], [], color='black', linestyle='--', linewidth=1, label=r'W/HCPB-Th$^{232}$'),
    ]
    ax.legend(handles=left_handles, fancybox=False, edgecolor='black', frameon=True, framealpha=0.75, loc='upper right')
    ax2.legend(handles=right_handles, fancybox=False, edgecolor='black', frameon=True, framealpha=0.75, loc='center right')
    for leg in [ax.get_legend(), ax2.get_legend()]:
        if leg is not None:
            leg.get_frame().set_linewidth(0.5)
    plt.tight_layout()
    os.makedirs("./Figures/pdf", exist_ok=True)
    os.makedirs("./Figures/png", exist_ok=True)
    base = "fig_tbr_prismC_vs_HCPB"
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
    collate_tbr_tallies(case=case)
    print(f"Collated TBR tallies for {case}")
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
plot_fissile_factor()
print("Plotted fissile factor (U238 + Th232)")
plot_tbr_factor()
print("Plotted TBR factor (U238 + Th232)")



print("All plots completed and saved to ./Figures/pdf and ./Figures/png.")
