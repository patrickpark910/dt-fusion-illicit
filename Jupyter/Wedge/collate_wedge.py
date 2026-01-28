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
#for L in loadings:
#    extract_tallies(case='C', fertile_isotope='U238', fertile_loading=L)
#    print(f"Extracted tallies for {L} kg/m³")  
# COLLATE TALLIES ###########################################
def collate_ng_tallies(case='C'):
    tally_folders = [x for x in os.listdir("./Jupyter/Wedge/OpenMC/") if (x.startswith(f"prism_{case}_U238_"))]
    print(tally_folders)
    ng_list = []
    fertile_isotope='U238'
    for folder in tally_folders:
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

    df_ng_collated.to_csv(f"./Figures/Data/prism_{case}_{fertile_isotope}_n-gamma_{900}K.csv",index=False)

# COLLATE FLUX TALLIES #######################################
def collate_flux_tallies(case='C'):
    tally_folders = [x for x in os.listdir("./Jupyter/Wedge/OpenMC/") if (x.startswith(f"prism_{case}_U238_"))]
    print(tally_folders)
    flux_list = []
    fertile_isotope = 'U238'
    for folder in tally_folders:
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
    tally_folders = [x for x in os.listdir("./Jupyter/Wedge/OpenMC/") if (x.startswith(f"prism_{case}_U238_"))]
    print(tally_folders)
    fission_list = []
    fertile_isotope='U238'
    for folder in tally_folders:
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
'''
def plot_cum_norm_histogram():
    print(f"\nPlotting cumulative, normalized fissile production vs. energy...")

    fig, ax = plt.subplots(figsize=(7,7)) # sharex='col', sharey='row'
    u_path = "./Figures/XSPlot/U238gamma.txt"
    u238_energy, u238_mxs = readtxtFile(u_path) # cross sections are returned as log(microxs)
    u238_mxs_shifted = (u238_mxs - np.min(u238_mxs)) * 0.8 / (np.max(u238_mxs) - np.min(u238_mxs)) + 0.1
    hcpb_u_ng_df = pd.read_csv('./Figures/Data/prism_C_U238_n-gamma_900K.csv') 

    title, df = r"HCPB-UO$_2$ Prism C", hcpb_u_ng_df

    ax.set_xlabel("Incident neutron energy [eV]")
    ax.set_ylabel(r"Cumulative fraction of fertile $($n,$\gamma)$ rxns")
    uranium, = ax.plot(u238_energy, u238_mxs_shifted, color='gray', linewidth=0.7, alpha=0.4, label=r'U238 (n, $\gamma$)')
    uranium_leg = ax.legend(handles=[uranium], loc='upper left', edgecolor='gray', frameon=True, framealpha=.75)
    ax.add_artist(uranium_leg)

    df = df[df["fertile_kg/m3"].isin([15,30,60,90,120,150,250,500,750,999.99])]
    df_mean = df.groupby("fertile_kg/m3")["mean"].sum().to_frame() # pd.DataFrame(, columns=["MT_fertile","sum"]) # 2-col df, colA = MT_fertile values, colB = sum of 'mean'
    df_mean.columns = ["sum"]   
    df = df.merge(df_mean, on="fertile_kg/m3", how="left") # adds sum as column to df
    df['norm_mean'] = df['mean'] / df['sum']
    
    edges = np.unique(np.concatenate([df['energy low [eV]'].values, df['energy high [eV]'].values]))
    bins  = np.sort(df['energy mid [eV]'].unique())
    lw=0.5

    print(df['fertile_kg/m3'].unique())

    sub   = df[df['fertile_kg/m3'] == 15]
    label = sub['fertile_kg/m3'].iloc[0]
    _, _, col15 = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='black', label=fr'{label} kg$/$m$^3$', lw=lw)

    sub   = df[df['fertile_kg/m3'] == 30]
    label = sub['fertile_kg/m3'].iloc[0]
    _, _, col30 = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='red', label=fr'{label} kg$/$m$^3$', lw=lw)

    sub   = df[df['fertile_kg/m3'] == 60]
    label = sub['fertile_kg/m3'].iloc[0]
    _, _, col60 = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='orange', label=fr'{label} kg$/$m$^3$', lw=lw)

    sub   = df[df['fertile_kg/m3'] == 90]
    label = sub['fertile_kg/m3'].iloc[0]
    _, _, col90 = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='yellow', label=fr'{label} kg$/$m$^3$',lw=lw)

    sub   = df[df['fertile_kg/m3'] == 120]
    label = sub['fertile_kg/m3'].iloc[0]
    _, _, col120 = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='green', label=fr'{label} kg$/$m$^3$',lw=lw)

    sub   = df[df['fertile_kg/m3'] == 150]
    label = sub['fertile_kg/m3'].iloc[0]
    _, _, col150 = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='blue', label=fr'{label} kg$/$m$^3$',lw=lw)

    sub   = df[df['fertile_kg/m3'] == 250]
    label = sub['fertile_kg/m3'].iloc[0]
    _, _, col250 = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='purple', label=fr'{label} kg$/$m$^3$',lw=lw)

    sub   = df[df['fertile_kg/m3'] == 500.0]
    label = sub['fertile_kg/m3'].iloc[0]
    _, _, col500 = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='deeppink', label=fr'{label} kg$/$m$^3$',lw=lw)

    sub   = df[df['fertile_kg/m3'] == 750.0]
    label = sub['fertile_kg/m3'].iloc[0]
    _, _, col750 = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='lightpink', label=fr'{label} kg$/$m$^3$',lw=lw)


    sub   = df[df['fertile_kg/m3'] == 999.99]
    label = 1000 # sub['fertile_kg/m3'].iloc[0]
    _, _, col1000 = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='aqua', label=fr'{label} kg$/$m$^3$',lw=lw)


    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
    ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

    ax.set_xscale('log')
    ax.set_xlim(0.67*1e0, 1.5*1e7)
    ax.set_ylim(-0.03, 1.03)
    fig.tight_layout()

    # black red orange green blue purple pink
    leg = ax.legend(handles=[col15[0],col30[0], col60[0], col90[0], col120[0], col150[0], col250[0], col500[0], col750[0], col1000[0]], 
                    title=title, fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=1, loc="lower right")
    leg.get_frame().set_linewidth(0.5) 

    plt.savefig(f'./Figures/pdf/fig_cum_norm_histogram_prism.pdf', bbox_inches='tight', format='pdf')
'''
def plot_cum_norm_histogram():
    print(f"\nPlotting cumulative, normalized fissile production vs. energy (Prism A, B, C)...")

# At the start of plot_cum_norm_histogram (or the loop over cases), replace the u_path block with:

    u_path = "./Figures/XSPlot/U238gamma.txt"
    uranium_handle = None
    if os.path.isfile(u_path):
        u238_energy, u238_mxs = readtxtFile(u_path)
        u238_mxs_shifted = (u238_mxs - np.min(u238_mxs)) * 0.8 / (np.max(u238_mxs) - np.min(u238_mxs)) + 0.1
    else:
        u238_energy = u238_mxs_shifted = None  # not used
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    loadings_wanted = [15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]
    colors = ['black', 'red', 'orange', 'yellow', 'green', 'blue', 'purple', 'deeppink', 'lightpink', 'aqua']
    lw = 0.5

    for i, case in enumerate(['A', 'B', 'C']):
        ax = axes[i]
        csv_path = f"./Figures/Data/prism_{case}_U238_n-gamma_900K.csv"
        if not os.path.isfile(csv_path):
            ax.set_title(f"Prism {case} (no data)")
            continue
        df = pd.read_csv(csv_path)
        title = f"Prism {case}"
        ax.set_xlabel("Incident neutron energy [eV]")
        ax.set_ylabel(r"Cumulative fraction of fertile $($n,$\gamma)$ rxns")
        if u238_energy is not None and u238_mxs_shifted is not None:
            uranium, = ax.plot(u238_energy, u238_mxs_shifted, color='gray', linewidth=0.7, alpha=0.4, label=r'U238 (n, $\gamma$)')
            ax.add_artist(ax.legend(handles=[uranium], loc='upper left', edgecolor='gray', frameon=True, framealpha=.75))

        df = df[df["fertile_kg/m3"].isin(loadings_wanted)]
        if df.empty:
            ax.set_title(title)
            continue
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
            _, _, patches = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color=color, label=fr'{label} kg$/$m$^3$', lw=lw)
            legend_handles.append(patches[0])

        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)
        ax.set_xscale('log')
        ax.set_xlim(0.67*1e0, 1.5*1e7)
        ax.set_ylim(-0.03, 1.03)
        ax.set_title(title)
        if legend_handles:
            leg = ax.legend(handles=legend_handles, title=title, fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=1, loc="lower right")
            leg.get_frame().set_linewidth(0.5)

    fig.tight_layout()
    os.makedirs("./Figures/pdf", exist_ok=True)
    plt.savefig(f'./Figures/pdf/fig_cum_norm_histogram_prism.pdf', bbox_inches='tight', format='pdf')
    plt.close()

def plot_flux_prism():
    """Log-log and lin-lin neutron flux spectrum for wedge Cases A, B, C U238. Saves to Figures/pdf and Figures/png."""
    print("\nPlotting wedge Cases neutron flux spectrum (log-log and lin-lin)...")

    def sci_notation(x, pos):
        return f"{x:.0e}"

    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    for i, case in enumerate(['A', 'B', 'C']):
        dst_csv = f"./Figures/Data/prism_{case}_U238_flux_900K.csv"
        if not os.path.isfile(dst_csv):
            print(f"Flux CSV not found: {dst_csv}. Run collate_flux_tallies() first.")
            plt.close()
            return
        df = pd.read_csv(dst_csv)
        loadings = sorted(df["fertile_kg/m3"].unique())
        if len(loadings) == 0:
            print(f"No loadings in flux CSV for {case}; skipping row.")
            continue
        ax_log = axes[i, 0]
        ax_lin = axes[i, 1]
        for loading in loadings:
            sub = df[df["fertile_kg/m3"] == loading].sort_values("energy mid [eV]")
            sub = sub[sub["mean"] > 0]
            if sub.empty:
                continue
            lbl = f"{int(loading) if loading == int(loading) else loading} kg/m³"
            ax_log.plot(sub["energy mid [eV]"], sub["mean"], label=lbl, linewidth=0.5, alpha=0.8)
            ax_lin.plot(sub["energy mid [eV]"], sub["mean"], label=lbl, linewidth=0.5, alpha=0.8)
        ax_log.set_xscale("log")
        ax_log.set_yscale("log")
        ax_log.set_xlim(1e0, 2e7)
        ax_log.set_ylim(1e-5, 1e2)
        ax_log.set_xlabel("Energy [eV]", fontsize=12)
        ax_log.set_ylabel("Flux [n/cm²-s]", fontsize=12)
        ax_log.set_title(f"Prism {case} flux (log-log)", fontsize=14)
        ax_log.legend(title="U238 [kg/m³]", fontsize=8, loc="best")
        ax_log.grid(True, which="both", alpha=0.3)

        ax_lin.set_xlim(1e5, 5e6)
        ax_lin.xaxis.set_major_formatter(FuncFormatter(sci_notation))
        ymax = df["mean"].max()
        ax_lin.set_ylim(0, max(0.1, ymax * 1.05) if ymax > 0 else 0.1)
        ax_lin.set_xlabel("Energy [eV]", fontsize=12)
        ax_lin.set_ylabel("Flux [n/cm²-s]", fontsize=12)
        ax_lin.set_title(rf"Prism {case} flux, lin-lin 0.1$-$5 MeV", fontsize=14)
        ax_lin.legend(title="U238 [kg/m³]", fontsize=8, loc="best")
        ax_lin.grid(True, alpha=0.3)

    plt.tight_layout()
    for d in ["./Figures/pdf", "./Figures/png"]:
        os.makedirs(d, exist_ok=True)
    plt.savefig("./Figures/pdf/fig_flux_prism.pdf", bbox_inches="tight", format="pdf")
    plt.savefig("./Figures/png/fig_flux_prism.png", bbox_inches="tight", format="png")
    print("Saved flux plots: ./Figures/pdf/fig_flux_prism.pdf, ./Figures/png/fig_flux_prism.png")
    plt.close()

def plot_wedge_c_fission_vs_loading(case='C'):
    """
    Plot total fertile fission reaction rate vs loading for wedge Prism {case}.

    Requires per-run files like:
      ./Jupyter/Wedge/OpenMC/prism_{case}_U238_030.00kgm3/U238_fission_Ebins.csv

    Uses total fission RR = sum(mean) over all cells and energy bins.
    """
    base = "./Jupyter/Wedge/OpenMC"
    if not os.path.isdir(base):
        print(f"Missing directory: {base} (run from repo root)")
        return

    # Find all prism {case} folders for this isotope
    folders = sorted(
        d for d in os.listdir(base)
        if d.startswith(f"prism_{case}_U238_") and d.endswith("kgm3")
    )

    xs, ys = [], []
    for folder in folders:
        try:
            loading = float(folder.split("_")[-1].replace("kgm3", ""))
        except Exception:
            continue

        fpath = os.path.join(base, folder, f"U238_fission_Ebins.csv")
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
    ax.set_title(f"Wedge Prism {case}: U238 fission RR vs loading")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=True)

    plt.tight_layout()
    os.makedirs("./Figures/pdf", exist_ok=True)
    os.makedirs("./Figures/png", exist_ok=True)
    plt.savefig(f"./Figures/pdf/fig_wedge_prism{case}_U238_fission_vs_loading.pdf",
                bbox_inches="tight", format="pdf")
    plt.savefig(f"./Figures/png/fig_wedge_prism{case}_U238_fission_vs_loading.png",
                bbox_inches="tight", format="png")
    plt.close("all")

    print("Saved fission-vs-loading plots to ./Figures/pdf and ./Figures/png")

for case in ['A', 'B', 'C']:
    collate_ng_tallies(case=case)
    print(f"Collated ng tallies for {case}")
    collate_flux_tallies(case=case)
    print(f"Collated flux tallies for {case}")
    collate_fission_tallies(case=case)
    print(f"Collated fission tallies for {case}")
    plot_wedge_c_fission_vs_loading(case=case)
    print("Plotted fission vs. loading")
plot_cum_norm_histogram()
print("Plotted cumulative, normalized fissile production vs. energy")
plot_flux_prism()
print("Plotted flux spectrum")


print("All plots completed and saved to ./Figures/pdf and ./Figures/png.")
