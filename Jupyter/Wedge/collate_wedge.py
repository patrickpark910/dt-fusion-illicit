import os, sys
import openmc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter, FuncFormatter
import numpy as np

current_dir = os.getcwd()
print(current_dir)

def extract_tallies(case='C', fertile_isotope='U238', fertile_loading=30): 
    """ Load tallies """
    path = os.path.join(current_dir, "Jupyter/Wedge/OpenMC", f"prism_{case}_{fertile_isotope}_{fertile_loading:06.2f}kgm3")
    os.chdir(path)

    sp_path = f'{path}/statepoint.100.h5'
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

    # Add new column for energy bin midpoint (for plotting)
    for df in [fertile_spec]:
        df['energy mid [eV]'] = (df['energy low [eV]'] + df['energy high [eV]'])/ 2

    U238_ng_Ebin_df  = fertile_spec[(fertile_spec['nuclide'] == fertile_isotope) & (fertile_spec['score'] == '(n,gamma)')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]
    U238_ng_Ebin_df.to_csv(f'{path}/{fertile_isotope}_n-gamma_Ebins.csv', index=False)

#############################################################
# 30, 60, 90, 120, 150
# extract_tallies(fertile_loading=0.1)


# COLLATE TALLIES ###########################################
def collate_ng_tallies():
    tally_folders = [x for x in os.listdir("./Jupyter/Wedge/OpenMC/") if (x.startswith(f"prism_C_U238_"))]
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

    df_ng_collated.to_csv(f"./Figures/Data/prism_C_{fertile_isotope}_n-gamma_{900}K.csv",index=False)

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

def plot_cum_norm_histogram():
    print(f"\nPlotting cumulative, normalized fissile production vs. energy...")

    fig, ax = plt.subplots(figsize=(7,7)) # sharex='col', sharey='row'
    u_path = "./Figures/XSPlot/U238gamma.txt"
    u238_energy, u238_mxs = readtxtFile(u_path) # cross sections are returned as log(microxs)
    u238_mxs_shifted = (u238_mxs - np.min(u238_mxs)) * 0.8 / (np.max(u238_mxs) - np.min(u238_mxs)) + 0.1
    hcpb_u_ng_df = pd.read_csv('./Figures/data/prism_C_U238_n-gamma_900K.csv') 

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

    sub   = df[df['fertile_kg/m3'] == 999.99]
    label = 1000 # sub['fertile_kg/m3'].iloc[0]
    _, _, col1000 = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='lightpink', label=fr'{label} kg$/$m$^3$',lw=lw)


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
    leg = ax.legend(handles=[col15[0],col30[0], col60[0], col90[0], col120[0], col150[0], col250[0], col500[0], col1000[0]], 
                    title=title, fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=1, loc="lower right")
    leg.get_frame().set_linewidth(0.5) 

    plt.savefig(f'./Figures/pdf/fig_cum_norm_histogram_prism.pdf', bbox_inches='tight', format='pdf')

plot_cum_norm_histogram()