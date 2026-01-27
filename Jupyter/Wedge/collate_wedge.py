import os, sys
import openmc
import pandas as pd
import matplotlib.pyplot as plt

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

# PLOT WEDGE ################################################
def plot_cum_norm_histogram(self):
    """
    Plots cumulative, normalized Pu production vs. energy for contours of MTU.
    """
    print(f"\nPlotting cumulative, normalized fissile production vs. energy...")

    fig, axes = plt.subplots(3, 2, figsize=(15,15)) # sharex='col', sharey='row', 

    # load in thorium and uranium background (n, gamma) cross section, shift to 0.1 and 0.9
    u_path = "./Figures/XSPlot/U238gamma.txt"
    th_path = "./Figures/XSPlot/Th232gamma.txt"
    u238_energy, u238_mxs = readtxtFile(u_path) # cross sections are returned as log(microxs)
    th232_energy, th232_mxs = readtxtFile(th_path)

    u238_mxs_shifted = (u238_mxs - np.min(u238_mxs)) * 0.8 / (np.max(u238_mxs) - np.min(u238_mxs)) + 0.1
    th232_mxs_shifted = (th232_mxs - np.min(th232_mxs)) * 0.8 / (np.max(th232_mxs) - np.min(th232_mxs)) + 0.1

    titles = [r"FLiBe-UF$_4$", r"FLiBe-ThF$_4$", r"HCPB-UO$_2$", r"HCPB-ThO$_2$", r"DCLL-UO$_2$", r"DCLL-ThO$_2$",]
    dfs    = [self.flibe_u_ng_df, self.flibe_th_ng_df, self.hcpb_u_ng_df, self.hcpb_th_ng_df, self.dcll_u_ng_df, self.dcll_th_ng_df,]

    for i in range(3):
        for j in range(2):
            ax = axes[i, j]
            if i == 2:
                ax.set_xlabel("Incident neutron energy [eV]")
            if j == 0:
                ax.set_ylabel(r"Cumulative fraction of fertile $($n,$\gamma)$ rxns")
                uranium,  = ax.plot(u238_energy, u238_mxs_shifted, color='gray', linewidth=0.7, alpha=0.4, label=r'U238 (n, $\gamma$)')
                uranium_leg = ax.legend(handles=[uranium], loc='upper left', edgecolor='gray', frameon=True, framealpha=.75)
                ax.add_artist(uranium_leg)
            if j == 1: 
                thorium,  = ax.plot(th232_energy, th232_mxs_shifted, color='gray', linewidth=0.7, alpha=0.4, label=r'Th232 (n, $\gamma$)')
                thorium_leg = ax.legend(handles=[thorium], loc='upper left', edgecolor='gray', frameon=True, framealpha=.75)
                ax.add_artist(thorium_leg)

    for ax, df, title in zip(axes.flatten(), dfs, titles):

        # Filter out MT_fertile loadings we want to plot 
        df = df[df["fertile_kg/m3"].isin([30,60,90,120,150])]

        # Compute sum of 'mean' for each MT_fertile. 
        df_mean = df.groupby("fertile_kg/m3")["mean"].sum().to_frame() # pd.DataFrame(, columns=["MT_fertile","sum"]) # 2-col df, colA = MT_fertile values, colB = sum of 'mean'
        df_mean.columns = ["sum"]   
        df = df.merge(df_mean, on="fertile_kg/m3", how="left") # adds sum as column to df

        df['norm_mean'] = df['mean'] / df['sum']

        # df.to_csv('test.csv',index=False)
        

        # Create bin edges from low and high energy boundaries
        edges = np.unique(np.concatenate([df['energy low [eV]'].values, df['energy high [eV]'].values]))
        bins  = np.sort(df['energy mid [eV]'].unique())

        sub   = df[df['fertile_kg/m3'] == 30]
        label = int(round(sub['fertile_kg/m3'].iloc[0], -1))
        _, _, red = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#ff1f5b', label=fr'{label} kg$/$m$^3$')

        sub = df[df['fertile_kg/m3'] == 60]
        label = int(round(sub['fertile_kg/m3'].iloc[0], -1))
        _, _, orange= ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#f48628', label=fr'{label} kg$/$m$^3$')

        sub = df[df['fertile_kg/m3'] == 90]
        label = int(round(sub['fertile_kg/m3'].iloc[0], -1))
        _, _, green= ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#04cc6c', label=fr'{label} kg$/$m$^3$')

        sub = df[df['fertile_kg/m3'] == 120]
        label = int(round(sub['fertile_kg/m3'].iloc[0], -1))
        _, _, blue = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#0c9edd', label=fr'{label} kg$/$m$^3$')

        sub = df[df['fertile_kg/m3'] == 150]
        label = int(round(sub['fertile_kg/m3'].iloc[0], -1))
        _, _, purple = ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#b367bd', label=fr'{label} kg$/$m$^3$')
        
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

        ax.set_xscale('log')
        ax.set_xlim(0.67*1e0, 1.5*1e7)
        ax.set_ylim(-0.03, 1.03)
        fig.tight_layout()

        leg = ax.legend(handles=[red[0], orange[0], green[0], blue[0], purple[0]], title=title, fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=1, loc="lower right")
        leg.get_frame().set_linewidth(0.5) 

    if self.save:
        plt.savefig(f'./Figures/pdf/fig_cum_norm_histogram_new.pdf', bbox_inches='tight', format='pdf')
        plt.savefig(f'./Figures/png/fig_cum_norm_histogram_new.png', bbox_inches='tight', format='png')
        print(f"Exported cumulative normalized histogram plots.")
    else:
        print(f"Did not export cumulative normalized histogram plot.")

    if self.show: plt.show()
    plt.close('all')