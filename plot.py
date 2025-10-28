import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from scipy.interpolate import Akima1DInterpolator
from scipy.interpolate import PchipInterpolator


""" import helper functions """
from Python.utilities import *
from Python.parameters import *

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


def plot_all():
    """ Read CSV data into pandas DataFrames """

    combined_plot = Plot(save=True, show=True)
    
    """ Pick which one to plot by uncommenting """
    # combined_plot.plot_tbr()
    # combined_plot.plot_pu_per_yr()
    combined_plot.plot_cum_norm_histogram()

    print("All plots completed and saved.")


class Plot:

    def __init__(self, save=False, show=True):

        # Dataframes of reaction rates (rr)
        self.flibe_u_rr_df  = pd.read_csv('./Figures/Data/FLiBe_U_rxns_900K.csv')
        self.flibe_th_rr_df = pd.read_csv('./Figures/Data/FLiBe_Th_rxns_900K.csv')
        self.pbli_u_rr_df   = pd.read_csv('./Figures/Data/LL_U_rxns_900K.csv')
        self.pbli_th_rr_df  = pd.read_csv('./Figures/Data/LL_Th_rxns_900K.csv')
        self.pb_u_rr_df     = pd.read_csv('./Figures/Data/PB_U_rxns_900K.csv')
        self.pb_th_rr_df    = pd.read_csv('./Figures/Data/PB_Th_rxns_900K.csv')

        self.flibe_u_eb_df  = pd.read_csv('./Figures/Data/FLiBe_U_n-gamma_900K.csv')
        self.flibe_th_eb_df = pd.read_csv('./Figures/Data/FLiBe_Th_n-gamma_900K.csv')
        self.pbli_u_eb_df   = pd.read_csv('./Figures/Data/LL_U_n-gamma_900K.csv')
        self.pbli_th_eb_df  = pd.read_csv('./Figures/Data/LL_Th_n-gamma_900K.csv')
        self.pb_u_eb_df     = pd.read_csv('./Figures/data/PB_U_n-gamma_900K.csv') 
        self.pb_th_eb_df    = pd.read_csv('./Figures/data/PB_Th_n-gamma_900K.csv') 

        self.save, self.show = save, show
        # self.name = 'All_Blankets'

        for sd in ['PDF','PNG','Data']:
            sd_path = f'./Figures/{sd}/'
            print(f"Ensuring directory exists: {sd_path}")
            os.makedirs(sd_path, exist_ok=True)


    def plot_tbr(self):
        print(f"\nPlotting tritium breeding ratio vs. fertile density...")

        # Setting up x, y separately here so you can remove the impurity/wppm-magnitude cases --ppark 2025-08-06
        x6, y6 =     self.pb_u_rr_df['fertile_kg/m3'],     self.pb_u_rr_df['Li6(n,t)'] +     self.pb_u_rr_df['Li7(n,Xt)']
        x5, y5 =    self.pb_th_rr_df['fertile_kg/m3'],    self.pb_th_rr_df['Li6(n,t)'] +    self.pb_th_rr_df['Li7(n,Xt)']     
        x4, y4 =   self.pbli_u_rr_df['fertile_kg/m3'],   self.pbli_u_rr_df['Li6(n,t)'] +   self.pbli_u_rr_df['Li7(n,Xt)']
        x3, y3 =  self.pbli_th_rr_df['fertile_kg/m3'],  self.pbli_th_rr_df['Li6(n,t)'] +  self.pbli_th_rr_df['Li7(n,Xt)']
        x2, y2 =  self.flibe_u_rr_df['fertile_kg/m3'],  self.flibe_u_rr_df['Li6(n,t)'] +  self.flibe_u_rr_df['Li7(n,Xt)']
        x1, y1 = self.flibe_th_rr_df['fertile_kg/m3'], self.flibe_th_rr_df['Li6(n,t)'] + self.flibe_th_rr_df['Li7(n,Xt)']

        # AkimaInterpolation ripped from my ELWR paper --ppark 2025-08-06
        x_fine = np.linspace(x1.min(), x1.max(), 300) # Evaluate on a fine grid
        y_akima1 = Akima1DInterpolator(x1, y1)(x_fine)
        y_akima2 = Akima1DInterpolator(x2, y2)(x_fine)
        y_akima3 = Akima1DInterpolator(x3, y3)(x_fine)
        y_akima4 = Akima1DInterpolator(x4, y4)(x_fine)
        y_akima5 = Akima1DInterpolator(x5, y5)(x_fine)
        y_akima6 = Akima1DInterpolator(x6, y6)(x_fine)

        for y in [y_akima1,y_akima2,y_akima3,y_akima4,y_akima5,y_akima6]: # y_akima5,
            for i in range(1, len(y)):
                if y[i] > y[i-1]:
                    y[i] = y[i-1] # adjust the current point to ensure it's not greater than the previous point


        plt.figure(figsize=(7.5,5))

        plt.scatter(x6, y6, marker='s', s=40, color='#b41f24') #  PB-UO2
        plt.scatter(x5, y5, marker='1', s=70, color='#b41f24') #  PB-ThO2
        plt.scatter(x4, y4, marker='^', s=40, color='#0047ba') #  LL-UO2
        plt.scatter(x3, y3, marker='x', s=50, color='#0047ba') #  LL-ThO2
        plt.scatter(x2, y2, marker='o', s=30, color='#66b420')   #  FLiBe-UF4
        plt.scatter(x1, y1, marker='+', s=60, color='#66b420')   #  FLiBe-ThF4 

        plt.plot(x_fine, y_akima6, '-',   linewidth=1, color='#b41f24',)   #  PB-UO2
        plt.plot(x_fine, y_akima5, '--',  linewidth=1, color='#b41f24',)   #  PB-ThO2
        plt.plot(x_fine, y_akima4, '-',   linewidth=1, color='#0047ba',)   #  LL-UO2
        plt.plot(x_fine, y_akima3, '--',  linewidth=1, color='#0047ba',)   #  LL-ThO2
        plt.plot(x_fine, y_akima2, '-',   linewidth=1, color='#66b420', )    #  FLiBe-UF4
        plt.plot(x_fine, y_akima1, '--',  linewidth=1, color='#66b420', )    #  FLiBe-ThF4 

        # Dummy plots for legend -- bit of a hack lmao -- ppark
        plt.plot([9e8,9e9], [9e8,9e9], '^-',  markersize=6, linewidth=1, color='#0047ba', label=r'LL-UO$_2$')     #  blue: #0047ba
        plt.plot([9e8,9e9], [9e8,9e9], 'x--', markersize=7, linewidth=1, color='#0047ba', label=r'LL-ThO$_2$')    #  
        plt.plot([9e8,9e9], [9e8,9e9], 's-',  markersize=5, linewidth=1, color='#b41f24', label=r'PB-UO$_2$')     #  red: #b41f24
        plt.plot([9e8,9e9], [9e8,9e9], '1--', markersize=9, linewidth=1, color='#b41f24', label=r'PB-ThO$_2$')    # 
        plt.plot([9e8,9e9], [9e8,9e9], 'o-',  markersize=5, linewidth=1, color='#66b420',   label=r'FLiBe-UF$_4$')  #  green: #66b420
        plt.plot([9e8,9e9], [9e8,9e9], '+--', markersize=8, linewidth=1, color='#66b420',   label=r'FLiBe-ThF$_4$') #  
        
        
        # plt.title(f'Tritium breeding ratio (All)') # Exclude title for production figs --ppark 2025-08-06
        plt.xlabel(r'Fertile isotope density in blanket [kg$/$m$^3$]')
        plt.ylabel('Tritium breeding ratio')

        plt.xlim(-5,165)
        plt.ylim(1.03,1.47) # plt.ylim(0.98,1.42)
        
        # Tick grid
        ax = plt.gca()
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)
        # ax.axhspan(.9,1, color='#e0ded8') # color='#e0ded8', # alpha=0.3)
        ax.axhspan(1,1.05, color='#e0ded8') # color='#f7f7f6' # , alpha=0.3)

        plt.tight_layout()

        leg = plt.legend(fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=3)
        leg.get_frame().set_linewidth(0.5) 

        if self.save:
            plt.savefig(f'./Figures/pdf/fig_tbr.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/fig_tbr.png', bbox_inches='tight', format='png')
            print(f"Exported tritium breeding ratio plots.")
        else:
            print(f"Did not export tritium breeding ratio plots.")

        if self.show: plt.show()
        plt.close('all')


        return 0
    

    def plot_pu_per_yr(self):
      
        print(f"\nPlotting Pu-239 production per year for all fuels...")

        # <utilities.py> Default NPS_FUS = 500 MJ/s * 3.546e17 n/MJ = 1.773e+20 n/s
        # seconds in a year = 3.156e+7 s/yr

        tot_n_per_yr = NPS_FUS * 3.156e+7

        # Gotta then convert fissile atoms/yr to kg/yr
        pu239_at_to_kg_per_yr = tot_n_per_yr / AVO * AMU_PU239 / 1e3
        u233_at_to_kg_per_yr  = tot_n_per_yr / AVO * AMU_U233  / 1e3
        
        # Setting up x, y separately here so you can remove the impurity/wppm-magnitude cases --ppark 2025-08-06
        x6, y6 =     self.pb_u_rr_df['fertile_kg/m3'], pu239_at_to_kg_per_yr *     self.pb_u_rr_df['U238(n,g)'] 
        x5, y5 =    self.pb_th_rr_df['fertile_kg/m3'],  u233_at_to_kg_per_yr *    self.pb_th_rr_df['Th232(n,g)']
        x4, y4 =   self.pbli_u_rr_df['fertile_kg/m3'], pu239_at_to_kg_per_yr *   self.pbli_u_rr_df['U238(n,g)'] 
        x3, y3 =  self.pbli_th_rr_df['fertile_kg/m3'],  u233_at_to_kg_per_yr *  self.pbli_th_rr_df['Th232(n,g)']
        x2, y2 =  self.flibe_u_rr_df['fertile_kg/m3'], pu239_at_to_kg_per_yr *  self.flibe_u_rr_df['U238(n,g)']
        x1, y1 = self.flibe_th_rr_df['fertile_kg/m3'],  u233_at_to_kg_per_yr * self.flibe_th_rr_df['Th232(n,g)']

        # AkimaInterpolation ripped from my ELWR paper --ppark 2025-08-06
        x_fine = np.linspace(x1.min(), x1.max(), 300) # Evaluate on a fine grid
        y_akima1 = Akima1DInterpolator(x1, y1)(x_fine)
        y_akima2 = Akima1DInterpolator(x2, y2)(x_fine)
        y_akima3 = Akima1DInterpolator(x3, y3)(x_fine)
        y_akima4 = Akima1DInterpolator(x4, y4)(x_fine)
        y_akima5 = Akima1DInterpolator(x5, y5)(x_fine)
        y_akima6 = Akima1DInterpolator(x6, y6)(x_fine)

        plt.figure(figsize=(7.5,5))

        plt.scatter(x6, y6, marker='s', s=40, color='#b41f24') # PB-UO2
        plt.scatter(x5, y5, marker='1', s=70, color='#b41f24') # PB-ThO2
        plt.scatter(x4, y4, marker='^', s=40, color='#0047ba') # LL-UO2
        plt.scatter(x3, y3, marker='x', s=50, color='#0047ba') # LL-ThO2
        plt.scatter(x2, y2, marker='o', s=30, color='#66b420')   # FLiBe-UF4
        plt.scatter(x1, y1, marker='+', s=60, color='#66b420')   # FLiBe-ThF4

        plt.plot(x_fine, y_akima6, '-',   linewidth=1, color='#b41f24',)   #  PB-UO2
        plt.plot(x_fine, y_akima5, '--',  linewidth=1, color='#b41f24',)   #  PB-ThO2
        plt.plot(x_fine, y_akima4, '-',   linewidth=1, color='#0047ba',)   #  LL-UO2
        plt.plot(x_fine, y_akima3, '--',  linewidth=1, color='#0047ba',)   #  LL-ThO2
        plt.plot(x_fine, y_akima2, '-',   linewidth=1, color='#66b420', )    #  FLiBe-UF4
        plt.plot(x_fine, y_akima1, '--',  linewidth=1, color='#66b420', )    #  FLiBe-ThF4 

        # Dummy plots for legend -- bit of a hack lmao -- ppark
        plt.plot([9e8,9e9], [9e8,9e9], '+--', markersize=8, linewidth=1, color='#66b420',   label=r'FLiBe-ThF$_4$') #  green: #66b420
        plt.plot([9e8,9e9], [9e8,9e9], 'o-',  markersize=5, linewidth=1, color='#66b420',   label=r'FLiBe-UF$_4$')  # 
        plt.plot([9e8,9e9], [9e8,9e9], '1--', markersize=9, linewidth=1, color='#b41f24', label=r'PB-ThO$_2$')    # 
        plt.plot([9e8,9e9], [9e8,9e9], 's-',  markersize=5, linewidth=1, color='#b41f24', label=r'PB-UO$_2$')     # 
        plt.plot([9e8,9e9], [9e8,9e9], 'x--', markersize=7, linewidth=1, color='#0047ba', label=r'LL-ThO$_2$')    #  red: #b41f24
        plt.plot([9e8,9e9], [9e8,9e9], '^-',  markersize=6, linewidth=1, color='#0047ba', label=r'LL-UO$_2$')     #  blue: #0047ba


        # plt.title(f'Tritium breeding ratio (All)') # Exclude title for production figs --ppark 2025-08-06
        plt.xlabel(r'Fertile isotope density in blanket [kg$/$m$^3$]')
        plt.ylabel('Initial fissile production rate [kg$/$yr]')

        plt.xlim(-5,165)
        plt.ylim(-15,2*(250+7.5))  # plt.ylim(-7.5,225+7.5) 

        # Tick grid
        ax = plt.gca()
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
 
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(MultipleLocator(50))   # ax.yaxis.set_major_locator(MultipleLocator(25))
        ax.yaxis.set_minor_locator(MultipleLocator(25)) # ax.yaxis.set_minor_locator(MultipleLocator(12.5))
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

        plt.tight_layout()

        leg = plt.legend(fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=1)
        leg.get_frame().set_linewidth(0.5) 
    
        if self.save:
            plt.savefig(f'./Figures/pdf/fig_fissile_per_yr.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/fig_fissile_per_yr.png', bbox_inches='tight', format='png')
            print("Exported fissile production per year plot for all blankets.")
        else:
            print("Did not export fissile production per year plot due to user setting.")

        if self.show: plt.show()
        plt.close('all')



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

        titles = [r"FLiBe-UF$_4$", r"FLiBe-ThF$_4$", r"PB-UO$_2$", r"PB-ThO$_2$", r"LL-UO$_2$", r"LL-ThO$_2$",]
        dfs    = [self.flibe_u_eb_df, self.flibe_th_eb_df, self.pb_u_eb_df, self.pb_th_eb_df, self.pbli_u_eb_df, self.pbli_th_eb_df,]

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
            plt.savefig(f'./Figures/pdf/fig_cum_norm_histogram.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/fig_cum_norm_histogram.png', bbox_inches='tight', format='png')
            print(f"Exported cumulative normalized histogram plots.")
        else:
            print(f"Did not export cumulative normalized histogram plot.")

        if self.show: plt.show()
        plt.close('all')




    
if __name__ == "__main__":
    plot_all()