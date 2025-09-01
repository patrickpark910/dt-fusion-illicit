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



def main():
    """ Read CSV data into pandas DataFrames """

    combined_plot = Plot(save=True, show=True)
    
    """ Pick which one to plot by uncommenting """
    # combined_plot.plot_tbr()
    combined_plot.plot_pu_per_yr()
    # combined_plot.plot_cum_norm_histogram()

    print("All plots completed and saved.")


class Plot:

    def __init__(self, save=False, show=True):

        self.flibe_u_rr_df  = pd.read_csv('./Figures/data/FLiBe_U_FW4cm_Li07.5_900K_2025-07-22_tot_rxn_rates.csv')
        self.flibe_th_rr_df = pd.read_csv('./Figures/data/FLiBe_Th_FW4cm_Li07.5_900K_2025-07-22_tot_rxn_rates.csv')
        self.pbli_u_rr_df   = pd.read_csv('./Figures/data/DCLL_FW4cm_Li90_900K_2025-07-22_tot_rxn_rates.csv')
        self.pebble_rr_df   = pd.read_csv('./Figures/data/HCPB_FW4cm_Li60_900K_2025-07-22_tot_rxn_rates.csv')

        self.flibe_u_eb_df  = pd.read_csv('./Figures/data/FLiBe_U_FW4cm_Li07.5_900K_2025-07-22_U238_n-gamma_Ebins.csv')
        self.flibe_th_eb_df = pd.read_csv('./Figures/data/FLiBe_Th_FW4cm_Li07.5_900K_2025-07-22_Th232_n-gamma_Ebins.csv')
        self.pbli_u_eb_df   = pd.read_csv('./Figures/data/DCLL_FW4cm_Li90_900K_2025-07-22_U238_n-gamma_Ebins.csv')
        self.pebble_eb_df   = pd.read_csv('./Figures/data/HCPB_FW4cm_Li60_900K_2025-07-22_U238_n-gamma_Ebins.csv') 

        self.save, self.show = save, show
        self.name = 'All_Blankets'

        for sd in ['pdf','png','data']:
            if sd == 'data': 
                sd_path = f'./Figures/{sd}/'
            else: 
                sd_path = f'./Figures/{sd}/{self.name}/'
            print(f"Ensuring directory exists: {sd_path}")
            os.makedirs(sd_path, exist_ok=True)


    def plot_tbr(self):
        print(f"\nPlotting tritium breeding ratio vs. fertile density...")

        # Setting up x, y separately here so you can remove the impurity/wppm-magnitude cases --ppark 2025-08-06
        x1, y1 =  self.flibe_u_rr_df['kg/m3_fertile'],  self.flibe_u_rr_df['Li-6(n,t)'] +  self.flibe_u_rr_df['Li-7(n,Xt)']
        x2, y2 = self.flibe_th_rr_df['kg/m3_fertile'], self.flibe_th_rr_df['Li-6(n,t)'] + self.flibe_th_rr_df['Li-7(n,Xt)']
        x3, y3 =   self.pbli_u_rr_df['kg/m3_fertile'],   self.pbli_u_rr_df['Li-6(n,t)'] +   self.pbli_u_rr_df['Li-7(n,Xt)']
        # x4, y4 =  self.pbli_th_rr_df['kg/m3_fertile'],  self.pbli_th_rr_df['Li-6(n,t)'] +  self.pbli_th_rr_df['Li-7(n,Xt)']
        x5, y5 =   self.pebble_rr_df['kg/m3_fertile'],   self.pebble_rr_df['Li-6(n,t)'] +   self.pebble_rr_df['Li-7(n,Xt)']

        # AkimaInterpolation ripped from my ELWR paper --ppark 2025-08-06
        x_fine = np.linspace(x1.min(), x1.max(), 300) # Evaluate on a fine grid
        y_akima1 = Akima1DInterpolator(x1, y1)(x_fine)
        y_akima2 = Akima1DInterpolator(x2, y2)(x_fine)
        y_akima3 = Akima1DInterpolator(x3, y3)(x_fine)
        # y_akima4 = akima4(x_fine)
        y_akima5 = Akima1DInterpolator(x5, y5)(x_fine)

        for y in [y_akima1,y_akima2,y_akima3,y_akima5]: # y_akima4,
            for i in range(1, len(y)):
                if y[i] > y[i-1]:
                    y[i] = y[i-1] # adjust the current point to ensure it's not greater than the previous point


        plt.figure(figsize=(7.5,5))

        plt.scatter(x5, y5, marker='o', s=40, color='#b41f24') # ZR clean
        plt.scatter(x3, y3, marker='^', s=50, color='#0047ba') # ZR clean
        plt.scatter(x1, y1, marker='s', s=40, color='black') # ZR clean
        plt.scatter(x2, y2, marker='x', s=60, color='#66b420') # ZR clean

        plt.plot(x_fine, y_akima5, linewidth=1, color='#b41f24',)   # 'o-', markersize=4,  label=r'Pebble-BISO'
        plt.plot(x_fine, y_akima3, linewidth=1, color='#0047ba', )     # '^-', markersize=5, label=r'PbLi-BISO'
        plt.plot(x_fine, y_akima1, linewidth=1, color='black', )    # 's-', markersize=4, label=r'FLiBe-UF$_4$'
        plt.plot(x_fine, y_akima2, linewidth=1, color='#66b420', ) # 'x-', markersize=6, label=r'FLiBe-ThF$_4$'

        # Dummy plots for legend -- bit of a hack lmao -- ppark
        plt.plot([9e98,9e99], [9e98,9e99], 'o-', markersize=4, linewidth=1, color='#b41f24', label=r'Pebble-BISO')   # 
        plt.plot([9e98,9e99], [9e98,9e99], '^-', markersize=5, linewidth=1, color='#0047ba', label=r'PbLi-BISO')     # 
        plt.plot([9e98,9e99], [9e98,9e99], 's-', markersize=4, linewidth=1, color='black',   label=r'FLiBe-UF$_4$')    # 
        plt.plot([9e98,9e99], [9e98,9e99], 'x-', markersize=6, linewidth=1, color='#66b420', label=r'FLiBe-ThF$_4$') # 
        
        
        # plt.title(f'Tritium breeding ratio (All)') # Exclude title for production figs --ppark 2025-08-06
        plt.xlabel(r'Fertile isotope density in blanket [kg$/$m$^3$]')
        plt.ylabel('Tritium breeding ratio')

        plt.xlim(-5,165)
        plt.ylim(1.03,1.42) # plt.ylim(0.98,1.42)
        
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

        leg = plt.legend(fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=1)
        leg.get_frame().set_linewidth(0.5) 

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_tbr.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_tbr.png', bbox_inches='tight', format='png')
            print(f"Exported tritium breeding ratio plots.")
        else:
            print(f"Did not export tritium breeding ratio plots.")

        if self.show: plt.show()
        plt.close('all')


    def plot_pu_per_yr(self):
      
        print(f"\nPlotting Pu-239 production per year for all fuels...")

        # <utilities.py> Default NPS_FUS = 500 MJ/s * 3.546e17 n/MJ = 1.773e+20 n/s
        # seconds in a year = 3.156e+7 s/yr

        tot_n_per_yr = NPS_FUS * 3.156e+7

        # Gotta then convert fissile atoms/yr to kg/yr
        pu239_at_to_kg_per_yr = tot_n_per_yr / AVO * AMU_PU239 / 1e3
        u233_at_to_kg_per_yr  = tot_n_per_yr / AVO * AMU_U233  / 1e3
        
        # Setting up x, y separately here so you can remove the impurity/wppm-magnitude cases --ppark 2025-08-06
        x1, y1 =  self.flibe_u_rr_df['kg/m3_fertile'], pu239_at_to_kg_per_yr *  self.flibe_u_rr_df['U-238(n,gamma)'] 
        x2, y2 = self.flibe_th_rr_df['kg/m3_fertile'], u233_at_to_kg_per_yr * self.flibe_th_rr_df['Th-232(n,gamma)']
        x3, y3 =   self.pbli_u_rr_df['kg/m3_fertile'], pu239_at_to_kg_per_yr *   self.pbli_u_rr_df['U-238(n,gamma)'] 
        # x4, y4 =  self.pbli_th_rr_df['kg/m3_fertile'], pu239_at_to_kg_per_yr * self.pbli_th_rr_df['U-238(n,gamma)']
        x5, y5 =   self.pebble_rr_df['kg/m3_fertile'], pu239_at_to_kg_per_yr *   self.pebble_rr_df['U-238(n,gamma)']

        # AkimaInterpolation ripped from my ELWR paper --ppark 2025-08-06
        x_fine = np.linspace(x1.min(), x1.max(), 300) # Evaluate on a fine grid
        y_akima1 = Akima1DInterpolator(x1, y1)(x_fine)
        y_akima2 = Akima1DInterpolator(x2, y2)(x_fine)
        y_akima3 = Akima1DInterpolator(x3, y3)(x_fine)
        # y_akima4 = akima4(x_fine)
        y_akima5 = Akima1DInterpolator(x5, y5)(x_fine)



        plt.figure(figsize=(7.5,5))

        plt.scatter(x5, y5, marker='o', s=40, color='#b41f24') # ZR clean
        plt.scatter(x3, y3, marker='^', s=50, color='#0047ba') # ZR clean
        plt.scatter(x1, y1, marker='s', s=40, color='black') # ZR clean
        plt.scatter(x2, y2, marker='x', s=60, color='#66b420') # ZR clean

        plt.plot(x_fine, y_akima5, linewidth=1, color='#b41f24',)   # 'o-', markersize=4,  label=r'Pebble-BISO'
        plt.plot(x_fine, y_akima3, linewidth=1, color='#0047ba', )     # '^-', markersize=5, label=r'PbLi-BISO'
        plt.plot(x_fine, y_akima1, linewidth=1, color='black', )    # 's-', markersize=4, label=r'FLiBe-UF$_4$'
        plt.plot(x_fine, y_akima2, linewidth=1, color='#66b420', ) # 'x-', markersize=6, label=r'FLiBe-ThF$_4$'

        # Dummy plots for legend -- bit of a hack lmao -- ppark
        plt.plot([9e98,9e99], [9e98,9e99], 'o-', markersize=4, linewidth=1, color='#b41f24', label=r'Pebble-BISO')   # 
        plt.plot([9e98,9e99], [9e98,9e99], '^-', markersize=5, linewidth=1, color='#0047ba', label=r'PbLi-BISO')     # 
        plt.plot([9e98,9e99], [9e98,9e99], 's-', markersize=4, linewidth=1, color='black',   label=r'FLiBe-UF$_4$')  # 
        plt.plot([9e98,9e99], [9e98,9e99], 'x-', markersize=6, linewidth=1, color='#66b420', label=r'FLiBe-ThF$_4$') # 
        
        
        # plt.title(f'Tritium breeding ratio (All)') # Exclude title for production figs --ppark 2025-08-06
        plt.xlabel(r'Fertile isotope density in blanket [kg$/$m$^3$]')
        plt.ylabel('Fissile production rate [kg$/$yr]')

        plt.xlim(-5,165)
        plt.ylim(-10,260) 

        # Tick grid
        ax = plt.gca()
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)

        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_minor_locator(MultipleLocator(25))
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

        plt.tight_layout()

        leg = plt.legend(fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=1)
        leg.get_frame().set_linewidth(0.5) 
    
        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_fissile_per_yr.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_fissile_per_yr.png', bbox_inches='tight', format='png')
            print("Exported Pu-239 production per year plot for all fuels.")
        else:
            print("Did not export Pu-239 production per year plot due to user setting.")

        if self.show: plt.show()
        plt.close('all')



    def plot_cum_norm_histogram(self):
        """
        Plots cumulative, normalized Pu production vs. energy for contours of MTU.
        """
        print(f"\nPlotting cumulative, normalized fissile production vs. energy...")

        fig, axes = plt.subplots(2, 2, figsize=(15,10)) # sharex='col', sharey='row', 

        titles = [r"FLiBe-UF$_4$", r"FLiBe-ThF$_4$", r"PbLi-BISO", r"Pebble-BISO"]
        dfs    = [self.flibe_u_eb_df, self.flibe_th_eb_df, self.pbli_u_eb_df, self.pebble_eb_df]

        for ax, df, title in zip(axes.flatten(), dfs, titles):

            # Filter out MT_fertile loadings we want to plot 
            df = df[df["MT_fertile"].isin([10, 20, 30, 40, 50])]

            # Compute sum of 'mean' for each MT_fertile. 
            df_mean = df.groupby("MT_fertile")["mean"].sum().to_frame() # pd.DataFrame(, columns=["MT_fertile","sum"]) # 2-col df, colA = MT_fertile values, colB = sum of 'mean'
            df_mean.columns = ["sum"]
            df = df.merge(df_mean, on="MT_fertile", how="left") # adds sum as column to df

            df['norm_mean'] = df['mean'] / df['sum']

            # df.to_csv('test.csv',index=False)
            

            # Create bin edges from low and high energy boundaries
            edges = np.unique(np.concatenate([df['energy low [eV]'].values, df['energy high [eV]'].values]))
            bins  = np.sort(df['energy mid [eV]'].unique())

            sub   = df[df['MT_fertile'] == 10]
            label = int(round(sub['kg/m3_fertile'].iloc[0], -1))
            ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#ff1f5b', label=fr'{label} kg$/$m$^3$')

            sub = df[df['MT_fertile'] == 20]
            label = int(round(sub['kg/m3_fertile'].iloc[0], -1))
            ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#f48628', label=fr'{label} kg$/$m$^3$')

            sub = df[df['MT_fertile'] == 30]
            label = int(round(sub['kg/m3_fertile'].iloc[0], -1))
            ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#04cc6c', label=fr'{label} kg$/$m$^3$')

            sub = df[df['MT_fertile'] == 40]
            label = int(round(sub['kg/m3_fertile'].iloc[0], -1))
            ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#0c9edd', label=fr'{label} kg$/$m$^3$')

            sub = df[df['MT_fertile'] == 50]
            label = int(round(sub['kg/m3_fertile'].iloc[0], -1))
            ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#b367bd', label=fr'{label} kg$/$m$^3$')
            

            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.yaxis.set_minor_locator(MultipleLocator(0.05))
            ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
            ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

            ax.set_xscale('log')
            ax.set_xlim(0.67*1e0, 1.5*1e7)
            ax.set_ylim(-0.03, 1.03)
            fig.tight_layout()

            leg = ax.legend(title=title, fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=1, loc="lower right")
            leg.get_frame().set_linewidth(0.5) 

        for i in range(2):
            for j in range(2):
                ax = axes[i, j]
                if i == 1:
                    ax.set_xlabel("Incident neutron energy [eV]")
                if j == 0:
                    ax.set_ylabel(r"Cumulative fraction of fertile $($n,$\gamma)$ rxns")

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_cum_norm_histogram.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_cum_norm_histogram.png', bbox_inches='tight', format='png')
            print(f"Exported cumulative normalized histogram plots.")
        else:
            print(f"Did not export cumulative normalized histogram plot.")

        if self.show: plt.show()
        plt.close('all')

            # print((grouped["normalized_mean"]))
        

            # ax.step(energy_cum.index, energy_cum.values, where="post")


                

        """
        # Compute cumulative sum
        cum_y = np.cumsum(y)

        # Normalize cumulative sum to max value
        cum_y_norm = cum_y / cum_y.iloc[-1] if cum_y.iloc[-1] != 0 else cum_y

        plt.plot(x, cum_y_norm, linewidth=0.75, label=f'{self.u_list[i]} MTU')

        plt.xlabel('Energy [eV]')
        plt.ylabel('Cumulative normalized reactions')
        plt.title(f'FLiBe Cumulative normalized U-238 (n,gamma) rxn rate')
        plt.xscale('log'), plt.yscale('linear')
        plt.xlim(1e1,1e3), plt.ylim(0,1.05)

        # Reposition legend
        leg = plt.legend(loc='lower right', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,)
        leg.get_frame().set_linewidth(1)

        # Export figure
        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_U238ng_cum_norm.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_U238ng_cum_norm.png', bbox_inches='tight', format='png')
            print(f"   Exported cumulative normalized Pu production vs. energy with MTU contours plot.")

            plt.xlim(1e0,2e7), plt.ylim(0,1.05)


            plt.savefig(f'./Figures/pdf/{self.name}/fig_U238ng_cum_norm_full.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_U238ng_cum_norm_full.png', bbox_inches='tight', format='png')

        if self.show:plt.show()
        plt.close('all')

        plt.figure(figsize=(9,6))

        for i, cell_id in enumerate(self.cell_ids):
            if self.u_list[i] in [10, 20, 30, 40, 50]:
                df = self.U238_ng_Ebin_rr_df[self.U238_ng_Ebin_rr_df['cell'] == cell_id]
                x = df['energy mid [eV]']
                y = df['mean']

                # Compute cumulative sum
                cum_y = np.cumsum(y)

                plt.plot(x, cum_y, linewidth=0.75, label=f'{self.u_list[i]} MTU')

        plt.xlabel('Energy [eV]')
        plt.ylabel('Cumulative normalized reactions')
        plt.title(f'FLiBe Cumulative U-238 (n,gamma) rxn rate')
        plt.xscale('log'), plt.yscale('linear')
        plt.xlim(1e1,1e3)

        # Reposition legend
        leg = plt.legend(loc='lower right', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,)
        leg.get_frame().set_linewidth(1)

        # Export figure
        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_cum_hist.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_cum_hist.png', bbox_inches='tight', format='png')
            print(f"   Exported cumulative normalized Pu production vs. energy with MTU contours plot.")

            plt.xlim(1e0,2e7)


            plt.savefig(f'./Figures/pdf/{self.name}/fig_U238ng_cum_full.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_U238ng_cum_full.png', bbox_inches='tight', format='png')

        if self.show:plt.show()
        plt.close('all')
        """



    
if __name__ == "__main__":
    main()