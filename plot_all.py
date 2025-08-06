import os, sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from scipy.interpolate import Akima1DInterpolator
from scipy.interpolate import PchipInterpolator


""" import helper functions """
from Python.utilities import *
from Python.parameters import *



def main():
    """ Read CSV data into pandas DataFrames """
    # Total reaction rates
     flibe_u_rr_df = pd.read_csv('./Figures/data/FLiBe_U_FW4cm_Li07.5_900K_2025-07-22_tot_rxn_rates.csv')
    flibe_th_rr_df = pd.read_csv('./Figures/data/FLiBe_Th_FW4cm_Li07.5_900K_2025-07-22_tot_rxn_rates.csv')
        pbli_rr_df = pd.read_csv('./Figures/data/DCLL_U_FW4cm_Li90_900K_2025-07-22_tot_rxn_rates.csv')
      pebble_rr_df = pd.read_csv('./Figures/data/HCPB_U_FW4cm_Li60_900K_2025-07-22_tot_rxn_rates.csv')

    # Fertile (n,gamma) per energy bin
     flibe_u_eb_df = pd.read_csv('./Figures/data/FLiBe_U_FW4cm_Li07.5_900K_2025-07-22_U238_n-gamma_Ebins.csv')
    flibe_th_eb_df = pd.read_csv('./Figures/data/FLiBe_Th_FW4cm_Li07.5_900K_2025-07-22_Th232_n-gamma_Ebins.csv')
        pbli_eb_df = pd.read_csv('./Figures/data/DCLL_U_FW4cm_Li90_900K_2025-07-22_U238_n-gamma_Ebins.csv')
      pebble_eb_df = pd.read_csv('./Figures/data/HCPB_U_FW4cm_Li60_900K_2025-07-22_U238_n-gamma_Ebins.csv') 

    df_list = [flibe_u_rr_df, flibe_th_rr_df, pbli_u_rr_df, pebble_rr_df, flibe_u_eb_df, flibe_th_eb_df, pbli_u_eb_df, pebble_eb_df]

    for df in df_list: 
        df['fertile_kg/m3'] = (df['MTU']*1e3) / (VOL_CC/1e6)
 
    
    combined_plot = Plot(df_list, save=True, show=True)
    
    combined_plot.plot_tbr()
    # combined_plot.plot_pu_per_yr()
    combined_plot.plot_histograms()

    print("All plots completed and saved.")


class Plot:

    def __init__(self, df_list, save=False, show=True):
        self.flibe_u_rr_df  = df_list[0]
        self.flibe_th_rr_df = df_list[1]
        self.pbli_u_rr_df   = df_list[2]
        self.pebble_rr_df   = df_list[3]

        self.flibe_u_eb_df  = df_list[4]
        self.flibe_th_eb_df = df_list[5]
        self.pbli_u_eb_df   = df_list[6]
        self.pebble_eb_df   = df_list[7]

        self.save, self.show = save, show
        self.name = 'All_Blankets'

        for sd in ['pdf','png',]:
            if sd == 'data': 
                sd_path = f'./Figures/{sd}/'
            else: 
                sd_path = f'./Figures/{sd}/{self.name}/'
            print(f"Ensuring directory exists: {sd_path}")
            os.makedirs(sd_path, exist_ok=True)


    def plot_tbr(self):
        print(f"\nPlotting tritium breeding ratio vs. fertile density...")

        # Setting up x, y separately here so you can remove the impurity/wppm-magnitude cases --ppark 2025-08-06
        x1, y1 =  self.flibe_u_rr_df['fertile_kg/m3'],  self.flibe_u_rr_df['Li-6(n,t)'] +  self.flibe_u_rr_df['Li-7(n,Xt)']
        x2, y2 = self.flibe_th_rr_df['fertile_kg/m3'], self.flibe_th_rr_df['Li-6(n,t)'] + self.flibe_th_rr_df['Li-7(n,Xt)']
        x3, y3 =   self.pbli_u_rr_df['fertile_kg/m3'],   self.pbli_u_rr_df['Li-6(n,t)'] +   self.pbli_u_rr_df['Li-7(n,Xt)']
        # x4, y4 =  self.pbli_th_rr_df['fertile_kg/m3'],  self.pbli_th_rr_df['Li-6(n,t)'] +  self.pbli_th_rr_df['Li-7(n,Xt)']
        x5, y5 =   self.pebble_rr_df['fertile_kg/m3'],   self.pebble_rr_df['Li-6(n,t)'] +   self.pebble_rr_df['Li-7(n,Xt)']

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

        # Dummy plots for legend -- bit of a hack lmao
        plt.plot([998,999], [998,999], 'o-', markersize=4, linewidth=1, color='#b41f24', label=r'Pebble-BISO')   # 
        plt.plot([998,999], [998,999], '^-', markersize=5, linewidth=1, color='#0047ba', label=r'PbLi-BISO')     # 
        plt.plot([998,999], [998,999], 's-', markersize=4, linewidth=1, color='black', label=r'FLiBe-UF$_4$')    # 
        plt.plot([998,999], [998,999], 'x-', markersize=6, linewidth=1, color='#66b420', label=r'FLiBe-ThF$_4$') # 
        
        
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

        """ Create new column for fissile rates """


        """ Initialize figure """
        plt.figure()
        plt.plot(self.flibe['MTU'], self.flibePu_per_yr_list,
                'o-', markersize=2, linewidth=0.75, color='#00FFFF', label='FLiBe UF4')
        plt.plot(self.flibeTh['MTU'], self.flibeU_per_yr_list,
                'o-', markersize=2, linewidth=0.75, color='#FF8000', label='FLiBe ThF4')
        plt.plot(self.pbli['MTU'][:-3], self.pbliPu_per_yr_list[:-3],
                'o-', markersize=2, linewidth=0.75, color='#FF00FF', label='PbLi')
        plt.plot(self.pebble['MTU'], self.pebblePu_per_yr_list,
                'o-', markersize=2, linewidth=0.75, color='#FFA500', label='Pebble')
        # The three MTU points for FLiBe annotation boxes
        mtu_fpoints = [0.0096, 5, 30]
        mtu_pbpoints = [2.65, 5, 50]
        mtu_ppoints = [0.076, 5, 30]

        # Function to get Pu value closest to given MTU for each fuel
        def get_closest_value(mtu_val, mtu_list, pu_list):
            idx = (abs(mtu_list - mtu_val)).idxmin()
            return pu_list[idx]

        # Gather Pu production values for all fuels at each MTU point
        # Annotation for FLiBe points
        textstr = '--- Pu-239 production summary ---\n'

        textstr += 'FLiBe UF4:\n'
        for mtu in mtu_fpoints:
            val = get_closest_value(mtu, self.flibe['MTU'], self.flibePu_per_yr_list)
            textstr += f"  {mtu:.4f} MTU: {val:.4e} kg/yr\n"

        textstr += 'FLiBe ThF4:\n'
        for mtu in mtu_fpoints:
            val = get_closest_value(mtu, self.flibeTh['MTU'], self.flibeU_per_yr_list)
            textstr += f"  {mtu:.4f} MTU: {val:.4e} kg/yr\n"

        textstr += 'PbLi:\n'
        for mtu in mtu_pbpoints:
            val = get_closest_value(mtu, self.pbli['MTU'], self.pbliPu_per_yr_list)
            textstr += f"  {mtu:.4f} MTU: {val:.4e} kg/yr\n"

        textstr += 'Pebble:\n'
        for mtu in mtu_ppoints:
            val = get_closest_value(mtu, self.pebble['MTU'], self.pebblePu_per_yr_list)
            textstr += f"  {mtu:.4f} MTU: {val:.4e} kg/yr\n"

        # Add the single annotation box on the side of the plot

        plt.title(f'Pu-239 / U-233 production per year (All Fuels)')
        plt.xlabel('Uranium / Thorium loaded [metric tons]')
        plt.ylabel(r'Pu-239 / U-233 produced [kg$/$yr]')
        plt.tight_layout()
        plt.legend()
    
        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_Pu_per_yr_all.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_Pu_per_yr_all.png', bbox_inches='tight', format='png')
            print("Exported Pu-239 production per year plot for all fuels.")
        else:
            print("Did not export Pu-239 production per year plot due to user setting.")

        if self.show: plt.show()
        plt.close('all')



    def plot_cum_norm_u_vs_energy(self):
        """
        Plots cumulative, normalized Pu production vs. energy for contours of MTU.
        """
        print(f"\nPlotting cumulative, normalized Pu production vs. energy with MTU contours ...")

        plt.figure(figsize=(9,6))

        for i, cell_id in enumerate(self.cell_ids):
            if self.u_list[i] in [10, 20, 30, 40, 50]:
                df = self.U238_ng_Ebin_rr_df[self.U238_ng_Ebin_rr_df['cell'] == cell_id]
                x = df['energy mid [eV]']
                y = df['mean']

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



    
if __name__ == "__main__":
    main()