import openmc
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
import imageio.v2 as iio # use v2 to avoid deprecation warnings --ppark


""" import helper functions """
sys.path.insert(0, f"{os.getcwd()}/Python")
from parameters import *
from utilities import *
from flibe_plot import PlotStatepoint
from pbli_plot import PlotStatepointPbLi
from pebble_plot import PlotStatepointPebble

def main():
    # Read CSV data into pandas DataFrames
    flibe = pd.read_csv('./Figures/data/FLiBe_Li07.5_7_22_tot_rxn_rates.csv')
    flibeTh = pd.read_csv('./Figures/data/FLiBe_Th_Li07.5_7_22_tot_rxn_rates.csv')
    pbli = pd.read_csv('./Figures/data/PbLi_Li90_7_22_tot_rxn_rates.csv')
    pebble = pd.read_csv('./Figures/data/Pebble_Li60_7_22_tot_rxn_rates.csv')
    
    combined_plot = PlotStatepointALL(flibe, flibeTh, pbli, pebble, save=True, show=False, to_csv=True)
    
    combined_plot.plot_tbr()
    combined_plot.plot_pu()
    combined_plot.plot_pu_per_yr()
    combined_plot.plot_pu_per_mtu()
    combined_plot. plot_fis()

    print("All plots completed and saved.")


class PlotStatepointALL:

    def __init__(self, flibe, flibeTh, pbli, pebble, save=False, show=True, to_csv=False):
        self.flibe = flibe
        self.flibeTh = flibeTh
        self.pbli = pbli
        self.pebble = pebble

        self.save = save
        self.show = show
        self.to_csv = to_csv
        self.name = 'All_Fuels'

        for sd in ['pdf','png','gif','data']:
            if sd == 'data': sd_path = f'./Figures/{sd}/'
            else: sd_path = f'./Figures/{sd}/{self.name}/'
            print(f"Ensuring directory exists: {sd_path}")
            os.makedirs(sd_path, exist_ok=True)
        self.flibePu_per_yr_list = []
        for FPu_per_srcn in self.flibe['U-238(n,gamma)']:
            self.flibePu_per_yr_list.append( FPu_per_srcn * NPS_FUS * SEC_PER_YR * AMU_PU239 / AVO / 1e3 )
        self.flibeU_per_yr_list = []
        for U_per_srcn in self.flibeTh['Th-232(n,gamma)']:
            self.flibeU_per_yr_list.append( U_per_srcn * NPS_FUS * SEC_PER_YR * AMU_U233 / AVO / 1e3 )
        self.pbliPu_per_yr_list = []
        for PBPu_per_srcn in self.pbli['U-238(n,gamma)']:
            self.pbliPu_per_yr_list.append( PBPu_per_srcn * NPS_FUS * SEC_PER_YR * AMU_PU239 / AVO / 1e3 )
        self.pebblePu_per_yr_list = []
        for PPu_per_srcn in self.pebble['U-238(n,gamma)']:
            self.pebblePu_per_yr_list.append( PPu_per_srcn * NPS_FUS * SEC_PER_YR * AMU_PU239 / AVO / 1e3 )



    def plot_tbr(self):
        print(f"\nPlotting tritium breeding ratio vs. uranium loading...")

        plt.figure()

        plt.plot(self.flibe['MTU'], self.flibe['Li-6(n,t)'],
                 'o-', markersize=2, linewidth=0.75, color='#00FFFF', label='FLiBe UF4')
        plt.plot(self.flibeTh['MTU'], self.flibeTh['Li-6(n,t)'],
                 'o-', markersize=2, linewidth=0.75, color='#FF8000', label='FLiBe ThF4')
        plt.plot(self.pbli['MTU'], self.pbli['Li-6(n,t)'],
                 'o-', markersize=2, linewidth=0.75, color='#FF00FF', label='PbLi')
        plt.plot(self.pebble['MTU'], self.pebble['Li-6(n,t)'],
                 'o-', markersize=2, linewidth=0.75, color='#FFA500', label='Pebble')

        plt.title(f'Tritium breeding ratio (All)')
        plt.xlabel('Uranium/Thorium loaded [metric tons]')
        plt.ylabel('Tritium breeding ratio')
        plt.tight_layout()
        plt.legend()

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_tbr.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_tbr.png', bbox_inches='tight', format='png')
            print(f"Exported tritium breeding ratio plots.")
        else:
            print(f"Did not export tritium breeding ratio plots.")

        if self.show: plt.show()
        plt.close('all')

    def plot_pu(self):
        """ Plot U-238 (n,gamma) reaction rate for all fuels on one plot """

        print(f"\nPlotting U-238/Th-232 (n,gamma) reaction rate vs. uranium loading for all fuels...")

        plt.figure()

        plt.plot(self.flibe['MTU'], self.flibe['U-238(n,gamma)'],
                 'o-', markersize=2, linewidth=0.75, color='#00FFFF', label='FLiBe UF4')
        plt.plot(self.flibeTh['MTU'], self.flibeTh['Th-232(n,gamma)'],
                 'o-', markersize=2, linewidth=0.75, color='#FF8000', label='FLiBe ThF4')
        plt.plot(self.pbli['MTU'], self.pbli['U-238(n,gamma)'],
                 'o-', markersize=2, linewidth=0.75, color='#FF00FF', label='PbLi')
        plt.plot(self.pebble['MTU'], self.pebble['U-238(n,gamma)'],
                 'o-', markersize=2, linewidth=0.75, color='#FFA500', label='Pebble')

        plt.title('U-238/Th-232 (n,gamma) reaction rate (All Fuels)')
        plt.xlabel('Uranium/Thorium loaded [metric tons]')
        plt.ylabel(r'U-238/Th-232 (n,$\gamma$) reaction rate')
        plt.legend()
        plt.tight_layout()

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_U238_ng_all.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_U238_ng_all.png', bbox_inches='tight', format='png')
            print("Exported U-238 (n,gamma) reaction rate plot for all fuels.")
        else:
            print("Did not export U-238 (n,gamma) reaction rate plot due to user setting.")

        if self.show:
            plt.show()

        plt.close('all')

    def plot_fis(self):
        """ Plot U-238 (n,fis) reaction rate for all fuels on one plot """

        print(f"\nPlotting U-238 (n,fis) reaction rate vs. uranium loading for all fuels...")

        plt.figure()

        plt.plot(self.flibe['MTU'], self.flibe['U-238(n,fis)'],
                 'o-', markersize=2, linewidth=0.75, color='#00FFFF', label='FLiBe UF4')
        plt.plot(self.flibeTh['MTU'], self.flibeTh['Th-232(n,fis)'],
                 'o-', markersize=2, linewidth=0.75, color='#FF8000', label='FLiBe ThF4')
        plt.plot(self.pbli['MTU'], self.pbli['U-238(n,fis)'],
                 'o-', markersize=2, linewidth=0.75, color='#FF00FF', label='PbLi')
        plt.plot(self.pebble['MTU'], self.pebble['U-238(n,fis)'],
                 'o-', markersize=2, linewidth=0.75, color='#FFA500', label='Pebble')

        plt.title('U-238/Th-232 (n,fis) reaction rate (All Fuels)')
        plt.xlabel('Uranium/Thorium loaded [metric tons]')
        plt.ylabel(r'U-238/Th-232 (n,fis) reaction rate')
        plt.legend()
        plt.tight_layout()

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_U238_fis_all.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_U238_fis_all.png', bbox_inches='tight', format='png')
            print("Exported U-238 (n,fis) reaction rate plot for all fuels.")
        else:
            print("Did not export U-238 (n,fis) reaction rate plot due to user setting.")

        if self.show:
            plt.show()

        plt.close('all')


    def plot_pu_per_mtu(self):
        
        print(f"\nPlotting U-238 (n,gamma) reaction rate per MTU vs. uranium loading for all fuels...")

        plt.figure()

        flibe_pu_per_mtu = self.flibe['U-238(n,gamma)'] / self.flibe['MTU']
        plt.plot(self.flibe['MTU'], flibe_pu_per_mtu,
                 'o-', markersize=2, linewidth=0.75, color='#00FFFF', label='FLiBe UF4')
        flibe_u_per_mtu = self.flibeTh['Th-232(n,gamma)'] / self.flibeTh['MTU']
        plt.plot(self.flibeTh['MTU'], flibe_u_per_mtu,
                 'o-', markersize=2, linewidth=0.75, color='#FF8000', label='FLiBe ThF4')
        pbli_pu_per_mtu = self.pbli['U-238(n,gamma)'] / self.pbli['MTU']
        plt.plot(self.pbli['MTU'], pbli_pu_per_mtu,
                 'o-', markersize=2, linewidth=0.75, color='#FF00FF', label='PbLi')
        pebble_pu_per_mtu = self.pebble['U-238(n,gamma)'] / self.pebble['MTU']
        plt.plot(self.pebble['MTU'], pebble_pu_per_mtu,
                 'o-', markersize=2, linewidth=0.75, color='#FFA500', label='Pebble')

        plt.title(f'U-238/Th-232 (n,gamma) reaction rate per MTU (All Fuels)')
        plt.xlabel('Uranium/Thorium loaded [metric tons]')
        plt.ylabel(r'U-238/Th-232 (n,$\gamma$) reaction rate per MTU')
        plt.tight_layout()
        plt.legend()

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_U238_ng_per_MTU_all.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_U238_ng_per_MTU_all.png', bbox_inches='tight', format='png')
            print("Exported U-238 (n,gamma) reaction rate per MTU plot for all fuels.")
        else:
            print("Did not export U-238 (n,gamma) reaction rate per MTU plot due to user setting.")

        if self.show: plt.show()
        plt.close('all')

    def plot_pu_per_yr(self):
      
        print(f"\nPlotting Pu-239 production per year for all fuels...")

        plt.figure()
        plt.plot(self.flibe['MTU'], self.flibePu_per_yr_list,
                'o-', markersize=2, linewidth=0.75, color='#00FFFF', label='FLiBe UF4')
        plt.plot(self.flibeTh['MTU'], self.flibeU_per_yr_list,
                'o-', markersize=2, linewidth=0.75, color='#FF8000', label='FLiBe ThF4')
        plt.plot(self.pbli['MTU'], self.pbliPu_per_yr_list,
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
        ax.text(1.05, 0.5, textstr, transform=ax.transAxes, fontsize=8,
                verticalalignment='center', bbox=box_props, fontfamily='monospace')

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


    
if __name__ == "__main__":
    main()