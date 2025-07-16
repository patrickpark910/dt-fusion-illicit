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


def main():

    # Read and plot tallies for each Li enrich case
    flibe = PlotStatepoint(enrich_li=7.5, temp_k=300, save=True, show=False, to_csv=True)
    pbli = PlotStatepointPbLi(enrich_li=90, u_list=MASS_U_LIST_PBLI, save=True, show=False, to_csv=True)
    pebble = PlotStatepointPebble(enrich_li=60, u_list=MASS_U_LIST_PBLI, save=True, show=False, to_csv=True)

    current_sp = PlotStatepointALL(flibe, pbli, pebble, save=False, show=True, to_csv=False)
    current_sp.plot_tbr()
    #current_sp.plot_pu()
    # current_sp.plot_pu_per_yr()
    #current_sp.plot_pu_vs_energy()
    #current_sp.plot_rel_pu_vs_energy()
    #current_sp.plot_flux_vs_energy()

    '''  
    current_sp.print_rxn_rates()
    current_sp.plot_tbr()
    current_sp.plot_pu()
    current_sp.plot_pu_per_mtu()
    
    current_sp.plot_rxn_rates()
    current_sp.plot_pu_vs_energy()
    current_sp.plot_rel_pu_vs_energy()
    
    '''
class PlotStatepointALL:

    def __init__(self, flibe, pbli, pebble, save=False, show=True, to_csv=False):
        self.flibe = flibe
        self.pbli = pbli
        self.pebble = pebble

        self.save, self.show, self.to_csv = save, show, to_csv
        self.name = f'All_Fuels'

    def plot_tbr(self):
        print(f"\nPlotting tritium breeding ratio vs. uranium loading...")

        plt.figure()
        plt.errorbar(self.flibe.u_list, self.flibe.tbr_list, yerr=self.flibe.tbr_err_list,
                     fmt='o-', markersize=2, capsize=0, linewidth=0.75, color='#00FFFF')
        plt.errorbar(self.pbli.u_list, self.pbli.tbr_list, yerr=self.pbli.tbr_err_list,
                     fmt='o-', markersize=2, capsize=0, linewidth=0.75, color='#FF00FF')
        plt.errorbar(self.pebble.u_list, self.pebble.tbr_list, yerr=self.pebble.tbr_err_list,
                     fmt='o-', markersize=2, capsize=0, linewidth=0.75, color='#00FF00')

        plt.title(f'Tritium breeding ratio (All)')
        plt.xlabel('Uranium loaded [metric tons]')
        plt.ylabel('Tritium breeding ratio')
        plt.tight_layout()

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_tbr.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_tbr.png', bbox_inches='tight', format='png')
            print(f"Exported tritium breeding ratio plots.")
        else:
            print(f"Did not export tritium breeding ratio plots.")

        if self.show:
            plt.show()
        plt.close('all')
if __name__ == "__main__":
    main()