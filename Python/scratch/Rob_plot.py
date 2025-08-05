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
    current_sp = PlotStatepointPbLi(enrich_li=90, u_list=N_TRISO_PBLI, save=True, show=False, to_csv=True)

    
    current_sp.plot_tbr()
    


    '''  
    current_sp.print_rxn_rates()
    current_sp.plot_tbr()
    current_sp.plot_pu()
    current_sp.plot_pu_per_mtu()
    
    current_sp.plot_rxn_rates()
    current_sp.plot_pu_vs_energy()
    current_sp.plot_rel_pu_vs_energy()
    
    '''


class PlotStatepointPbLi:

    def __init__(self, enrich_li=90, u_list=N_TRISO_PBLI, save=False, show=True, to_csv=False):

        self.e = enrich_li
        self.u_list = u_list
        self.temp = 900
        self.save, self.show, self.to_csv = save, show, to_csv
        self.name = f'PbLi_Li{self.e}_Rob'

        """ Load tallies """
        sp_path = f'./OpenMC/{self.name}/statepoint.100.h5'
        print("\n\n")
        print("="*42)
        print(f"Loading statepoint: {sp_path}")
        
        try:
            sp = openmc.StatePoint(sp_path) # _Li6-7.5wt  _Li6-20wt
        except Exception as e:
            print(f"\n{e}\n")
            sys.exit(f"oopsie woopsie fucky wucky can't read the sp")

        # Make directories for Figures 
        for sd in ['pdf','png','gif','data']:
            if sd == 'data': sd_path = f'./Figures/{sd}/'
            else: sd_path = f'./Figures/{sd}/{self.name}/'
            print(f"Ensuring directory exists: {sd_path}")
            os.makedirs(sd_path, exist_ok=True)

        """ Convert tallies into usable forms """
        # Read tallies
        print(f"\nReading tallies...")
        flux = sp.get_tally(name='flux')
        U    = sp.get_tally(name='uranium rxn rates')
        Li   = sp.get_tally(name='lithium rxn rates')

        # Convert Tally objects into pandas dataframes
        self.flux_df = flux.get_pandas_dataframe() # Convert Tally object! Not 'XXX_rr'! --ppark
        U_df  = U.get_pandas_dataframe()
        Li_df = Li.get_pandas_dataframe()
        # F_df = F.get_pandas_dataframe() # unused

        # Add new column for energy bin midpoint (for plotting)
        for df in [self.flux_df, U_df, Li_df]: # , F_df]:
            df['energy mid [eV]'] = (df['energy low [eV]'] + df['energy high [eV]'])/ 2

        """ Reaction rates for every mtu loading and for every energy bin """
        self.U238_ng_Ebin_df = U_df[(U_df['nuclide'] == 'U238') & (U_df['score'] == '(n,gamma)')][['energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        self.U238_fis_Ebin_df = U_df[(U_df['nuclide'] == 'U238') & (U_df['score'] == 'fission')][['energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        self.Li6_nt_Ebin_df   = Li_df[(Li_df['nuclide'] == 'Li6') & (Li_df['score'] == '(n,Xt)')][['energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        self.Li7_nt_Ebin_df   = Li_df[(Li_df['nuclide'] == 'Li7') & (Li_df['score'] == '(n,Xt)')][['energy mid [eV]', 'cell', 'mean', 'std. dev.']]

        """ Reaction rate for every mtu loading summed over all energy bins"""
        Li6_rr  = Li.summation(filter_type=openmc.EnergyFilter, nuclides=['Li6'], remove_filter=True)
        Li7_rr  = Li.summation(filter_type=openmc.EnergyFilter, nuclides=['Li7'], remove_filter=True)
        U238_rr = U.summation(filter_type=openmc.EnergyFilter, nuclides=['U238'], remove_filter=True)

        # Lithium reaction rates for each mtu loading summed over all energies
        Li6_df = Li6_rr.get_pandas_dataframe()
        Li7_df = Li7_rr.get_pandas_dataframe()
        self.Li6_nt_list = Li6_df[Li6_df['score'] == '(n,Xt)'][['cell', 'mean', 'std. dev.']]['mean'].tolist()
        self.Li6_nt_err_list = Li6_df[Li6_df['score'] == '(n,Xt)'][['cell', 'mean', 'std. dev.']]['std. dev.'].tolist()
        self.Li7_nt_list = Li7_df[Li7_df['score'] == '(n,Xt)'][['cell', 'mean', 'std. dev.']]['mean'].tolist()
        self.Li7_nt_err_list = Li7_df[Li7_df['score'] == '(n,Xt)'][['cell', 'mean', 'std. dev.']]['std. dev.'].tolist()
        self.tbr_list = [x + y for x, y in zip(self.Li6_nt_list, self.Li7_nt_list)]
        self.tbr_err_list = [x + y for x, y in zip(self.Li6_nt_err_list, self.Li7_nt_err_list)]

        # U-238 reaction rates for each mtu loading summed over all energies
        U238_df = U238_rr.get_pandas_dataframe()
        self.U238_fis_list     = U238_df[U238_df['score'] == 'fission'][['cell', 'mean', 'std. dev.']]['mean'].tolist()
        self.U238_fis_err_list = U238_df[U238_df['score'] == 'fission'][['cell', 'mean', 'std. dev.']]['std. dev.'].tolist()
        self.U238_ng_list     = U238_df[U238_df['score'] == '(n,gamma)'][['cell', 'mean', 'std. dev.']]['mean'].tolist()
        self.U238_ng_err_list = U238_df[U238_df['score'] == '(n,gamma)'][['cell', 'mean', 'std. dev.']]['std. dev.'].tolist()

        # Plutonium kg per year
        self.Pu_per_yr_list = []
        for Pu_per_srcn in self.U238_ng_list:
            self.Pu_per_yr_list.append( Pu_per_srcn * NPS_FUS * SEC_PER_YR * AMU_PU239 / AVO / 1e3 )

        """Create list of cell IDs randomly assigned by OpenMC to match mtu loading to cell ID"""
        self.cell_ids = U_df['cell'].unique().tolist()


    def print_rxn_rates(self):
        """
        Prints and exports as CSV total reaction rates in each mtu loading summed over energies
        """
        df = pd.DataFrame({'MTU':self.u_list,
                           'Li-6(n,t)':self.Li6_nt_list,
                           'Li-7(n,Xt)':self.Li7_nt_list,
                           'U-238(n,gamma)':self.U238_ng_list,
                           'U-238(n,fis)':self.U238_fis_list,})

        print(f"\nTotal reaction rates in PbLi with {self.e}wt% Li-6 are:")
        print(f"{df.to_string(index=False)}\n") # ensures the whole df gets printed

        if self.to_csv:
            path = f"./Figures/data/{self.name}_tot_rxn_rates.csv"
            df.to_csv(path, index=False) # keep as 'self.pu_path'!
            print(f"Exported total reaction rates CSV to: {path}")


    def plot_tbr(self):
        """ Plot tritium breeding ratio """
        print(f"\nPlotting tritium breeding ratio vs. BISO loading...")

        plt.figure()
        plt.errorbar(self.u_list, self.tbr_list, yerr=self.tbr_err_list, fmt='o-', markersize=2, capsize=0, linewidth=0.75, color='#000000',) # turn capsize > 0 to show error bars, but they're super smol

        plt.xlim(-1.5,450)
        
        plt.title(f'Tritium breeding ratio (PbLi, {self.e}wt% Li-6)')
        plt.xlabel('BISO/cc')
        plt.ylabel('Tritium breeding ratio')
        plt.tight_layout()

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_tbr.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_tbr.png', bbox_inches='tight', format='png')
            print(f"Exported tritium breeding ratio plots.")
        else:
            print(f"Did not export tritium breeding ratio plots.")

        if self.show: plt.show()
        plt.close('all')

if __name__ == "__main__":
    main()
