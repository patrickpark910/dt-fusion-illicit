import openmc
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
import imageio.v2 as iio # use v2 to avoid deprecation warnings --ppark

# Import helper functions
from Python.parameters import *
from Python.utilities import *
from Python.flibe_plots_extra import *


def main():
    current_sp = PlotStatepoint(enrich_li=7.5, temp_k=900, save=True, show=False, to_csv=True)

    current_sp.print_rxn_rates()
    current_sp.plot_tbr()
    current_sp.plot_u()
    current_sp.plot_u_per_mtu()
    
    #current_sp.plot_rxn_rates()
    current_sp.plot_u_vs_energy()
    current_sp.plot_rel_u_vs_energy()
    current_sp.plot_flux_vs_energy()
    current_sp.plot_u_per_yr()
    current_sp.plot_cum_norm_th_vs_energy()

    '''
    
    current_sp.print_rxn_rates()
    current_sp.plot_tbr()
    current_sp.plot_pu()
    current_sp.plot_pu_per_mtu()
    
    current_sp.plot_rxn_rates()
    current_sp.plot_pu_vs_energy()
    current_sp.plot_rel_pu_vs_energy()
    current_sp.plot_flux_vs_energy()
    '''


class PlotStatepoint:

    def __init__(self, enrich_li=7.5, temp_k=900, save=False, show=True, to_csv=False):
        self.e = enrich_li
        self.temp = temp_k
        self.name = f"FLiBe_Th_Li{self.e:04.1f}_2025-07-22"
        self.save, self.show, self.to_csv = save, show, to_csv
        self.u_list = MASS_U_LIST
        plt.rcParams.update({
            'axes.titlesize': 16,       # Title font size
            'axes.labelsize': 14,       # Axis label font size
            'xtick.labelsize': 12,      # X-axis tick label size
            'ytick.labelsize': 12,      # Y-axis tick label size
            'legend.fontsize': 12,      # Legend font size
            'figure.titlesize': 16,     # Figure title font size (if using suptitle)
        })


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
        Th    = sp.get_tally(name='thorium rxn rates')
        Li   = sp.get_tally(name='lithium rxn rates')
        F    = sp.get_tally(name='fluorine rxn rates') # unused
        # Be   = sp.get(name='beryllium rxn rates')

        # Convert Tally objects into pandas dataframes
        self.flux_df = flux.get_pandas_dataframe() # Convert Tally object! Not 'XXX_rr'! --ppark
        Th_df  = Th.get_pandas_dataframe()
        Li_df = Li.get_pandas_dataframe()
        # F_df = F.get_pandas_dataframe() # unused

        # Add new column for energy bin midpoint (for plotting)
        for df in [self.flux_df, Th_df, Li_df]: # , F_df]:
            df['energy mid [eV]'] = (df['energy low [eV]'] + df['energy high [eV]'])/ 2

        """ Reaction rates for every mtu loading and for every energy bin """
        self.Th232_ng_Ebin_df = Th_df[(Th_df['nuclide'] == 'Th232') & (Th_df['score'] == '(n,gamma)')][['energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        self.Th232_fis_Ebin_df = Th_df[(Th_df['nuclide'] == 'Th232') & (Th_df['score'] == 'fission')][['energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        self.Li6_nt_Ebin_df   = Li_df[(Li_df['nuclide'] == 'Li6') & (Li_df['score'] == '(n,Xt)')][['energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        self.Li7_nt_Ebin_df   = Li_df[(Li_df['nuclide'] == 'Li7') & (Li_df['score'] == '(n,Xt)')][['energy mid [eV]', 'cell', 'mean', 'std. dev.']]

        """ Reaction rate for every mtu loading summed over all energy bins"""
        Li6_rr  = Li.summation(filter_type=openmc.EnergyFilter, nuclides=['Li6'], remove_filter=True)
        Li7_rr  = Li.summation(filter_type=openmc.EnergyFilter, nuclides=['Li7'], remove_filter=True)
        F19_rr  = F.summation(filter_type=openmc.EnergyFilter, nuclides=['F19'], remove_filter=True)
        Th232_rr = Th.summation(filter_type=openmc.EnergyFilter, nuclides=['Th232'], remove_filter=True)

        # Lithium reaction rates for each mtu loading summed over all energies
        Li6_df = Li6_rr.get_pandas_dataframe()
        Li7_df = Li7_rr.get_pandas_dataframe()
        self.Li6_nt_list = Li6_df[Li6_df['score'] == '(n,Xt)'][['cell', 'mean', 'std. dev.']]['mean'].tolist()
        self.Li6_nt_err_list = Li6_df[Li6_df['score'] == '(n,Xt)'][['cell', 'mean', 'std. dev.']]['std. dev.'].tolist()
        self.Li7_nt_list = Li7_df[Li7_df['score'] == '(n,Xt)'][['cell', 'mean', 'std. dev.']]['mean'].tolist()
        self.Li7_nt_err_list = Li7_df[Li7_df['score'] == '(n,Xt)'][['cell', 'mean', 'std. dev.']]['std. dev.'].tolist()
        self.tbr_list = [x + y for x, y in zip(self.Li6_nt_list, self.Li7_nt_list)]
        self.tbr_err_list = [x + y for x, y in zip(self.Li6_nt_err_list, self.Li7_nt_err_list)]

        # Th232 reaction rates for each mtu loading summed over all energies
        Th232_df = Th232_rr.get_pandas_dataframe()
        self.Th232_fis_list     = Th232_df[Th232_df['score'] == 'fission'][['cell', 'mean', 'std. dev.']]['mean'].tolist()
        self.Th232_fis_err_list = Th232_df[Th232_df['score'] == 'fission'][['cell', 'mean', 'std. dev.']]['std. dev.'].tolist()
        self.Th232_ng_list     = Th232_df[Th232_df['score'] == '(n,gamma)'][['cell', 'mean', 'std. dev.']]['mean'].tolist()
        self.Th232_ng_err_list = Th232_df[Th232_df['score'] == '(n,gamma)'][['cell', 'mean', 'std. dev.']]['std. dev.'].tolist()

        # Uranium kg per year
        self.U_per_yr_list = []
        for U_per_srcn in self.Th232_ng_list:
            self.U_per_yr_list.append( U_per_srcn * NPS_FUS * SEC_PER_YR * AMU_U233 / AVO / 1e3 )

        """Create list of cell IDs randomly assigned by OpenMC to match mtu loading to cell ID"""
        self.cell_ids = Th_df['cell'].unique().tolist()


    def print_rxn_rates(self):
        """
        Prints and exports as CSV total reaction rates in each mtu loading summed over energies
        """
        df = pd.DataFrame({'MTU':MASS_U_LIST,
                           'Li-6(n,t)':self.Li6_nt_list,
                           'Li-7(n,Xt)':self.Li7_nt_list,
                           'Th-232(n,gamma)':self.Th232_ng_list,
                           'Th-232(n,fis)':self.Th232_fis_list,})

        print(f"\nTotal reaction rates in FLiBe with {self.e}wt% Li-6 are:")
        print(f"{df.to_string(index=False)}\n") # ensures the whole df gets printed

        if self.to_csv:
            path = f"./Figures/data/{self.name}_tot_rxn_rates.csv"
            df.to_csv(path, index=False) # keep as 'self.pu_path'!
            print(f"Exported total reaction rates CSV to: {path}")


    def plot_tbr(self):
        """ Plot tritium breeding ratio """
        print(f"\nPlotting tritium breeding ratio vs. thorium loading...")

        plt.figure()
        plt.errorbar(MASS_U_LIST, self.tbr_list, yerr=self.tbr_err_list, fmt='o-', markersize=2, capsize=0, linewidth=0.75, color='#000000',) # turn capsize > 0 to show error bars, but they're super smol

        plt.xlim(-1.5,51.5)
    

        plt.title(f'Tritium breeding ratio (Li-6 {self.e}wt%, {self.temp} K)')
        plt.xlabel('Thorium loaded [metric tons]')
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


    def plot_u(self):
        """ Plot U-238 (n,gamma) reaction rate """

        print(f"\nPlotting Th-232 (n,gamma) reaction rate vs. thorium loading...")

        plt.figure()
        plt.errorbar(MASS_U_LIST, self.Th232_ng_list, yerr=self.Th232_ng_err_list, fmt='o-', markersize=2, capsize=0, linewidth=0.75, color='#000000',) # turn capsize > 0 to show error bars, but they're super smol
        plt.xlim(-1.5,51.5)
        
        plt.title(f'U-233 production rate (Li-6 {self.e}wt%, {self.temp} K)')
        plt.xlabel('Thorium loaded [metric tons]')
        plt.ylabel(r'U-233 produced [atoms$/$source neutron]')
        plt.tight_layout()

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_U.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_U.png', bbox_inches='tight', format='png')
            print(f"Exported Th-232 (n,gamma) plots.")
        else:
            print(f"Did not export Th-232 (n,gamma) reaction rate vs. throium loading plots due to user setting.")

        if self.show: plt.show()
        plt.close('all')


    def plot_u_per_mtu(self):
        """
        Plot of Th-232 (n,gamma) reaction rate divided by MTU
        """
        print(f"\nPlotting Th-232 (n,gamma) reaction rate per MTU vs. thorium loading...")

        U_per_mtu_list = []
        for i, m in enumerate(MASS_U_LIST[1:]):
            U_per_mtu_list.append(self.Th232_ng_list[i+1]/m) # fixed to [i+1] -ezoccoli 2025-07-11

        plt.figure()
        plt.plot(MASS_U_LIST[1:], U_per_mtu_list, markersize=2, linewidth=0.75, color='#000000',) # turn capsize > 0 to show error bars, but they're super smol

        plt.xlim(-1.5,51.5)
        # plt.ylim(-0.005,0.105)
        plt.title(f'U-233 production per MTU (Li-6 {self.e}wt%, {self.temp} K)')
        plt.xlabel('Thorium loaded [metric tons]')
        plt.ylabel(r'U-233 produced per MTU [atoms$/$source neutron$/$MTU]')
        plt.tight_layout()

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_U_per_MTU.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_U_per_MTU.png', bbox_inches='tight', format='png')
            print(f"Exported Th-232 (n,gamma) per MTU plots.")
        else:
            print(f"Did not export Th-232 (n,gamma) per MTU plots due to user setting.")

        if self.show: plt.show()
        plt.close('all')


    def plot_u_per_yr(self):
        """
        Plot of Th-232 (n,gamma) reaction rate per year
        """
        print(f"\nPlotting kg of U produced per year...")

        try: i = MASS_U_LIST.index(0.0096)
        except Exception as e:
            print(f"\n{e}\n")
            print(f"Fatal: This function requires you to have '0.0096' in your `MASS_U_LIST`.")

        print(MASS_U_LIST)
        U_rate_per_mtu = self.U_per_yr_list[i] / 0.0096 # = Pu kg/yr/mtu  # 0.3 / 0.0096 #
        y = [m * U_rate_per_mtu for m in MASS_U_LIST]
        print(f'\nU-233 production per year per 12.1 kg U: {U_rate_per_mtu*0.01209475} | {NPS_FUS} nps')
        print(y)

        plt.figure()
        plt.plot(MASS_U_LIST, self.U_per_yr_list, markersize=2, linewidth=0.75, color='#000000', label=f'OpenMC') # turn capsize > 0 to show error bars, but they're super smol
        #plt.plot(MASS_U_LIST, y, 'r--', linewidth=0.75, label=f'200 wppm U')
        mtu_points = [0.0096, 5, 50]
        def get_closest_value(mtu_val, mtu_list, pu_list):
            idx = (abs(mtu_list - mtu_val)).argmin()
            return pu_list[idx]

        # Gather Pu production values for all fuels at each MTU point
        annotations = []
        for mtu in mtu_points:
            pbli_val = get_closest_value(mtu, np.array(self.u_list), self.U_per_yr_list)
            annotations.append((mtu, pbli_val))
        ax = plt.gca()
        box_props = dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="black", alpha=0.85)
        box_positions = {
            50: (30, max(self.U_per_yr_list) * 0.5),
            5: (15, max(self.U_per_yr_list) * 0.7),
            0.0096: (0.01, max(self.U_per_yr_list) * 0.4),
        }

        # Create one box for each MTU value showing all three fuels' Pu production
        for mtu, flibe_val in annotations:
            x, y = box_positions[mtu]
            text = (f"MTU: {mtu}\n"
                    f"FLiBe: {flibe_val:.3f} kg/yr")
            ax.text(x, y, text, fontsize=12, bbox=box_props)

        

        plt.xlim(-1.5,51.5)
        # plt.ylim(-0.005,0.105)
        plt.title(f'FLiBe U-233 production per year vs. Th-232 Loading')
        plt.xlabel('Thorium loaded [metric tons]')
        plt.ylabel(r'U-233 produced [kg$/$yr]')
        plt.tight_layout()

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_U_per_yr.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_U_per_yr.png', bbox_inches='tight', format='png')
            print(f"Exported kg of U produced per year plots.")
        else:
            print(f"Did not export kg of U produced per year plots due to user setting.")

        if self.show: plt.show()
        plt.close('all')


    def plot_rxn_rates(self):
        """
        Plot all reaction rates as function of energies, for every mtu loading, in both log-log and lin-log scales
        """

        print(f"\nPlotting reaction rates vs. energy...")

        for i, cell_id in enumerate(self.cell_ids):

            print(f"...for {MASS_U_LIST[i]} MTU")
            x = self.flux_df[(self.flux_df['cell'] == cell_id)]['energy mid [eV]']
            y_flux     = self.flux_df[self.flux_df['cell'] == cell_id]['mean']
            y_Th232_ng  = self.Th232_ng_Ebin_df[self.Th232_ng_Ebin_df['cell'] == cell_id]['mean']
            y_Th232_fis = self.Th232_fis_Ebin_df[self.Th232_fis_Ebin_df['cell'] == cell_id]['mean']
            y_Li6_nt   = self.Li6_nt_Ebin_df[self.Li6_nt_Ebin_df['cell'] == cell_id]['mean']
            y_Li7_nt   = self.Li7_nt_Ebin_df[self.Li7_nt_Ebin_df['cell'] == cell_id]['mean']

            plt.figure()
            flux_line, = plt.plot(x, y_flux, linewidth=0.75, color='#000000', label=r'Flux') # black / flux_line, = is to remove it in lin scale bc it doesn't fit
            plt.plot(x, y_Li6_nt,   linewidth=0.75, color='#009ade', label=r'Li-6 (n,t)') # blue
            plt.plot(x, y_Th232_ng,  linewidth=0.75, color='#ff1f5b', label=r'Th-232 (n,g)') # red
            plt.plot(x, y_Th232_fis, linewidth=0.75, color='#f28522', label=r'Th-232 (n,fis)') # yellow / #00cd6c = green
            plt.plot(x, y_Li7_nt,   linewidth=0.75, color='#00cd6c', label=r'Li-7 (n,t)') # green
            plt.xlabel('Energy [eV]')
            plt.ylabel('Counts $/$ source neutron')

            """ Export figure in log-log scale """
            plt.title(f'Reaction rates, {MASS_U_LIST[i]} MTU (Li-6 {self.e}wt%, {self.temp} K, log-log)')
            plt.xscale('log'), plt.yscale('log')
            plt.xlim(1e0,2e7), plt.ylim(1e-14,1e2)

            # Reposition legend
            leg = plt.legend(loc='lower right', ncols=2, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1,
            leg.get_frame().set_linewidth(1)

            # Export figure
            if self.save:
                plt.savefig(f'./Figures/pdf/{self.name}/fig_rxnrates_log_{MASS_U_LIST[i]}mtu.pdf', bbox_inches='tight', format='pdf')
                plt.savefig(f'./Figures/png/{self.name}/fig_rxnrates_log_{MASS_U_LIST[i]}mtu.png', bbox_inches='tight', format='png') # you want 'log' before 'mtu' so you can flip thru them in File Explorer
                print(f"   Exported LOG-LOG reaction rates plots.")
            if self.show: plt.show()


            """ Export figure in lin-log scale """

            flux_line.remove() # remove flux plot in lin scale bc it doesn't fit

            plt.title(f'Reaction rates, {MASS_U_LIST[i]} MTU (Li-6 {self.e}wt%, {self.temp} K, lin-log)')
            plt.xscale('log'), plt.yscale('linear')
            plt.xlim(1e0,2e7), plt.ylim(0,0.03)
            plt.tight_layout()


            # Force scientific notation for y-axis
            plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            plt.gca().ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            # plt.gca().yaxis.set_major_locator(MultipleLocator(0.001)) # force y ticks at every integer

            # Reposition legend
            leg = plt.legend(loc='upper left', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1,
            leg.get_frame().set_linewidth(1)

            if self.save:
                plt.savefig(f'./Figures/pdf/{self.name}/fig_rxnrates_lin_{MASS_U_LIST[i]}mtu.pdf', bbox_inches='tight', format='pdf')
                plt.savefig(f'./Figures/png/{self.name}/fig_rxnrates_lin_{MASS_U_LIST[i]}mtu.png', bbox_inches='tight', format='png')
                print(f"   Exported LIN-LOG reaction rates plots.")
            if self.show: plt.show()


        """ Make GIF of reaction rates! """

        if self.save:
            filepaths_lin, filepaths_log = [], []
            for i in MASS_U_LIST:
                filepaths_log.append(f"./Figures/png/{self.name}/fig_rxnrates_log_{i}mtu.png")
                filepaths_lin.append(f"./Figures/png/{self.name}/fig_rxnrates_lin_{i}mtu.png")

            frames_log = [iio.imread(path) for path in filepaths_log]
            frames_lin = [iio.imread(path) for path in filepaths_lin]

            iio.mimsave(f"./Figures/gif/{self.name}/fig_rxnrates_log.gif", frames_log, fps=1, loop=0) # loop=0 : infinite loop
            iio.mimsave(f"./Figures/gif/{self.name}/fig_rxnrates_lin.gif", frames_lin, fps=1, loop=0) # loop=0 : infinite loop
            print(f"Exported GIFs of reaction rates for varying MTU.")

        plt.close('all')


    # noinspection DuplicatedCode
    def plot_u_vs_energy(self):
        """
        Plots Pu production vs. energy, for contours of MTU
        """
        print(f"\nPlotting U production vs. energy with MTU contours ...")

        plt.figure(figsize=(18,4))

        for i, cell_id in enumerate(self.cell_ids):
            if MASS_U_LIST[i] in [10, 20, 30, 40, 50]:
                x = self.Th232_ng_Ebin_df[(self.Th232_ng_Ebin_df['cell'] == cell_id)]['energy mid [eV]']
                y  = self.Th232_ng_Ebin_df[self.Th232_ng_Ebin_df['cell'] == cell_id]['mean']
                plt.plot(x, y,   linewidth=0.75, label=f'{MASS_U_LIST[i]} MTU') # green color='#00cd6c',

        plt.xlabel('Energy [eV]')
        plt.ylabel('Reactions $/$ source neutron')
        plt.title(f'Th-232 (n,gamma) rxn rate (Li-6 {self.e}wt%, {self.temp} K)')
        plt.xscale('log'), plt.yscale('log')
        # plt.xlim(1e0,2e7), plt.ylim(0,0.007)
        plt.xlim(1e1,1e3), plt.ylim(1e-7,1e-2)

        # Force scientific notation for y-axis
        # plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        # plt.gca().ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        # Reposition legend
        leg = plt.legend(loc='upper left', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1, 
        leg.get_frame().set_linewidth(1)

        # Export figure
        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_Th232ng.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_Th232ng.png', bbox_inches='tight', format='png')
            
            plt.xlim(1e0,2e7), plt.ylim(1e-7,1e-2)
            plt.savefig(f'./Figures/pdf/{self.name}/fig_Th232ng_full.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_Th232ng_full.png', bbox_inches='tight', format='png')

            print(f"   Exported U-233 production vs. energy with MTU contours plot.")
            
        if self.show: plt.show()
        plt.close('all')


    def plot_rel_u_vs_energy(self):
        """ Plots factor change in Pu rel 1 MTU vs. energy, for contours of MTU """
        print(f"\nPlotting factor change in U rel 1 MTU vs. energy with MTU contours ...")

        plt.figure(figsize=(18,4))

        for i, cell_id in enumerate(self.cell_ids):
            if MASS_U_LIST[i] in [1]:
                ref_Th232_ng = self.Th232_ng_Ebin_df[self.Th232_ng_Ebin_df['cell'] == cell_id]['mean']

            elif MASS_U_LIST[i] in [10, 20, 30, 40, 50]:
                x     = self.Th232_ng_Ebin_df[(self.Th232_ng_Ebin_df['cell'] == cell_id)]['energy mid [eV]']
                y     = self.Th232_ng_Ebin_df[self.Th232_ng_Ebin_df['cell'] == cell_id]['mean']
                y_rel = [a / b if b != 0.0 else 0.0 for a, b in zip(y, ref_Th232_ng)]
                plt.plot(x, y_rel,   linewidth=0.75, label=f'{MASS_U_LIST[i]} MTU') # green color='#00cd6c',

        plt.xlabel('Energy [eV]')
        plt.ylabel('Factor change relative to rxn rate in 1 MTU')
        plt.title(f'Relative Th-232 (n,gamma) rxn rate w.r.t. 1 MTU (Li-6 {self.e}wt%, {self.temp} K)')
        plt.xscale('log'), plt.yscale('linear')
        plt.xlim(1e1,1e3), plt.ylim(0,60)

        # Reposition legend
        leg = plt.legend(loc='upper left', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1, 
        leg.get_frame().set_linewidth(1)

        # Export figure
        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_Th232ng_rel.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_Th232ng_rel.png', bbox_inches='tight', format='png')
            
            plt.xlim(1e0,2e7)

            plt.savefig(f'./Figures/pdf/{self.name}/fig_Th232ng_rel_full.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_Th232ng_rel_full.png', bbox_inches='tight', format='png')

            print(f"   Exported factor change in Pu rel 1 MTU vs. energy with MTU contours plot.")
            
        if self.show: plt.show()
        plt.close('all')


    # noinspection DuplicatedCode
    def plot_flux_vs_energy(self):
        """
        Plots flux vs. energy, for contours of MTU
        """
        print(f"\nPlotting flux vs. energy with MTU contours...")

        plt.figure(figsize=(14,3))
        for i, cell_id in enumerate(self.cell_ids):
            if MASS_U_LIST[i] in [10, 20, 30, 40, 50]:
                x = self.flux_df[(self.flux_df['cell'] == cell_id)]['energy mid [eV]']
                y = self.flux_df[self.flux_df['cell'] == cell_id]['mean']
                plt.plot(x, y, linewidth=0.75, label=f'{MASS_U_LIST[i]} MTU') # green color='#00cd6c',

        plt.xlabel('Energy [eV]')
        plt.ylabel('Neutrons $/$ source neutron')
        plt.title(f'FLiBe Neutron flux')
        plt.xscale('log'), plt.yscale('log')
        plt.xlim(1e1,1e3), plt.ylim(1e-3,1e0)

        # Reposition legend
        leg = plt.legend(loc='upper left', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1, 
        leg.get_frame().set_linewidth(1)

        # Export figure
        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_flux.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_flux.png', bbox_inches='tight', format='png') # you want 'log' before 'mtu' so you can flip thru them in File Explorer
            
            plt.xlim(1e0,2e7), plt.ylim(1e-3,1e2)

            plt.savefig(f'./Figures/pdf/{self.name}/fig_flux_full.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_flux_full.png', bbox_inches='tight', format='png')

            print(f"   Exported flux vs. energy with MTU contours plot.")
            
        if self.show: plt.show()
        plt.close('all')

    def plot_cum_norm_th_vs_energy(self):
        """
        Plots cumulative, normalized Pu production vs. energy for contours of MTU.
        """
        print(f"\nPlotting cumulative, normalized Pu production vs. energy with MTU contours ...")

        plt.figure(figsize=(9,6))

        for i, cell_id in enumerate(self.cell_ids):
            if self.u_list[i] in [10, 20, 30, 40, 50]:
                df = self.Th232_ng_Ebin_df[self.Th232_ng_Ebin_df['cell'] == cell_id]
                x = df['energy mid [eV]']
                y = df['mean']

                # Compute cumulative sum
                cum_y = np.cumsum(y)

                # Normalize cumulative sum to max value
                cum_y_norm = cum_y / cum_y.iloc[-1] if cum_y.iloc[-1] != 0 else cum_y

                plt.plot(x, cum_y_norm, linewidth=0.75, label=f'{self.u_list[i]} MTU')

        plt.xlabel('Energy [eV]')
        plt.ylabel('Cumulative normalized reactions')
        plt.title(f'FLiBe Cumulative normalized Th-232 (n,gamma) rxn rate')
        plt.xscale('log'), plt.yscale('linear')
        plt.xlim(1e1,1e3), plt.ylim(0,1.05)

        # Reposition legend
        leg = plt.legend(loc='lower right', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,)
        leg.get_frame().set_linewidth(1)

        # Export figure
        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_Th232ng_cum_norm.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_Th232ng_cum_norm.png', bbox_inches='tight', format='png')
            print(f"   Exported cumulative normalized U-233 production vs. energy with MTU contours plot.")

            plt.xlim(1e0,2e7), plt.ylim(0,1.05)


            plt.savefig(f'./Figures/pdf/{self.name}/fig_Th232ng_cum_norm_full.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_Th232ng_cum_norm_full.png', bbox_inches='tight', format='png')

        if self.show:plt.show()
        plt.close('all')

        plt.figure(figsize=(9,6))

        for i, cell_id in enumerate(self.cell_ids):
            if self.u_list[i] in [10, 20, 30, 40, 50]:
                df = self.Th232_ng_Ebin_df[self.Th232_ng_Ebin_df['cell'] == cell_id]
                x = df['energy mid [eV]']
                y = df['mean']

                # Compute cumulative sum
                cum_y = np.cumsum(y)


                plt.plot(x, cum_y, linewidth=0.75, label=f'{self.u_list[i]} MTU')

        plt.xlabel('Energy [eV]')
        plt.ylabel('Cumulative normalized reactions')
        plt.title(f'FliBe CumulativeTh-232 (n,gamma) rxn rate')
        plt.xscale('log'), plt.yscale('linear')
        plt.xlim(1e1,1e3)

        # Reposition legend
        leg = plt.legend(loc='lower right', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,)
        leg.get_frame().set_linewidth(1)

        # Export figure
        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_Th232ng_cum.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_Th232ng_cum.png', bbox_inches='tight', format='png')
            print(f"   Exported cumulative normalized U-233 production vs. energy with MTU contours plot.")

            plt.xlim(1e0,2e7)


            plt.savefig(f'./Figures/pdf/{self.name}/fig_Th232ng_cum_full.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_Th232ng_cum_full.png', bbox_inches='tight', format='png')

        if self.show:plt.show()
        plt.close('all')



if __name__ == "__main__":
    main()
