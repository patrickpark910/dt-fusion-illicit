import openmc
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
import imageio.v2 as iio # use v2 to avoid deprecation warnings --ppark

""" import helper functions """
from Python.parameters import *
from Python.utilities import *


def main():

    # Read and plot tallies for each Li enrich case
    current_sp = PlotStatepointPebble(enrich_li=60, u_list=MASS_U_LIST_HCPB, save=True, show=False, to_csv=True)

    current_sp.print_tallies()
    # current_sp.plot_tbr()
    # current_sp.plot_pu()
    # current_sp.plot_pu_per_yr()
    # current_sp.plot_pu_vs_energy()
    # current_sp.plot_rel_pu_vs_energy()
    # current_sp.plot_flux_vs_energy()
    # current_sp.plot_pu_per_mtu()
    # current_sp.plot_cum_norm_u_vs_energy()



class PlotStatepointPebble:

    def __init__(self, enrich_li=60, u_list=MASS_U_LIST_HCPB, save=False, show=True, to_csv=False):

        self.e = enrich_li
        self.fertile_mt_list = u_list
        self.temp = 900
        self.save, self.show, self.to_csv = save, show, to_csv
        self.name = f'HCPB_FW4cm_Li{self.e}_900K_2025-07-22'


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

        # In each of those Ebin_df, add a new column translating the cell # into the metric tons of fertile mass loaded """
        for df in [self.U238_ng_Ebin_df, self.U238_fis_Ebin_df, self.Li6_nt_Ebin_df, self.Li7_nt_Ebin_df]:
            df['fertile_MT'] = df['cell'].apply(lambda c: self.fertile_mt_list[c - 1])

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


    def print_tallies(self):
        """
        Prints and exports as CSV total reaction rates in each mtu loading summed over energies
        """
        df = pd.DataFrame({'MTU':MASS_U_LIST_HCPB,
                           'Li-6(n,t)':self.Li6_nt_list,
                           'Li-7(n,Xt)':self.Li7_nt_list,
                           'U-238(n,gamma)':self.U238_ng_list,
                           'U-238(n,fis)':self.U238_fis_list,})

        print(f"\nTotal reaction rates in PbLi with {self.e}wt% Li-6 are:")
        print(f"{df.to_string(index=False)}\n") # ensures the whole df gets printed

        if self.to_csv:

            csv1 = f"./Figures/data/{self.name}_tot_rxn_rates.csv"
            csv2 = f"./Figures/data/{self.name}_U238_n-gamma_Ebins.csv"
            df.to_csv(csv1, index=False) 
            self.U238_ng_Ebin_df.to_csv(csv2, index=False) 

            print(f"Exported total reaction rates CSV to:\n  {csv1}")
            print(f"Exported U-238 (n,gamma) tallies per energy bin CSV to:\n  {csv2}")


    def plot_tbr(self):
        """ Plot tritium breeding ratio """
        print(f"\nPlotting tritium breeding ratio vs. uranium loading...")

        plt.figure()
        plt.errorbar(self.fertile_mt_list, self.tbr_list, yerr=self.tbr_err_list, fmt='o-', markersize=2, capsize=0, linewidth=0.75, color='#000000',) # turn capsize > 0 to show error bars, but they're super smol

        plt.xlim(-1.5,51.5)
        
        plt.title(f'Tritium breeding ratio (Pebble, {self.e}wt% Li-6)')
        plt.xlabel('Uranium loaded [metric tons]')
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


    def plot_pu(self):
        """ Plot U-238 (n,gamma) reaction rate """

        print(f"\nPlotting U-238 (n,gamma) reaction rate vs. uranium loading...")

        plt.figure()
        plt.errorbar(self.fertile_mt_list, self.U238_ng_list, yerr=self.U238_ng_err_list, fmt='o-', markersize=2, capsize=0, linewidth=0.75, color='#000000',) # turn capsize > 0 to show error bars, but they're super smol
        plt.xlim(-1.5,51.5)
        plt.ylim(-0.005,0.105)

        plt.title(f'Pu-239 production rate (Pebble, {self.e}wt% Li-6)')
        plt.xlabel('Uranium loaded [metric tons]')
        plt.ylabel(r'Pu-239 produced [atoms$/$source neutron]')
        plt.tight_layout()

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_Pu.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_Pu.png', bbox_inches='tight', format='png')
            print(f"Exported U-238 (n,gamma) plots.")
        else:
            print(f"Did not export U-238 (n,gamma) reaction rate vs. uranium loading plots due to user setting.")

        if self.show: plt.show()
        plt.close('all')


    def plot_pu_per_mtu(self):
        """
        Plot of U-238 (n,gamma) reaction rate divided by MTU
        """
        print(f"\nPlotting U-238 (n,gamma) reaction rate per MTU vs. uranium loading...")

        Pu_per_mtu_list = []
        for i, m in enumerate(self.fertile_mt_list[1:]):
            Pu_per_mtu_list.append(self.U238_ng_list[i+1]/m) # fixed to [i+1] -ezoccoli 2025-07-11

        plt.figure()
        plt.plot(self.fertile_mt_list[1:], Pu_per_mtu_list, markersize=2, linewidth=0.75, color='#000000',) # turn capsize > 0 to show error bars, but they're super smol

        plt.xlim(-1.5,51.5)
        # plt.ylim(-0.005,0.105)
        plt.title(f'Pu-239 production per MTU (Pebble, {self.e}wt% Li-6)')
        plt.xlabel('Uranium loaded [metric tons]')
        plt.ylabel(r'Pu-239 produced per MTU [atoms$/$source neutron$/$MTU]')
        plt.tight_layout()

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_Pu_per_MTU.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_Pu_per_MTU.png', bbox_inches='tight', format='png')
            print(f"Exported U-238 (n,gamma) per MTU plots.")
        else:
            print(f"Did not export U-238 (n,gamma) per MTU plots due to user setting.")

        if self.show: plt.show()
        plt.close('all')


    def plot_pu_per_yr(self):
        """
        Plot of U-238 (n,gamma) reaction rate per year
        """
        print(f"\nPlotting kg of Pu produced per year...")

        try: i = self.fertile_mt_list.index(0.076)
        except Exception as e:
            print(f"\n{e}\n")
            print(f"Fatal: This function requires you to have '0.076' in your `MASS_U_LIST_HCPB`.")

        print(self.fertile_mt_list)
        Pu_rate_per_mtu = self.Pu_per_yr_list[i] / 0.076 # = Pu kg/yr/mtu  # 0.3 / 0.0096 #
        y = [m * Pu_rate_per_mtu for m in self.fertile_mt_list]
        print(f'\nPu-239 production per year per 12.1 kg U: {Pu_rate_per_mtu*0.01209475} | {NPS_FUS} nps')
        print(y)

        plt.figure()
        plt.plot(self.fertile_mt_list, self.Pu_per_yr_list, markersize=2, linewidth=0.75, color='#000000', label=f'OpenMC') # turn capsize > 0 to show error bars, but they're super smol
        plt.plot(self.fertile_mt_list, y, 'r--', linewidth=0.75, label=f'200 wppm U')

        mtu_points = [0.076, 5, 30]
        def get_closest_value(mtu_val, mtu_list, pu_list):
            idx = (abs(mtu_list - mtu_val)).argmin()
            return pu_list[idx]

        # Gather Pu production values for all fuels at each MTU point
        annotations = []
        for mtu in mtu_points:
            pebble_val = get_closest_value(mtu, np.array(self.fertile_mt_list), self.Pu_per_yr_list)
            annotations.append((mtu, pebble_val))
        ax = plt.gca()
        box_props = dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="black", alpha=0.85)
        box_positions = {
            30: (30, max(self.Pu_per_yr_list) * 0.6),
            5: (15, max(self.Pu_per_yr_list) * 0.5),
            0.076: (0.01, max(self.Pu_per_yr_list) * 0.2),
        }

        # Create one box for each MTU value showing all three fuels' Pu production
        for mtu, pebble_val in annotations:
            x, y = box_positions[mtu]
            text = (f"MTU: {mtu}\n"
                    f"Pebble: {pebble_val:.3f} kg/yr")
            ax.text(x, y, text, fontsize=9, bbox=box_props)

        plt.legend()

        plt.xlim(-1.5,51.5)
        # plt.ylim(-0.005,0.105)
        plt.title(f'Pu-239 production per year (Pebble, {self.e}wt% Li-6, {P_FUS_MW} MW = {NPS_FUS:.2e} n/s)')
        plt.xlabel('Uranium loaded [metric tons]')
        plt.ylabel(r'Pu-239 produced [kg$/$yr]')
        plt.tight_layout()

        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_Pu_per_yr.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/{self.name}/fig_Pu_per_yr.png', bbox_inches='tight', format='png')
            print(f"Exported kg of Pu produced per year plots.")
        else:
            print(f"Did not export kg of Pu produced per year plots due to user setting.")

        if self.show: plt.show()
        plt.close('all')


    def plot_rxn_rates(self):
        """
        Plot all reaction rates as function of energies, for every mtu loading, in both log-log and lin-log scales
        """

        print(f"\nPlotting reaction rates vs. energy...")

        for i, cell_id in enumerate(self.cell_ids):

            print(f"...for {self.fertile_mt_list[i]} MTU")
            x = self.flux_df[(self.flux_df['cell'] == cell_id)]['energy mid [eV]']
            y_flux     = self.flux_df[self.flux_df['cell'] == cell_id]['mean']
            y_U238_ng  = self.U238_ng_Ebin_df[self.U238_ng_Ebin_df['cell'] == cell_id]['mean']
            y_U238_fis = self.U238_fis_Ebin_df[self.U238_fis_Ebin_df['cell'] == cell_id]['mean']
            y_Li6_nt   = self.Li6_nt_Ebin_df[self.Li6_nt_Ebin_df['cell'] == cell_id]['mean']
            y_Li7_nt   = self.Li7_nt_Ebin_df[self.Li7_nt_Ebin_df['cell'] == cell_id]['mean']

            plt.figure()
            flux_line, = plt.plot(x, y_flux, linewidth=0.75, color='#000000', label=r'Flux') # black / flux_line, = is to remove it in lin scale bc it doesn't fit
            plt.plot(x, y_Li6_nt,   linewidth=0.75, color='#009ade', label=r'Li-6 (n,t)') # blue
            plt.plot(x, y_U238_ng,  linewidth=0.75, color='#ff1f5b', label=r'U-238 (n,g)') # red
            plt.plot(x, y_U238_fis, linewidth=0.75, color='#f28522', label=r'U-238 (n,fis)') # yellow / #00cd6c = green
            plt.plot(x, y_Li7_nt,   linewidth=0.75, color='#00cd6c', label=r'Li-7 (n,t)') # green
            plt.xlabel('Energy [eV]')
            plt.ylabel('Counts $/$ source neutron')

            """ Export figure in log-log scale """
            plt.title(f'Reaction rates, {self.fertile_mt_list[i]} MTU (Pebble, {self.e}wt% Li-6, log-log)')
            plt.xscale('log'), plt.yscale('log')
            plt.xlim(1e0,2e7), plt.ylim(1e-14,1e2)

            # Reposition legend
            leg = plt.legend(loc='lower right', ncols=2, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1,
            leg.get_frame().set_linewidth(1)

            # Export figure
            if self.save:
                plt.savefig(f'./Figures/pdf/{self.name}/fig_rxnrates_log_{self.fertile_mt_list[i]}mtu.pdf', bbox_inches='tight', format='pdf')
                plt.savefig(f'./Figures/png/{self.name}/fig_rxnrates_log_{self.fertile_mt_list[i]}mtu.png', bbox_inches='tight', format='png') # you want 'log' before 'mtu' so you can flip thru them in File Explorer
                print(f"   Exported LOG-LOG reaction rates plots.")
            if self.show: plt.show()


            """ Export figure in lin-log scale """

            flux_line.remove() # remove flux plot in lin scale bc it doesn't fit

            plt.title(f'Reaction rates, {self.fertile_mt_list[i]} MTU (Pebble, {self.e}wt% Li-6, lin-log)')
            plt.xscale('log'), plt.yscale('linear')
            plt.xlim(1e0,2e7), plt.ylim(0,0.03)

            # Force scientific notation for y-axis
            plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            plt.gca().ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            # plt.gca().yaxis.set_major_locator(MultipleLocator(0.001)) # force y ticks at every integer

            # Reposition legend
            leg = plt.legend(loc='upper left', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1,
            leg.get_frame().set_linewidth(1)

            if self.save:
                plt.savefig(f'./Figures/pdf/{self.name}/fig_rxnrates_lin_{self.fertile_mt_list[i]}mtu.pdf', bbox_inches='tight', format='pdf')
                plt.savefig(f'./Figures/png/{self.name}/fig_rxnrates_lin_{self.fertile_mt_list[i]}mtu.png', bbox_inches='tight', format='png')
                print(f"   Exported LIN-LOG reaction rates plots.")
            if self.show: plt.show()


        """ Make GIF of reaction rates! """

        if self.save:
            filepaths_lin, filepaths_log = [], []
            for i in self.fertile_mt_list:
                filepaths_log.append(f"./Figures/png/{self.name}/fig_rxnrates_log_{i}mtu.png")
                filepaths_lin.append(f"./Figures/png/{self.name}/fig_rxnrates_lin_{i}mtu.png")

            frames_log = [iio.imread(path) for path in filepaths_log]
            frames_lin = [iio.imread(path) for path in filepaths_lin]

            iio.mimsave(f"./Figures/gif/{self.name}/fig_rxnrates_log.gif", frames_log, fps=1, loop=0) # loop=0 : infinite loop
            iio.mimsave(f"./Figures/gif/{self.name}/fig_rxnrates_lin.gif", frames_lin, fps=1, loop=0) # loop=0 : infinite loop
            print(f"Exported GIFs of reaction rates for varying MTU.")

        plt.close('all')


    # noinspection DuplicatedCode
    def plot_pu_vs_energy(self):
        """
        Plots Pu production vs. energy, for contours of MTU
        """
        print(f"\nPlotting Pu production vs. energy with MTU contours ...")

        plt.figure(figsize=(18,4))

        for i, cell_id in enumerate(self.cell_ids):
            if self.fertile_mt_list[i] in [10, 20, 30, 40, 50]:
                x = self.U238_ng_Ebin_df[(self.U238_ng_Ebin_df['cell'] == cell_id)]['energy mid [eV]']
                y  = self.U238_ng_Ebin_df[self.U238_ng_Ebin_df['cell'] == cell_id]['mean']
                plt.plot(x, y,   linewidth=0.75, label=f'{self.fertile_mt_list[i]} MTU') # green color='#00cd6c',

        plt.xlabel('Energy [eV]')
        plt.ylabel('Reactions $/$ source neutron')
        plt.title(f'U-238 (n,gamma) rxn rate (Pebble, {self.e}wt% Li-6)')
        plt.xscale('log'), plt.yscale('linear')
        # plt.xlim(1e0,2e7), plt.ylim(0,0.007)
        plt.xlim(1e1,1e3)# , plt.ylim(1e-7,1e-2)

        # Force scientific notation for y-axis
        # plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        # plt.gca().ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        # Reposition legend
        leg = plt.legend(loc='upper left', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1, 
        leg.get_frame().set_linewidth(1)

        # Export figure
        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_U238ng.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_U238ng.png', bbox_inches='tight', format='png')
            
            plt.xlim(1e0,2e7) # , plt.ylim(1e-7,1e-2)
            plt.savefig(f'./Figures/pdf/{self.name}/fig_U238ng_full.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_U238ng_full.png', bbox_inches='tight', format='png')

            print(f"   Exported Pu production vs. energy with MTU contours plot.")
            
        if self.show: plt.show()
        plt.close('all')


    def plot_rel_pu_vs_energy(self):
        """ Plots factor change in Pu rel 1 MTU vs. energy, for contours of MTU """
        print(f"\nPlotting factor change in Pu rel 1 MTU vs. energy with MTU contours ...")

        plt.figure(figsize=(18,4))

        for i, cell_id in enumerate(self.cell_ids):
            if self.fertile_mt_list[i] in [1]:
                ref_U238_ng = self.U238_ng_Ebin_df[self.U238_ng_Ebin_df['cell'] == cell_id]['mean']

            elif self.fertile_mt_list[i] in [10, 20, 30, 40, 50]:
                x     = self.U238_ng_Ebin_df[(self.U238_ng_Ebin_df['cell'] == cell_id)]['energy mid [eV]']
                y     = self.U238_ng_Ebin_df[self.U238_ng_Ebin_df['cell'] == cell_id]['mean']
                y_rel = [a / b if b != 0.0 else 0.0 for a, b in zip(y, ref_U238_ng)]
                plt.plot(x, y_rel,   linewidth=0.75, label=f'{self.fertile_mt_list[i]} MTU') # green color='#00cd6c',

        plt.xlabel('Energy [eV]')
        plt.ylabel('Factor change relative to rxn rate in 1 MTU')
        plt.title(f'Relative U-238 (n,gamma) rxn rate w.r.t. 1 MTU (Pebble, {self.e}wt% Li-6)')
        plt.xscale('log'), plt.yscale('linear')
        plt.xlim(1e1,1e3), plt.ylim(0,100)

        # Reposition legend
        leg = plt.legend(loc='upper left', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1, 
        leg.get_frame().set_linewidth(1)

        # Export figure
        if self.save:
            plt.savefig(f'./Figures/pdf/{self.name}/fig_U238ng_rel.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_U238ng_rel.png', bbox_inches='tight', format='png')
            
            plt.xlim(1e0,2e7)

            plt.savefig(f'./Figures/pdf/{self.name}/fig_U238ng_rel_full.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./Figures/png/{self.name}/fig_U238ng_rel_full.png', bbox_inches='tight', format='png')

            print(f"   Exported factor change in Pu rel 1 MTU vs. energy with MTU contours plot.")
            
        if self.show: plt.show()
        plt.close('all')


    # noinspection DuplicatedCode
    def plot_flux_vs_energy(self):
        """
        Plots flux vs. energy, for contours of MTU
        """
        print(f"\nPlotting flux vs. energy with MTU contours...")

        plt.figure(figsize=(18,4))
        for i, cell_id in enumerate(self.cell_ids):
            if self.fertile_mt_list[i] in [10, 20, 30, 40, 50]:
                x = self.flux_df[(self.flux_df['cell'] == cell_id)]['energy mid [eV]']
                y = self.flux_df[self.flux_df['cell'] == cell_id]['mean']
                plt.plot(x, y, linewidth=0.75, label=f'{self.fertile_mt_list[i]} MTU') # green color='#00cd6c',

        plt.xlabel('Energy [eV]')
        plt.ylabel('Neutrons $/$ source neutron')
        plt.title(f'Neutron flux (Pebble, {self.e}wt% Li-6)')
        plt.xscale('log'), plt.yscale('log')
        plt.xlim(1e1,1e3) 

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

    def plot_cum_norm_u_vs_energy(self):
        """
        Plots cumulative, normalized Pu production vs. energy for contours of MTU.
        """
        print(f"\nPlotting cumulative, normalized Pu production vs. energy with MTU contours ...")

        plt.figure(figsize=(18,4))

        for i, cell_id in enumerate(self.cell_ids):
            if self.fertile_mt_list[i] in [10, 20, 30, 40, 50]:
                df = self.U238_ng_Ebin_df[self.U238_ng_Ebin_df['cell'] == cell_id]
                x = df['energy mid [eV]']
                y = df['mean']

                # Compute cumulative sum
                cum_y = np.cumsum(y)

                # Normalize cumulative sum to max value
                cum_y_norm = cum_y / cum_y.iloc[-1] if cum_y.iloc[-1] != 0 else cum_y

                plt.plot(x, cum_y_norm, linewidth=0.75, label=f'{self.fertile_mt_list[i]} MTU')

        plt.xlabel('Energy [eV]')
        plt.ylabel('Cumulative normalized reactions')
        plt.title(f'Cumulative normalized U-238 (n,gamma) rxn rate (Li-6 {self.e}wt%, {self.temp} K)')
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



if __name__ == "__main__":
    main()
