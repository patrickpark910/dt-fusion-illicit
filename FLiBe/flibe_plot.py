import openmc
import os, sys, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from scipy.interpolate import make_interp_spline
import imageio.v2 as iio # use v2 to avoid deprecation warnings --ppark

""" Helper functions """
from parameters import *
sys.path.insert(0, f"{os.getcwd()}/helper")
from utilities import *
from flibe_plots_extra import *


def main():
    tbr_paths, pu_paths = [], []
    for e in ENRICH_LI_LIST:
        # Read and plot tallies for each Li enrich case
        current_sp = PlotStatepoint(enrich_li=e, save=True, show=False, to_csv=True)
        
        # Appends path to TBR and Pu data CSVs to plot them all together
        # tbr_paths.append(current_sp.tbr_path)
        # pu_paths.append(current_sp.Pu_path)

    # plot_all_invert(tbr_paths, pu_paths, save=True, show=False, to_csv=True)
    # plot_all_tbr(tbr_paths, save=True, show=False, to_csv=True)
    # plot_all_pu(pu_paths, save=True, show=False, to_csv=True)


class PlotStatepoint:

    def __init__(self, enrich_li=7.5, save=False, show=True, to_csv=False):
        self.e = enrich_li
        self.name = f"FLiBe_Li{self.e:04.1f}"
        self.save, self.show, self.to_csv = save, show, to_csv


        """ Load tallies """
        sp_path = f'./openmc/{self.name}/statepoint.100.h5'
        print("\n\n")
        print("="*42)
        print(f"Loading statepoint: {sp_path}")
        
        try:
            sp = openmc.StatePoint(sp_path) # _Li6-7.5wt  _Li6-20wt
        except:
            sys.exit(f"oopsie woopsie fucky wucky")

        # Make directories for figures 
        for sd in ['pdf','png','gif','data']:
            sd_path = f'./figures/{sd}/{self.name}/'
            print(f"Ensuring directory exists: {sd_path}")
            os.makedirs(f'./figures/{sd}/{self.name}/', exist_ok=True)


        """ Convert tallies into usable forms """

        # Read tallies
        print(f"Reading tallies...")
        flux = sp.get_tally(name='flux')
        U    = sp.get_tally(name='uranium rxn rates')
        Li   = sp.get_tally(name='lithium rxn rates')
        F    = sp.get_tally(name='fluorine rxn rates') # unused
        # Be   = sp.get(name='beryllium rxn rates')

        # Sum over energy bins for total reaction rates
        Li6_rr  = Li.summation(filter_type=openmc.EnergyFilter, nuclides=['Li6'], remove_filter=True)
        Li7_rr  = Li.summation(filter_type=openmc.EnergyFilter, nuclides=['Li7'], remove_filter=True)
        F19_rr  = F.summation(filter_type=openmc.EnergyFilter, nuclides=['F19'], remove_filter=True)
        U235_rr = U.summation(filter_type=openmc.EnergyFilter, nuclides=['U238'], remove_filter=True)
        U238_rr = U.summation(filter_type=openmc.EnergyFilter, nuclides=['U238'], remove_filter=True)

        # Convert Tally objects into pandas dataframes
        flux_df = flux.get_pandas_dataframe() # Convert Tally object! Not 'XXX_rr'! --ppark
        U_df    = U.get_pandas_dataframe()
        Li_df   = Li.get_pandas_dataframe()
        # F_df    = F.get_pandas_dataframe() # unused

        # Add new column for energy bin midpoint (for plotting) 
        for df in [flux_df, U_df, Li_df]: # , F_df]:
            df['energy mid [eV]'] = (df['energy low [eV]'] + df['energy high [eV]'])/ 2

        cell_ids = U_df['cell'].unique().tolist()


        # """ Plot tritium breeding ratio """

        # print(f"Plotting tritium breeding ratio vs. uranium loading...")
        # Li6_df = Li6_rr.get_pandas_dataframe()
        # Li7_df = Li7_rr.get_pandas_dataframe()
        # Li6_nt_df = Li6_df[Li6_df['score'] == '(n,Xt)'][['cell', 'mean', 'std. dev.']]
        # Li7_nt_df = Li7_df[Li7_df['score'] == '(n,Xt)'][['cell', 'mean', 'std. dev.']]

        # self.tbr_list = [x + y for x, y in zip(Li6_nt_df['mean'].tolist(), Li7_nt_df['mean'].tolist())]
        # tbr_err_list = [x + y for x, y in zip(Li6_nt_df['std. dev.'].tolist(), Li7_nt_df['std. dev.'].tolist())]

        # plt.figure()
        # plt.errorbar(MASS_U_LIST, self.tbr_list, yerr=tbr_err_list, fmt='o-', markersize=2, capsize=0, linewidth=0.75, color='#000000',) # turn capsize > 0 to show error bars, but they're super smol
        
        # plt.xlim(-1.5,51.5)
        # plt.ylim(1.196,1.264)
        
        # plt.title(f'Tritium breeding ratio (Li-6 {self.e}wt%)')
        # plt.xlabel('Uranium loaded [metric tons]')
        # plt.ylabel('Tritium breeding ratio')
        # plt.tight_layout()

        # if self.save:
            # plt.savefig(f'./figures/pdf/{self.name}/fig_tbr.pdf', bbox_inches='tight', format='pdf')
            # plt.savefig(f'./figures/png/{self.name}/fig_tbr.png', bbox_inches='tight', format='png')
            # print(f"Exported tritium breeding ratio plots.")
        # else:
            # print(f"Did not export tritium breeding ratio plots.")

        # if self.show:
            # plt.show()

        # # Export dataframe to CSV to plot later
        # if self.to_csv:
            # tbr_df = pd.DataFrame({'mtu':MASS_U_LIST, 'li6_tbr':Li6_nt_df['mean'].tolist(), 'li7_tbr':Li7_nt_df['mean'].tolist(), 'tbr':self.tbr_list, 'tbr_err':tbr_err_list})
            # self.tbr_path = f"./figures/data/{self.name}_tbr.csv"
            # tbr_df.to_csv(self.tbr_path,index=False)
            # print(f"Exported tritium breeding ratio CSV.")

 
        # """ Plot U-238 (n,gamma) reaction rate """

        # print(f"Tabulating U-238 (n,gamma) reaction rate vs. uranium loading...")
        # U238_df = U238_rr.get_pandas_dataframe()
        # U238_ng_df = U238_df[U238_df['score'] == '(n,gamma)'][['cell', 'mean', 'std. dev.']]
        # self.Pu_list, Pu_err_list = U238_ng_df['mean'].tolist(), U238_ng_df['std. dev.'].tolist()
        
        # plt.figure()
        # plt.errorbar(MASS_U_LIST, self.Pu_list, yerr=Pu_err_list, fmt='o-', markersize=2, capsize=0, linewidth=0.75, color='#000000',) # turn capsize > 0 to show error bars, but they're super smol
        # plt.xlim(-1.5,51.5)
        # plt.ylim(-0.005,0.105)

        # plt.title(f'Pu-239 production rate (Li-6 {self.e}wt%)')
        # plt.xlabel('Uranium loaded [metric tons]')
        # plt.ylabel(r'Pu-239 produced [atoms$/$source neutron]')
        # plt.tight_layout()

        # if self.save:
            # plt.savefig(f'./figures/pdf/{self.name}/fig_Pu.pdf', bbox_inches='tight', format='pdf')
            # plt.savefig(f'./figures/png/{self.name}/fig_Pu.png', bbox_inches='tight', format='png')
            # print(f"Exported U-238 (n,gamma) plots.")
        # else:
            # print(f"Did not export U-238 (n,gamma) reaction rate vs. uranium loading plots...")

        # if self.show:
            # plt.show()

        # # Export dataframe to CSV to plot later
        # if self.to_csv:
            # U238_fis_list     = U238_df[U238_df['score'] == 'fission'][['mean']]['mean'].tolist()
            # U238_fis_err_list = U238_df[U238_df['score'] == 'fission'][['std. dev.']]['std. dev.'].tolist()

            # PU_df = pd.DataFrame({'mtu':MASS_U_LIST, 
                                  # 'U-238(n,gamma)':self.Pu_list, 
                                  # 'U-238(n,gamma)_err':Pu_err_list, 
                                  # 'U-238(n,fis)':U238_fis_list,
                                  # 'U-238(n,fis)_err':U238_fis_err_list})
            # self.Pu_path = f"./figures/data/{self.name}_u238.csv"
            # PU_df.to_csv(self.Pu_path, index=False)
            # print(f"Exported U-238 (n,gamma) and (n,fis) CSV.")

 
        # """ Plot of U-238 (n,gamma) reaction rate divided by MTU """        

        # print(f"Tabulating U-238 (n,gamma) reaction rate per MTU vs. uranium loading...")

        # Pu_per_mtu_list = []
        # for i, m in enumerate(MASS_U_LIST[1:]):
            # Pu_per_mtu_list.append(self.Pu_list[i]/m)

        # plt.figure()
        # plt.plot(MASS_U_LIST[1:], Pu_per_mtu_list, markersize=2, linewidth=0.75, color='#000000',) # turn capsize > 0 to show error bars, but they're super smol
        
        # plt.xlim(-1.5,51.5)
        # # plt.ylim(-0.005,0.105)
        # plt.title(f'Pu-239 production per MTU (Li-6 {self.e}wt%)')
        # plt.xlabel('Uranium loaded [metric tons]')
        # plt.ylabel(r'Pu-239 produced per MTU [atoms$/$source neutron$/$MTU]')
        # plt.tight_layout()

        # if self.save:
            # plt.savefig(f'./figures/pdf/{self.name}/fig_Pu_per_MTU.pdf', bbox_inches='tight', format='pdf')
            # plt.savefig(f'./figures/png/{self.name}/fig_Pu_per_MTU.png', bbox_inches='tight', format='png')
            # print(f"Exported U-238 (n,gamma) per MTU plots.")
        # else:
            # print(f"Did not save U-238 (n,gamma) per MTU plots.")
        
        # if self.show:
            # plt.show()

 
        # """ Plot all reaction rates as function of energies, 
        # for every uranium loading, in both log-log and lin-log scales """

        # print(f"Tabulating reaction rates vs. energy...")

        U238_ng_Ebin_df = U_df[(U_df['nuclide'] == 'U238') & (U_df['score'] == '(n,gamma)')][['energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        # U238_fis_Ebin_df = U_df[(U_df['nuclide'] == 'U238') & (U_df['score'] == 'fission')][['energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        # Li6_nt_Ebin_df   = Li_df[(Li_df['nuclide'] == 'Li6') & (Li_df['score'] == '(n,Xt)')][['energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        # Li7_nt_Ebin_df   = Li_df[(Li_df['nuclide'] == 'Li7') & (Li_df['score'] == '(n,Xt)')][['energy mid [eV]', 'cell', 'mean', 'std. dev.']]

        # for i, cell_id in enumerate(cell_ids):

            # print(f"...for {MASS_U_LIST[i]} MTU")

            # x = flux_df[(flux_df['cell'] == cell_id)]['energy mid [eV]']
            # y_flux     = flux_df[flux_df['cell'] == cell_id]['mean']
            # y_U238_ng  = U238_ng_Ebin_df[U238_ng_Ebin_df['cell'] == cell_id]['mean']
            # y_U238_fis = U238_fis_Ebin_df[U238_fis_Ebin_df['cell'] == cell_id]['mean']
            # y_Li6_nt   = Li6_nt_Ebin_df[Li6_nt_Ebin_df['cell'] == cell_id]['mean'] 
            # y_Li7_nt   = Li7_nt_Ebin_df[Li7_nt_Ebin_df['cell'] == cell_id]['mean']

            # plt.figure()

            # flux_line, = plt.plot(x, y_flux, linewidth=0.75, color='#000000', label=r'Flux') # black / flux_line, = is to remove it in lin scale bc it doesn't fit
            # plt.plot(x, y_Li6_nt,   linewidth=0.75, color='#009ade', label=r'Li-6 (n,t)') # blue
            # plt.plot(x, y_U238_ng,  linewidth=0.75, color='#ff1f5b', label=r'U-238 (n,g)') # red
            # plt.plot(x, y_U238_fis, linewidth=0.75, color='#f28522', label=r'U-238 (n,fis)') # yellow / #00cd6c = green
            # plt.plot(x, y_Li7_nt,   linewidth=0.75, color='#00cd6c', label=r'Li-7 (n,t)') # green

            # plt.xlabel('Energy [eV]')
            # plt.ylabel('Counts $/$ source neutron')
            
            # 
            # ''' Export figure in log-log scale ''' 

            # plt.title(f'Reaction rates, {MASS_U_LIST[i]} MTU (Li-6 {self.e}wt%, log-log)')
            # plt.xscale('log'), plt.yscale('log')
            # plt.xlim(1e0,2e7), plt.ylim(1e-14,1e2)

            # # Reposition legend
            # leg = plt.legend(loc='lower right', ncols=2, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1, 
            # leg.get_frame().set_linewidth(1)

            # # Export figure
            # if self.save:
                # plt.savefig(f'./figures/pdf/{self.name}/fig_rxnrates_log_{MASS_U_LIST[i]}mtu.pdf', bbox_inches='tight', format='pdf') 
                # plt.savefig(f'./figures/png/{self.name}/fig_rxnrates_log_{MASS_U_LIST[i]}mtu.png', bbox_inches='tight', format='png') # you want 'log' before 'mtu' so you can flip thru them in File Explorer
                # print(f"   Exported LOG-LOG reaction rates plots.")
                
            # if self.show:
                # plt.show()

            # 
            # ''' Export figure in lin-log scale '''

            # flux_line.remove() # remove flux plot in lin scale bc it doesn't fit
            
            # plt.title(f'Reaction rates, {MASS_U_LIST[i]} MTU (Li-6 {self.e}wt%, lin-log)')
            # plt.xscale('log'), plt.yscale('linear')
            # plt.xlim(1e0,2e7), plt.ylim(0,0.03)
            
            # # Force scientific notation for y-axis
            # plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            # plt.gca().ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            # # plt.gca().yaxis.set_major_locator(MultipleLocator(0.001)) # force y ticks at every integer

            # # Reposition legend
            # leg = plt.legend(loc='upper left', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1, 
            # leg.get_frame().set_linewidth(1)

            # if self.save:
                # plt.savefig(f'./figures/pdf/{self.name}/fig_rxnrates_lin_{MASS_U_LIST[i]}mtu.pdf', bbox_inches='tight', format='pdf')
                # plt.savefig(f'./figures/png/{self.name}/fig_rxnrates_lin_{MASS_U_LIST[i]}mtu.png', bbox_inches='tight', format='png')
                # print(f"   Exported LIN-LOG reaction rates plots.") 
            # if self.show:
                # plt.show()

 
        # """ Make GIF of reaction rates! """

        # if self.save:

            # filepaths_lin, filepaths_log = [], []
            # for i in MASS_U_LIST:
                # filepaths_log.append(f"./figures/png/{self.name}/fig_rxnrates_log_{i}mtu.png")
                # filepaths_lin.append(f"./figures/png/{self.name}/fig_rxnrates_lin_{i}mtu.png")

            # frames_log = [iio.imread(path) for path in filepaths_log]
            # frames_lin = [iio.imread(path) for path in filepaths_lin]

            # iio.mimsave(f"./figures/gif/{self.name}/fig_rxnrates_log.gif", frames_log, fps=1, loop=0) # loop=0 : infinite loop
            # iio.mimsave(f"./figures/gif/{self.name}/fig_rxnrates_lin.gif", frames_lin, fps=1, loop=0) # loop=0 : infinite loop
            # print(f"Exported GIFs of reaction rates for varying MTU.")


        """ Plots Pu production vs. energy, for contours of MTU """

        print(f"Tabulating Pu production vs. energy with MTU contours ...")

        plt.figure(figsize=(18,9))

        for i, cell_id in enumerate(cell_ids):
            if MASS_U_LIST[i] in [0, 10, 20, 30, 40, 50]:
                x_U238_ng = U238_ng_Ebin_df[(U238_ng_Ebin_df['cell'] == cell_id)]['energy mid [eV]']
                y_U238_ng  = U238_ng_Ebin_df[U238_ng_Ebin_df['cell'] == cell_id]['mean']
                plt.plot(x_U238_ng, y_U238_ng,   linewidth=0.75, label=f'{MASS_U_LIST[i]} MTU') # green color='#00cd6c',

        plt.xlabel('Energy [eV]')
        plt.ylabel('Reactions $/$ source neutron')
        plt.title(f'U-238 (n,gamma) rxn rate (Li-6 {self.e}wt%)')
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
            plt.savefig(f'./figures/pdf/{self.name}/fig_U238ng.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./figures/png/{self.name}/fig_U238ng.png', bbox_inches='tight', format='png')
            
            plt.xlim(1e0,2e7), plt.ylim(1e-7,1e-2)

            plt.savefig(f'./figures/pdf/{self.name}/fig_U238ng_full.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./figures/png/{self.name}/fig_U238ng_full.png', bbox_inches='tight', format='png')

            print(f"   Exported Pu production vs. energy with MTU contours plot.")
            
        if self.show:
            plt.show()


        """ Plots factor change in Pu rel 1 MTU vs. energy, for contours of MTU """

        print(f"Tabulating factor change in Pu rel 1 MTU vs. energy with MTU contours ...")

        plt.figure(figsize=(18,9))

        for i, cell_id in enumerate(cell_ids):
            if MASS_U_LIST[i] in [1]:
                ref_U238_ng = U238_ng_Ebin_df[U238_ng_Ebin_df['cell'] == cell_id]['mean'] 

            elif MASS_U_LIST[i] in [10, 20, 30, 40, 50]:
                x_U238_ng = U238_ng_Ebin_df[(U238_ng_Ebin_df['cell'] == cell_id)]['energy mid [eV]']
                y_U238_ng  = U238_ng_Ebin_df[U238_ng_Ebin_df['cell'] == cell_id]['mean']
                y_U238_ng_rel = [a / b if b != 0.0 else 0.0 for a, b in zip(y_U238_ng, ref_U238_ng)]
                plt.plot(x_U238_ng, y_U238_ng_rel,   linewidth=0.75, label=f'{MASS_U_LIST[i]} MTU') # green color='#00cd6c',

        plt.xlabel('Energy [eV]')
        plt.ylabel('Factor change relative to rxn rate in 1 MTU')
        plt.title(f'Relative U-238 (n,gamma) rxn rate w.r.t. 1 MTU (Li-6 {self.e}wt%)')
        plt.xscale('log'), plt.yscale('linear')
        plt.xlim(1e1,1e3), plt.ylim(0,60)

        # Reposition legend
        leg = plt.legend(loc='upper left', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1, 
        leg.get_frame().set_linewidth(1)

        # Export figure
        if self.save:
            plt.savefig(f'./figures/pdf/{self.name}/fig_U238ng_rel.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./figures/png/{self.name}/fig_U238ng_rel.png', bbox_inches='tight', format='png')
            
            plt.xlim(1e0,2e7)

            plt.savefig(f'./figures/pdf/{self.name}/fig_U238ng_rel_full.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./figures/png/{self.name}/fig_U238ng_rel_full.png', bbox_inches='tight', format='png')

            print(f"   Exported factor change in Pu rel 1 MTU vs. energy with MTU contours plot.")
            
        if self.show:
            plt.show()





        """ Plots flux vs. energy, for contours of MTU """

        print(f"Tabulating flux vs. energy with MTU contours...")

        plt.figure(figsize=(18,9))

        for i, cell_id in enumerate(cell_ids):
            if MASS_U_LIST[i] in [0, 10, 20, 30, 40, 50]:
                x_flux = flux_df[(flux_df['cell'] == cell_id)]['energy mid [eV]']
                y_flux = flux_df[flux_df['cell'] == cell_id]['mean']
                plt.plot(x_flux, y_flux,   linewidth=0.75, label=f'{MASS_U_LIST[i]} MTU') # green color='#00cd6c', 

        plt.xlabel('Energy [eV]')
        plt.ylabel('Neutrons $/$ source neutron')
        plt.title(f'Neutron flux (Li-6 {self.e}wt%)')
        plt.xscale('log'), plt.yscale('log')
        # plt.xlim(1e0,2e7), plt.ylim(1e-6,1e2)
        plt.xlim(1e1,1e3), plt.ylim(1e-3,1e0)
        # plt.xlim(1e1,1e3), plt.ylim(0,0.4) # x log y lin

        # Reposition legend
        leg = plt.legend(loc='upper left', ncols=1, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1, 
        leg.get_frame().set_linewidth(1)

        # Export figure
        if self.save:
            plt.savefig(f'./figures/pdf/{self.name}/fig_flux.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./figures/png/{self.name}/fig_flux.png', bbox_inches='tight', format='png') # you want 'log' before 'mtu' so you can flip thru them in File Explorer
            
            plt.xlim(1e0,2e7), plt.ylim(1e-3,1e2)

            plt.savefig(f'./figures/pdf/{self.name}/fig_flux_full.pdf', bbox_inches='tight', format='pdf') 
            plt.savefig(f'./figures/png/{self.name}/fig_flux_full.png', bbox_inches='tight', format='png')

            print(f"   Exported flux vs. energy with MTU contours plot.")
            
        if self.show:
            plt.show()



        """ End loop """
        plt.close('all')



if __name__ == "__main__":
    main()