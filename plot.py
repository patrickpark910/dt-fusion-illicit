import os, sys
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter
from scipy.interpolate import Akima1DInterpolator
from scipy.optimize import curve_fit


""" import helper functions """
from Python.utilities import *
from Python.parameters import *

class Plot:

    def __init__(self, show=True, save=False, ):

        # Dataframes of reaction rates (rr)
        self.flibe_u_rr_df  = read_and_sort("./Figures/Data/FLiBe_900K_Li07.5_U238_rxns.csv")
        self.flibe_th_rr_df = read_and_sort("./Figures/Data/FLiBe_900K_Li07.5_Th232_rxns.csv")
        self.dcll_u_rr_df   = read_and_sort("./Figures/Data/DCLL_900K_Li90.0_U238_rxns.csv")
        self.dcll_th_rr_df  = read_and_sort("./Figures/Data/DCLL_900K_Li90.0_Th232_rxns.csv")
        self.hcpb_u_rr_df   = read_and_sort("./Figures/Data/HCPB_900K_Li60.0_U238_rxns.csv")
        self.hcpb_th_rr_df  = read_and_sort("./Figures/Data/HCPB_900K_Li60.0_Th232_rxns.csv")

        self.flibe_u_ebins_df  = read_and_sort("./Figures/Data/FLiBe_900K_Li07.5_U238_Ebins_rxns.csv")
        self.flibe_th_ebins_df = read_and_sort("./Figures/Data/FLiBe_900K_Li07.5_Th232_Ebins_rxns.csv")
        self.dcll_u_ebins_df   = read_and_sort("./Figures/Data/DCLL_900K_Li90.0_U238_Ebins_rxns.csv")
        self.dcll_th_ebins_df  = read_and_sort("./Figures/Data/DCLL_900K_Li90.0_Th232_Ebins_rxns.csv")
        self.hcpb_u_ebins_df   = read_and_sort("./Figures/Data/HCPB_900K_Li60.0_U238_Ebins_rxns.csv")
        self.hcpb_th_ebins_df  = read_and_sort("./Figures/Data/HCPB_900K_Li60.0_Th232_Ebins_rxns.csv")

        # Dataframes of energy-binned flux spectra (for flux ratio plot)
        self.flibe_u_flux_df  = read_and_sort("./Figures/Data/FLiBe_900K_Li07.5_U238_flux.csv")
        self.flibe_th_flux_df = read_and_sort("./Figures/Data/FLiBe_900K_Li07.5_Th232_flux.csv")
        self.dcll_u_flux_df   = read_and_sort("./Figures/Data/DCLL_900K_Li90.0_U238_flux.csv")
        self.dcll_th_flux_df  = read_and_sort("./Figures/Data/DCLL_900K_Li90.0_Th232_flux.csv")
        self.hcpb_u_flux_df   = read_and_sort("./Figures/Data/HCPB_900K_Li60.0_U238_flux.csv")
        self.hcpb_th_flux_df  = read_and_sort("./Figures/Data/HCPB_900K_Li60.0_Th232_flux.csv")

        self.flibe_u_rr_df = self.flibe_u_rr_df[self.flibe_u_rr_df['fertile_kg/m3'] <= LIMIT_UF4_FLIBE]
        self.flibe_th_rr_df = self.flibe_th_rr_df[self.flibe_th_rr_df['fertile_kg/m3'] <= LIMIT_THF4_FLIBE]
        self.flibe_u_ebins_df = self.flibe_u_ebins_df[self.flibe_u_ebins_df['fertile_kg/m3'] <= LIMIT_UF4_FLIBE]
        self.flibe_th_ebins_df = self.flibe_th_ebins_df[self.flibe_th_ebins_df['fertile_kg/m3'] <= LIMIT_THF4_FLIBE]

        # ------------------------------------------------------------------
        # Weigh TBRs and FPRs by correction factor derived from wedge models
        # -------------------------------------------------------------------

        def get_ratio(x, array):
            """
            Bc you might be trying to plot a fertile_kg/m3 not listed in XXXX_CONV_X_TBR <parameters.py>
            (e.g., 350 kg/m³) we piecewise linearly interpolate the correction ratio for any fertile_kg/m3 x

            Usage: np.interp(x_value_to_calculate, original_x, original_y)
            """
            df_fertile = pd.DataFrame(array, columns=['fertile_kgm3', 'ratio'])
            return np.interp(x, df_fertile['fertile_kgm3'], df_fertile['ratio'])

        # Calculate corrected TBR values by multiplying the appropriate correction factor for each fertile_kg/m3
        self.hcpb_u_rr_df_corr  = self.hcpb_u_rr_df.copy()
        self.hcpb_th_rr_df_corr = self.hcpb_th_rr_df.copy()
        self.hcpb_u_rr_df_corr['U238(n,g)']   *= get_ratio(self.hcpb_u_rr_df_corr['fertile_kg/m3'],  HCPB_CONV_U_FPR)
        self.hcpb_th_rr_df_corr['Th232(n,g)'] *= get_ratio(self.hcpb_th_rr_df_corr['fertile_kg/m3'], HCPB_CONV_TH_FPR)        

        self.show, self.save = show, save


    def plot_tbr(self):

        print(f"\nComment. <plot.py/plot_tbr()> Plotting tritium breeding ratio vs. fertile density...")

        # -------------------------------------------------------------
        # Plot them!
        # -------------------------------------------------------------

        # Load the dataframes to plot, label, color, marker, markersize, linestyle, and polynomial fit
        LONG_DASH = (0, (10, 2))
        datasets = [ (self.flibe_u_rr_df,  r'FLiBe-UF$_4$',  '#66b420', 'o',  5,  '-',  'makima'),
                     (self.hcpb_u_rr_df,   r'HCPB-UO$_2$',   '#b41f24', 's',  6,  '-',  2),
                     (self.dcll_u_rr_df,   r'DCLL-UO$_2$',   '#0047ba', '^',  8,  '-',  2),
                     (self.flibe_th_rr_df, r'FLiBe-ThF$_4$', '#66b420', '+', 12, LONG_DASH,  'makima'),
                     (self.hcpb_th_rr_df,  r'HCPB-ThO$_2$',  '#b41f24', '1', 13, LONG_DASH,  2),
                     (self.dcll_th_rr_df,  r'DCLL-ThO$_2$',  '#0047ba', 'x', 10, LONG_DASH,  2),
                   ]

        # Select a subset of fertile kg/m3 to SHOW to avoid cluttering the low end with markers
        selected = [0, 25, 50, 75, 100, 150, 250, 500, 750, 999.99]  # [0.0, 0.10, 0.50, 1, 10, 25, 50, 75, 100, 150, 250, 500, 750, 999.99] 

        # Settings

        v = {'xlim': (-25, 1025), 'xmax':1000, 'x_major': 100, 'x_minor': 50, 
             'ylim': (0.73, 1.47), 'suffix': '1000kgm3', 'fertile_kgm3': selected,
             'leg_loc': 'upper right'}
  
        plt.figure(figsize=(3.5, 3.0))  # 3.5 x 2.75
        ax = plt.gca()
        ax.axhspan(0, 1.00, color='#F0F0F0')            

        for df, label, color, marker, markersize, linestyle, degree in datasets:

            # selected = v['fertile_kgm3']
            mask = df['fertile_kg/m3'].apply(lambda d: any(np.isclose(d, s) for s in selected))
            df_filtered = df[mask]

            x_filtered = df_filtered['fertile_kg/m3']
            y_filtered = BLANKET_COVERAGE * df_filtered['tbr']

            # Plot filtered data
            # plt.scatter(x_filtered, y_filtered, marker=marker, s=markersize*5, color=color, zorder=3)

            # But still fit the polynomial to all of the data you did compute
            x_fit = df['fertile_kg/m3']
            y_fit = BLANKET_COVERAGE * df['tbr']

            if degree == 'makima':
                x_fine = np.linspace(0, v['xmax'], 500)
                y_fine = Akima1DInterpolator(x_fit, y_fit, method='makima')(x_fine) # poly_func(x_fine)

            else:
                coeffs = np.polyfit(x_fit, y_fit, degree)
                y_fine = np.poly1d(coeffs)(x_fine)

            # Plot interpolation
            plt.plot(x_fine, y_fine, linestyle=linestyle, linewidth=0.75, color=color)

            # Dummy plots for legend -- bit of a hack lmao -- ppark
            # plt.plot([9e8,9e9], [9e8,9e9], marker+linestyle, markersize=markersize, linewidth=0.75, color=color, label=label)

        # Labeling
        plt.xlabel(r'Fertile isotope density [kg$/$m³]')  # specifically use unicode superscript m³ and not m$^3$
        plt.ylabel('Tritium breeding ratio')

        # Set Dynamic Limits and Locators
        plt.xlim(v['xlim'])
        plt.ylim(v['ylim']) 
        
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_major_locator(MultipleLocator(v['x_major']))
        ax.xaxis.set_minor_locator(MultipleLocator(v['x_minor']))
        
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(MultipleLocator(0.025))
        
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

        plt.tight_layout()
        # leg = plt.legend(loc=v['leg_loc'], fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=3)
        # leg = plt.legend(loc=v['leg_loc'], fancybox=False, edgecolor='none', frameon=True, framealpha=1.0, ncol=3)
        # leg.get_frame().set_linewidth(0.5) 

        if self.save:
            plt.savefig(f'./Figures/PDF/fig_tbr_{v["suffix"]}.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig(f'./Figures/PNG/fig_tbr_{v["suffix"]}.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print(f"Comment. <plot.py/plot_tbr()> Exported TBR plot: fig_tbr_{v['suffix']}")
        else:
            print(f"{C.YELLOW}Comment.{C.END} <plot.py/plot_tbr()> Did NOT export TBR plot due to user setting: fig_tbr_{v['suffix']}")

        if self.show: 
            plt.show()
        else:
            print(f"{C.YELLOW}Comment.{C.END} <plot.py/plot_tbr()> Did NOT show TBR plot due to user setting: fig_tbr_{v['suffix']}")

        plt.close('all')


    def plot_fpr(self):
      
        print(f"\nPlotting Pu-239 production per year for all blankets...")

        # Constants and conversion factors
        # <utilities.py> default NPS_FUS = 500 MJ/s * 3.546e17 n/MJ = 1.773e+20 n/s
        tot_n_per_yr = NPS_FUS * 3.156e+7
        pu239_conv   = (tot_n_per_yr / AVO * AMU_PU239 / 1e3) * BLANKET_COVERAGE
        u233_conv    = (tot_n_per_yr / AVO * AMU_U233  / 1e3) * BLANKET_COVERAGE

        # -----------------------------------------------------------------------------------------
        # Dataset Configuration: 
        # (dataframe, label, color, marker, markersize, linestyle, reaction_key, conversion_factor)
        # -----------------------------------------------------------------------------------------
        
        datasets = [ (self.flibe_th_rr_df,      r'FLiBe-ThF$_4$', '#66b420', '+', 12, LONG_DASH, 'makima', 'Th232(n,g)', u233_conv),
                     (self.hcpb_th_rr_df_corr,  r'HCPB-ThO$_2$',  '#b41f24', '1', 13, LONG_DASH, 'makima', 'Th232(n,g)', u233_conv),
                     (self.dcll_th_rr_df,       r'DCLL-ThO$_2$',  '#0047ba', 'x', 10, LONG_DASH, 'makima', 'Th232(n,g)', u233_conv),
                     (self.flibe_u_rr_df,       r'FLiBe-UF$_4$',  '#66b420', 'o',  5,       '-', 'makima',  'U238(n,g)', pu239_conv),
                     (self.hcpb_u_rr_df_corr,   r'HCPB-UO$_2$',   '#b41f24', 's',  6,       '-', 'makima',  'U238(n,g)', pu239_conv),
                     (self.dcll_u_rr_df,        r'DCLL-UO$_2$',   '#0047ba', '^',  8,       '-', 'makima',  'U238(n,g)', pu239_conv),      ]

        # Select a subset of fertile kg/m3 to SHOW to avoid cluttering the low end with markers
        selected1 = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150,]  # [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 250, 500, 750, 999.99]
        selected2 = [0, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]  # [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 250, 500, 750, 999.99]

        # Set view limits for each of the plots
        view_limits = [ # Zoomed in on (0, 155) fertile kg/m3
                        # {'xlim': (-2.5, 152.5),   'xmax': 150, 'x_major': 15,  'x_minor': 5,  
                        #  'ylim': (-15, 565), 'y_major': 50, 'y_minor':25 ,
                        #  'suffix': '0150kgm3', 'fertile_kgm3': selected1,
                        #  'leg_loc': 'upper left'},
                        # Zoomed out to (0, 1000) fertile kg/m3
                        {'xlim': (-25, 1025), 'xmax':1000, 'x_major': 100, 'x_minor': 50, 
                         'ylim': (-50, 2250), 'y_major': 200, 'y_minor':100,
                         'suffix': '1000kgm3', 'fertile_kgm3': selected2,
                         'leg_loc': 'upper left'} ]


        for v in view_limits:

            plt.figure(figsize=(3.5, 3.0))  # 3.5 x 2.75
            ax = plt.gca()
            x_fine = np.linspace(0, 1000, 500)

            for df, label, color, marker, markersize, linestyle, degree, rxn, conv, in datasets:

                selected = v['fertile_kgm3']
                mask = df['fertile_kg/m3'].apply(lambda d: any(np.isclose(d, s) for s in selected))
                df_filtered = df[mask]

                x_filtered = df_filtered['fertile_kg/m3']
                y_filtered = df_filtered[rxn] * conv

                # Plot filtered data
                # plt.scatter(x_filtered, y_filtered, marker=marker, s=markersize*5, color=color, zorder=3)

                # But still fit the polynomial to all of the data you did compute
                x_fit = df['fertile_kg/m3']
                y_fit = df[rxn] * conv

                if degree == 'makima':
                    x_fine = np.linspace(0, v['xmax'], 500)
                    y_fine = Akima1DInterpolator(x_fit, y_fit, method='makima')(x_fine) # poly_func(x_fine)

                else:
                    coeffs = np.polyfit(x_fit, y_fit, degree)
                    y_fine = np.poly1d(coeffs)(x_fine)

                # Plot interpolation
                plt.plot(x_fine, y_fine, linestyle=linestyle, linewidth=0.75, color=color)

                # Dummy plots for legend -- bit of a hack lmao -- ppark
                # plt.plot([9e8,9e9], [9e8,9e9], marker+linestyle, markersize=markersize, linewidth=0.75, color=color, label=label)

            # Labeling and Titles
            plt.xlabel(r"Fertile isotope density [kg$/$m³]")  # specifically use unicode superscript m³ and not m$^3$
            plt.ylabel(r"Initial fissile production rate [kg$/$yr]")

            # Set Dynamic limits and locators
            plt.xlim(v['xlim'])
            plt.ylim(v['ylim']) 
            
            ax.xaxis.set_ticks_position('both')
            ax.xaxis.set_major_locator(MultipleLocator(v['x_major']))
            ax.xaxis.set_minor_locator(MultipleLocator(v['x_minor']))
            
            ax.yaxis.set_ticks_position('both')
            ax.yaxis.set_major_locator(MultipleLocator(v['y_major']))
            ax.yaxis.set_minor_locator(MultipleLocator(v['y_minor']))
            
            ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
            ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

            plt.tight_layout()
            # leg = plt.legend(loc=v['leg_loc'], fancybox=False, edgecolor='black', frameon=False, framealpha=.75, ncol=2)
            # leg.get_frame().set_linewidth(0.5) 

            if self.save:
                plt.savefig(f'./Figures/PDF/fig_fpr_{v["suffix"]}.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
                plt.savefig(f'./Figures/PNG/fig_fpr_{v["suffix"]}.png', bbox_inches='tight', pad_inches=0.01, format='png')
                print("Exported fissile production per year plot for all blankets.")
            else:
                print("Did not export fissile production per year plot due to user setting.")

            if self.show: 
                plt.show()
                
            plt.close('all')



    def plot_hist(self):
        """
        Plots cumulative, normalized Pu production vs. energy for contours of MTU.
        """
        print(f"\nPlotting cumulative, normalized fissile production vs. energy...")

        fig, axes = plt.subplots(3, 2, figsize=(7,9), sharex=True, sharey=True)

        # load in thorium and uranium background (n, gamma) cross section, shift to 0.1 and 0.9
        u_path = "./Figures/XSPlot/U238gamma.txt"
        th_path = "./Figures/XSPlot/Th232gamma.txt"
        u238_energy, u238_mxs = readtxtFile(u_path)
        th232_energy, th232_mxs = readtxtFile(th_path)

        u238_mxs_shifted = (u238_mxs - np.min(u238_mxs)) * 0.8 / (np.max(u238_mxs) - np.min(u238_mxs)) + 0.1
        th232_mxs_shifted = (th232_mxs - np.min(th232_mxs)) * 0.8 / (np.max(th232_mxs) - np.min(th232_mxs)) + 0.1

        for ax in axes.flatten():
            if ax.get_subplotspec().is_last_row():
                ax.set_xlabel("Incident neutron energy [eV]")
            if ax.get_subplotspec().is_first_col():
                ax.set_ylabel(r"Cumulative fraction of fertile $($n,$\gamma)$ rxns")
                ax.plot(u238_energy, u238_mxs_shifted, linewidth=1.0, color='gray', alpha=0.25, label=fr'U-238 $($n,$\gamma)$')
            else:
                ax.plot(th232_energy, th232_mxs_shifted, linewidth=1.0, color='gray', alpha=0.25, label=fr'Th-232 $($n,$\gamma)$')

        titles  = [r"FLiBe-UF$_4$", r"FLiBe-ThF$_4$", r"HCPB-UO$_2$", r"HCPB-ThO$_2$", r"DCLL-UO$_2$", r"DCLL-ThO$_2$",]
        dfs     = [self.flibe_u_ebins_df, self.flibe_th_ebins_df, self.hcpb_u_ebins_df, self.hcpb_th_ebins_df, self.dcll_u_ebins_df, self.dcll_th_ebins_df,]
        ng_cols = ['U238_ng', 'Th232_ng', 'U238_ng', 'Th232_ng', 'U238_ng', 'Th232_ng']

        densities_to_plot = [0.1, 10, 100, 999.99]
        colors = {0.1: '#ff1f5b', 10: '#f48628', 100: '#04cc6c', 999.99: '#0c9edd'}
        labels = {0.1: r'0.1 kg$/$m³', 10: r'10 kg$/$m³', 100: r'100 kg$/$m³', 999.99: r'1000 kg$/$m³'}

        print(f"\n{'Panel':<20} {'Density [kg/m³]':>16} {'Mean [eV]':>12} {'Median [eV]':>13} {'Mode [eV]':>12}")
        print("-" * 75)

        for ax, df, title, ng_col in zip(axes.flatten(), dfs, titles, ng_cols):

            df = df[df["fertile_kg/m3"].isin(densities_to_plot)]

            df_sum = df.groupby("fertile_kg/m3")[ng_col].sum().to_frame()
            df_sum.columns = ["sum"]   
            df = df.merge(df_sum, on="fertile_kg/m3", how="left")
            df['norm_mean'] = df[ng_col] / df['sum']

            bins = np.sort(df['energy mid [eV]'].unique())

            for d in densities_to_plot:
                sub = df[df['fertile_kg/m3'] == d]
                ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', linewidth=0.75, color=colors[d], label=labels[d])

                # Compute mean, median, mode
                energies = sub['energy mid [eV]'].values
                weights = sub['norm_mean'].values
                idx = np.argsort(energies)
                energies, weights = energies[idx], weights[idx]

                w_mean = np.average(energies, weights=weights)
                cumw = np.cumsum(weights)
                w_median = energies[np.searchsorted(cumw, 0.5 * cumw[-1])]
                w_mode = energies[np.argmax(weights)]

                # Strip LaTeX for clean printing
                clean_title = title.replace(r'$_4$', '4').replace(r'$_2$', '2')
                d_label = f"{d:.1f}" if d < 999 else "1000.0"
                print(f"{clean_title:<20} {d_label:>16} {w_mean:>12.2e} {w_median:>13.2e} {w_mode:>12.2e}")

            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.yaxis.set_minor_locator(MultipleLocator(0.05))
            ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
            ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'1e{int(np.log10(x))}' if x > 0 else ''))
            ax.set_xlim(0.5*1e1, 1.5*1e6)
            ax.set_ylim(-0.03, 1.03)
            fig.tight_layout()

            leg = ax.legend(title=title, fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=1, loc="lower right")
            leg.get_frame().set_linewidth(0.5) 

        if self.save:
            plt.savefig(f'./Figures/PDF/fig_hist.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig(f'./Figures/PNG/fig_hist.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print(f"\nExported cumulative normalized histogram plots.")
        else:
            print(f"\nDid not export cumulative normalized histogram plot.")

        if self.show: plt.show()
        plt.close('all')


    def plot_flux_ratio(self):
        """
        Flux spectral ratio plot: phi(rho) / phi(rho_ref) vs. energy.

        3x2 grid (rows = FLiBe / HCPB / DCLL, cols = U238 / Th232).
        Each panel overlays the ratio at several fertile loadings relative
        to the clean (zero fertile) reference case.  A ratio > 1
        means more flux at that energy; < 1 means less.

        In FLiBe and HCPB, the ratio dips below 1 in the resonance region
        (fertile absorber depleting resonance-energy neutrons) and rises
        above 1 at higher energies.  In DCLL, the ratio dips above ~50 keV
        (fast flux depletion by the growing absorber + its inelastic
        scattering) and rises in the 1--50 keV range (fertile
        self-moderation feeding neutrons into the resonance region).
        """
        print(f"\nPlotting flux spectral ratios vs. energy...")

        panels = [
            (self.flibe_u_flux_df,  r'FLiBe-UF$_4$',  0, 0),
            (self.flibe_th_flux_df, r'FLiBe-ThF$_4$', 0, 1),
            (self.hcpb_u_flux_df,   r'HCPB-UO$_2$',   1, 0),
            (self.hcpb_th_flux_df,  r'HCPB-ThO$_2$',  1, 1),
            (self.dcll_u_flux_df,   r'DCLL-UO$_2$',   2, 0),
            (self.dcll_th_flux_df,  r'DCLL-ThO$_2$',  2, 1),
        ]

        rho_ref = 0                         # reference density: clean (no fertile) case
        rho_compare = [10, 100, 999.99]     # densities to compare # 0.1, 
        colors = {0.1: '#ff1f5b', 10: '#f48628', 100: '#04cc6c', 999.99: '#0c9edd'}
        labels = {0.1: r'0.1 kg$/$m³', 10: r'10 kg$/$m³', 100: r'100 kg$/$m³', 999.99: r'1000 kg$/$m³'}

        fig, axes = plt.subplots(3, 2, figsize=(7, 9), sharex=True, sharey=False)

        # Per-row y-axis limits and tick spacing:
        #   rows 0-1 (FLiBe, HCPB): wider range to capture resonance depletion
        #   row  2   (DCLL):         tighter range around unity
        row_ylim  = {0: (0.15*0.97, 1.3*1.03), 1: (0.15*0.97, 1.3*1.03), 2: (0.50*0.97, 5.5*1.03)}
        row_major = {0: 0.1,          1: 0.1,           2: 0.5}
        row_minor = {0: 0.05,         1: 0.05,          2: 0.25}

        for df_flux, title, row, col in panels:
            ax = axes[row, col]

            # Match closest available density to rho_ref
            available = df_flux['fertile_kg/m3'].unique()
            ref_rho = min(available, key=lambda x: abs(x - rho_ref))
            ref = df_flux[df_flux['fertile_kg/m3'] == ref_rho].sort_values('energy mid [eV]')
            E_mid_ref = ref['energy mid [eV]'].values
            phi_ref   = ref['mean'].values

            for rho in rho_compare:
                comp_rho = min(available, key=lambda x: abs(x - rho))
                comp = df_flux[df_flux['fertile_kg/m3'] == comp_rho].sort_values('energy mid [eV]')
                phi_comp = comp['mean'].values

                # Only compute ratio where reference flux is non-negligible
                # (avoids 0/0 in thermal bins for DCLL where flux is essentially zero)
                ratio = np.ones_like(phi_ref, dtype=float)
                mask = phi_ref > 1e-20
                ratio[mask] = phi_comp[mask] / phi_ref[mask]
                ratio[~mask] = np.nan

                ax.semilogx(E_mid_ref, ratio,
                            color=colors[rho], linewidth=0.75, label=labels[rho])

            # Reference line at ratio = 1
            ax.axhline(1.0, color='gray', linewidth=0.5, linestyle='-', zorder=0)

            # Axis formatting
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.grid(axis='both', which='major', linestyle='-', linewidth=0.5)
            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'1e{int(np.log10(x))}' if x > 0 else ''))
            ax.set_xlim(log_buffer(1e1,1e7))
            ax.set_ylim(row_ylim[row])
            ax.yaxis.set_major_locator(MultipleLocator(row_major[row]))
            ax.yaxis.set_minor_locator(MultipleLocator(row_minor[row]))

            if ax.get_subplotspec().is_last_row():
                ax.set_xlabel('Neutron energy [eV]')
            if ax.get_subplotspec().is_first_col():
                ax.set_ylabel(r'Flux ratio $\phi(\rho)\,/\,\phi_0$')

            # Panel title
            # ax.text(0.02, 0.97, title, transform=ax.transAxes, va='top', ha='left')

            # Legend: top for FLiBe/HCPB (rows 0-1), bottom-right for DCLL (row 2)
            leg_loc = 'lower right' if row <= 1 else 'upper right'
            leg = ax.legend(title=title, loc=leg_loc, fontsize=7, ncol=1,
                            fancybox=False, edgecolor='black',
                            frameon=True, framealpha=0.75)
            leg.get_frame().set_linewidth(0.5)

        fig.tight_layout()

        if self.save:
            plt.savefig('./Figures/PDF/fig_flux_ratio.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig('./Figures/PNG/fig_flux_ratio.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print("Exported flux spectral ratio plot: fig_flux_ratio")
        else:
            print("Did not export flux ratio plot due to user setting.")

        if self.show:
            plt.show()
        plt.close('all')


    def plot_flux_spectrum(self):
        """
        Flux spectrum plot: phi(E) vs. energy.

        3x2 grid (rows = FLiBe / HCPB / DCLL, cols = U238 / Th232).
        Each panel overlays the flux at several fertile loadings
        (0, 10, 100, 1000 kg/m³).
        """
        print(f"\nPlotting flux spectra vs. energy...")

        panels = [
            (self.flibe_u_flux_df,  r'FLiBe-UF$_4$',  0, 0),
            (self.flibe_th_flux_df, r'FLiBe-ThF$_4$', 0, 1),
            (self.hcpb_u_flux_df,   r'HCPB-UO$_2$',   1, 0),
            (self.hcpb_th_flux_df,  r'HCPB-ThO$_2$',  1, 1),
            (self.dcll_u_flux_df,   r'DCLL-UO$_2$',   2, 0),
            (self.dcll_th_flux_df,  r'DCLL-ThO$_2$',  2, 1),
        ]

        densities_to_plot = [0, 10, 100, 999.99]
        colors = {0: '#ff1f5b', 10: '#f48628', 100: '#04cc6c', 999.99: '#0c9edd'}
        labels = {0: r'0 kg$/$m³', 10: r'10 kg$/$m³', 100: r'100 kg$/$m³', 999.99: r'1000 kg$/$m³'}
        row_ylim = {0: log_buffer(1e-4, 1e1), 1: log_buffer(1e-4, 1e1), 2: log_buffer(1e-6, 1e1)}

        fig, axes = plt.subplots(3, 2, figsize=(7, 9), sharex=True, sharey=False)
        

        for df_flux, title, row, col in panels:
            ax = axes[row, col]

            available = df_flux['fertile_kg/m3'].unique()

            for rho in densities_to_plot:
                closest = min(available, key=lambda x: abs(x - rho))
                spec = df_flux[df_flux['fertile_kg/m3'] == closest].sort_values('energy mid [eV]')

                E_mid = spec['energy mid [eV]'].values
                phi   = spec['mean'].values.copy()

                # Mask zeros / negatives to avoid log-scale artifacts
                phi[phi <= 0] = np.nan

                ax.loglog(E_mid, phi,
                          color=colors[rho], linewidth=0.75, label=labels[rho])

            # Axis formatting
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.grid(axis='both', which='major', linestyle='-', linewidth=0.5)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'1e{int(np.log10(x))}' if x > 0 else ''))
            ax.set_xlim(log_buffer(1e1,1e7))
            ax.set_ylim(row_ylim[row])

            from matplotlib.ticker import LogLocator
            ax.yaxis.set_major_locator(LogLocator(base=10, numticks=20))
            ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'1e{int(np.log10(x))}' if x > 0 else ''))

            if ax.get_subplotspec().is_last_row():
                ax.set_xlabel('Neutron energy [eV]')
            if ax.get_subplotspec().is_first_col():
                ax.set_ylabel(r'Flux [n$/$cm²$/$src-n]')

            # Legend with panel title
            leg = ax.legend(title=title, loc='lower right', fontsize=7, ncol=1,
                            fancybox=False, edgecolor='black',
                            frameon=True, framealpha=1)
            leg.get_frame().set_linewidth(0.5)

        fig.tight_layout()

        if self.save:
            plt.savefig('./Figures/PDF/fig_flux_spectrum.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig('./Figures/PNG/fig_flux_spectrum.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print("Exported flux spectrum plot: fig_flux_spectrum")
        else:
            print("Did not export flux spectrum plot due to user setting.")

        if self.show:
            plt.show()
        plt.close('all')


    def plot_hist_diff(self):
        """
        Differential cumulative capture distributions:
        CDF(rho) - CDF(rho_ref) vs. energy.

        3x2 grid (rows = FLiBe / HCPB / DCLL, cols = U238 / Th232).
        Each panel shows the difference in cumulative fractional (n,gamma)
        distribution relative to the lowest-density case (0.1 kg/m3).

        A positive region means that density's CDF runs ahead of the
        reference (more captures have occurred by that energy), i.e.
        captures shifted *downward* in energy.  A negative region means
        the CDF lags (captures shifted *upward*).

        Expected signatures:
          FLiBe / HCPB: negative in the resonance region, positive above
                        -- resonance saturation pushes captures upward.
          DCLL:         positive in the 1--50 keV region, negative above
                        -- fertile self-moderation pushes captures downward.
        """
        print(f"\nPlotting differential cumulative capture distributions...")

        fig, axes = plt.subplots(3, 2, figsize=(7, 9), sharex=True, sharey=False)

        # Per-row y-axis limits and tick spacing
        row_ylim  = {0: (-0.30, 0.10), 1: (-0.20, 0.10), 2: (-0.05, 0.12)}
        row_major = {0: 0.05,          1: 0.05,           2: 0.02}
        row_minor = {0: 0.025,         1: 0.025,          2: 0.01}

        panels = [
            (self.flibe_u_ebins_df,  'U238_ng',  r'FLiBe-UF$_4$',  0, 0),
            (self.flibe_th_ebins_df, 'Th232_ng', r'FLiBe-ThF$_4$', 0, 1),
            (self.hcpb_u_ebins_df,   'U238_ng',  r'HCPB-UO$_2$',   1, 0),
            (self.hcpb_th_ebins_df,  'Th232_ng', r'HCPB-ThO$_2$',  1, 1),
            (self.dcll_u_ebins_df,   'U238_ng',  r'DCLL-UO$_2$',   2, 0),
            (self.dcll_th_ebins_df,  'Th232_ng', r'DCLL-ThO$_2$',  2, 1),
        ]

        rho_ref = 0.1                            # reference density for CDF difference
        rho_compare = [10, 100, 999.99]
        colors = {10: '#f48628', 100: '#04cc6c', 999.99: '#0c9edd'}
        labels = {10: r'10 kg$/$m³', 100: r'100 kg$/$m³', 999.99: r'1000 kg$/$m³'}

        for df_ebin, ng_col, title, row, col in panels:
            ax = axes[row, col]

            # --- Build the reference CDF at rho_ref ---
            available = df_ebin['fertile_kg/m3'].unique()
            ref_rho = min(available, key=lambda x: abs(x - rho_ref))
            ref = df_ebin[df_ebin['fertile_kg/m3'] == ref_rho].sort_values('energy mid [eV]')
            E_mid = ref['energy mid [eV]'].values
            ref_weights = ref[ng_col].values
            ref_total = ref_weights.sum()
            cdf_ref = np.cumsum(ref_weights) / ref_total if ref_total > 0 else np.zeros_like(ref_weights)

            for rho in rho_compare:
                comp_rho = min(available, key=lambda x: abs(x - rho))
                comp = df_ebin[df_ebin['fertile_kg/m3'] == comp_rho].sort_values('energy mid [eV]')
                comp_weights = comp[ng_col].values
                comp_total = comp_weights.sum()
                cdf_comp = np.cumsum(comp_weights) / comp_total if comp_total > 0 else np.zeros_like(comp_weights)

                delta_cdf = cdf_comp - cdf_ref

                ax.semilogx(E_mid, delta_cdf,
                            color=colors[rho], linewidth=0.75, label=labels[rho])

            # Reference line at zero
            ax.axhline(0.0, color='gray', linewidth=0.5, linestyle='-', zorder=0)

            # Axis formatting
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.grid(axis='both', which='major', linestyle='-', linewidth=0.5)
            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'1e{int(np.log10(x))}' if x > 0 else ''))
            ax.set_xlim(0.5e1, 1.5e6)
            ax.set_ylim(row_ylim[row])
            ax.yaxis.set_major_locator(MultipleLocator(row_major[row]))
            ax.yaxis.set_minor_locator(MultipleLocator(row_minor[row]))

            if ax.get_subplotspec().is_last_row():
                ax.set_xlabel('Neutron energy [eV]')
            if ax.get_subplotspec().is_first_col():
                ax.set_ylabel(r'$\Delta$ CDF$(E)$  vs. 0.1 kg$/$m³')

            # Legend with panel title
            leg_loc = 'lower left' if row <= 1 else 'upper left'
            leg = ax.legend(title=title, loc=leg_loc, fontsize=7, ncol=1,
                            fancybox=False, edgecolor='black',
                            frameon=True, framealpha=0.75)
            leg.get_frame().set_linewidth(0.5)

        fig.tight_layout()

        if self.save:
            plt.savefig('./Figures/PDF/fig_hist_diff.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig('./Figures/PNG/fig_hist_diff.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print("Exported differential cumulative capture distribution plot: fig_hist_diff")
        else:
            print("Did not export differential CDF plot due to user setting.")

        if self.show:
            plt.show()
        plt.close('all')


    def plot_dRdn(self):
        print(f"\nPlotting rate of change of fissile material with respect to fertile density vs. fertile density...")

        tot_n_per_yr = NPS_FUS * 3.156e+7

        # Gotta then convert fissile atoms/yr to kg/yr
        pu239_at_to_kg_per_yr = tot_n_per_yr / AVO * AMU_PU239 / 1e3
        u233_at_to_kg_per_yr  = tot_n_per_yr / AVO * AMU_U233  / 1e3

        # Setting up x, y separately here so you can remove the impurity/wppm-magnitude cases --ppark 2025-08-06
        # pebble bed
        x6, y6 =     self.hcpb_u_rr_df['fertile_kg/m3'], pu239_at_to_kg_per_yr * self.hcpb_u_rr_df['U238(n,g)']
        x5, y5 =    self.hcpb_th_rr_df['fertile_kg/m3'], u233_at_to_kg_per_yr * self.hcpb_th_rr_df['Th232(n,g)']

        # lead lithium     
        x4, y4 =   self.dcll_u_rr_df['fertile_kg/m3'],  pu239_at_to_kg_per_yr *  self.dcll_u_rr_df['U238(n,g)']
        x3, y3 =  self.dcll_th_rr_df['fertile_kg/m3'],  u233_at_to_kg_per_yr * self.dcll_th_rr_df['Th232(n,g)']

        # flibe
        x2, y2 =  self.flibe_u_rr_df['fertile_kg/m3'],  pu239_at_to_kg_per_yr *  self.flibe_u_rr_df['U238(n,g)']
        x1, y1 = self.flibe_th_rr_df['fertile_kg/m3'],  u233_at_to_kg_per_yr * self.flibe_th_rr_df['Th232(n,g)']

        # find derivatives 
        quartic_deriv1 = fitquartic(x1[:17], y1[:17])[1] # fit quartic up until 250 kg enrichment
        quartic_deriv2 = fitquartic(x2[:17], y2[:17])[1]
        quartic_deriv3 = fitquartic(x3[:17], y3[:17])[1]
        quartic_deriv4 = fitquartic(x4[:17], y4[:17])[1]
        quartic_deriv5 = fitquartic(x5[:17], y5[:17])[1]
        quartic_deriv6 = fitquartic(x6[:17], y6[:17])[1]

        # evaluate
        x_fine = np.linspace(x1.min(), x1.max(), 300)
        y1eval, y1fine = quartic_deriv1(x1), quartic_deriv1(x_fine)
        y2eval, y2fine = quartic_deriv2(x2), quartic_deriv2(x_fine)
        y3eval, y3fine = quartic_deriv3(x3), quartic_deriv3(x_fine)
        y4eval, y4fine = quartic_deriv4(x4), quartic_deriv4(x_fine)
        y5eval, y5fine = quartic_deriv5(x5), quartic_deriv5(x_fine)
        y6eval, y6fine = quartic_deriv6(x6), quartic_deriv6(x_fine)

        # plot
        plt.figure(figsize=(7.5,5))

        # plt.scatter(x6, y6eval, marker='s', s=40, color='#b41f24')   # HCPB-UO2
        # plt.scatter(x5, y5eval, marker='1', s=70, color='#b41f24')   # HCPB-ThO2
        # plt.scatter(x4, y4eval, marker='^', s=40, color='#0047ba')   # DCLL-UO2
        # plt.scatter(x3, y3eval, marker='x', s=50, color='#0047ba')   # DCLL-ThO2
        # plt.scatter(x2, y2eval, marker='o', s=30, color='#66b420')   # FLiBe-UF4
        # plt.scatter(x1, y1eval, marker='+', s=60, color='#66b420')   # FLiBe-ThF4

        plt.plot(x_fine, y6fine, '-',   linewidth=1, color='#b41f24',)   #  HCPB-UO2
        plt.plot(x_fine, y5fine, '--',  linewidth=1, color='#b41f24',)   #  HCPB-ThO2
        plt.plot(x_fine, y4fine, '-',   linewidth=1, color='#0047ba',)   #  DCLL-UO2
        plt.plot(x_fine, y3fine, '--',  linewidth=1, color='#0047ba',)   #  DCLL-ThO2
        plt.plot(x_fine, y2fine, '-',   linewidth=1, color='#66b420', )  #  FLiBe-UF4
        plt.plot(x_fine, y1fine, '--',  linewidth=1, color='#66b420', )  #  FLiBe-ThF4 

        # Dummy plots for legend -- ppark
        # plt.plot([9e8,9e9], [9e8,9e9], '+--', markersize=8, linewidth=1, color='#66b420', label=r'FLiBe-ThF$_4$')   #  green: #66b420
        # plt.plot([9e8,9e9], [9e8,9e9], 'o-',  markersize=5, linewidth=1, color='#66b420', label=r'FLiBe-UF$_4$')    #  green: #66b420
        # plt.plot([9e8,9e9], [9e8,9e9], '1--', markersize=9, linewidth=1, color='#b41f24', label=r'HCPB-ThO$_2$')    #  red:   #b41f24
        # plt.plot([9e8,9e9], [9e8,9e9], 's-',  markersize=5, linewidth=1, color='#b41f24', label=r'HCPB-UO$_2$')     #  red:   #b41f24
        # plt.plot([9e8,9e9], [9e8,9e9], 'x--', markersize=7, linewidth=1, color='#0047ba', label=r'DCLL-ThO$_2$')    #  blue:  #0047ba
        # plt.plot([9e8,9e9], [9e8,9e9], '^-',  markersize=6, linewidth=1, color='#0047ba', label=r'DCLL-UO$_2$')     #  blue:  #0047ba

        # plt.title(f'Tritium breeding ratio (All)') # Exclude title for production figs --ppark 2025-08-06
        plt.xlabel(r'Fertile isotope density [kg$/$m³]')  # specifically use unicode superscript m³ and not m$^3$
        plt.ylabel(r'$dR_{\text{fis}}/dn$ [kg$/\text{yr}^2$]')

        plt.xlim(-25,1025)
        plt.ylim(-0.33,6.33)  # plt.ylim(-7.5,225+7.5) REMEMBER TO CHANGE 

        # Tick grid
        ax = plt.gca()
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
 
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(MultipleLocator(1))   # ax.yaxis.set_major_locator(MultipleLocator(25))
        ax.yaxis.set_minor_locator(MultipleLocator(0.25)) # ax.yaxis.set_minor_locator(MultipleLocator(12.5))
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

        plt.tight_layout()

        leg = plt.legend(fancybox=False, edgecolor='black', frameon=False, framealpha=.75, ncol=1)
        leg.get_frame().set_linewidth(0.5) 
    
        if self.save:
            plt.savefig(f'./Figures/PDF/fig_dRdn.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig(f'./Figures/PNG/fig_dRdn.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print("Exported rate of change for fissile material production with respect to fertile material plot for all blankets.")
        else:
            print("Did not export dR / dn plot due to user setting.")

        if self.show: plt.show()
        plt.close('all')


    def plot_Rovern(self):
        print(f"\nPlotting fissile material divided by fertile density vs. fertile density...")

        tot_n_per_yr = NPS_FUS * 3.156e+7

        # Gotta then convert fissile atoms/yr to kg/yr
        pu239_at_to_kg_per_yr = tot_n_per_yr / AVO * AMU_PU239 / 1e3
        u233_at_to_kg_per_yr  = tot_n_per_yr / AVO * AMU_U233  / 1e3

        # Setting up x, y separately here so you can remove the impurity/wppm-magnitude cases --ppark 2025-08-06
        # cut out the first point since dividing by zero.
        # helium cooled pebble bed
        x6, y6 =     self.hcpb_u_rr_df['fertile_kg/m3'][1:], pu239_at_to_kg_per_yr * self.hcpb_u_rr_df['U238(n,g)'][1:]
        x5, y5 =    self.hcpb_th_rr_df['fertile_kg/m3'][1:], u233_at_to_kg_per_yr * self.hcpb_th_rr_df['Th232(n,g)'][1:]

        # dual coolant lead lithium     
        x4, y4 =   self.dcll_u_rr_df['fertile_kg/m3'][1:],  pu239_at_to_kg_per_yr *  self.dcll_u_rr_df['U238(n,g)'][1:]
        x3, y3 =  self.dcll_th_rr_df['fertile_kg/m3'][1:],  u233_at_to_kg_per_yr * self.dcll_th_rr_df['Th232(n,g)'][1:]

        # flibe
        x2, y2 =  self.flibe_u_rr_df['fertile_kg/m3'][1:],  pu239_at_to_kg_per_yr *  self.flibe_u_rr_df['U238(n,g)'][1:]
        x1, y1 = self.flibe_th_rr_df['fertile_kg/m3'][1:],  u233_at_to_kg_per_yr * self.flibe_th_rr_df['Th232(n,g)'][1:]

        x_fine = np.linspace(x1.min(), x1.max(), 300) # Evaluate on a fine grid
        y_akima1 = Akima1DInterpolator(x1, y1/x1)(x_fine)
        y_akima2 = Akima1DInterpolator(x2, y2/x2)(x_fine)
        y_akima3 = Akima1DInterpolator(x3, y3/x3)(x_fine)
        y_akima4 = Akima1DInterpolator(x4, y4/x4)(x_fine)
        y_akima5 = Akima1DInterpolator(x5, y5/x5)(x_fine)
        y_akima6 = Akima1DInterpolator(x6, y6/x6)(x_fine)

        # plot
        plt.figure(figsize=(3.5, 3.0))

        # plt.scatter(x6, y6/x6, marker='s', s=40/3, color='#b41f24') # HCPB-UO2
        # plt.scatter(x5, y5/x5, marker='1', s=70/3, color='#b41f24') # HCPB-ThO2
        # plt.scatter(x4, y4/x4, marker='^', s=40/3, color='#0047ba') # DCLL-UO2
        # plt.scatter(x3, y3/x3, marker='x', s=50/3, color='#0047ba') # DCLL-ThO2
        # plt.scatter(x2, y2/x2, marker='o', s=30/3, color='#66b420')   # FLiBe-UF4
        # plt.scatter(x1, y1/x1, marker='+', s=60/3, color='#66b420')   # FLiBe-ThF4

        plt.plot(x_fine, y_akima2, linestyle='-',       linewidth=0.75, color='#66b420',label=r'FLiBe-UF$_4$' )    #  FLiBe-UF4
        plt.plot(x_fine, y_akima1, linestyle=LONG_DASH, linewidth=0.75, color='#66b420',label=r'FLiBe-ThF$_4$')    #  FLiBe-ThF4 
        plt.plot(x_fine, y_akima6, linestyle='-',       linewidth=0.75, color='#b41f24',label=r'HCPB-UO$_2$'  )   #  HCPB-UO2
        plt.plot(x_fine, y_akima5, linestyle=LONG_DASH, linewidth=0.75, color='#b41f24',label=r'HCPB-ThO$_2$' )   #  HCPB-ThO2       
        plt.plot(x_fine, y_akima3, linestyle=LONG_DASH, linewidth=0.75, color='#0047ba',label=r'DCLL-ThO$_2$' )   #  DCLL-ThO2
        plt.plot(x_fine, y_akima4, linestyle='-',       linewidth=0.75, color='#0047ba',label=r'DCLL-UO$_2$'  )   #  DCLL-UO2


        # Dummy plots for legend -- ppark
        # plt.plot([9e8,9e9], [9e8,9e9], '+--', markersize=8, linewidth=1, color='#66b420', label=r'FLiBe-ThF$_4$') #  green: #66b420
        # plt.plot([9e8,9e9], [9e8,9e9], 'o-',  markersize=5, linewidth=1, color='#66b420', label=r'FLiBe-UF$_4$')  # 
        # plt.plot([9e8,9e9], [9e8,9e9], '1--', markersize=9, linewidth=1, color='#b41f24', label=r'HCPB-ThO$_2$')    # 
        # plt.plot([9e8,9e9], [9e8,9e9], 's-',  markersize=5, linewidth=1, color='#b41f24', label=r'HCPB-UO$_2$')     # 
        # plt.plot([9e8,9e9], [9e8,9e9], 'x--', markersize=7, linewidth=1, color='#0047ba', label=r'DCLL-ThO$_2$')    #  red: #b41f24
        # plt.plot([9e8,9e9], [9e8,9e9], '^-',  markersize=6, linewidth=1, color='#0047ba', label=r'DCLL-UO$_2$')     #  blue: #0047ba

        # plt.title(f'Tritium breeding ratio (All)') # Exclude title for production figs --ppark 2025-08-06
        plt.xlabel(r'Fertile isotope density [kg$/$m³]')  # specifically use unicode superscript m³ and not m$^3$
        plt.ylabel(r'Fissile prod rate [kg$/$yr] per fertile density')

        plt.xlim(-25,1025)
        plt.ylim(-0.25,8.25) 

        # Tick grid
        ax = plt.gca()
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_major_locator(MultipleLocator(100)) 
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
 
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(MultipleLocator(1)) 
        ax.yaxis.set_minor_locator(MultipleLocator(0.5)) 
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

        plt.tight_layout()

        leg = plt.legend(fancybox=False, edgecolor='black', frameon=False, framealpha=.75, ncol=1)
        leg.get_frame().set_linewidth(0.5) 
    
        if self.save:
            plt.savefig(f'./Figures/PDF/fig_Rn.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig(f'./Figures/PNG/fig_Rn.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print("Exported fissile material production divided by fertile material for all blankets.")
        else:
            print("Did not export R/n plot due to user setting.")

        if self.show: 
            plt.show()
        plt.close('all')


    def plot_leakage_spectra(self):
        """
        Neutron leakage energy spectrum for four FLiBe configurations.
        Prints mean, median, and mode energy for each case.
        """
        print(f"\nComment. <plot.py/plot_leakage_spectra()> Plotting leakage energy spectra...")

        cases = [
            {"file": "./Figures/Data/FLiBe_900K_Li07.5_U238_leak.csv",  "rho": 0,    "color": "black",      "label": r"7.5at% Li6 / Clean"},
            {"file": "./Figures/Data/FLiBe_900K_Li00.0_U238_leak.csv",  "rho": 0,    "color": "#0047ba",  "label": r"Li7 / Clean"},
            {"file": "./Figures/Data/FLiBe_900K_Li00.0_U238_leak.csv",  "rho": 1000, "color": "#b41f24",  "label": r"Li7 / 1000 kg(U238)/m³"},
            {"file": "./Figures/Data/FLiBe_900K_Li00.0_Th232_leak.csv", "rho": 1000, "color": "#66b420",  "label": r"Li7 / 1000 kg(Th232)/m³"},
        ]

        plt.figure(figsize=(7, 3.0))
        ax = plt.gca()

        print(f"\n  {'Case':<45s} {'Mean [eV]':>14s} {'Median [eV]':>14s} {'Mode [eV]':>14s}")
        print(f"  {'─'*45} {'─'*14} {'─'*14} {'─'*14}")

        for case in cases:
            try:
                df = pd.read_csv(case["file"])
            except FileNotFoundError:
                print(f"{C.YELLOW}Warning.{C.END} Leakage file not found: {case['file']}")
                continue

            # Match target density (handles 999.99 ≈ 1000)
            available = df["fertile_kg/m3"].unique()
            closest = min(available, key=lambda x: abs(x - case["rho"]))
            spec = df[df["fertile_kg/m3"] == closest].sort_values("energy mid [eV]")

            E_mid  = spec["energy mid [eV]"].values
            E_low  = spec["energy low [eV]"].values
            E_high = spec["energy high [eV]"].values
            dE     = E_high - E_low
            counts = spec["mean"].values                # leakage current per bin [n/src-n]
            flux_per_lethargy = counts / (dE / E_mid)

            # ── Spectral statistics ──────────────────────────
            total = counts.sum()
            if total > 0:
                mean_E   = np.sum(E_mid * counts) / total
                cdf      = np.cumsum(counts) / total
                median_E = np.interp(0.5, cdf, E_mid)
                mode_E   = E_mid[np.argmax(flux_per_lethargy)]   # peak of per-lethargy spectrum
            else:
                mean_E = median_E = mode_E = 0.0

            short_label = case["label"].replace(r"$^6$", "⁶")   # for terminal printing
            print(f"  {short_label:<45s} {mean_E:14.4e} {median_E:14.4e} {mode_E:14.4e}")

            ax.loglog(E_mid, flux_per_lethargy,
                      color=case["color"], linewidth=0.75, label=case["label"])

        print()

        ax.set_xlabel("Energy [eV]")
        ax.set_ylabel("Leakage current / lethargy [n/src-n]")
        ax.set_xlim(1e-3, 20e6)

        # "1eX" tick labels on both axes, major ticks every decade, minor ticks at 2–9
        from matplotlib.ticker import LogLocator, LogFormatterSciNotation
        efmt = LogFormatterSciNotation(base=10, labelOnlyBase=True)
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_locator(LogLocator(base=10, numticks=20))
            axis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
            axis.set_major_formatter(FuncFormatter(lambda x, _: f'1e{int(np.log10(x))}' if x > 0 else ''))

        ax.grid(axis='both', which='major', linestyle='-', linewidth=0.5)
        ax.tick_params(which='minor', length=3)
        ax.tick_params(which='major', length=5)

        plt.tight_layout()

        leg = plt.legend(fontsize=6, fancybox=False, edgecolor='black', frameon=False, framealpha=0.75, ncol=1)
        leg.get_frame().set_linewidth(0.5)

        if self.save:
            plt.savefig('./Figures/PDF/fig_leakage_spectra.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig('./Figures/PNG/fig_leakage_spectra.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print(f"Comment. <plot.py/plot_leakage_spectra()> Exported leakage spectra plot.")
        else:
            print(f"{C.YELLOW}Comment.{C.END} <plot.py/plot_leakage_spectra()> Did NOT export leakage spectra plot due to user setting.")

        if self.show: plt.show()
        plt.close('all')


    def plot_neutron_balance(self):
        """
        Cumulative neutron balance vs. energy.
 
        3x2 grid: rows = FLiBe / HCPB / DCLL, columns = U / Th.
        Each panel shows cumulative neutron gain/loss curves for four
        reaction channels at two fertile loadings (0.1 and 1000 kg/m3).
 
        Channels and neutron weights:
          Be9 or Pb (n,2n)                         -> +1
          U235 + fertile isotope (n,2n)             -> +1
          U235 + fertile isotope fis               -> +(nu_bar - 1)   [approximate]
          Fertile isotope (n,gamma)                -> -1
          Li6 (n,Xt)                               -> -1
        """
        print(f"\nPlotting cumulative neutron balance vs energy...")
 
        # Approximate nu-bar values
        # DONE // TODO: replace with nu-fission tally for energy-dependent accuracy
        # NU_U235  = 2.5
        # NU_U238  = 2.8
        # NU_TH232 = 2.3
 
        # Panel layout: (dataframe, title, fertile_iso, row, col)
        panels = [
            (self.flibe_u_ebins_df,  r'FLiBe-UF$_4$',  'U238',  0, 0),
            (self.flibe_th_ebins_df, r'FLiBe-ThF$_4$', 'Th232', 0, 1),
            (self.hcpb_u_ebins_df,   r'HCPB-UO$_2$',   'U238',  1, 0),
            (self.hcpb_th_ebins_df,  r'HCPB-ThO$_2$',  'Th232', 1, 1),
            (self.dcll_u_ebins_df,   r'DCLL-UO$_2$',   'U238',  2, 0),
            (self.dcll_th_ebins_df,  r'DCLL-ThO$_2$',  'Th232', 2, 1),
        ]
 
        densities_to_plot = [0.1, 999.99]
        linestyles = {0.1: '-', 999.99: '--'}
        linewidths = {0.1: 0.75, 999.99: 1.0}
 
        # Channel colors
        colors = {
            'n2n':      '#1D9E75',   # teal   -- neutron multiplier (Be/Pb)
            'n2n_act':  '#7B2D8E',   # purple -- neutron multiplier (actinide)
            'fis':      '#b41f24',   # red    -- fission
            'ng':       '#0047ba',   # blue   -- fertile capture
            'li6':      '#EF9F27',   # amber  -- tritium breeding
        }

        from matplotlib.lines import Line2D
        legend_handles = [
            Line2D([0], [0], color=colors['n2n'],     lw=1.5, label=r'Be$/$Pb $($n,2n$)$: +1'),
            # Line2D([0], [0], color=colors['n2n_act'], lw=1.5, label=r'U235$/$fert $($n,2n$)$: +1'),
            Line2D([0], [0], color=colors['fis'], lw=1.5, label=r'Fission: +$(\bar\nu - 1)$'),
            Line2D([0], [0], color='black', ls='-',  lw=1.0, label=r'0.1 kg$/$m³'),
            Line2D([0], [0], color=colors['ng'],  lw=1.5, label=r'Fertile $($n,$\gamma)$: $-$1'),
            Line2D([0], [0], color=colors['li6'], lw=1.5, label=r'Li-6 $($n,t$)$: $-$1'),
            Line2D([0], [0], color='black', ls='--', lw=1.0,  label=r'1000 kg$/$m³'),
        ]
 
        fig, axes = plt.subplots(3, 2, figsize=(7, 9), sharex=True, sharey=True)
 
        for df_Erxn, title, fert_iso, row, col in panels:
            ax = axes[row, col]
 
            for rho in densities_to_plot:
 
                # Match closest available density
                available = df_Erxn['fertile_kg/m3'].unique()
                closest = min(available, key=lambda x: abs(x - rho))
                spec = df_Erxn[df_Erxn['fertile_kg/m3'] == closest].copy()
                spec = spec.sort_values('energy mid [eV]')
 
                E_mid = spec['energy mid [eV]'].values
                ls = linestyles[rho]
                lw = linewidths[rho]
 
                # Neutron-weighted channels

                # Multiplier (n,2n): +1 per reaction
                n2n = (+1) * (spec['Be9_n2n'].values + spec['Pb_n2n'].values)

                # Actinide (n,2n): +1 per reaction
                n2n_act = (+1) * (spec['U235_n2n'].values + spec[f'{fert_iso}_n2n'].values)
 
                # Fission: +(nu-1) for fertile + +(nu_U235 - 1) for U235
                fis = (spec[f'{fert_iso}_nufis'].values + spec['U235_nufis'].values)
 
                # Fertile (n,gamma): -1 per reaction
                ng = (-1) * spec[f'{fert_iso}_ng'].values
 
                # Li-6 (n,t): -1 per reaction
                li6 = (-1) * spec['Li6_nt'].values
 
                # Cumulative sum from low to high energy
                for key, vals in [('n2n', n2n), ('fis', fis), ('ng', ng), ('li6', li6)]: # ('n2n_act', n2n_act), 
                    ax.semilogx(E_mid, np.cumsum(vals),
                                color=colors[key], linestyle=ls, linewidth=lw)
 
            # Axis formatting
            # ax.axhline(0, color='gray', linewidth=0.5, zorder=0)
 
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.grid(axis='both', which='major', linestyle='-', linewidth=0.5)
            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'1e{int(np.log10(x))}' if x > 0 else ''))
            ax.set_xlim(log_buffer(1e1,1e7))
            ax.set_yscale('linear')
            ax.set_ylim(-1.5*1.03, 0.75*1.03)
            ax.yaxis.set_major_locator(MultipleLocator(0.25))
 
            if ax.get_subplotspec().is_last_row():
                ax.set_xlabel('Neutron energy [eV]')
            if ax.get_subplotspec().is_first_col():
                ax.set_ylabel(r'Cumulative $\Delta n$ / src-n')
 
            # Panel title in upper-left corner
            ax.text(0.02, 0.97, title, transform=ax.transAxes,
                    va='top', ha='left') # fontsize=10, fontweight='bold', -- use default font size --ppark 2026-05-21

            # Legend keys right beneath the title
            if row in [0, 1, 2]:
                leg = ax.legend(handles=legend_handles, loc='upper left',
                                bbox_to_anchor=(0.00, 0.935),
                                fontsize=7, ncol=2,
                                fancybox=False, edgecolor='black',
                                frameon=False, framealpha=0.75)
                # leg.get_frame().set_linewidth(0.5)
 
        fig.tight_layout()
 
        if self.save:
            plt.savefig('./Figures/PDF/fig_nbal.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig('./Figures/PNG/fig_nbal.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print("Exported cumulative neutron balance plot: fig_nbal")
        else:
            print("Did not export neutron balance plot due to user setting.")
 
        if self.show:
            plt.show()
        plt.close('all')


    def plot_rxn_spectra(self):
        """
        Per-lethargy reaction rate spectra vs. energy.

        3x2 grid (rows = FLiBe / HCPB / DCLL, cols = U238 / Th232).
        Each panel shows reaction rate per unit lethargy [rxns / src-n]
        for seven channels at two fertile loadings (0.1 and 1000 kg/m3):

          Be/Pb (n,2n)          -- multiplier material
          Be/Pb (n,elastic)     -- multiplier elastic scattering
          Actinide (n,2n)       -- fertile + U235
          Fertile (n,inelastic) -- fertile MT 4
          Fission               -- total (fertile + U235), raw rate
          Fertile (n,gamma)     -- capture / fissile breeding
          Li-6 (n,t)            -- tritium breeding

        No neutron-weight coefficients are applied.  Per-lethargy
        normalization (rate / d_lnE) puts all energy decades on equal
        visual footing, making resonance peaks and threshold reactions
        directly comparable on a log-x axis.
        """
        print(f"\nPlotting per-lethargy reaction rate spectra vs energy...")

        panels = [
            (self.flibe_u_ebins_df,  r'FLiBe-UF$_4$',  'U238',  0, 0),
            (self.flibe_th_ebins_df, r'FLiBe-ThF$_4$', 'Th232', 0, 1),
            (self.hcpb_u_ebins_df,   r'HCPB-UO$_2$',   'U238',  1, 0),
            (self.hcpb_th_ebins_df,  r'HCPB-ThO$_2$',  'Th232', 1, 1),
            (self.dcll_u_ebins_df,   r'DCLL-UO$_2$',   'U238',  2, 0),
            (self.dcll_th_ebins_df,  r'DCLL-ThO$_2$',  'Th232', 2, 1),
        ]

        densities_to_plot = [0.1, 999.99]
        linestyles = {0.1: '-', 999.99: '--'}
        linewidths = {0.1: 0.75, 999.99: 1.0}

        # Channel colors
        colors = {
            'mult_n2n':  '#1D9E75',   # teal        -- Be/Pb (n,2n)
            'mult_elas': '#85C77E',   # light green -- Be/Pb elastic
            'fert_n2n':  '#7B2D8E',   # purple      -- actinide (n,2n)
            'fert_inel': '#C77DBA',   # lilac       -- fertile inelastic
            'fis':       '#b41f24',   # red         -- fission
            'ng':        '#0047ba',   # blue        -- fertile capture
            'li':        '#EF9F27',   # amber       -- Li-6 tritium breeding
        }

        from matplotlib.lines import Line2D
        legend_handles = [
            Line2D([0], [0], color=colors['mult_n2n'],  lw=1.5, label=r'Be$/$Pb $($n,2n$)$'),
            Line2D([0], [0], color=colors['mult_elas'], lw=1.5, label=r'Be$/$Pb elastic'),
            Line2D([0], [0], color=colors['fert_n2n'],  lw=1.5, label=r'Actinide $($n,2n$)$'),
            Line2D([0], [0], color=colors['fert_inel'], lw=1.5, label=r'Fertile inelastic'),
            Line2D([0], [0], color='black', ls='-',     lw=1.0, label=r'0.1 kg$/$m³'),
            Line2D([0], [0], color=colors['fis'],       lw=1.5, label=r'Fission'),
            Line2D([0], [0], color=colors['ng'],        lw=1.5, label=r'Fertile $($n,$\gamma)$'),
            Line2D([0], [0], color=colors['li'],        lw=1.5, label=r'Li $($n,t$)$'),
            Line2D([0], [0], color='black', ls='--',    lw=1.0,  label=r'1000 kg$/$m³'),
        ]

        fig, axes = plt.subplots(3, 2, figsize=(7, 9), sharex=True, sharey=True)

        for df_Erxn, title, fert_iso, row, col in panels:
            ax = axes[row, col]

            for rho in densities_to_plot:

                # Match closest available density
                available = df_Erxn['fertile_kg/m3'].unique()
                closest = min(available, key=lambda x: abs(x - rho))
                spec = df_Erxn[df_Erxn['fertile_kg/m3'] == closest].copy()
                spec = spec.sort_values('energy mid [eV]')

                E_mid  = spec['energy mid [eV]'].values
                E_low  = spec['energy low [eV]'].values
                E_high = spec['energy high [eV]'].values
                d_lnE  = (E_high - E_low) / E_mid   # lethargy bin width
                ls = linestyles[rho]
                lw = linewidths[rho]

                # Raw reaction-rate channels (no neutron weights)
                channels = {
                    'mult_n2n':  spec['Be9_n2n'].values + spec['Pb_n2n'].values,
                    'mult_elas': spec['Be9_elas'].values + spec['Pb_elas'].values,
                    'fert_n2n':  spec[f'{fert_iso}_n2n'].values  + spec['U235_n2n'].values,
                    'fert_inel': spec[f'{fert_iso}_inel'].values + spec['U235_inel'].values,
                    'fis':       spec[f'{fert_iso}_fis'].values  + spec['U235_fis'].values,
                    'ng':        spec[f'{fert_iso}_ng'].values,
                    'li':        spec['Li6_nt'].values + spec['Li7_nt'].values,
                }

                for key, vals in channels.items():
                    per_lethargy = vals / d_lnE
                    # Mask zeros to avoid log-scale artifacts
                    per_lethargy[per_lethargy <= 0] = np.nan
                    ax.loglog(E_mid, per_lethargy,
                              color=colors[key], linestyle=ls, linewidth=lw)

            # Axis formatting
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.grid(axis='both', which='major', linestyle='-', linewidth=0.5)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'1e{int(np.log10(x))}' if x > 0 else ''))
            ax.set_xlim(log_buffer(1e1,1e7))
            # ax.set_ylim(0.5e-11, 1.5e1)

            from matplotlib.ticker import LogLocator
            ax.yaxis.set_major_locator(LogLocator(base=10, numticks=20))
            ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'1e{int(np.log10(x))}' if x > 0 else ''))

            if ax.get_subplotspec().is_last_row():
                ax.set_xlabel('Neutron energy [eV]')
            if ax.get_subplotspec().is_first_col():
                ax.set_ylabel(r'Reaction rate / lethargy [rxns$/$src-n]')

            # Panel title
            ax.text(0.02, 0.97, title, transform=ax.transAxes,
                    va='top', ha='left')

            # Legend
            leg = ax.legend(handles=legend_handles, loc='upper left',
                            bbox_to_anchor=(0.00, 0.935),
                            fontsize=6, ncol=3,
                            fancybox=False, edgecolor='black',
                            frameon=False, framealpha=0.75)

        fig.tight_layout()

        if self.save:
            plt.savefig('./Figures/PDF/fig_rxn_spectra.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig('./Figures/PNG/fig_rxn_spectra.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print("Exported per-lethargy reaction rate spectra: fig_rxn_spectra")
        else:
            print("Did not export reaction rate spectra plot due to user setting.")

        if self.show:
            plt.show()
        plt.close('all')


    def plot_rxn_ratio(self):
        """
        Per-lethargy reaction-rate ratio: R(1000) / R(0.1) vs. energy.

        3x2 grid (rows = FLiBe / HCPB / DCLL, cols = U238 / Th232).
        For each of the seven reaction channels used in plot_rxn_spectra,
        the ratio of the per-lethargy reaction rate at 1000 kg/m³ to that
        at 0.1 kg/m³ is plotted.  A ratio > 1 means the channel increased
        with fertile loading; < 1 means it decreased.

        Uses the same channel color scheme as plot_rxn_spectra.
        """
        print(f"\nPlotting per-lethargy reaction rate ratios (1000 / 0.1 kg/m³) vs energy...")

        panels = [
            (self.flibe_u_ebins_df,  r'FLiBe-UF$_4$',  'U238',  0, 0),
            (self.flibe_th_ebins_df, r'FLiBe-ThF$_4$', 'Th232', 0, 1),
            (self.hcpb_u_ebins_df,   r'HCPB-UO$_2$',   'U238',  1, 0),
            (self.hcpb_th_ebins_df,  r'HCPB-ThO$_2$',  'Th232', 1, 1),
            (self.dcll_u_ebins_df,   r'DCLL-UO$_2$',   'U238',  2, 0),
            (self.dcll_th_ebins_df,  r'DCLL-ThO$_2$',  'Th232', 2, 1),
        ]

        rho_num = 999.99   # numerator density   (1000 kg/m³)
        rho_den = 0.1      # denominator density  (0.1  kg/m³)

        # Channel colors — identical to plot_rxn_spectra
        colors = {
            'mult_n2n':  '#1D9E75',   # teal        -- Be/Pb (n,2n)
            'mult_elas': '#85C77E',   # light green -- Be/Pb elastic
            'fert_n2n':  '#7B2D8E',   # purple      -- actinide (n,2n)
            'fert_inel': '#C77DBA',   # lilac       -- fertile inelastic
            'fis':       '#b41f24',   # red         -- fission
            'ng':        '#0047ba',   # blue        -- fertile capture
            'li':        '#EF9F27',   # amber       -- Li-6 tritium breeding
        }

        # Per-row y-axis limits and tick spacing
        row_ylim  = {0: (0.0, 1.5e4),  1: (0.0, 1.5e4),  2: (0.0, 4e4)}
        # row_major = {0: 0.25,        1: 0.25,         2: 0.5}
        # row_minor = {0: 0.125,       1: 0.125,        2: 0.25}

        from matplotlib.lines import Line2D
        legend_handles = [ Line2D([0], [0], color=colors['mult_n2n'],  lw=1.5, label=r'Be$/$Pb $($n,2n$)$'),
                           Line2D([0], [0], color=colors['mult_elas'], lw=1.5, label=r'Be$/$Pb elastic'),
                           Line2D([0], [0], color=colors['fert_n2n'],  lw=1.5, label=r'Actinide $($n,2n$)$'),
                           Line2D([0], [0], color=colors['fert_inel'], lw=1.5, label=r'Fertile inelastic'),
                           Line2D([0], [0], color=colors['fis'],       lw=1.5, label=r'Fission'),
                           Line2D([0], [0], color=colors['ng'],        lw=1.5, label=r'Fertile $($n,$\gamma)$'),
                           Line2D([0], [0], color=colors['li'],        lw=1.5, label=r'Li $($n,t$)$'),           ]

        fig, axes = plt.subplots(3, 2, figsize=(7, 9), sharex=True, sharey=False)

        for df_Erxn, title, fert_iso, row, col in panels:
            ax = axes[row, col]

            available = df_Erxn['fertile_kg/m3'].unique()

            # ---- denominator (0.1 kg/m³) ----
            den_rho = min(available, key=lambda x: abs(x - rho_den))
            den = df_Erxn[df_Erxn['fertile_kg/m3'] == den_rho].copy().sort_values('energy mid [eV]')

            E_mid_den  = den['energy mid [eV]'].values
            E_low_den  = den['energy low [eV]'].values
            E_high_den = den['energy high [eV]'].values
            d_lnE_den  = (E_high_den - E_low_den) / E_mid_den

            den_channels = {
                'mult_n2n':  den['Be9_n2n'].values  + den['Pb_n2n'].values,
                'mult_elas': den['Be9_elas'].values  + den['Pb_elas'].values,
                'fert_n2n':  den[f'{fert_iso}_n2n'].values  + den['U235_n2n'].values,
                'fert_inel': den[f'{fert_iso}_inel'].values + den['U235_inel'].values,
                'fis':       den[f'{fert_iso}_fis'].values  + den['U235_fis'].values,
                'ng':        den[f'{fert_iso}_ng'].values,
                'li':        den['Li6_nt'].values + den['Li7_nt'].values,
            }

            # ---- numerator (1000 kg/m³) ----
            num_rho = min(available, key=lambda x: abs(x - rho_num))
            num = df_Erxn[df_Erxn['fertile_kg/m3'] == num_rho].copy().sort_values('energy mid [eV]')

            E_mid_num  = num['energy mid [eV]'].values
            E_low_num  = num['energy low [eV]'].values
            E_high_num = num['energy high [eV]'].values
            d_lnE_num  = (E_high_num - E_low_num) / E_mid_num

            num_channels = {
                'mult_n2n':  num['Be9_n2n'].values  + num['Pb_n2n'].values,
                'mult_elas': num['Be9_elas'].values  + num['Pb_elas'].values,
                'fert_n2n':  num[f'{fert_iso}_n2n'].values  + num['U235_n2n'].values,
                'fert_inel': num[f'{fert_iso}_inel'].values + num['U235_inel'].values,
                'fis':       num[f'{fert_iso}_fis'].values  + num['U235_fis'].values,
                'ng':        num[f'{fert_iso}_ng'].values,
                'li':        num['Li6_nt'].values + num['Li7_nt'].values,
            }

            for key in colors:
                num_pl = num_channels[key] / d_lnE_num
                den_pl = den_channels[key] / d_lnE_den

                ratio = np.ones_like(den_pl, dtype=float)
                valid = (den_pl > 0) & (num_pl > 0)
                ratio[valid]  = num_pl[valid] / den_pl[valid]
                ratio[~valid] = np.nan

                ax.semilogx(E_mid_den, ratio, color=colors[key], linewidth=0.75)

            # Reference line at ratio = 1
            ax.axhline(1.0, color='gray', linewidth=0.5, linestyle='-', zorder=0)

            # Axis formatting
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.grid(axis='both', which='major', linestyle='-', linewidth=0.5)
            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'1e{int(np.log10(x))}' if x > 0 else ''))
            ax.set_xlim(log_buffer(1e1, 1e7))
            ax.set_yscale('log')
            ax.set_ylim(row_ylim[row])
            # ax.yaxis.set_major_locator(MultipleLocator(row_major[row]))
            # ax.yaxis.set_minor_locator(MultipleLocator(row_minor[row]))
            # ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

            if ax.get_subplotspec().is_last_row():
                ax.set_xlabel('Neutron energy [eV]')
            if ax.get_subplotspec().is_first_col():
                ax.set_ylabel(r'Reaction rate ratio $R(1000)\,/\,R(0.1)$')

            # Panel title
            ax.text(0.02, 0.97, title, transform=ax.transAxes,
                    va='top', ha='left')

            # Legend
            leg = ax.legend(handles=legend_handles, loc='upper left',
                            bbox_to_anchor=(0.00, 0.935),
                            fontsize=6, ncol=3,
                            fancybox=False, edgecolor='black',
                            frameon=False, framealpha=0.75)

        fig.tight_layout()

        if self.save:
            plt.savefig('./Figures/PDF/fig_rxn_ratio.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig('./Figures/PNG/fig_rxn_ratio.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print("Exported per-lethargy reaction rate ratio plot: fig_rxn_ratio")
        else:
            print("Did not export reaction rate ratio plot due to user setting.")

        if self.show:
            plt.show()
        plt.close('all')

    
    def plot_rxn_cumulative(self):
        """
        Cumulative reaction rate vs. energy.

        3x2 grid (rows = FLiBe / HCPB / DCLL, cols = U238 / Th232).
        Each panel shows the running integral of the reaction rate
        [rxns / src-n] from low energy upward for seven channels at
        two fertile loadings (0.1 and 1000 kg/m³).

        Uses the same channel color scheme as plot_rxn_spectra.
        """
        print(f"\nPlotting cumulative reaction rate vs energy...")

        panels = [
            (self.flibe_u_ebins_df,  r'FLiBe-UF$_4$',  'U238',  0, 0),
            (self.flibe_th_ebins_df, r'FLiBe-ThF$_4$', 'Th232', 0, 1),
            (self.hcpb_u_ebins_df,   r'HCPB-UO$_2$',   'U238',  1, 0),
            (self.hcpb_th_ebins_df,  r'HCPB-ThO$_2$',  'Th232', 1, 1),
            (self.dcll_u_ebins_df,   r'DCLL-UO$_2$',   'U238',  2, 0),
            (self.dcll_th_ebins_df,  r'DCLL-ThO$_2$',  'Th232', 2, 1),
        ]

        densities_to_plot = [0.1, 999.99]
        linestyles = {0.1: '-', 999.99: '--'}
        linewidths = {0.1: 0.75, 999.99: 1.0}

        # Channel colors — identical to plot_rxn_spectra
        colors = {
            'mult_n2n':  '#1D9E75',   # teal        -- Be/Pb (n,2n)
            'mult_elas': '#85C77E',   # light green -- Be/Pb elastic
            'fert_n2n':  '#7B2D8E',   # purple      -- actinide (n,2n)
            'fert_inel': '#C77DBA',   # lilac       -- fertile inelastic
            'fis':       '#b41f24',   # red         -- fission
            'ng':        '#0047ba',   # blue        -- fertile capture
            'li':        '#EF9F27',   # amber       -- Li-6 tritium breeding
        }

        from matplotlib.lines import Line2D
        legend_handles = [
            Line2D([0], [0], color=colors['mult_n2n'],  lw=1.5, label=r'Be$/$Pb $($n,2n$)$'),
            # Line2D([0], [0], color=colors['mult_elas'], lw=1.5, label=r'Be$/$Pb elastic'),
            Line2D([0], [0], color=colors['fert_n2n'],  lw=1.5, label=r'Actinide $($n,2n$)$'),
            Line2D([0], [0], color=colors['fert_inel'], lw=1.5, label=r'Fertile inelastic'),
            Line2D([0], [0], color='black', ls='-',     lw=1.0, label=r'0.1 kg$/$m³'),
            Line2D([0], [0], color=colors['fis'],       lw=1.5, label=r'Fission'),
            Line2D([0], [0], color=colors['ng'],        lw=1.5, label=r'Fertile $($n,$\gamma)$'),
            Line2D([0], [0], color=colors['li'],        lw=1.5, label=r'Li $($n,t$)$'),
            Line2D([0], [0], color='black', ls='--',    lw=1.0, label=r'1000 kg$/$m³'),
        ]

        fig, axes = plt.subplots(3, 2, figsize=(7, 9), sharex=True, sharey=True)

        for df_Erxn, title, fert_iso, row, col in panels:
            ax = axes[row, col]

            for rho in densities_to_plot:

                # Match closest available density
                available = df_Erxn['fertile_kg/m3'].unique()
                closest = min(available, key=lambda x: abs(x - rho))
                spec = df_Erxn[df_Erxn['fertile_kg/m3'] == closest].copy()
                spec = spec.sort_values('energy mid [eV]')

                E_mid = spec['energy mid [eV]'].values
                ls = linestyles[rho]
                lw = linewidths[rho]

                # Raw reaction-rate channels (no neutron weights, no per-lethargy)
                channels = {
                    'mult_n2n':  spec['Be9_n2n'].values  + spec['Pb_n2n'].values,
                    # 'mult_elas': spec['Be9_elas'].values + spec['Pb_elas'].values,
                    'fert_n2n':  spec[f'{fert_iso}_n2n'].values  + spec['U235_n2n'].values,
                    'fert_inel': spec[f'{fert_iso}_inel'].values + spec['U235_inel'].values,
                    'fis':       spec[f'{fert_iso}_fis'].values  + spec['U235_fis'].values,
                    'ng':        spec[f'{fert_iso}_ng'].values,
                    'li':        spec['Li6_nt'].values + spec['Li7_nt'].values,
                }

                for key, vals in channels.items():
                    cumulative = np.cumsum(vals)
                    # Mask leading zeros so the line starts where reactions begin
                    cumulative_masked = np.where(cumulative > 0, cumulative, np.nan)
                    ax.semilogx(E_mid, cumulative_masked,
                                color=colors[key], linestyle=ls, linewidth=lw)

            # Axis formatting
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.grid(axis='both', which='major', linestyle='-', linewidth=0.5)
            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'1e{int(np.log10(x))}' if x > 0 else ''))
            ax.set_xlim(log_buffer(1e1, 1e7))

            if ax.get_subplotspec().is_last_row():
                ax.set_xlabel('Neutron energy [eV]')
            if ax.get_subplotspec().is_first_col():
                ax.set_ylabel(r'Cumulative reaction rate [rxns$/$src-n]')

            # Panel title
            ax.text(0.02, 0.97, title, transform=ax.transAxes,
                    va='top', ha='left')

            # Legend
            leg = ax.legend(handles=legend_handles, loc='upper left',
                            bbox_to_anchor=(0.00, 0.935),
                            fontsize=6, ncol=3,
                            fancybox=False, edgecolor='black',
                            frameon=False, framealpha=0.75)

        fig.tight_layout()

        if self.save:
            plt.savefig('./Figures/PDF/fig_rxn_cumulative.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig('./Figures/PNG/fig_rxn_cumulative.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print("Exported cumulative reaction rate plot: fig_rxn_cumulative")
        else:
            print("Did not export cumulative reaction rate plot due to user setting.")

        if self.show:
            plt.show()
        plt.close('all')




if __name__ == "__main__":

    for sd in ['PNG','PDF','Data']:
        sd_path = f'./Figures/{sd}/'
        print(f"Comment. <plot.py/plot_all()> Ensuring directory exists: {sd_path}")
        os.makedirs(sd_path, exist_ok=True)


    plot_types = ['tbr','fpr','hist','histd','fluxr','fluxs','rxns','rxnr','dfis','fisn']

    parser = argparse.ArgumentParser(description=f"Choose plots with -p flag, multiple separated by spaces: {plot_types}")

    parser.add_argument("-p", "--plot_type", 
                        type=str, nargs="+", default=plot_types, 
                        help=f"Specify plot, multiple separated by spaces: {plot_types}. Defaults to all plots." )
      
    parser.add_argument("--no_show", 
                        dest="show", action="store_false",
                        help="Disable showing the plots.")
    
    parser.add_argument("--no_save", 
                        dest="save", action="store_false",
                        help="Disable saving the plots. ")

    parser.set_defaults(show=True, save=True)

    plot_type  = [p.lower() if isinstance(p, str) else p for p in parser.parse_args().plot_type]  # parser.parse_args().plot_type 
    plot_show  = parser.parse_args().show   
    plot_save  = parser.parse_args().save

    combined_plot = Plot(show=plot_show, save=plot_save)

    for p in plot_type:
        if   p == 'tbr':
            combined_plot.plot_tbr()
        elif p == 'fpr':
            combined_plot.plot_fpr()
        elif p == 'hist':
            combined_plot.plot_hist()
        elif p == 'histd':
            combined_plot.plot_hist_diff()
        elif p == 'fluxr':
            combined_plot.plot_flux_ratio()
        elif p == 'fluxs':
            combined_plot.plot_flux_spectrum()
        elif p == 'drdn':
            combined_plot.plot_dRdn()
        elif p == 'rn':
            combined_plot.plot_Rovern()
        elif p == 'leak':
            combined_plot.plot_leakage_spectra()
        elif p == 'bal':
            combined_plot.plot_neutron_balance()
        elif p == 'rxns':
            combined_plot.plot_rxn_spectra()
        elif p == 'rxnr':
            combined_plot.plot_rxn_ratio()
        elif p == 'rxnc':
            combined_plot.plot_rxn_cumulative()

    print("\nComment. <plot.py/plot_all()> All plotting commands completed.")