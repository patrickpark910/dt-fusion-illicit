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

        self.flibe_u_ng_df  = read_and_sort("./Figures/Data/FLiBe_900K_Li07.5_U238_n-gamma.csv")
        self.flibe_th_ng_df = read_and_sort("./Figures/Data/FLiBe_900K_Li07.5_Th232_n-gamma.csv")
        self.dcll_u_ng_df   = read_and_sort("./Figures/Data/DCLL_900K_Li90.0_U238_n-gamma.csv")
        self.dcll_th_ng_df  = read_and_sort("./Figures/Data/DCLL_900K_Li90.0_Th232_n-gamma.csv")
        self.hcpb_u_ng_df   = read_and_sort("./Figures/Data/HCPB_900K_Li60.0_U238_n-gamma.csv")
        self.hcpb_th_ng_df  = read_and_sort("./Figures/Data/HCPB_900K_Li60.0_Th232_n-gamma.csv")

        self.flibe_u_rr_df = self.flibe_u_rr_df[self.flibe_u_rr_df['fertile_kg/m3'] <= LIMIT_UF4_FLIBE]
        self.flibe_th_rr_df = self.flibe_th_rr_df[self.flibe_th_rr_df['fertile_kg/m3'] <= LIMIT_THF4_FLIBE]
        self.flibe_u_ng_df = self.flibe_u_ng_df[self.flibe_u_ng_df['fertile_kg/m3'] <= LIMIT_UF4_FLIBE]
        self.flibe_th_ng_df = self.flibe_th_ng_df[self.flibe_th_ng_df['fertile_kg/m3'] <= LIMIT_THF4_FLIBE]

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

        titles = [r"FLiBe-UF$_4$", r"FLiBe-ThF$_4$", r"HCPB-UO$_2$", r"HCPB-ThO$_2$", r"DCLL-UO$_2$", r"DCLL-ThO$_2$",]
        dfs    = [self.flibe_u_ng_df, self.flibe_th_ng_df, self.hcpb_u_ng_df, self.hcpb_th_ng_df, self.dcll_u_ng_df, self.dcll_th_ng_df,]

        densities_to_plot = [0.1, 10, 100, 999.99]
        colors = {0.1: '#ff1f5b', 10: '#f48628', 100: '#04cc6c', 999.99: '#0c9edd'}
        labels = {0.1: r'0.1 kg$/$m³', 10: r'10 kg$/$m³', 100: r'100 kg$/$m³', 999.99: r'1000 kg$/$m³'}

        print(f"\n{'Panel':<20} {'Density [kg/m³]':>16} {'Mean [eV]':>12} {'Median [eV]':>13} {'Mode [eV]':>12}")
        print("-" * 75)

        for ax, df, title in zip(axes.flatten(), dfs, titles):

            df = df[df["fertile_kg/m3"].isin(densities_to_plot)]

            df_mean = df.groupby("fertile_kg/m3")["mean"].sum().to_frame()
            df_mean.columns = ["sum"]   
            df = df.merge(df_mean, on="fertile_kg/m3", how="left")
            df['norm_mean'] = df['mean'] / df['sum']

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




    
if __name__ == "__main__":

    for sd in ['PNG','PDF','Data']:
        sd_path = f'./Figures/{sd}/'
        print(f"Comment. <plot.py/plot_all()> Ensuring directory exists: {sd_path}")
        os.makedirs(sd_path, exist_ok=True)


    plot_types = ['tbr','fpr','hist','dfis','fisn']

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
        elif p == 'drdn':
            combined_plot.plot_dRdn()
        elif p == 'rn':
            combined_plot.plot_Rovern()
        elif p == 'leak':
            combined_plot.plot_leakage_spectra()

    print("\nComment. <plot.py/plot_all()> All plotting commands completed.")