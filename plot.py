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

    def __init__(self, show=True, save=False, src_n=None):

        def _load(path):
            df = read_and_sort(path)
            if src_n is not None and 'src-n' in df.columns:
                df = df[df['src-n'] == src_n]
            return df

        # Dataframes of reaction rates (rr)
        self.flibe_u_rr_df  = _load("./Figures/Data/FLiBe_900K_Li07.5_U238_summary.csv")
        self.flibe_th_rr_df = _load("./Figures/Data/FLiBe_900K_Li07.5_Th232_summary.csv")
        self.dcll_u_rr_df   = _load("./Figures/Data/DCLL_900K_Li90.0_U238_summary.csv")
        self.dcll_th_rr_df  = _load("./Figures/Data/DCLL_900K_Li90.0_Th232_summary.csv")
        self.hcpb_u_rr_df   = _load("./Figures/Data/HCPB_900K_Li60.0_U238_summary.csv")
        self.hcpb_th_rr_df  = _load("./Figures/Data/HCPB_900K_Li60.0_Th232_summary.csv")

        self.flibe_u_ebins_df  = _load("./Figures/Data/FLiBe_900K_Li07.5_U238_rxns.csv")
        self.flibe_th_ebins_df = _load("./Figures/Data/FLiBe_900K_Li07.5_Th232_rxns.csv")
        self.dcll_u_ebins_df   = _load("./Figures/Data/DCLL_900K_Li90.0_U238_rxns.csv")
        self.dcll_th_ebins_df  = _load("./Figures/Data/DCLL_900K_Li90.0_Th232_rxns.csv")
        self.hcpb_u_ebins_df   = _load("./Figures/Data/HCPB_900K_Li60.0_U238_rxns.csv")
        self.hcpb_th_ebins_df  = _load("./Figures/Data/HCPB_900K_Li60.0_Th232_rxns.csv")

        # Dataframes of energy-binned flux spectra (for flux ratio plot)
        self.flibe_u_flux_df  = _load("./Figures/Data/FLiBe_900K_Li07.5_U238_flux.csv")
        self.flibe_th_flux_df = _load("./Figures/Data/FLiBe_900K_Li07.5_Th232_flux.csv")
        self.dcll_u_flux_df   = _load("./Figures/Data/DCLL_900K_Li90.0_U238_flux.csv")
        self.dcll_th_flux_df  = _load("./Figures/Data/DCLL_900K_Li90.0_Th232_flux.csv")
        self.hcpb_u_flux_df   = _load("./Figures/Data/HCPB_900K_Li60.0_U238_flux.csv")
        self.hcpb_th_flux_df  = _load("./Figures/Data/HCPB_900K_Li60.0_Th232_flux.csv")

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

        datasets = [ (self.flibe_th_rr_df, r'FLiBe-ThF$_4$', '#66b420', LONG_DASH, 'makima'),
                     (self.hcpb_th_rr_df,  r'HCPB-ThO$_2$',  '#b41f24', LONG_DASH,  2),
                     (self.dcll_th_rr_df,  r'DCLL-ThO$_2$',  '#0047ba', LONG_DASH,  2),
                     (self.flibe_u_rr_df,  r'FLiBe-UF$_4$',  '#66b420', '-',  'makima'),
                     (self.hcpb_u_rr_df,   r'HCPB-UO$_2$',   '#b41f24', '-',  2),
                     (self.dcll_u_rr_df,   r'DCLL-UO$_2$',   '#0047ba', '-',  2), ]

        plt.figure(figsize=(3.5, 3.0))
        ax = plt.gca()
        ax.axhspan(0, 1.00, color='#F0F0F0')

        for df, label, color, linestyle, degree in datasets:

            x_fit = df['fertile_kg/m3']
            y_fit = BLANKET_COVERAGE * df['tbr']

            x_fine = np.linspace(0, 1000, 500)
            if degree == 'makima':
                y_fine = Akima1DInterpolator(x_fit, y_fit, method='makima')(x_fine)
            else:
                y_fine = np.poly1d(np.polyfit(x_fit, y_fit, degree))(x_fine)

            plt.plot(x_fine, y_fine, linestyle=linestyle, linewidth=0.75, color=color, label=label)

        plt.xlabel(r'Fertile isotope density [kg$/$m³]')
        plt.ylabel('Tritium breeding ratio')
        plt.xlim(-25, 1025)
        plt.ylim(0.73, 1.47)

        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_major_locator(MultipleLocator(100))
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(MultipleLocator(0.025))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

        plt.tight_layout()
        leg = plt.legend(loc='lower left', fancybox=False, edgecolor='none',
                         frameon=True, framealpha=0.0, ncol=2, fontsize=6.5)
        leg.get_frame().set_linewidth(0.5)

        if self.save:
            plt.savefig('./Figures/PDF/fig_tbr.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig('./Figures/PNG/fig_tbr.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print(f"Comment. <plot.py/plot_tbr()> Exported TBR plot.")
        else:
            print(f"{C.YELLOW}Comment.{C.END} <plot.py/plot_tbr()> Did NOT export TBR plot due to user setting.")

        if self.show:
            plt.show()

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
                plt.plot(x_fine, y_fine, linestyle=linestyle, linewidth=0.75, color=color, label=label)

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
            leg = plt.legend(loc='upper left', fancybox=False, edgecolor='none',
                            frameon=True, framealpha=0.0, ncol=2, fontsize=6.5)
            leg.get_frame().set_linewidth(0.5)

            if self.save:
                plt.savefig(f'./Figures/PDF/fig_fpr_{v["suffix"]}.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
                plt.savefig(f'./Figures/PNG/fig_fpr_{v["suffix"]}.png', bbox_inches='tight', pad_inches=0.01, format='png')
                print("Exported fissile production per year plot for all blankets.")
            else:
                print("Did not export fissile production per year plot due to user setting.")

            if self.show:
                plt.show()

            plt.close('all')



    def stats(self):
        """Print spectral statistics for flux and fertile (n,gamma) capture, plus TBR and FPR."""

        # Conversion factors: rxns/src-n  →  kg/yr
        #   1. NPS_FUS [n/s] × 3.156e7 [s/yr]  =  source neutrons per year
        #   2. ÷ AVO [atoms/mol]                =  moles of fissile atoms per year
        #   3. × AMU [g/mol]                    =  grams per year
        #   4. ÷ 1e3                            =  kilograms per year
        #   5. × BLANKET_COVERAGE               =  geometric correction
        tot_n_per_yr = NPS_FUS * 3.156e+7
        pu239_conv = (tot_n_per_yr / AVO * AMU_PU239 / 1e3) * BLANKET_COVERAGE
        u233_conv  = (tot_n_per_yr / AVO * AMU_U233  / 1e3) * BLANKET_COVERAGE

        titles  = [r"FLiBe-UF$_4$", r"FLiBe-ThF$_4$", r"HCPB-UO$_2$", r"HCPB-ThO$_2$", r"DCLL-UO$_2$", r"DCLL-ThO$_2$"]
        ebins   = [self.flibe_u_ebins_df, self.flibe_th_ebins_df, self.hcpb_u_ebins_df, self.hcpb_th_ebins_df, self.dcll_u_ebins_df, self.dcll_th_ebins_df]
        ng_cols = ['U238(n,g)', 'Th232(n,g)', 'U238(n,g)', 'Th232(n,g)', 'U238(n,g)', 'Th232(n,g)']
        fis_amu = [AMU_PU239, AMU_U233, AMU_PU239, AMU_U233, AMU_PU239, AMU_U233]
        convs   = [pu239_conv, u233_conv, pu239_conv, u233_conv, pu239_conv, u233_conv]
        fluxes  = [self.flibe_u_flux_df, self.flibe_th_flux_df, self.hcpb_u_flux_df, self.hcpb_th_flux_df, self.dcll_u_flux_df, self.dcll_th_flux_df]
        rr_dfs  = [self.flibe_u_rr_df, self.flibe_th_rr_df, self.hcpb_u_rr_df_corr, self.hcpb_th_rr_df_corr, self.dcll_u_rr_df, self.dcll_th_rr_df]

        densities = [0, 150, 999.99]

        print(f"\n{'Panel':<12} {'rho':>8}  "
              f"{'TBR':>15} {'FPR [rxn/src]':>15} {'FPR [kg/yr]':>15}  "
              f"{'Flux Mean':>15} {'Flux Med':>15} {'Flux Mode':>15}  "
              f"{'(n,g) Mean':>15} {'(n,g) Med':>15} {'(n,g) Mode':>15}  "
              f"{'(n,Xt) Mean':>15} {'(n,Xt) Med':>15} {'(n,Xt) Mode':>15}")
        print("-" * 210)

        for title, ebin_df, ng_col, conv, flux_df, rr_df in zip(titles, ebins, ng_cols, convs, fluxes, rr_dfs):
            for d in densities:
                clean_title = strip_latex(title)
                d_label = f"{d:.1f}" if d < 999 else "1000.0"

                # TBR and FPR from summary dataframe
                rr_avail = rr_df['fertile_kg/m3'].unique()
                rr_rho = min(rr_avail, key=lambda x: abs(x - d))
                rr_row = rr_df[rr_df['fertile_kg/m3'] == rr_rho].iloc[0]
                tbr = BLANKET_COVERAGE * rr_row['tbr']
                fpr_rxn = rr_row[ng_col]
                fpr_kg = fpr_rxn * conv

                # Flux spectrum stats
                flux_avail = flux_df['fertile_kg/m3'].unique()
                flux_rho = min(flux_avail, key=lambda x: abs(x - d))
                flux_sub = flux_df[flux_df['fertile_kg/m3'] == flux_rho].sort_values('energy mid [eV]')
                f_mean, f_med, f_mode, _ = spectrum_stats(
                    flux_sub['energy mid [eV]'].values, flux_sub['mean'].values)

                # Fertile (n,gamma) capture spectrum stats
                ebin_avail = ebin_df['fertile_kg/m3'].unique()
                ebin_rho = min(ebin_avail, key=lambda x: abs(x - d))
                ebin_sub = ebin_df[ebin_df['fertile_kg/m3'] == ebin_rho].sort_values('energy mid [eV]')
                ng_mean, ng_med, ng_mode, _ = spectrum_stats(
                    ebin_sub['energy mid [eV]'].values, ebin_sub[ng_col].values)

                # Li(n,Xt) tritium production spectrum stats (Li-6 + Li-7)
                li_vals = ebin_sub['Li6(n,t)'].values + ebin_sub['Li7(n,t)'].values
                li_mean, li_med, li_mode, _ = spectrum_stats(
                    ebin_sub['energy mid [eV]'].values, li_vals)

                print(f"{clean_title:<12} {d_label:>8}  "
                      f"{tbr:>15.4e} {fpr_rxn:>15.4e} {fpr_kg:>15.4e}  "
                      f"{f_mean:>15.4e} {f_med:>15.4e} {f_mode:>15.4e}  "
                      f"{ng_mean:>15.4e} {ng_med:>15.4e} {ng_mode:>15.4e}  "
                      f"{li_mean:>15.4e} {li_med:>15.4e} {li_mode:>15.4e}")


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
        ng_cols = ['U238(n,g)', 'Th232(n,g)', 'U238(n,g)', 'Th232(n,g)', 'U238(n,g)', 'Th232(n,g)']

        densities_to_plot = [0.01, 100, 1000, 4000] # 10,
        colors = {0.01: '#ff1f5b',     100: '#04cc6c',   1000: '#0c9edd', 4000: '#f48628',  }
        labels = {0.01: r'0.01 kg$/$m³', 100: r'100 kg$/$m³', 1000: r'1000 kg$/$m³', 4000: r'4000 kg$/$m³'}

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
            plt.savefig(f'./Figures/PDF/fig_hist_4000kgm3.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig(f'./Figures/PNG/fig_hist_4000kgm3.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print(f"\nExported cumulative normalized histogram plots.")
        else:
            print(f"\nDid not export cumulative normalized histogram plot.")

        if self.show: plt.show()
        plt.close('all')


    def plot_flux_spectrum(self):
        """
        Flux spectrum plot: phi(E) vs. energy.

        3x2 grid (rows = FLiBe / HCPB / DCLL, cols = U238 / Th232).
        Each panel overlays the flux at several fertile loadings
        (0, 10, 100, 1000 kg/m³).

        Also prints the characteristic energies of every spectrum drawn:
        the flux-weighted MEAN energy, the cumulative-flux MEDIAN energy,
        the MODE (energy of the peak flux bin, i.e. the 14.07 MeV source),
        and the 2ND MODE (tallest peak below the source), all in eV.
        The second mode is also marked on each curve with a triangle.
        """
        print(f"\nPlotting flux spectra vs. energy...")

        from matplotlib.lines import Line2D

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
        labels_plain = {0: '0 kg/m3', 10: '10 kg/m3', 100: '100 kg/m3', 999.99: '1000 kg/m3'}
        row_ylim = {0: log_buffer(1e-3, 1e1), 1: log_buffer(1e-3, 1e1), 2: log_buffer(1e-3, 1e1)}  # 2: log_buffer(1e-5, 1e1)}

        # ---- printed-statistics table header ----
        print("\nFlux-spectrum characteristic energies [eV]")
        print("  (mean = flux-weighted; 5/25/median/75/95 % = energy at that fraction of")
        print("   cumulative flux; mode = peak bin; 2nd mode = tallest peak below 14.07 MeV,")
        print("   also marked on the plot)\n")
        print(f"  {'Panel':<12s} {'Loading':<11s} {'Mean E':>11s} {'5%':>11s} {'25%':>11s} {'Median E':>11s} {'75%':>11s} {'95%':>11s} {'Mode E':>11s} {'2nd Mode':>11s}")
        print("  " + "-" * 120)

        fig, axes = plt.subplots(3, 2, figsize=(7, 9), sharex=True, sharey=False)

        for df_flux, title, row, col in panels:
            ax = axes[row, col]

            available = df_flux['fertile_kg/m3'].unique()

            for rho in densities_to_plot:
                closest = min(available, key=lambda x: abs(x - rho))
                spec = df_flux[df_flux['fertile_kg/m3'] == closest].sort_values('energy mid [eV]')

                E_mid = spec['energy mid [eV]'].values
                phi   = spec['mean'].values.copy()

                # --- characteristic energies of this spectrum ---
                mean_E, median_E, mode_E, q = spectrum_stats(E_mid, phi, quantiles=(0.05, 0.25, 0.75, 0.95))
                mode2_E, mode2_phi = second_mode(E_mid, phi)   # tallest peak below 14.07 MeV
                print(f"  {strip_latex(title):<12s} {labels_plain[rho]:<11s} "
                      f"{mean_E:>11.3e} {q[0.05]:>11.3e} {q[0.25]:>11.3e} {median_E:>11.3e} "
                      f"{q[0.75]:>11.3e} {q[0.95]:>11.3e} {mode_E:>11.3e} {mode2_E:>11.3e}")

                # Mask zeros / negatives to avoid log-scale artifacts
                phi[phi <= 0] = np.nan

                ax.loglog(E_mid, phi, color=colors[rho], linewidth=0.75, label=labels[rho])

                # Mark the second mode (peak below the 14.07 MeV source peak)
                # if np.isfinite(mode2_E):
                #     ax.plot(mode2_E, mode2_phi, marker='v', markersize=5,
                #             color=colors[rho], markeredgecolor='black',
                #             markeredgewidth=0.4, linestyle='None', zorder=5,
                #             label='_nolegend_')

            print("  " + "-" * 120)  # separate panels in the table

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

            # Legend with panel title (+ one proxy entry for the 2nd-mode markers)
            handles, lbls = ax.get_legend_handles_labels()
            handles.append(Line2D([], [], marker='v', linestyle='None', markersize=5,
                                   color='0.4', markeredgecolor='black', markeredgewidth=0.4))
            # lbls.append('2nd mode')
            leg = ax.legend(handles, lbls, title=title, loc='lower right', fontsize=7, ncol=1,
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


    def flux_nominal(self):
        """
        Nominal (0 kg/m³) flux spectra for FLiBe, HCPB, and DCLL.
        Uses the U-series data since the clean blanket is identical for U and Th.

        Produces two figures:
          fig_flux_nominal       — clean spectra
          fig_flux_nominal_stats — same spectra annotated with mean, median,
                                   2nd mode, and quantiles from spectrum_stats
        """
        print(f"\nPlotting nominal (0 kg/m³) flux spectra...")

        from matplotlib.ticker import LogLocator
        from matplotlib.lines import Line2D

        datasets = [
            (self.flibe_u_flux_df, 'FLiBe', '#66b420'),
            (self.hcpb_u_flux_df,  'HCPB',  '#b41f24'),
            (self.dcll_u_flux_df,  'DCLL',  '#0047ba'),
        ]

        spectra = []
        for df_flux, label, color in datasets:
            available = df_flux['fertile_kg/m3'].unique()
            closest = min(available, key=lambda x: abs(x - 0))
            spec = df_flux[df_flux['fertile_kg/m3'] == closest].sort_values('energy mid [eV]')
            spectra.append((spec['energy mid [eV]'].values,
                            spec['mean'].values.copy(), label, color))

        def setup_flux_axes(ax):
            ax.set_xlabel('Neutron energy [eV]')
            ax.set_ylabel(r'Flux [n$/$cm²$/$src-n]')
            ax.set_xlim(log_buffer(1e0, 20e6))
            ax.set_ylim(log_buffer(1e-3, 1e1))
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.grid(axis='both', which='major', linestyle='-', linewidth=0.5)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'1e{int(np.log10(x))}' if x > 0 else ''))
            ax.yaxis.set_major_locator(LogLocator(base=10, numticks=20))
            ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'1e{int(np.log10(x))}' if x > 0 else ''))

        def plot_spectra(ax):
            for E_mid, phi, label, color in spectra:
                phi_plot = phi.copy()
                phi_plot[phi_plot <= 0] = np.nan
                ax.loglog(E_mid, phi_plot, color=color, linewidth=0.75, label=label)

        # ---- Figure 1: Clean spectra ----
        plt.figure(figsize=(3.5, 3.0))
        ax = plt.gca()
        plot_spectra(ax)
        setup_flux_axes(ax)
        plt.tight_layout()
        leg = plt.legend(loc='lower right', fancybox=False, edgecolor='none',
                         frameon=True, framealpha=0.0, ncol=1, fontsize=6.5)
        leg.get_frame().set_linewidth(0.5)

        if self.save:
            plt.savefig('./Figures/PDF/fig_flux_nominal.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig('./Figures/PNG/fig_flux_nominal.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print("Exported nominal flux spectrum plot: fig_flux_nominal")
        else:
            print("Did not export nominal flux spectrum plot due to user setting.")
        if self.show:
            plt.show()
        plt.close('all')

        # ---- Figure 2: Annotated with spectral statistics ----
        print(f"\n  {'Blanket':<8s} {'Mean':>11s} {'q10':>11s} {'q25':>11s} "
              f"{'Median':>11s} {'q75':>11s} {'q90':>11s} {'2nd Mode':>11s}")
        print("  " + "-" * 88)

        plt.figure(figsize=(3.5, 3.0))
        ax = plt.gca()
        plot_spectra(ax)

        for E_mid, phi, label, color in spectra:
            mean_E, median_E, mode_E, q = spectrum_stats(
                E_mid, phi, quantiles=(0.10, 0.25, 0.75, 0.90))
            mode2_E, mode2_phi = second_mode(E_mid, phi)

            print(f"  {label:<8s} {mean_E:>11.3e} {q[0.10]:>11.3e} {q[0.25]:>11.3e} "
                  f"{median_E:>11.3e} {q[0.75]:>11.3e} {q[0.90]:>11.3e} {mode2_E:>11.3e}")

            ax.axvline(mean_E,   color=color, linestyle='-',  linewidth=0.5, alpha=0.6)
            ax.axvline(median_E, color=color, linestyle='--', linewidth=0.5, alpha=0.6)
            for p_val in q.values():
                ax.axvline(p_val, color=color, linestyle=':', linewidth=0.4, alpha=0.4)
            if np.isfinite(mode2_E):
                ax.plot(mode2_E, mode2_phi, marker='v', markersize=5,
                        color=color, markeredgecolor='black',
                        markeredgewidth=0.4, linestyle='None', zorder=5)

        setup_flux_axes(ax)
        plt.tight_layout()

        blanket_handles = [Line2D([0], [0], color=c, lw=0.75, label=l)
                           for _, _, l, c in spectra]
        stat_handles = [
            Line2D([0], [0], color='gray', ls='-',  lw=0.75, label='Mean'),
            Line2D([0], [0], color='gray', ls='--', lw=0.75, label='Median'),
            Line2D([], [], marker='v', linestyle='None', markersize=5,
                   color='gray', markeredgecolor='black',
                   markeredgewidth=0.4, label='2nd mode'),
            Line2D([0], [0], color='gray', ls=':', lw=0.75, label='Quantiles'),
        ]
        leg = ax.legend(handles=blanket_handles + stat_handles,
                        loc='lower right', fancybox=False, edgecolor='none',
                        frameon=True, framealpha=0.0, ncol=2, fontsize=6.5)
        leg.get_frame().set_linewidth(0.5)

        if self.save:
            plt.savefig('./Figures/PDF/fig_flux_nominal_stats.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
            plt.savefig('./Figures/PNG/fig_flux_nominal_stats.png', bbox_inches='tight', pad_inches=0.01, format='png')
            print("Exported annotated nominal flux spectrum plot: fig_flux_nominal_stats")
        else:
            print("Did not export annotated nominal flux spectrum plot due to user setting.")
        if self.show:
            plt.show()
        plt.close('all')


    def plot_RoverN(self):
        print(f"\nPlotting fissile material divided by fertile density vs. fertile density...")

        tot_n_per_yr = NPS_FUS * 3.156e+7

        # Gotta then convert fissile atoms/yr to kg/yr
        pu239_at_to_kg_per_yr = (tot_n_per_yr / AVO * AMU_PU239 / 1e3) * BLANKET_COVERAGE
        u233_at_to_kg_per_yr  = (tot_n_per_yr / AVO * AMU_U233  / 1e3) * BLANKET_COVERAGE

        # Setting up x, y separately here so you can remove the impurity/wppm-magnitude cases --ppark 2025-08-06
        # cut out the first point since dividing by zero.
        # helium cooled pebble bed
        x6, y6 =     self.hcpb_u_rr_df_corr['fertile_kg/m3'][1:], pu239_at_to_kg_per_yr * self.hcpb_u_rr_df_corr['U238(n,g)'][1:]
        x5, y5 =    self.hcpb_th_rr_df_corr['fertile_kg/m3'][1:], u233_at_to_kg_per_yr * self.hcpb_th_rr_df_corr['Th232(n,g)'][1:]

        # dual coolant lead lithium
        x4, y4 =   self.dcll_u_rr_df['fertile_kg/m3'][1:],  pu239_at_to_kg_per_yr *  self.dcll_u_rr_df['U238(n,g)'][1:]
        x3, y3 =  self.dcll_th_rr_df['fertile_kg/m3'][1:],  u233_at_to_kg_per_yr * self.dcll_th_rr_df['Th232(n,g)'][1:]

        # flibe
        x2, y2 =  self.flibe_u_rr_df['fertile_kg/m3'][1:],  pu239_at_to_kg_per_yr *  self.flibe_u_rr_df['U238(n,g)'][1:]
        x1, y1 = self.flibe_th_rr_df['fertile_kg/m3'][1:],  u233_at_to_kg_per_yr * self.flibe_th_rr_df['Th232(n,g)'][1:]

        akima1 = Akima1DInterpolator(x1, y1/x1)
        akima2 = Akima1DInterpolator(x2, y2/x2)
        akima3 = Akima1DInterpolator(x3, y3/x3)
        akima4 = Akima1DInterpolator(x4, y4/x4)
        akima5 = Akima1DInterpolator(x5, y5/x5)
        akima6 = Akima1DInterpolator(x6, y6/x6)

        from matplotlib.ticker import LogLocator

        for scale in ['lin', 'log']:
            if scale == 'lin':
                x_fine = np.linspace(x1.min(), x1.max(), 300)
            else:
                x_fine = np.geomspace(x1.min(), x1.max(), 300)

            plt.figure(figsize=(3.5, 3.0))

            if scale == 'lin':
                plt.plot(x_fine, akima2(x_fine), linestyle='-',       linewidth=0.75, color='#66b420',label=r'FLiBe-UF$_4$' )
                plt.plot(x_fine, akima1(x_fine), linestyle=LONG_DASH, linewidth=0.75, color='#66b420',label=r'FLiBe-ThF$_4$')
                plt.plot(x_fine, akima6(x_fine), linestyle='-',       linewidth=0.75, color='#b41f24',label=r'HCPB-UO$_2$'  )
                plt.plot(x_fine, akima5(x_fine), linestyle=LONG_DASH, linewidth=0.75, color='#b41f24',label=r'HCPB-ThO$_2$' )
                plt.plot(x_fine, akima3(x_fine), linestyle=LONG_DASH, linewidth=0.75, color='#0047ba',label=r'DCLL-ThO$_2$' )
                plt.plot(x_fine, akima4(x_fine), linestyle='-',       linewidth=0.75, color='#0047ba',label=r'DCLL-UO$_2$'  )

            else:
                plt.plot(x_fine, akima1(x_fine), linestyle=LONG_DASH, linewidth=0.75, color='#66b420',label=r'FLiBe-ThF$_4$')
                plt.plot(x_fine, akima5(x_fine), linestyle=LONG_DASH, linewidth=0.75, color='#b41f24',label=r'HCPB-ThO$_2$' )
                plt.plot(x_fine, akima3(x_fine), linestyle=LONG_DASH, linewidth=0.75, color='#0047ba',label=r'DCLL-ThO$_2$' )
                plt.plot(x_fine, akima2(x_fine), linestyle='-',       linewidth=0.75, color='#66b420',label=r'FLiBe-UF$_4$' )
                plt.plot(x_fine, akima6(x_fine), linestyle='-',       linewidth=0.75, color='#b41f24',label=r'HCPB-UO$_2$'  )
                plt.plot(x_fine, akima4(x_fine), linestyle='-',       linewidth=0.75, color='#0047ba',label=r'DCLL-UO$_2$'  )

            plt.xlabel(r'Fertile isotope density [kg$/$m³]')
            plt.ylabel(r'Fissile prod rate [kg$/$yr] per fertile density')
            plt.ylim(-0.25,8.25)

            ax = plt.gca()
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

            if scale == 'lin':
                plt.xlim(-25, 1025)
                ax.xaxis.set_major_locator(MultipleLocator(100))
                ax.xaxis.set_minor_locator(MultipleLocator(50))
            else:
                plt.xlim(10**(-2.15), 10**(3.15))
                ax.set_xscale('log')
                ax.xaxis.set_major_locator(LogLocator(base=10))
                ax.xaxis.set_minor_locator(LogLocator(base=10, subs='auto', numticks=10))

            ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)

            plt.tight_layout()

            if scale == 'lin':
                leg = plt.legend(fancybox=False, edgecolor='black', frameon=False, framealpha=.75, ncol=1)
            else:
                leg = plt.legend(loc='center left', fancybox=False, edgecolor='black', frameon=False, framealpha=.75, ncol=2)
            leg.get_frame().set_linewidth(0.5)

            if self.save:
                plt.savefig(f'./Figures/PDF/fig_RoverN_{scale}.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
                plt.savefig(f'./Figures/PNG/fig_RoverN_{scale}.png', bbox_inches='tight', pad_inches=0.01, format='png')
                print(f"Exported R/n plot ({scale}).")
            else:
                print(f"Did not export R/n plot ({scale}) due to user setting.")

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

            # ── Spectral statistics (per-lethargy peak defines the mode) ──
            mean_E, median_E, mode_E, _ = spectrum_stats(
                E_mid, counts, mode_weight=flux_per_lethargy)

            short_label = strip_latex(case["label"])   # for terminal printing
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


    plot_types = ['stat','tbr','fpr','hist','fluxs','fluxn','rovern','leak']

    parser = argparse.ArgumentParser(description=f"Choose plots with -t flag, multiple separated by spaces: {plot_types}")

    parser.add_argument("-t", "--plot_type",
                        type=str, nargs="+", default=plot_types,
                        help=f"Specify plot, multiple separated by spaces: {plot_types}. Defaults to all plots." )

    parser.add_argument("-p", "--n_particles",
                        type=lambda x: int(float(x)), default=N_PARTICLES,
                        help=f"Specify number of particles in integer scientific notation, e.g. 4e5" )

    parser.add_argument("-c", "--n_cycles",
                        type=lambda x: int(float(x)), default=N_CYCLES,
                        help=f"Specify number of cycles in integers, e.g. 25" )

    parser.add_argument("--no_show",
                        dest="show", action="store_false",
                        help="Disable showing the plots.")

    parser.add_argument("--no_save",
                        dest="save", action="store_false",
                        help="Disable saving the plots. ")

    parser.set_defaults(show=True, save=True)

    args = parser.parse_args()
    plot_type  = [p.lower() if isinstance(p, str) else p for p in args.plot_type]
    plot_show  = args.show
    plot_save  = args.save
    src_n      = args.n_particles * args.n_cycles

    combined_plot = Plot(show=plot_show, save=plot_save, src_n=src_n)

    for p in plot_type:
        if   p == 'stat':
            combined_plot.stats()
        elif p == 'tbr':
            combined_plot.plot_tbr()
        elif p == 'fpr':
            combined_plot.plot_fpr()
        elif p == 'hist':
            combined_plot.plot_hist()
        elif p == 'fluxs':
            combined_plot.plot_flux_spectrum()
        elif p == 'fluxn':
            combined_plot.flux_nominal()
        elif p == 'rovern':
            combined_plot.plot_RoverN()
        elif p == 'leak':
            combined_plot.plot_leakage_spectra()


    print("\nComment. <plot.py/plot_all()> All plotting commands completed.")
