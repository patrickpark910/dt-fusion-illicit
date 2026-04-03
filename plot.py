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


def plot_all():

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
        if p == 'tbr':
            combined_plot.plot_tbr()
        elif p == 'fpr':
            combined_plot.plot_fpr()
        elif p == 'hist':
            combined_plot.plot_cum_norm_histogram()
        elif p == 'dfis':
            combined_plot.plot_dfisdfer()
        elif p == 'fisn':
            combined_plot.plot_fisovern()

    print("\nComment. <plot.py/plot_all()> All plotting commands completed.")


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

        self.show, self.save = show, save


    def plot_flux(self):
        
        # Create figure with 2x2 subplots
        fig, axes = plt.subplots(2, 2, figsize=(18, 14))
        
        # Define formatter for scientific notation on each tick
        def sci_notation(x, pos):
            return f'{x:.0e}'
    
        # Row 0: FLiBe-U
        # Plot 0,0: FLiBe-U Log-log plot (full range)
        df1 = self.flibe_u_flux_df
        df2 = self.dcll_u_flux_df

        loadings = [0.0, 15.0, 150.0, 1000.0]
        for loading in loadings:
            df_subset = df1[df1['fertile_kg/m3'] == loading]
            # df_zoom = df_subset[(df_subset['energy mid [eV]'] >= energy_min) & 
            #                    (df_subset['energy mid [eV]'] <= energy_max)]
            axes[0,0].plot(df_subset['energy mid [eV]'], df_subset['mean'], 
                    label=f'{int(loading)} kg/m³', linewidth=2, alpha=0.8)
        
        axes[0,0].set_xlim(1e0, 2e7)
        axes[0,0].set_ylim(1e-5, 1e2)
        axes[0,0].set_xscale('log')
        axes[0,0].set_yscale('log')
        axes[0,0].set_xlabel('Energy [eV]', fontsize=12)
        axes[0,0].set_ylabel('Flux [n/cm²-s]', fontsize=12)
        axes[0,0].set_title('Neutron flux spectrum, log-log (FLiBe-U)', fontsize=16)
        axes[0,0].legend(title='U238/breeder [kg/m³]', fontsize=12, loc='best')
        axes[0,0].grid(True, which='both', alpha=0.3)
        
        # Plot 0,1: FLiBe-U Linear-linear plot (0-10 MeV range)
        for loading in loadings:
            df_subset = df1[df1['fertile_kg/m3'] == loading]

            # df_zoom = df_subset[(df_subset['energy mid [eV]'] >= energy_min) & 
            #                     (df_subset['energy mid [eV]'] <= energy_max)]

            axes[0,1].plot(df_subset['energy mid [eV]'], df_subset['mean'], 
                    label=f'{int(loading)} kg/m³', linewidth=2, alpha=0.8)
        
        axes[0,1].set_xlim(1e5, 5e6)
        axes[0,1].xaxis.set_major_formatter(FuncFormatter(sci_notation))
        axes[0,1].set_ylim(0,0.6)
        axes[0,1].set_xlabel('Energy [eV]', fontsize=12)
        axes[0,1].set_ylabel('Flux [n/cm²-s]', fontsize=12)
        axes[0,1].set_title(r'Neutron flux spectrum, lin-lin, 0.1$-$5 MeV (FLiBe-U)', fontsize=16)
        axes[0,1].legend(title='U238/breeder [kg/m³]', fontsize=12, loc='best')
        axes[0,1].grid(True, alpha=0.3)
        
        # Row 1: DCLL-U
        # Plot 1,0: DCLL-U Log-log plot (full range)

        for loading in loadings:
            df_subset = df2[df2['fertile_kg/m3'] == loading]

            axes[1,0].plot(df_subset['energy mid [eV]'], df_subset['mean'], 
                    label=f'{int(loading)} kg/m³', linewidth=2, alpha=0.8)
        
        axes[1,0].set_xlim(1e0,2e7)
        axes[1,0].set_ylim(1e-5,1e2)
        axes[1,0].set_xscale('log')
        axes[1,0].set_yscale('log')
        axes[1,0].set_xlabel('Energy [eV]', fontsize=12)
        axes[1,0].set_ylabel('Flux [n/cm²-s]', fontsize=12)
        axes[1,0].set_title('Neutron flux spectrum, log-log (DCLL-U)', fontsize=16)
        axes[1,0].legend(title='U238/breeder [kg/m³]', fontsize=12, loc='best')
        axes[1,0].grid(True, which='both', alpha=0.3)
        
        # Plot 1,1: DCLL-U Linear-linear plot (0-10 MeV range)
        for loading in loadings:
            df_subset = df2[df2['fertile_kg/m3'] == loading]

            # df_zoom = df_subset[(df_subset['energy mid [eV]'] >= energy_min) & 
            #                    (df_subset['energy mid [eV]'] <= energy_max)]

            axes[1,1].plot(df_subset['energy mid [eV]'], df_subset['mean'], 
                    label=f'{int(loading)} kg/m³', linewidth=2, alpha=0.8)
        
        axes[1,1].set_xlim(1e3, 5e6)
        axes[1,1].xaxis.set_major_formatter(FuncFormatter(sci_notation))
        axes[1,1].set_ylim(0,3.5)
        axes[1,1].set_xlabel('Energy [eV]', fontsize=12)
        axes[1,1].set_ylabel('Flux [n/cm²-s]', fontsize=12)
        axes[1,1].set_title(r'Neutron flux spectrum, lin-lin, 0.1$-$5 MeV (DCLL-U)', fontsize=16)
        axes[1,1].legend(title='U238/breeder [kg/m³]', fontsize=12, loc='best')
        axes[1,1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if self.save:
            plt.savefig(f'./Figures/PDF/fig_flux_comparison.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/PNG/fig_flux_comparison.png', bbox_inches='tight', format='png')
            print(f"Exported flux spectrum plots.")
            print(f"Plot saved with {len(loadings)} different fertile loadings")
            print(f"Fertile loadings: {loadings} kg/m³")
        else:
            print(f"Did not export flux spectrum plots.")
        
        if self.show: 
            plt.show()
        
        plt.close('all')


    def plot_tbr(self):

        print(f"\nComment. <plot.py/plot_tbr()> Plotting tritium breeding ratio vs. fertile density...")

        # -------------------------------------------------------------
        # Weigh HCPB TBR by correction factor derived from wedge modelinestyle
        # -------------------------------------------------------------

        df_u    = pd.DataFrame(HCPB_CONV_U_TBR,  columns=['fertile_kgm3', 'ratio'])
        df_th   = pd.DataFrame(HCPB_CONV_TH_TBR, columns=['fertile_kgm3', 'ratio'])

        # NB. np.interp(x_values_to_calculate, original_x, original_y)
        def get_ratio(x, df_fertile):
            return np.interp(x, df_fertile['fertile_kgm3'], df_fertile['ratio'])

        # Apply interpolation to the Uranium dataframe
        hcpb_u_rr_df_corr = self.hcpb_u_rr_df.copy()
        hcpb_u_rr_df_corr['tbr'] = hcpb_u_rr_df_corr['tbr'] * get_ratio(hcpb_u_rr_df_corr['fertile_kg/m3'], df_u)

        # Apply interpolation to the Thorium dataframe
        hcpb_th_rr_df_corr = self.hcpb_th_rr_df.copy()
        hcpb_th_rr_df_corr['tbr'] = hcpb_th_rr_df_corr['tbr'] * get_ratio(hcpb_th_rr_df_corr['fertile_kg/m3'], df_th)


        # -------------------------------------------------------------
        # Plot them!
        # -------------------------------------------------------------

        # Load the dataframes to plot, label, color, marker, markersize, linestyle, and polynomial fit
        datasets = [ (self.flibe_u_rr_df,  r'FLiBe-UF$_4$',  '#66b420', 'o', 5, '-',  'makima'),
                     (self.flibe_th_rr_df, r'FLiBe-ThF$_4$', '#66b420', '+', 12, '--', 'makima'),
                     (hcpb_u_rr_df_corr,   r'HCPB-UO$_2$',   '#b41f24', 's', 6, '-',  2),
                     (hcpb_th_rr_df_corr,  r'HCPB-ThO$_2$',  '#b41f24', '1', 13, '--', 2),
                     (self.dcll_u_rr_df,   r'DCLL-UO$_2$',   '#0047ba', '^', 8, '-',  2),
                     (self.dcll_th_rr_df,  r'DCLL-ThO$_2$',  '#0047ba', 'x', 10, '--', 2),
                   ]

        # Select a subset of fertile kg/m3 to SHOW to avoid cluttering the low end with markers
        selected1 = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150,]  # [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 250, 500, 750, 999.99]
        selected2 = [0, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]  # [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 250, 500, 750, 999.99]

        # Set view limits for each of the plots
        view_limits = [ # Zoomed in on (0, 155) fertile kg/m3
                        {'xlim': (-5, 155),   'xmax': 150, 'x_major': 15,  'x_minor': 5,  
                         'ylim': (1.04, 1.46), 'suffix': '0150kgm3', 'fertile_kgm3': selected1,   # 'ylim': (0.94, 1.26)
                         'leg_loc': 'center right',}, # leg_loc options: best or upper/center/lower left/right/center
                        # Zoomed out to (0, 1000) fertile kg/m3
                        {'xlim': (-25, 1025), 'xmax':1000, 'x_major': 100, 'x_minor': 50, 
                         'ylim': (0.84, 1.46), 'suffix': '1000kgm3', 'fertile_kgm3': selected2,
                         'leg_loc': 'upper right'} ] 
  

        for v in view_limits:

            plt.figure(figsize=(8, 6))
            ax = plt.gca()
            ax.axhspan(0, 1.00, color='#e0ded8')            

            for df, label, color, marker, markersize, linestyle, degree in datasets:

                selected = v['fertile_kgm3']
                mask = df['fertile_kg/m3'].apply(lambda d: any(np.isclose(d, s) for s in selected))
                df_filtered = df[mask]

                x_filtered = df_filtered['fertile_kg/m3']
                y_filtered = BLANKET_COVERAGE * df_filtered['tbr']

                # Plot filtered data
                plt.scatter(x_filtered, y_filtered, marker=marker, s=markersize*5, color=color, zorder=3)

                # But still fit the polynomial to all of the data you did compute
                x_fit = df['fertile_kg/m3']
                y_fit = BLANKET_COVERAGE * df['tbr']
                # print(x_fit)

                if degree == 'makima':
                    x_fine = np.linspace(0, v['xmax'], 500)
                    y_fine    = Akima1DInterpolator(x_fit, y_fit, method='makima')(x_fine) # poly_func(x_fine)

                else:
                    coeffs    = np.polyfit(x_fit, y_fit, degree)
                    y_fine    = np.poly1d(coeffs)(x_fine)

                # Plot interpolation
                plt.plot(x_fine, y_fine, linestyle, linewidth=1, color=color)

                # Dummy plots for legend -- bit of a hack lmao -- ppark
                plt.plot([9e8,9e9], [9e8,9e9], marker+linestyle, markersize=markersize, linewidth=1, color=color, label=label)

            # Labeling
            plt.xlabel(r'Fertile isotope density in breeder [kg$/$m$^3$]')
            plt.ylabel('Tritium breeding ratio')

            # Set Dynamic Limits and Locators
            plt.xlim(v['xlim'])
            plt.ylim(v['ylim']) 
            
            ax.xaxis.set_ticks_position('both')
            ax.xaxis.set_major_locator(MultipleLocator(v['x_major']))
            ax.xaxis.set_minor_locator(MultipleLocator(v['x_minor']))
            
            ax.yaxis.set_ticks_position('both')
            ax.yaxis.set_major_locator(MultipleLocator(0.05))
            ax.yaxis.set_minor_locator(MultipleLocator(0.01))
            
            ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
            ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

            plt.tight_layout()
            leg = plt.legend(loc=v['leg_loc'], fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=3)
            leg.get_frame().set_linewidth(0.5) 

            if self.save:
                plt.savefig(f'./Figures/PDF/fig_tbr_{v["suffix"]}.pdf', bbox_inches='tight', format='pdf')
                plt.savefig(f'./Figures/PNG/fig_tbr_{v["suffix"]}.png', bbox_inches='tight', format='png')
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

        # -------------------------------------------------------------
        # Weigh HCPB TBR by correction factor derived from wedge modelinestyle
        # -------------------------------------------------------------

        df_u    = pd.DataFrame(HCPB_CONV_U_FPR,  columns=['fertile_kgm3', 'ratio'])
        df_th   = pd.DataFrame(HCPB_CONV_TH_FPR, columns=['fertile_kgm3', 'ratio'])

        # NB. np.interp(x_values_to_calculate, original_x, original_y)
        def get_u_ratio(x):
            return np.interp(x, df_u['fertile_kgm3'], df_u['ratio'])

        def get_th_ratio(x):
            return np.interp(x, df_th['fertile_kgm3'], df_th['ratio'])

        # Apply interpolation to the Uranium dataframe
        hcpb_u_rr_df_corr = self.hcpb_u_rr_df.copy()
        hcpb_u_rr_df_corr['U238(n,g)'] = hcpb_u_rr_df_corr['U238(n,g)'] * get_u_ratio(hcpb_u_rr_df_corr['fertile_kg/m3'])

        # Apply interpolation to the Thorium dataframe
        hcpb_th_rr_df_corr = self.hcpb_th_rr_df.copy()
        hcpb_th_rr_df_corr['Th232(n,g)'] = hcpb_th_rr_df_corr['Th232(n,g)'] * get_th_ratio(hcpb_th_rr_df_corr['fertile_kg/m3'])

        # Constants and conversion factors
        # <utilities.py> default NPS_FUS = 500 MJ/s * 3.546e17 n/MJ = 1.773e+20 n/s
        tot_n_per_yr = NPS_FUS * 3.156e+7
        pu239_conv   = (tot_n_per_yr / AVO * AMU_PU239 / 1e3) * BLANKET_COVERAGE
        u233_conv    = (tot_n_per_yr / AVO * AMU_U233  / 1e3) * BLANKET_COVERAGE

        # -----------------------------------------------------------------------------------------
        # Dataset Configuration: 
        # (dataframe, label, color, marker, markersize, linestyle, reaction_key, conversion_factor)
        # -----------------------------------------------------------------------------------------
        datasets = [
            (self.flibe_th_rr_df, r'FLiBe-ThF$_4$', '#66b420', '+', 12, '--', 'makima', 'Th232(n,g)', u233_conv),
            (hcpb_th_rr_df_corr,  r'HCPB-ThO$_2$',  '#b41f24', '1', 13, '--', 'makima', 'Th232(n,g)', u233_conv),
            (self.dcll_th_rr_df,  r'DCLL-ThO$_2$',  '#0047ba', 'x', 10, '--', 'makima', 'Th232(n,g)', u233_conv),
            (self.flibe_u_rr_df,  r'FLiBe-UF$_4$',  '#66b420', 'o', 5,  '-',  'makima',  'U238(n,g)',  pu239_conv),
            (hcpb_u_rr_df_corr,   r'HCPB-UO$_2$',   '#b41f24', 's', 6,  '-',  'makima',  'U238(n,g)',  pu239_conv),
            (self.dcll_u_rr_df,   r'DCLL-UO$_2$',   '#0047ba', '^', 8,  '-',  'makima',  'U238(n,g)',  pu239_conv),
        ]

        # Select a subset of fertile kg/m3 to SHOW to avoid cluttering the low end with markers
        selected1 = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150,]  # [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 250, 500, 750, 999.99]
        selected2 = [0, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]  # [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 250, 500, 750, 999.99]

        # Set view limits for each of the plots
        view_limits = [ # Zoomed in on (0, 155) fertile kg/m3
                        {'xlim': (-2.5, 152.5),   'xmax': 150, 'x_major': 15,  'x_minor': 5,  
                         'ylim': (-15, 565), 'y_major': 50, 'y_minor':25 ,
                         'suffix': '0150kgm3', 'fertile_kgm3': selected1,
                         'leg_loc': 'upper left'},
                        # Zoomed out to (0, 1000) fertile kg/m3
                        {'xlim': (-25, 1025), 'xmax':1000, 'x_major': 100, 'x_minor': 50, 
                         'ylim': (-50, 1850), 'y_major': 200, 'y_minor':100,
                         'suffix': '1000kgm3', 'fertile_kgm3': selected2,
                         'leg_loc': 'upper left'} ]


        for v in view_limits:

            plt.figure(figsize=(8, 6))
            ax = plt.gca()
            x_fine = np.linspace(0, 1000, 500)

            for df, label, color, marker, markersize, linestyle, degree, rxn, conv, in datasets:

                selected = v['fertile_kgm3']
                mask = df['fertile_kg/m3'].apply(lambda d: any(np.isclose(d, s) for s in selected))
                df_filtered = df[mask]

                x_filtered = df_filtered['fertile_kg/m3']
                y_filtered = df_filtered[rxn] * conv

                # Plot filtered data
                plt.scatter(x_filtered, y_filtered, marker=marker, s=markersize*5, color=color, zorder=3)

                # But still fit the polynomial to all of the data you did compute
                x_fit = df['fertile_kg/m3']
                y_fit = df[rxn] * conv

                if degree == 'makima':
                    x_fine = np.linspace(0, v['xmax'], 500)
                    y_fine    = Akima1DInterpolator(x_fit, y_fit, method='makima')(x_fine) # poly_func(x_fine)

                else:
                    coeffs    = np.polyfit(x_fit, y_fit, degree)
                    y_fine    = np.poly1d(coeffs)(x_fine)

                # Plot interpolation
                plt.plot(x_fine, y_fine, linestyle, linewidth=1, color=color)

                # Dummy plots for legend -- bit of a hack lmao -- ppark
                plt.plot([9e8,9e9], [9e8,9e9], marker+linestyle, markersize=markersize, linewidth=1, color=color, label=label)

            # Labeling and Titles
            # plt.title(r'Fissile prod rates up to 1000 kg$/$m$^3$ (2025-10-28)')
            plt.xlabel(r'Fertile isotope density in breeder [kg$/$m$^3$]')
            plt.ylabel('Initial fissile production rate [kg$/$yr]')

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
            leg = plt.legend(loc=v['leg_loc'], fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=2)
            leg.get_frame().set_linewidth(0.5) 

            if self.save:
                plt.savefig(f'./Figures/PDF/fig_fpr_{v["suffix"]}.pdf', bbox_inches='tight', format='pdf')
                plt.savefig(f'./Figures/PNG/fig_fpr_{v["suffix"]}.png', bbox_inches='tight', format='png')
                print("Exported fissile production per year plot for all blankets.")
            else:
                print("Did not export fissile production per year plot due to user setting.")

            if self.show: 
                plt.show()
                
            plt.close('all')



    def plot_cum_norm_histogram(self):
        """
        Plots cumulative, normalized Pu production vs. energy for contours of MTU.
        """
        print(f"\nPlotting cumulative, normalized fissile production vs. energy...")

        fig, axes = plt.subplots(3, 2, figsize=(15,15)) # sharex='col', sharey='row', 

        titles = [r"FLiBe-UF$_4$", r"FLiBe-ThF$_4$", r"PB-UO$_2$", r"PB-ThO$_2$", r"LL-UO$_2$", r"LL-ThO$_2$",]
        dfs    = [self.flibe_u_ng_df, self.flibe_th_ng_df, self.hcpb_u_ng_df, self.hcpb_th_ng_df, self.dcll_u_ng_df, self.dcll_th_ng_df,]

        for ax, df, title in zip(axes.flatten(), dfs, titles):

            # Filter out MT_fertile loadings we want to plot 
            df = df[df["fertile_kg/m3"].isin([30,60,90,120,150])]

            # Compute sum of 'mean' for each MT_fertile. 
            df_mean = df.groupby("fertile_kg/m3")["mean"].sum().to_frame() # pd.DataFrame(, columns=["MT_fertile","sum"]) # 2-col df, colA = MT_fertile values, colB = sum of 'mean'
            df_mean.columns = ["sum"]   
            df = df.merge(df_mean, on="fertile_kg/m3", how="left") # adds sum as column to df

            df['norm_mean'] = df['mean'] / df['sum']

            # df.to_csv('test.csv',index=False)
            

            # Create bin edges from low and high energy boundaries
            edges = np.unique(np.concatenate([df['energy low [eV]'].values, df['energy high [eV]'].values]))
            bins  = np.sort(df['energy mid [eV]'].unique())

            sub   = df[df['fertile_kg/m3'] == 30]
            label = int(round(sub['fertile_kg/m3'].iloc[0], -1))
            ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#ff1f5b', label=fr'{label} kg$/$m$^3$')

            sub = df[df['fertile_kg/m3'] == 60]
            label = int(round(sub['fertile_kg/m3'].iloc[0], -1))
            ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#f48628', label=fr'{label} kg$/$m$^3$')

            sub = df[df['fertile_kg/m3'] == 90]
            label = int(round(sub['fertile_kg/m3'].iloc[0], -1))
            ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#04cc6c', label=fr'{label} kg$/$m$^3$')

            sub = df[df['fertile_kg/m3'] == 120]
            label = int(round(sub['fertile_kg/m3'].iloc[0], -1))
            ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#0c9edd', label=fr'{label} kg$/$m$^3$')

            sub = df[df['fertile_kg/m3'] == 150]
            label = int(round(sub['fertile_kg/m3'].iloc[0], -1))
            ax.hist(sub['energy mid [eV]'], bins=bins, weights=sub['norm_mean'], cumulative=True, histtype='step', color='#b367bd', label=fr'{label} kg$/$m$^3$')
            

            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.yaxis.set_minor_locator(MultipleLocator(0.05))
            ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
            ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

            ax.set_xscale('log')
            ax.set_xlim(0.67*1e0, 1.5*1e7)
            ax.set_ylim(-0.03, 1.03)
            fig.tight_layout()

            leg = ax.legend(title=title, fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=1, loc="lower right")
            leg.get_frame().set_linewidth(0.5) 

        for i in range(3):
            for j in range(2):
                ax = axes[i, j]
                if i == 2:
                    ax.set_xlabel("Incident neutron energy [eV]")
                if j == 0:
                    ax.set_ylabel(r"Cumulative fraction of fertile $($n,$\gamma)$ rxns")

        if self.save:
            plt.savefig(f'./Figures/PDF/fig_cum_norm_histogram.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/PNG/fig_cum_norm_histogram.png', bbox_inches='tight', format='png')
            print(f"Exported cumulative normalized histogram plots.")
        else:
            print(f"Did not export cumulative normalized histogram plot.")

        if self.show: plt.show()
        plt.close('all')


    def plot_dfisdfer(self):
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

        plt.scatter(x6, y6eval, marker='s', s=40, color='#b41f24')   # HCPB-UO2
        plt.scatter(x5, y5eval, marker='1', s=70, color='#b41f24')   # HCPB-ThO2
        plt.scatter(x4, y4eval, marker='^', s=40, color='#0047ba')   # DCLL-UO2
        plt.scatter(x3, y3eval, marker='x', s=50, color='#0047ba')   # DCLL-ThO2
        plt.scatter(x2, y2eval, marker='o', s=30, color='#66b420')   # FLiBe-UF4
        plt.scatter(x1, y1eval, marker='+', s=60, color='#66b420')   # FLiBe-ThF4

        plt.plot(x_fine, y6fine, '-',   linewidth=1, color='#b41f24',)   #  HCPB-UO2
        plt.plot(x_fine, y5fine, '--',  linewidth=1, color='#b41f24',)   #  HCPB-ThO2
        plt.plot(x_fine, y4fine, '-',   linewidth=1, color='#0047ba',)   #  DCLL-UO2
        plt.plot(x_fine, y3fine, '--',  linewidth=1, color='#0047ba',)   #  DCLL-ThO2
        plt.plot(x_fine, y2fine, '-',   linewidth=1, color='#66b420', )  #  FLiBe-UF4
        plt.plot(x_fine, y1fine, '--',  linewidth=1, color='#66b420', )  #  FLiBe-ThF4 

        # Dummy plots for legend -- ppark
        plt.plot([9e8,9e9], [9e8,9e9], '+--', markersize=8, linewidth=1, color='#66b420', label=r'FLiBe-ThF$_4$')   #  green: #66b420
        plt.plot([9e8,9e9], [9e8,9e9], 'o-',  markersize=5, linewidth=1, color='#66b420', label=r'FLiBe-UF$_4$')    #  green: #66b420
        plt.plot([9e8,9e9], [9e8,9e9], '1--', markersize=9, linewidth=1, color='#b41f24', label=r'HCPB-ThO$_2$')    #  red:   #b41f24
        plt.plot([9e8,9e9], [9e8,9e9], 's-',  markersize=5, linewidth=1, color='#b41f24', label=r'HCPB-UO$_2$')     #  red:   #b41f24
        plt.plot([9e8,9e9], [9e8,9e9], 'x--', markersize=7, linewidth=1, color='#0047ba', label=r'DCLL-ThO$_2$')    #  blue:  #0047ba
        plt.plot([9e8,9e9], [9e8,9e9], '^-',  markersize=6, linewidth=1, color='#0047ba', label=r'DCLL-UO$_2$')     #  blue:  #0047ba

        # plt.title(f'Tritium breeding ratio (All)') # Exclude title for production figs --ppark 2025-08-06
        plt.xlabel(r'Fertile isotope density $n$ in blanket [kg$/$m$^3$]')
        plt.ylabel(r'$dR_{\text{fis}}/dn$ [kg$/\text{yr}^2$]')

        plt.xlim(-5,160)
        plt.ylim(-0.25,6)  # plt.ylim(-7.5,225+7.5) REMEMBER TO CHANGE 

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

        leg = plt.legend(fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=1)
        leg.get_frame().set_linewidth(0.5) 
    
        if self.save:
            plt.savefig(f'./Figures/PDF/fig_dRdn.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/PNG/fig_dRdn.png', bbox_inches='tight', format='png')
            print("Exported rate of change for fissile material production with respect to fertile material plot for all blankets.")
        else:
            print("Did not export dR / dn plot due to user setting.")

        if self.show: plt.show()
        plt.close('all')


    def plot_fisovern(self):
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
        plt.figure(figsize=(7.5,5))

        plt.scatter(x6, y6/x6, marker='s', s=40, color='#b41f24') # HCPB-UO2
        plt.scatter(x5, y5/x5, marker='1', s=70, color='#b41f24') # HCPB-ThO2
        plt.scatter(x4, y4/x4, marker='^', s=40, color='#0047ba') # DCLL-UO2
        plt.scatter(x3, y3/x3, marker='x', s=50, color='#0047ba') # DCLL-ThO2
        plt.scatter(x2, y2/x2, marker='o', s=30, color='#66b420')   # FLiBe-UF4
        plt.scatter(x1, y1/x1, marker='+', s=60, color='#66b420')   # FLiBe-ThF4

        plt.plot(x_fine, y_akima6, '-',   linewidth=1, color='#b41f24',)   #  HCPB-UO2
        plt.plot(x_fine, y_akima5, '--',  linewidth=1, color='#b41f24',)   #  HCPB-ThO2
        plt.plot(x_fine, y_akima4, '-',   linewidth=1, color='#0047ba',)   #  DCLL-UO2
        plt.plot(x_fine, y_akima3, '--',  linewidth=1, color='#0047ba',)   #  DCLL-ThO2
        plt.plot(x_fine, y_akima2, '-',   linewidth=1, color='#66b420', )    #  FLiBe-UF4
        plt.plot(x_fine, y_akima1, '--',  linewidth=1, color='#66b420', )    #  FLiBe-ThF4 

        # Dummy plots for legend -- ppark
        plt.plot([9e8,9e9], [9e8,9e9], '+--', markersize=8, linewidth=1, color='#66b420',   label=r'FLiBe-ThF$_4$') #  green: #66b420
        plt.plot([9e8,9e9], [9e8,9e9], 'o-',  markersize=5, linewidth=1, color='#66b420',   label=r'FLiBe-UF$_4$')  # 
        plt.plot([9e8,9e9], [9e8,9e9], '1--', markersize=9, linewidth=1, color='#b41f24', label=r'HCPB-ThO$_2$')    # 
        plt.plot([9e8,9e9], [9e8,9e9], 's-',  markersize=5, linewidth=1, color='#b41f24', label=r'HCPB-UO$_2$')     # 
        plt.plot([9e8,9e9], [9e8,9e9], 'x--', markersize=7, linewidth=1, color='#0047ba', label=r'DCLL-ThO$_2$')    #  red: #b41f24
        plt.plot([9e8,9e9], [9e8,9e9], '^-',  markersize=6, linewidth=1, color='#0047ba', label=r'DCLL-UO$_2$')     #  blue: #0047ba

        # plt.title(f'Tritium breeding ratio (All)') # Exclude title for production figs --ppark 2025-08-06
        plt.xlabel(r'Fertile isotope density $n$ in blanket [kg$/$m$^3$]')
        plt.ylabel(r'$R_{\text{fis}}/n$ [kg$/$yr]')

        plt.xlim(-5,160)
        plt.ylim(-0.25,6) 

        # Tick grid
        ax = plt.gca()
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
 
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(MultipleLocator(1)) 
        ax.yaxis.set_minor_locator(MultipleLocator(0.25)) 
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

        plt.tight_layout()

        leg = plt.legend(fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=1)
        leg.get_frame().set_linewidth(0.5) 
    
        if self.save:
            plt.savefig(f'./Figures/PDF/fig_Rn.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/PNG/fig_Rn.png', bbox_inches='tight', format='png')
            print("Exported fissile material production divided by fertile material for all blankets.")
        else:
            print("Did not export R/n plot due to user setting.")

        if self.show: plt.show()
        plt.close('all')




    
if __name__ == "__main__":
    plot_all()