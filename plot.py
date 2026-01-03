import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter, FuncFormatter
from scipy.interpolate import Akima1DInterpolator
from scipy.interpolate import PchipInterpolator


""" import helper functions """
from Python.utilities import *
from Python.parameters import *



def plot_all():
    """ Read CSV data into pandas DataFrames """

    combined_plot = Plot(save=True, show=True)
    
    """ Pick which one to plot by uncommenting """
    # combined_plot.plot_tbr()
    # combined_plot.plot_pu_per_yr()
    # combined_plot.plot_cum_norm_histogram()
    combined_plot.plot_flux()

    print("All plots completed and saved.")


class Plot:

    def __init__(self, save=False, show=True):

        # Dataframes of reaction rates (rr)
        self.flibe_u_rr_df  = pd.read_csv('./Figures/Data/FLiBe_U238_rxns_900K.csv')
        self.flibe_th_rr_df = pd.read_csv('./Figures/Data/FLiBe_Th232_rxns_900K.csv')
        self.dcll_u_rr_df   = pd.read_csv('./Figures/Data/DCLL_U238_rxns_900K.csv')
        self.dcll_th_rr_df  = pd.read_csv('./Figures/Data/DCLL_Th232_rxns_900K.csv')
        self.hcpb_u_rr_df   = pd.read_csv('./Figures/Data/HCPB_U238_rxns_900K.csv')
        self.hcpb_th_rr_df  = pd.read_csv('./Figures/Data/HCPB_Th232_rxns_900K.csv')

        self.flibe_u_ng_df  = pd.read_csv('./Figures/Data/FLiBe_U238_n-gamma_900K.csv')
        self.flibe_th_ng_df = pd.read_csv('./Figures/Data/FLiBe_Th232_n-gamma_900K.csv')
        self.dcll_u_ng_df   = pd.read_csv('./Figures/Data/DCLL_U238_n-gamma_900K.csv')
        self.dcll_th_ng_df  = pd.read_csv('./Figures/Data/DCLL_Th232_n-gamma_900K.csv')
        self.hcpb_u_ng_df   = pd.read_csv('./Figures/data/HCPB_U238_n-gamma_900K.csv') 
        self.hcpb_th_ng_df  = pd.read_csv('./Figures/data/HCPB_Th232_n-gamma_900K.csv') 

        self.flibe_u_flux_df  = pd.read_csv('./Figures/Data/FLiBe_U238_flux_900K.csv')
        self.flibe_th_flux_df = pd.read_csv('./Figures/Data/FLiBe_Th232_flux_900K.csv')
        self.dcll_u_flux_df   = pd.read_csv('./Figures/Data/DCLL_U238_flux_900K.csv')
        self.dcll_th_flux_df  = pd.read_csv('./Figures/Data/DCLL_Th232_flux_900K.csv')

        self.save, self.show = save, show
        # self.name = 'All_Blankets'

        for sd in ['PDF','PNG','Data']:
            sd_path = f'./Figures/{sd}/'
            print(f"Ensuring directory exists: {sd_path}")
            os.makedirs(sd_path, exist_ok=True)


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
        print(f"\nPlotting tritium breeding ratio vs. fertile density...")

        # Setting up x, y separately here so you can remove the impurity/wppm-magnitude cases --ppark 2025-08-06
        x6, y6 =   self.hcpb_u_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * (  self.hcpb_u_rr_df['Li6(n,t)'] +   self.hcpb_u_rr_df['Li7(n,Xt)'] )
        x5, y5 =  self.hcpb_th_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * ( self.hcpb_th_rr_df['Li6(n,t)'] +  self.hcpb_th_rr_df['Li7(n,Xt)'] )
        x4, y4 =   self.dcll_u_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * (  self.dcll_u_rr_df['Li6(n,t)'] +   self.dcll_u_rr_df['Li7(n,Xt)'] )
        x3, y3 =  self.dcll_th_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * ( self.dcll_th_rr_df['Li6(n,t)'] +  self.dcll_th_rr_df['Li7(n,Xt)'] )
        x2, y2 =  self.flibe_u_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * ( self.flibe_u_rr_df['Li6(n,t)'] +  self.flibe_u_rr_df['Li7(n,Xt)'] )
        x1, y1 = self.flibe_th_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * (self.flibe_th_rr_df['Li6(n,t)'] + self.flibe_th_rr_df['Li7(n,Xt)'] )

        # AkimaInterpolation ripped from my ELWR paper --ppark 2025-08-06
        x_fine = np.linspace(x1.min(), x1.max(), 300) # Evaluate on a fine grid
        y_akima1 = Akima1DInterpolator(x1, y1)(x_fine)
        y_akima2 = Akima1DInterpolator(x2, y2)(x_fine)
        y_akima3 = Akima1DInterpolator(x3, y3)(x_fine)
        y_akima4 = Akima1DInterpolator(x4, y4)(x_fine)
        y_akima5 = Akima1DInterpolator(x5, y5)(x_fine)
        y_akima6 = Akima1DInterpolator(x6, y6)(x_fine)

        for y in [y_akima1,y_akima2,y_akima3,y_akima4,y_akima5,y_akima6]: # y_akima5,
            for i in range(1, len(y)):
                if y[i] > y[i-1]:
                    y[i] = y[i-1] # adjust the current point to ensure it's not greater than the previous point


        plt.figure(figsize=(7.5,5))

        # Tick grid
        ax = plt.gca()

        # ax.axhspan(.9,1, color='#e0ded8') # color='#e0ded8', # alpha=0.3)
        ax.axhspan(0,1.00, color='#e0ded8') # color='#f7f7f6' # , alpha=0.3)

        plt.scatter(x6, y6, marker='s', s=40, color='#b41f24') #  PB-UO2
        plt.scatter(x5, y5, marker='1', s=70, color='#b41f24') #  PB-ThO2
        plt.scatter(x4, y4, marker='^', s=40, color='#0047ba') #  LL-UO2
        plt.scatter(x3, y3, marker='x', s=50, color='#0047ba') #  LL-ThO2
        plt.scatter(x2, y2, marker='o', s=30, color='#66b420')   #  FLiBe-UF4
        plt.scatter(x1, y1, marker='+', s=60, color='#66b420')   #  FLiBe-ThF4 

        plt.plot(x_fine, y_akima6, '-',   linewidth=1, color='#b41f24',)   #  PB-UO2
        plt.plot(x_fine, y_akima5, '--',  linewidth=1, color='#b41f24',)   #  PB-ThO2
        plt.plot(x_fine, y_akima4, '-',   linewidth=1, color='#0047ba',)   #  LL-UO2
        plt.plot(x_fine, y_akima3, '--',  linewidth=1, color='#0047ba',)   #  LL-ThO2
        plt.plot(x_fine, y_akima2, '-',   linewidth=1, color='#66b420', )    #  FLiBe-UF4
        plt.plot(x_fine, y_akima1, '--',  linewidth=1, color='#66b420', )    #  FLiBe-ThF4 

        # Dummy plots for legend -- bit of a hack lmao -- ppark
        plt.plot([9e8,9e9], [9e8,9e9], '^-',  markersize=6, linewidth=1, color='#0047ba', label=r'LL-UO$_2$')     #  blue: #0047ba
        plt.plot([9e8,9e9], [9e8,9e9], 'x--', markersize=7, linewidth=1, color='#0047ba', label=r'LL-ThO$_2$')    #  
        plt.plot([9e8,9e9], [9e8,9e9], 's-',  markersize=5, linewidth=1, color='#b41f24', label=r'PB-UO$_2$')     #  red: #b41f24
        plt.plot([9e8,9e9], [9e8,9e9], '1--', markersize=9, linewidth=1, color='#b41f24', label=r'PB-ThO$_2$')    # 
        plt.plot([9e8,9e9], [9e8,9e9], 'o-',  markersize=5, linewidth=1, color='#66b420',   label=r'FLiBe-UF$_4$')  #  green: #66b420
        plt.plot([9e8,9e9], [9e8,9e9], '+--', markersize=8, linewidth=1, color='#66b420',   label=r'FLiBe-ThF$_4$') #  
        
        
        plt.title(r'TBRs going up to 1000 kg$/$m$^3$ (2025-10-28)') # Exclude title for production figs --ppark 2025-08-06
        plt.xlabel(r'Fertile isotope density in blanket [kg$/$m$^3$]')
        plt.ylabel('Tritium breeding ratio')

        plt.xlim(-25,1025)  # plt.xlim(-5,165)
        plt.ylim(.89,1.16)  # plt.ylim(0.98,1.42)
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_major_locator(MultipleLocator(100))
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

        plt.tight_layout()

        leg = plt.legend(fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=3)
        leg.get_frame().set_linewidth(0.5) 

        if self.save:
            plt.savefig(f'./Figures/pdf/fig_tbr_1000kgm3_ALX.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/fig_tbr_1000kgm3_ALX.png', bbox_inches='tight', format='png')
            print(f"Exported tritium breeding ratio plots.")
        else:
            print(f"Did not export tritium breeding ratio plots.")

        if self.show: plt.show()
        plt.close('all')


    def plot_pu_per_yr(self):
      
        print(f"\nPlotting Pu-239 production per year for all fuels...")

        # <utilities.py> Default NPS_FUS = 500 MJ/s * 3.546e17 n/MJ = 1.773e+20 n/s
        # seconds in a year = 3.156e+7 s/yr

        tot_n_per_yr = NPS_FUS * 3.156e+7

        # Gotta then convert fissile atoms/yr to kg/yr
        pu239_at_to_kg_per_yr = tot_n_per_yr / AVO * AMU_PU239 / 1e3
        u233_at_to_kg_per_yr  = tot_n_per_yr / AVO * AMU_U233  / 1e3
        
        # Setting up x, y separately here so you can remove the impurity/wppm-magnitude cases --ppark 2025-08-06
        x6, y6 =   self.hcpb_u_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * pu239_at_to_kg_per_yr *   self.hcpb_u_rr_df['U238(n,g)'] 
        x5, y5 =  self.hcpb_th_rr_df['fertile_kg/m3'], BLANKET_COVERAGE *  u233_at_to_kg_per_yr *  self.hcpb_th_rr_df['Th232(n,g)']
        x4, y4 =   self.dcll_u_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * pu239_at_to_kg_per_yr *   self.dcll_u_rr_df['U238(n,g)'] 
        x3, y3 =  self.dcll_th_rr_df['fertile_kg/m3'], BLANKET_COVERAGE *  u233_at_to_kg_per_yr *  self.dcll_th_rr_df['Th232(n,g)']
        x2, y2 =  self.flibe_u_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * pu239_at_to_kg_per_yr *  self.flibe_u_rr_df['U238(n,g)']
        x1, y1 = self.flibe_th_rr_df['fertile_kg/m3'], BLANKET_COVERAGE *  u233_at_to_kg_per_yr * self.flibe_th_rr_df['Th232(n,g)']

        # AkimaInterpolation ripped from my ELWR paper --ppark 2025-08-06
        x_fine = np.linspace(x1.min(), x1.max(), 300) # Evaluate on a fine grid
        y_akima1 = Akima1DInterpolator(x1, y1)(x_fine)
        y_akima2 = Akima1DInterpolator(x2, y2)(x_fine)
        y_akima3 = Akima1DInterpolator(x3, y3)(x_fine)
        y_akima4 = Akima1DInterpolator(x4, y4)(x_fine)
        y_akima5 = Akima1DInterpolator(x5, y5)(x_fine)
        y_akima6 = Akima1DInterpolator(x6, y6)(x_fine)

        plt.figure(figsize=(7.5,5))

        plt.scatter(x6, y6, marker='s', s=40, color='#b41f24') # PB-UO2
        plt.scatter(x5, y5, marker='1', s=70, color='#b41f24') # PB-ThO2
        plt.scatter(x4, y4, marker='^', s=40, color='#0047ba') # LL-UO2
        plt.scatter(x3, y3, marker='x', s=50, color='#0047ba') # LL-ThO2
        plt.scatter(x2, y2, marker='o', s=30, color='#66b420')   # FLiBe-UF4
        plt.scatter(x1, y1, marker='+', s=60, color='#66b420')   # FLiBe-ThF4

        plt.plot(x_fine, y_akima6, '-',   linewidth=1, color='#b41f24',)   #  PB-UO2
        plt.plot(x_fine, y_akima5, '--',  linewidth=1, color='#b41f24',)   #  PB-ThO2
        plt.plot(x_fine, y_akima4, '-',   linewidth=1, color='#0047ba',)   #  LL-UO2
        plt.plot(x_fine, y_akima3, '--',  linewidth=1, color='#0047ba',)   #  LL-ThO2
        plt.plot(x_fine, y_akima2, '-',   linewidth=1, color='#66b420', )    #  FLiBe-UF4
        plt.plot(x_fine, y_akima1, '--',  linewidth=1, color='#66b420', )    #  FLiBe-ThF4 

        # Dummy plots for legend -- bit of a hack lmao -- ppark
        plt.plot([9e8,9e9], [9e8,9e9], '+--', markersize=8, linewidth=1, color='#66b420',   label=r'FLiBe-ThF$_4$') #  green: #66b420
        plt.plot([9e8,9e9], [9e8,9e9], 'o-',  markersize=5, linewidth=1, color='#66b420',   label=r'FLiBe-UF$_4$')  # 
        plt.plot([9e8,9e9], [9e8,9e9], '1--', markersize=9, linewidth=1, color='#b41f24', label=r'PB-ThO$_2$')    # 
        plt.plot([9e8,9e9], [9e8,9e9], 's-',  markersize=5, linewidth=1, color='#b41f24', label=r'PB-UO$_2$')     # 
        plt.plot([9e8,9e9], [9e8,9e9], 'x--', markersize=7, linewidth=1, color='#0047ba', label=r'LL-ThO$_2$')    #  red: #b41f24
        plt.plot([9e8,9e9], [9e8,9e9], '^-',  markersize=6, linewidth=1, color='#0047ba', label=r'LL-UO$_2$')     #  blue: #0047ba


        plt.title(r'Fissile prod rates up to 1000 kg$/$m$^3$ (2025-10-28)')  # Exclude title for production figs --ppark 2025-08-06
        plt.xlabel(r'Fertile isotope density in blanket [kg$/$m$^3$]')
        plt.ylabel('Initial fissile production rate [kg$/$yr]')

        plt.xlim(-25,1025) # plt.xlim(-5,165)
        plt.ylim(-25,1425) # plt.ylim(-15,2*(250+7.5))  # plt.ylim(-7.5,225+7.5) 

        # Tick grid
        ax = plt.gca()
        ax.xaxis.set_ticks_position('both')
        ax.xaxis.set_minor_locator(MultipleLocator(100)) # ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
 
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(MultipleLocator(200)) # ax.yaxis.set_major_locator(MultipleLocator(50))   # ax.yaxis.set_major_locator(MultipleLocator(25))
        ax.yaxis.set_minor_locator(MultipleLocator(50)) # ax.yaxis.set_minor_locator(MultipleLocator(25))   # ax.yaxis.set_minor_locator(MultipleLocator(12.5))
        ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

        plt.tight_layout()

        leg = plt.legend(fancybox=False, edgecolor='black', frameon=True, framealpha=.75, ncol=1)
        leg.get_frame().set_linewidth(0.5) 
    
        if self.save:
            plt.savefig(f'./Figures/pdf/fig_fissile_per_yr_1000kgm3.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/fig_fissile_per_yr_1000kgm3.png', bbox_inches='tight', format='png')
            print("Exported fissile production per year plot for all blankets.")
        else:
            print("Did not export fissile production per year plot due to user setting.")

        if self.show: plt.show()
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
            plt.savefig(f'./Figures/pdf/fig_cum_norm_histogram.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/png/fig_cum_norm_histogram.png', bbox_inches='tight', format='png')
            print(f"Exported cumulative normalized histogram plots.")
        else:
            print(f"Did not export cumulative normalized histogram plot.")

        if self.show: plt.show()
        plt.close('all')




    
if __name__ == "__main__":
    plot_all()