import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from scipy.interpolate import Akima1DInterpolator
from scipy.interpolate import PchipInterpolator


""" import helper functions """
from Python.utilities import *
from Python.parameters import *



def plot_all():
    """ Read CSV data into pandas DataFrames """

    combined_plot = Plot(save=False, show=True)
    
    """ Pick which one to plot by uncommenting """
    combined_plot.plot_tbr()
    # combined_plot.plot_pu_per_yr()
    # combined_plot.plot_cum_norm_histogram()
    combined_plot.plot_pb_coupon()

    print("All plots completed and saved.")


class Plot:

    def __init__(self, save=False, show=True):

        # Dataframes of reaction rates (rr)
        self.flibe_u_rr_df  = pd.read_csv('./Figures/Data/FLiBe_U_rxns_900K.csv')
        self.flibe_th_rr_df = pd.read_csv('./Figures/Data/FLiBe_Th_rxns_900K.csv')
        self.pbli_u_rr_df   = pd.read_csv('./Figures/Data/LL_U_rxns_900K.csv')
        self.pbli_th_rr_df  = pd.read_csv('./Figures/Data/LL_Th_rxns_900K.csv')
        self.pb_u_rr_df     = pd.read_csv('./Figures/Data/PB_U_rxns_900K.csv')
        self.pb_th_rr_df    = pd.read_csv('./Figures/Data/PB_Th_rxns_900K.csv')
        self.pb_het_u_rr_df = pd.read_csv('./Figures/Data/PBHet_U_rxns_900K.csv')
        self.pb_het_th_rr_df = pd.read_csv('./Figures/Data/PBHet_Th_rxns_900K.csv')
        self.pb_hom_u_rr_df = pd.read_csv('./Figures/Data/PBHom_U_rxns_900K.csv')
        self.pb_hom_th_rr_df = pd.read_csv('./Figures/Data/PBHom_Th_rxns_900K.csv')

        self.flibe_u_eb_df  = pd.read_csv('./Figures/Data/FLiBe_U_n-gamma_900K.csv')
        self.flibe_th_eb_df = pd.read_csv('./Figures/Data/FLiBe_Th_n-gamma_900K.csv')
        self.pbli_u_eb_df   = pd.read_csv('./Figures/Data/LL_U_n-gamma_900K.csv')
        self.pbli_th_eb_df  = pd.read_csv('./Figures/Data/LL_Th_n-gamma_900K.csv')
        self.pb_u_eb_df     = pd.read_csv('./Figures/data/PB_U_n-gamma_900K.csv') 
        self.pb_th_eb_df    = pd.read_csv('./Figures/data/PB_Th_n-gamma_900K.csv') 
        self.pb_het_u_eb_df = pd.read_csv('./Figures/Data/PBHet_U_n-gamma_900K.csv')
        self.pb_het_th_eb_df = pd.read_csv('./Figures/Data/PBHet_Th_n-gamma_900K.csv')
        self.pb_hom_u_eb_df = pd.read_csv('./Figures/Data/PBHom_U_n-gamma_900K.csv')
        self.pb_hom_th_eb_df = pd.read_csv('./Figures/Data/PBHom_Th_n-gamma_900K.csv')

        self.save, self.show = save, show
        # self.name = 'All_Blankets'

        for sd in ['PDF','PNG','Data']:
            sd_path = f'./Figures/{sd}/'
            print(f"Ensuring directory exists: {sd_path}")
            os.makedirs(sd_path, exist_ok=True)

        
    def plot_tbr(self):
        print(f"\nPlotting tritium breeding ratio vs. fertile density...")

        # Setting up x, y separately here so you can remove the impurity/wppm-magnitude cases --ppark 2025-08-06
        x6, y6 =     self.pb_u_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * (    self.pb_u_rr_df['Li6(n,t)'] +     self.pb_u_rr_df['Li7(n,Xt)'] )
        x5, y5 =    self.pb_th_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * (   self.pb_th_rr_df['Li6(n,t)'] +    self.pb_th_rr_df['Li7(n,Xt)'] )
        x4, y4 =   self.pbli_u_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * (  self.pbli_u_rr_df['Li6(n,t)'] +   self.pbli_u_rr_df['Li7(n,Xt)'] )
        x3, y3 =  self.pbli_th_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * ( self.pbli_th_rr_df['Li6(n,t)'] +  self.pbli_th_rr_df['Li7(n,Xt)'] )
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
        x6, y6 =     self.pb_u_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * pu239_at_to_kg_per_yr *     self.pb_u_rr_df['U238(n,g)'] 
        x5, y5 =    self.pb_th_rr_df['fertile_kg/m3'], BLANKET_COVERAGE *  u233_at_to_kg_per_yr *    self.pb_th_rr_df['Th232(n,g)']
        x4, y4 =   self.pbli_u_rr_df['fertile_kg/m3'], BLANKET_COVERAGE * pu239_at_to_kg_per_yr *   self.pbli_u_rr_df['U238(n,g)'] 
        x3, y3 =  self.pbli_th_rr_df['fertile_kg/m3'], BLANKET_COVERAGE *  u233_at_to_kg_per_yr *  self.pbli_th_rr_df['Th232(n,g)']
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
        dfs    = [self.flibe_u_eb_df, self.flibe_th_eb_df, self.pb_u_eb_df, self.pb_th_eb_df, self.pbli_u_eb_df, self.pbli_th_eb_df,]

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


    def plot_pb_coupon(self):
        """
        Plots TBR and kg/yr comparisons for PBHet (pebble coupon heterogeneous) vs PBHom (pebble coupon homogeneous)
        """
        print(f"\nPlotting PB coupon comparison: PBHet vs PBHom...")
        
        # Plot material comparison table
        self.plot_pb_coupon_material_comparison()
        
        # Plot TBR comparison
        self.plot_pb_coupon_tbr()
        
        # Plot kg/yr comparison
        self.plot_pb_coupon_kgyr()


    def plot_pb_coupon_tbr(self):
        """
        Plots TBR comparison for PBHet vs PBHom
        """
        print(f"\nPlotting TBR comparison: PBHet vs PBHom...")

        # Create figure with subplots for U and Th
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))

        for idx, (fertile_element, ax) in enumerate(zip(['U', 'Th'], axes)):
            # Prepare data for PBHet
            if fertile_element == 'U':
                pbhet_df = self.pb_het_u_rr_df
            else:
                pbhet_df = self.pb_het_th_rr_df

            # Prepare data for PBHom
            if fertile_element == 'U':
                pbhom_df = self.pb_hom_u_rr_df
            else:
                pbhom_df = self.pb_hom_th_rr_df

            # Calculate TBR for PBHet
            x_het = pbhet_df['fertile_kg/m3']
            y_het = BLANKET_COVERAGE * (pbhet_df['Li6(n,t)'] + pbhet_df['Li7(n,t)'])

            # Calculate TBR for PBHom
            x_hom = pbhom_df['fertile_kg/m3']
            y_hom = BLANKET_COVERAGE * (pbhom_df['Li6(n,t)'] + pbhom_df['Li7(n,t)'])

            # Plot scatter points
            ax.scatter(x_het, y_het, marker='s', s=50, color='#b41f24', label='PBHet (Heterogeneous)', zorder=3)
            ax.scatter(x_hom, y_hom, marker='o', s=50, color='#0047ba', label='PBHom (Homogeneous)', zorder=3)

            # Interpolation for smooth lines
            if len(x_het) > 3:
                x_fine = np.linspace(x_het.min(), x_het.max(), 300)
                y_akima_het = Akima1DInterpolator(x_het, y_het)(x_fine)
                # Ensure monotonic decrease
                for i in range(1, len(y_akima_het)):
                    if y_akima_het[i] > y_akima_het[i-1]:
                        y_akima_het[i] = y_akima_het[i-1]
                ax.plot(x_fine, y_akima_het, '-', linewidth=2, color='#b41f24', alpha=0.7, zorder=2)

            if len(x_hom) > 3:
                x_fine = np.linspace(x_hom.min(), x_hom.max(), 300)
                y_akima_hom = Akima1DInterpolator(x_hom, y_hom)(x_fine)
                # Ensure monotonic decrease
                for i in range(1, len(y_akima_hom)):
                    if y_akima_hom[i] > y_akima_hom[i-1]:
                        y_akima_hom[i] = y_akima_hom[i-1]
                ax.plot(x_fine, y_akima_hom, '--', linewidth=2, color='#0047ba', alpha=0.7, zorder=2)

            # Formatting
            ax.axhspan(0, 1.00, color='#e0ded8', alpha=0.3)
            ax.set_xlabel(r'Fertile isotope density in blanket [kg$/$m$^3$]')
            ax.set_ylabel('Tritium breeding ratio')
            ax.set_title(fr'PBHet vs PBHom - {fertile_element}O$_2$')
            
            # Set limits based on data
            all_x = list(x_het.values) + list(x_hom.values)
            all_y = list(y_het.values) + list(y_hom.values)
            
            ax.set_xlim(-25, max(1025, max(all_x) + 25))
            ax.set_ylim(0.89, max(1.16, max(all_y) + 0.05))

            ax.xaxis.set_ticks_position('both')
            ax.xaxis.set_major_locator(MultipleLocator(100))
            ax.xaxis.set_minor_locator(MultipleLocator(50))
            ax.yaxis.set_ticks_position('both')
            ax.yaxis.set_major_locator(MultipleLocator(0.05))
            ax.yaxis.set_minor_locator(MultipleLocator(0.01))
            ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
            ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

            leg = ax.legend(fancybox=False, edgecolor='black', frameon=True, framealpha=.75, loc='best')
            leg.get_frame().set_linewidth(0.5)

        plt.tight_layout()

        if self.save:
            plt.savefig(f'./Figures/PDF/fig_pb_coupon_tbr_comparison.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/PNG/fig_pb_coupon_tbr_comparison.png', bbox_inches='tight', format='png')
            print(f"Exported PB coupon TBR comparison plots.")
        else:
            print(f"Did not export PB coupon TBR comparison plots.")

        if self.show: plt.show()
        plt.close('all')


    def plot_pb_coupon_kgyr(self):
        """
        Plots kg/yr comparison for PBHet vs PBHom
        """
        print(f"\nPlotting kg/yr comparison: PBHet vs PBHom...")

        # Calculate conversion factors (same as in plot_pu_per_yr)
        tot_n_per_yr = NPS_FUS * 3.156e+7
        pu239_at_to_kg_per_yr = tot_n_per_yr / AVO * AMU_PU239 / 1e3
        u233_at_to_kg_per_yr = tot_n_per_yr / AVO * AMU_U233 / 1e3

        # Create figure with subplots for U and Th
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))

        for idx, (fertile_element, ax) in enumerate(zip(['U', 'Th'], axes)):
            fertile_isotope = 'U238' if fertile_element == 'U' else 'Th232'
            conversion_factor = pu239_at_to_kg_per_yr if fertile_element == 'U' else u233_at_to_kg_per_yr
            reaction_col = f'{fertile_isotope}(n,g)'
            
            # Prepare data for PBHet
            if fertile_element == 'U':
                pbhet_df = self.pb_het_u_rr_df
            else:
                pbhet_df = self.pb_het_th_rr_df

            # Prepare data for PBHom
            if fertile_element == 'U':
                pbhom_df = self.pb_hom_u_rr_df
            else:
                pbhom_df = self.pb_hom_th_rr_df

            # Calculate kg/yr for PBHet
            x_het = pbhet_df['fertile_kg/m3']
            y_het = BLANKET_COVERAGE * conversion_factor * pbhet_df[reaction_col]

            # Calculate kg/yr for PBHom
            x_hom = pbhom_df['fertile_kg/m3']
            y_hom = BLANKET_COVERAGE * conversion_factor * pbhom_df[reaction_col]

            # Plot scatter points
            ax.scatter(x_het, y_het, marker='s', s=50, color='#b41f24', label='PBHet (Heterogeneous)', zorder=3)
            ax.scatter(x_hom, y_hom, marker='o', s=50, color='#0047ba', label='PBHom (Homogeneous)', zorder=3)

            # Interpolation for smooth lines
            if len(x_het) > 3:
                x_fine = np.linspace(x_het.min(), x_het.max(), 300)
                y_akima_het = Akima1DInterpolator(x_het, y_het)(x_fine)
                ax.plot(x_fine, y_akima_het, '-', linewidth=2, color='#b41f24', alpha=0.7, zorder=2)

            if len(x_hom) > 3:
                x_fine = np.linspace(x_hom.min(), x_hom.max(), 300)
                y_akima_hom = Akima1DInterpolator(x_hom, y_hom)(x_fine)
                ax.plot(x_fine, y_akima_hom, '--', linewidth=2, color='#0047ba', alpha=0.7, zorder=2)

            # Formatting
            ax.set_xlabel(r'Fertile isotope density in blanket [kg$/$m$^3$]')
            ax.set_ylabel('Initial fissile production rate [kg$/$yr]')
            ax.set_title(fr'PBHet vs PBHom - {fertile_element}O$_2$')
            
            # Set limits based on data
            all_x = list(x_het.values) + list(x_hom.values)
            all_y = list(y_het.values) + list(y_hom.values)
            
            ax.set_xlim(-25, max(1025, max(all_x) + 25))
            ax.set_ylim(-25, max(1425, max(all_y) + 25))

            ax.xaxis.set_ticks_position('both')
            ax.xaxis.set_major_locator(MultipleLocator(100))
            ax.xaxis.set_minor_locator(MultipleLocator(50))
            ax.yaxis.set_ticks_position('both')
            ax.yaxis.set_major_locator(MultipleLocator(200))
            ax.yaxis.set_minor_locator(MultipleLocator(50))
            ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
            ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

            leg = ax.legend(fancybox=False, edgecolor='black', frameon=True, framealpha=.75, loc='best')
            leg.get_frame().set_linewidth(0.5)

        plt.tight_layout()

        if self.save:
            plt.savefig(f'./Figures/PDF/fig_pb_coupon_kgyr_comparison.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/PNG/fig_pb_coupon_kgyr_comparison.png', bbox_inches='tight', format='png')
            print(f"Exported PB coupon kg/yr comparison plots.")
        else:
            print(f"Did not export PB coupon kg/yr comparison plots.")

        if self.show: plt.show()
        plt.close('all')


    def plot_pb_coupon_material_comparison(self):
        """
        Creates a comparison table/chart showing material quantities for PBHet vs PBHom
        to verify both models have the same amount of material
        """
        print(f"\nCreating material comparison chart: PBHet vs PBHom...")

        # Get dataframes
        pbhet_u_df = self.pb_het_u_rr_df
        pbhet_th_df = self.pb_het_th_rr_df
        pbhom_u_df = self.pb_hom_u_rr_df
        pbhom_th_df = self.pb_hom_th_rr_df

        # Blanket volumes (same for both models)
        V_in_m3 = 151.6878656 * 1e-6  # m³
        V_out_m3 = 368.4721344 * 1e-6  # m³
        V_total_m3 = V_in_m3 + V_out_m3
        
        # BISO particle parameters
        V_k_cm3 = (4.0/3.0) * np.pi * (BISO_KERNEL_RADIUS**3)
        V_p_cm3 = (4.0/3.0) * np.pi * (BISO_RADIUS**3)
        
        # Prepare data for comparison
        comparison_data = []
        
        for fertile_element in ['U', 'Th']:
            if fertile_element == 'U':
                pbhet_df = pbhet_u_df
                pbhom_df = pbhom_u_df
                rho_k = DENSITY_UO2  # g/cm³
                fert_enrich = ENRICH_U * 0.01
                M_U = (1-fert_enrich)*AMU_U238 + fert_enrich*AMU_U235
                M_UO2 = M_U + 2.0*AMU_O
                w_fertile_in_oxide = M_U / M_UO2
                AMU_fertile = M_U
                AMU_oxide = M_UO2
            else:
                pbhet_df = pbhet_th_df
                pbhom_df = pbhom_th_df
                rho_k = DENSITY_ThO2  # g/cm³
                M_ThO2 = AMU_Th + 2.0*AMU_O
                w_fertile_in_oxide = AMU_Th / M_ThO2
                AMU_fertile = AMU_Th
                AMU_oxide = M_ThO2
            
            # Calculate mass per particle for PBHet
            m_k_g = rho_k * V_k_cm3
            m_fertile_per_particle_g = m_k_g * w_fertile_in_oxide
            
            for idx, row in pbhet_df.iterrows():
                fertile_density_kgm3 = row['fertile_kg/m3']
                
                if fertile_density_kgm3 == 0:
                    continue
                
                # PBHet calculations
                # Note: PBHet uses total blanket cuboid volume (V_in + V_out)
                m_in_g = fertile_density_kgm3 * (V_in_m3 *  1e3)  # kg/m³ * m³ = kg, then *1000 = g
                m_out_g = fertile_density_kgm3 * (V_out_m3 * 1e3)
                m_total_het_g = m_in_g + m_out_g
                
                if m_fertile_per_particle_g > 0:
                    N_in = int(np.rint(m_in_g / m_fertile_per_particle_g))
                    N_out = int(np.rint(m_out_g / m_fertile_per_particle_g))
                    N_total = N_in + N_out
                else:
                    N_in = N_out = N_total = 0
                
                pf_in = (N_in * V_p_cm3 * 1e-6) / V_in_m3 if N_in > 0 else 0.0
                pf_out = (N_out * V_p_cm3 * 1e-6) / V_out_m3 if N_out > 0 else 0.0
                pf_avg = (N_total * V_p_cm3 * 1e-6) / V_total_m3 if N_total > 0 else 0.0
                
                # PBHom calculations (from calc_biso_blanket_vol_fracs)
                # Note: PBHom uses breeding material volume only (Li+Be, not structure)
                vol_li4sio4_be = PBCoupon_BR_VOL * (0.1304 + 0.3790)  # m³ of Li+Be breeding material
                vf_breeder, vf_biso = calc_biso_blanket_vol_fracs(fertile_density_kgm3, vol_li4sio4_be, fertile_element=fertile_element)
                
                # Mass calculation for PBHom
                # fertile_density_kgm3 is per m³ of breeding material (Li+Be volume)
                m_total_hom_g = fertile_density_kgm3 * 1e3 * vol_li4sio4_be  # kg/m³ * 1000 * m³ = g
                
                # For fair comparison, calculate what PBHet mass would be if using same volume basis
                # (i.e., if fertile_density was per breeding material volume)
                m_total_het_equiv_g = fertile_density_kgm3 * 1e3 * vol_li4sio4_be
                
                # Calculate moles (using PBHet mass for now)
                moles_fertile = (m_total_het_g / 1e3) / (AMU_fertile / 1000)  # kg / (g/mol / 1000) = mol
                moles_per_m3 = moles_fertile / V_total_m3
                
                # Equivalent number of particles (if we were to pack them)
                if m_fertile_per_particle_g > 0:
                    N_equiv = int(np.rint(m_total_hom_g / m_fertile_per_particle_g))
                else:
                    N_equiv = 0
                
                comparison_data.append({
                    'Fertile Element': fertile_element,
                    'Density [kg/m³]': fertile_density_kgm3,
                    'PBHet: N_in': N_in,
                    'PBHet: N_out': N_out,
                    'PBHet: N_total': N_total,
                    'PBHet: pf_in': f'{pf_in:.4f}',
                    'PBHet: pf_out': f'{pf_out:.4f}',
                    'PBHet: pf_avg': f'{pf_avg:.4f}',
                    'PBHet: Mass [g]': f'{m_total_het_g:.2f}',
                    'PBHet: Mass (equiv) [g]': f'{m_total_het_equiv_g:.2f}',
                    'PBHom: vf_BISO': f'{vf_biso:.6f}',
                    'PBHom: N_equiv': N_equiv,
                    'PBHom: Mass [g]': f'{m_total_hom_g:.2f}',
                    'Moles [mol]': f'{moles_fertile:.4f}',
                    'Moles/m³ [mol/m³]': f'{moles_per_m3:.2f}',
                    'Mass Match (equiv)': '✓' if abs(m_total_het_equiv_g - m_total_hom_g) < 0.01 else '✗'
                })
        
        # Create DataFrame
        df_comp = pd.DataFrame(comparison_data)
        
        # Print to console
        print("\n" + "="*120)
        print("MATERIAL COMPARISON: PBHet vs PBHom")
        print("="*120)
        print(df_comp.to_string(index=False))
        print("="*120)
        
        # Create a nice table plot
        fig, ax = plt.subplots(figsize=(16, max(8, len(df_comp) * 0.5)))
        ax.axis('tight')
        ax.axis('off')
        
        # Create table
        table_data = []
        headers = ['Fertile', 'ρ [kg/m³]', 'PBHet N_in', 'PBHet N_out', 'PBHet N_tot', 
                   'PBHet pf_avg', 'PBHet Mass [g]', 'PBHet Mass(equiv) [g]',
                   'PBHom vf_BISO', 'PBHom N_equiv', 'PBHom Mass [g]', 'Moles [mol]', 'Match']
        
        for _, row in df_comp.iterrows():
            table_data.append([
                row['Fertile Element'],
                f"{row['Density [kg/m³]']:.2f}",
                f"{row['PBHet: N_in']}",
                f"{row['PBHet: N_out']}",
                f"{row['PBHet: N_total']}",
                row['PBHet: pf_avg'],
                row['PBHet: Mass [g]'],
                row['PBHet: Mass (equiv) [g]'],
                row['PBHom: vf_BISO'],
                f"{row['PBHom: N_equiv']}",
                row['PBHom: Mass [g]'],
                row['Moles [mol]'],
                row['Mass Match (equiv)']
            ])
        
        table = ax.table(cellText=table_data, colLabels=headers, 
                        cellLoc='center', loc='center',
                        colWidths=[0.06, 0.07, 0.07, 0.07, 0.07, 0.08, 0.09, 0.09, 0.08, 0.08, 0.09, 0.08, 0.05])
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 2)
        
        # Style the header
        for i in range(len(headers)):
            table[(0, i)].set_facecolor('#4CAF50')
            table[(0, i)].set_text_props(weight='bold', color='white')
        
        # Highlight matching rows
        for i, row in enumerate(df_comp.iterrows(), 1):
            if row[1]['Mass Match'] == '✓':
                for j in range(len(headers)):
                    table[(i, j)].set_facecolor('#E8F5E9')
            else:
                for j in range(len(headers)):
                    table[(i, j)].set_facecolor('#FFEBEE')
        
        plt.title('Material Comparison: PBHet (Heterogeneous) vs PBHom (Homogeneous)\n' +
                 'Verification that both models contain the same amount of fertile material',
                 fontsize=12, fontweight='bold', pad=20)
        
        plt.tight_layout()
        
        if self.save:
            plt.savefig(f'./Figures/PDF/fig_pb_coupon_material_comparison.pdf', bbox_inches='tight', format='pdf')
            plt.savefig(f'./Figures/PNG/fig_pb_coupon_material_comparison.png', bbox_inches='tight', format='png', dpi=300)
            print(f"Exported material comparison table.")
        else:
            print(f"Did not export material comparison table.")
        
        # Also save as CSV for easy reference
        df_comp.to_csv('./Figures/Data/PBHet_vs_PBHom_material_comparison.csv', index=False)
        print(f"Saved material comparison data to CSV: ./Figures/Data/PBHet_vs_PBHom_material_comparison.csv")
        
        if self.show: plt.show()
        plt.close('all')




    
if __name__ == "__main__":
    plot_all()
