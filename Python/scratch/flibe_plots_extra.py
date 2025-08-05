import openmc
import os, sys, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from scipy.interpolate import make_interp_spline
import imageio.v2 as iio # use v2 to avoid deprecation warnings --ppark

"""
Parking here some extra functions I wrote for flibe_plot.py while checking data but aren't actively used 
"""


def plot_all_tbr(paths, save=True, show=False, to_csv=True):
    """
    Plots TBR per MTU for all Li-6 enrichments on one plot

    Args:
        paths : list of str : list of directories

    Returns:
        saves a plot
    """

    plt.figure()

    plt.xlim(-1.5,51.5)
    plt.ylim(1.196,1.264)
    
    plt.title(f'Tritium breeding ratio')
    plt.xlabel('Uranium loaded [metric tons]')
    plt.ylabel('Tritium breeding ratio')
    plt.tight_layout()

    for path in paths:
        df = pd.read_csv(path)
        x = df.iloc[:, 0]
        y = df.iloc[:, 3]
        plt.plot(x, y, label=os.path.basename(path))

    if save:
        plt.savefig(f'./figures/pdf/fig_tbr_all.pdf', bbox_inches='tight', format='pdf')
        plt.savefig(f'./figures/png/fig_tbr_all.png', bbox_inches='tight', format='png')
        print(f"Exported tritium breeding ratio plots.")
    if show:
        plt.show()


def plot_all_u238(paths, save=True, show=False, to_csv=True):
    """
    Plots U-238 (n,gamma) and (n,fis) per MTU for all Li-6 enrichments on one plot

    Args:
        paths : list of str : list of directories

    Returns:
        saves a plot
    """

    plt.figure()

    plt.xlim(-1.5,51.5)
    # plt.ylim(1.196,1.264)
    
    plt.title(f'Uranium-238 total reaction rates')
    plt.xlabel('Uranium loaded [metric tons]')
    plt.ylabel('Reactions $/$ soruce neutron')
    plt.tight_layout()

    for path in paths:
        df = pd.read_csv(path)
        x = df.iloc[:, 0]
        y = df.iloc[:, 3]
        plt.plot(x, y, label=os.path.basename(path))

    if save:
        plt.savefig(f'./figures/pdf/fig_u238_all.pdf', bbox_inches='tight', format='pdf')
        plt.savefig(f'./figures/png/fig_u238_all.png', bbox_inches='tight', format='png')
        print(f"Exported tritium breeding ratio plots.")
    if show:
        plt.show()


def plot_all(tbr_paths, u238_paths, save=True, show=False, to_csv=True):

    fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(15,6))
    ax = axes.flatten()

    ''' lithium-6, lithium-7 tbr '''

    for path in tbr_paths:
        df = pd.read_csv(path)        
        x = df.iloc[:, 0]
        yt = df.iloc[:,3]
        y6 = df.iloc[:,1]
        y7 = df.iloc[:,2]
        
        ax[0].plot(x, yt, label=os.path.basename(path))
        ax[1].plot(x, y6, label=os.path.basename(path))
        ax[2].plot(x, y7, label=os.path.basename(path))


    ''' uranium-238 gamma, fis '''

    for path in u238_paths:
        df = pd.read_csv(path)        
        x = df.iloc[:, 0]
        yg = df.iloc[:,1]
        yf = df.iloc[:,3]
        
        ax[3].plot(x, yg, label=os.path.basename(path))
        ax[4].plot(x, yf, label=os.path.basename(path))

    ax[0].set_title('Total TBR')
    ax[1].set_title('Li-6 TBR')
    ax[2].set_title('Li-7 TBR')
    ax[3].set_title('U-238 (n,gamma) tot rxn rate')
    ax[4].set_title('U-238 (n,fis) tot rxn rate')

    # leg = ax[0].legend(loc='lower center', ncols=6, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1, 
    # leg.get_frame().set_linewidth(1)

    fig.tight_layout() # you need this here
    fig.savefig(f'test.png', bbox_inches='tight', format='png')
    fig.show()


def plot_all_invert(tbr_paths, u238_paths, save=True, show=False, to_csv=True):

    fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(15,3))
    ax = axes.flatten()

    ''' lithium-6, lithium-7 tbr '''

    tbr_df = pd.DataFrame(columns=['mtu', 'li6_tbr', 'li7_tbr', 'tbr', 'tbr_err', 'li_enrich'])
    for path in tbr_paths:
        df = pd.read_csv(path)
        df['li_enrich'] = extract_lie(path) 
        tbr_df = pd.concat([tbr_df, df], ignore_index=True)

    for lie, group in tbr_df.groupby('mtu'):
        ax[0].plot(group['li_enrich'],group['tbr'],label=f"{lie}")
        ax[1].plot(group['li_enrich'],group['li6_tbr'],label=f"{lie}")
        ax[2].plot(group['li_enrich'],group['li7_tbr'],label=f"{lie}")


    ''' uranium-238 gamma, fis '''

    u238_df = pd.DataFrame(columns=['mtu', 'U-238(n,gamma)', 'U-238(n,gamma)_err', 'U-238(n,fis)', 'U-238(n,fis)_err', 'li_enrich'])
    for path in u238_paths:
        df = pd.read_csv(path)
        df['li_enrich'] = extract_lie(path) 
        u238_df = pd.concat([u238_df, df], ignore_index=True)

    for lie, group in u238_df.groupby('mtu'):
        ax[3].plot(group['li_enrich'],group['U-238(n,gamma)'],label=f"{lie}")
        ax[4].plot(group['li_enrich'],group['U-238(n,fis)'],label=f"{lie}")

    ax[0].set_title('Total TBR')
    ax[1].set_title('Li-6 TBR')
    ax[2].set_title('Li-7 TBR')
    ax[3].set_title('U-238 (n,gamma) tot rxn rate')
    ax[4].set_title('U-238 (n,fis) tot rxn rate')

    # leg = ax[0].legend(loc='lower center', ncols=6, frameon=True, fancybox=False, edgecolor='black', framealpha=.75,) # fontsize='small', ncols=1, 
    # leg.get_frame().set_linewidth(1)

    fig.tight_layout() # you need this here
    fig.savefig(f'test_inv.png', bbox_inches='tight', format='png')
    fig.show()



