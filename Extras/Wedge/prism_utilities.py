"""
Set of helper functions
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


# GLOBAL PARAMETERS
TEMP_K = 900
BISO_VOLUME = 4/3*np.pi*(0.05)**3
KERNEL_VOLUME = 4/3*np.pi*(0.04)**3
BISO_KERNEL_VOL_FRAC = KERNEL_VOLUME / BISO_VOLUME
BISO_COAT_VOL_FRAC   = 1.0 - BISO_KERNEL_VOL_FRAC

DENSITY_UO2  = 10.50 # [g/cm³]
DENSITY_ThO2 = 10.00 # [g/cm³]

ENRICH_U = 0.71 # [wt%]
AMU_U235 = 235.0439299
AMU_U238 = 238.05078826
AMU_U = (ENRICH_U * AMU_U235 + (100-ENRICH_U) * AMU_U238)/100
AMU_O = 15.999
AMU_UO2 = AMU_U + 2 * AMU_O
AMU_Th232 = 232.0381
AMU_ThO2 = AMU_Th232 + 2 * AMU_O

# Lead-lithium (DCLL, 83 at% Pb - 17 at% Li)
DENSITY_DCLL    =  9.40  # [g/cm³]
ENRICH_DCLL     = 90.00  # [at%] enrich of Li-6

# Volume fractions in DCLL blanket (Glaser & Goldston 2012, Tb.1)
DCLL_VF_FS_NOM = 0.019 
DCLL_VF_LL_NOM = 0.808  
DCLL_VF_SI_NOM = 0.076
DCLL_VF_HE_NOM = 0.097

# Volume fractions in HCPB blanket
HCPB_VF_LI_NOM = 0.1304
HCPB_VF_BE_NOM = 0.3790
HCPB_VF_EU_NOM = 0.1176
HCPB_VF_HE_NOM = 0.3730

HCPB_VF_BREEDER_NOM = HCPB_VF_LI_NOM + HCPB_VF_BE_NOM
HCPB_VF_LI_IN_BREEDER_NOM = HCPB_VF_LI_NOM / HCPB_VF_BREEDER_NOM
HCPB_VF_BE_IN_BREEDER_NOM = HCPB_VF_BE_NOM / HCPB_VF_BREEDER_NOM


font_path = './fonts/arial.ttf'
fm.fontManager.addfont(font_path)
prop = fm.FontProperties(fname=font_path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams['mathtext.default'] = 'regular'

plt.rcParams['xtick.direction']   = 'in'
plt.rcParams['ytick.direction']   = 'in'
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['xtick.major.size']  = 6
plt.rcParams['ytick.major.size']  = 6
plt.rcParams['xtick.minor.size']  = 3
plt.rcParams['ytick.minor.size']  = 3
plt.rcParams['axes.linewidth']    = 0.5
plt.rcParams['grid.color']        = '#DBDBDB'

plt.rcParams['font.size']         = 18
plt.rcParams['axes.titlesize']    = 16
plt.rcParams['axes.labelsize']    = 18
plt.rcParams['xtick.labelsize']   = 14
plt.rcParams['ytick.labelsize']   = 14
plt.rcParams['legend.fontsize']   = 14
plt.rcParams['figure.titlesize']  = 16


def fertile_kgm3_to_biso_per_cc(fertile_kgm3, fertile_isotope='U238'):
    """
    Converts given kg/m³ of fertile material to the number of BISO particles per cm³

    BISO/cc = (g_isotope / cm³_ref) → (g_compound / cm³_ref) → (kernels / cm³_ref)

    Args:
        fertile_kgm3  (float): kg of fertile isotope per m³ of nominal breeder
        fertile_isotope (str): one of ['U238', 'Th232']

    Returns:
        biso_per_cc (float): number of biso particles per m³ of nominal breeder
    """
    if fertile_isotope == 'U238':
        x28 = (100 - ENRICH_U) / 100 # = 0.9929
        biso_per_cc = fertile_kgm3 * AMU_UO2 / (x28 * AMU_U238) / KERNEL_VOLUME / DENSITY_UO2 / 100**3 * 1000
    elif fertile_isotope == 'Th232':
        biso_per_cc = fertile_kgm3 * AMU_ThO2 / AMU_Th232 / KERNEL_VOLUME / DENSITY_ThO2 / 100**3 * 1000
    return biso_per_cc


def calc_biso_vol_fracs(fertile_kgm3, fertile_isotope='U238'):
    """
    Calculate volume fractions of BISO particles and breeder material.
    Per Glaser & Goldston (2012), we assume the BISO/TRISO particles are homogenized
    throughout the blanket. The kernel/coating geometry is used only to set the
    relative mass or volume fractions of fertile vs. coating material.

    Args:
        fertile_kgm3  (float): kg of fertile isotope per m³ of nominal breeder
        fertile_isotope (str): one of ['U238', 'Th232']

    Returns:
        vf_biso     (float): vol frac of BISO relative to nominal breeder volume
        vf_breeder  (float): vol frac of background material (everything else) relative to nominal breeder volume
        biso_per_cc (float): number of BISO spheres per cm³ of nominal breeder
    """

    # Number of BISO spheres per cc of nominal breeder
    biso_per_cc  = fertile_kgm3_to_biso_per_cc(fertile_kgm3, fertile_isotope=fertile_isotope)
    vf_biso      = biso_per_cc * BISO_VOLUME
    vf_breeder   = 1 - vf_biso

    if vf_biso > 1.0:
        print(f"Fatal. Your fertile kg/m³ exceeds the nominal breeder volume!")
        sys.exit()

    return vf_biso, vf_breeder, biso_per_cc


def has_statepoint(directory_path):
    """
    Check if any file starting with 'statepoint' exists in the given directory.
    
    Args:
        directory_path (str): Path to the directory to search
    
    Returns:
        bool: True if a file starting with 'statepoint' is found, False otherwise
    """
    found = False
    for filename in os.listdir(directory_path):
        if filename.startswith("statepoint"):
            found = True
    return found


def logspace_per_decade(start, stop, pts_per_decade):
    """
    Returns values from 'start' to 'stop' so that each factor-of-10
    interval contains 'pts_per_decade' points (including its first endpoint).
    Might be a little off if 'stop' isn't precisely at a decade, ex. 20e6 eV

    example: 10 points per decade from 1e-5 → 2e7
    grid = logspace_per_decade(1e-5, 20e6, pts_per_decade=10)
    for i in grid:
        print(np.log10(i))
    # print(np.log10(i) for i in grid)
    """
    log_start = np.log10(start)
    log_stop  = np.log10(stop)
    total_decades = log_stop - log_start
    # how many intervals of size 1 decade we need, as a float
    total_steps = total_decades * pts_per_decade
    # +1 so that we include the very last point at `stop`
    npts = int(np.ceil(total_steps)) + 1
    return np.logspace(log_start, log_stop, num=npts)


