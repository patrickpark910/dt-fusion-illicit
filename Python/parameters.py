
""" Universal parameters for ALL models """

# Run settings
#FERTILE_BULK_DENSITY_KGM3 = [0, 0.03, 0.3, 0.6, 1.5, 3, 7.5, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 250, 500, 750, 1000]
#BREEDERS = ['ARC','ARCBall','FLiBe','LL','PB']
FERTILE_BULK_DENSITY_KGM3 = [0, 0.03, 1.5, 15, 90, 150]
BREEDERS = ['PB']
BLANKET_COVERAGE = 0.8 # Assume blanket covers 80% of plasma surface, rest for divertor
TEMP_K = 900 # [K] ENDF data generally has: [250, 294, 600, 900, 1200, 2500 K]

#N_PARTICLES, N_CYCLES = int(1e6), 10
N_PARTICLES, N_CYCLES = int(1e2), 4

""" Material parameters of breeders """

# FLiBe (F4Li2Be)
DENSITY_FLIBE   =  1.94  # [g/cm³]
ENRICH_FLIBE    =  7.50  # wt% enrich of Li-6

# Lead-lithium (LL, 83 at% Pb - 17 at% Li)
DENSITY_LL      =  9.40  # [g/cm³]
ENRICH_LL       = 90.00  # at% enrich of Li-6

# Pebble bed (PB, Li4SiO4-Be)
DENSITY_LI4SIO4 =  2.17  # [g/cm³]
DENSITY_BE      =  1.80  # [g/cm³]
ENRICH_PB       = 60.00  # at% enrich of Li-6 


""" Material parameters of fertile material """

# BISO particle parameters
ENRICH_U     =  0.71 # wt% 
DENSITY_UO2  = 10.50 # [g/cm³]
DENSITY_ThO2 = 10.00 # [g/cm³]
DENSITY_UF4  =  6.70 # [g/cm³]
DENSITY_ThF4 =  6.30 # [g/cm³]
DENSITY_SiC  =  3.20 # [g/cm³]

BISO_KERNEL_RADIUS   = 0.04 # [cm]  (400 μm / 800 μm wide)
BISO_RADIUS          = 0.05 # [cm]  (500 μm / 100 μm thick)
BISO_VOLUME          = (4/3)*3.1416*(BISO_RADIUS)**3         # volume of single BISO particle
KERNEL_VOLUME        = (4/3)*3.1416*(BISO_KERNEL_RADIUS)**3  # volume of UO2/ThO2 kernel in single BISO particle
BISO_KERNEL_VOL_FRAC = KERNEL_VOLUME / BISO_VOLUME           # (= 0.512) volume fraction of the kernel in BISO
BISO_COAT_VOL_FRAC   = 1.0 - BISO_KERNEL_VOL_FRAC            # (= 0.488) volume fraction of the SiC coating


""" Geometric parameters for tokamaks """

# ARC-class tokamak with FLiBe Breeder
ARC_R0    = 400     # [cm] major radius
ARC_A     = 100     # [cm] minor radius
ARC_KAPPA = 1.6     #      elongation
ARC_DELTA = 0.5     #      triangularity
ARC_BR_VOL = 291.1  # [m³] from ./OpenMC/volume_ARC_900K_Li7.5_U000.00kgm3/volume_1.csv

# ARC-class tokamak in J. L. Ball et al. (25)
BALL_R0     = 400   # [cm] major radius
BALL_A      = 120   # [cm] minor radius
BALL_KAPPA  = 1.5   #      elongation
BALL_DELTA  = 0.5   #      triangularity
BALL_BR_VOL = 320   # [m³] from ./OpenMC/volume_ARCBall_900K_Li7.5_U000.00kgm3/volume_1.csv

ARC_FW_CM  =   0.3  # [cm] first wall
ARC_ST1_CM =   1.0  # [cm] structural region 1
ARC_BR1_CM =   2.0  # [cm]   breeding region 1
ARC_ST2_CM =   3.0  # [cm] structural region 2
ARC_BR2_CM = 100.0  # [cm]   breeding region 2
ARC_ST3_CM =   3.0  # [cm] structural region 3

# Our FLiBe blanket
FLIBE_R0     = 600    # [cm] major radius
FLIBE_A      = 200    # [cm] minor radius
FLIBE_KAPPA  = 1.72   #      elongation
FLIBE_DELTA  = 0.4    #      triangularity

FLIBE_FW_CM  =   0.2  # [cm] first wall
FLIBE_ST1_CM =   1.0  # [cm] structural region 1
FLIBE_BR1_CM =   2.0  # [cm]   breeding region 1
FLIBE_ST2_CM =   3.0  # [cm] structural region 2
FLIBE_BR2_CM = 100.0  # [cm]   breeding region 2
FLIBE_ST3_CM =   3.0  # [cm] structural region 3

FLIBE_BR_VOL = 771.8  # [m³] from ./OpenMC/volume_FLiBe_900K_Li7.5_U000.00kgm3/volume_1.csv


# Our LL blanket
LL_R0       = 620     # [cm] major radius
LL_A        = 207     # [cm] minor radius
LL_KAPPA    = 1.72    #      elongation
LL_DELTA    = 0.4     #      triangularity

LL_FW_CM    =   0.2  # [cm]          first wall
LL_FWF_CM   =   0.4  # [cm]    first wall front
LL_FWC_CM   =   2.0  # [cm]  first wall cooling
LL_FWB_CM   =   0.4  # [cm]     first wall back
LL_BR1_CM   =  22.5  # [cm]   breeding region 1
LL_D1_CM    =   3.2  # [cm]           divider 1
LL_BR2_CM   =  21.0  # [cm]   breeding region 2
LL_D2_O_CM  =   3.2  # [cm]           divider 2 (only on outboard)
LL_BR3_O_CM =  21.0  # [cm]   breeding region 3 (only on outboard)
LL_IM_CM    =   8.0  # [cm]      inner manifold 

LL_BR_VOL   =  427.3 # [m³] from ./OpenMC/volume_LL_900K_Li90.0_U000.00kgm3/volume_1.csv


# Our PB blanket
PB_R0       = 620    # [cm] major radius
PB_A        = 207    # [cm] minor radius
PB_KAPPA    = 1.72   #      elongation
PB_DELTA    = 0.4    #      triangularity
             
PB_FW_CM    =   0.2  # [cm] first wall
PB_ST1_CM   =   2.5  # [cm] structural region 1
PB_BR1_O_CM =  82.0  # [cm]   breeding region outboard
PB_BR1_I_CM =  45.0  # [cm]   breeding region inboard
PB_ST2_CM   =   3.0  # [cm] structural region 2

PB_BR_VOL   = 520.6 # [m³] from ./OpenMC/volume_PB_900K_Li60.0_U000.00kgm3/volume_1.csv
PBCoupon_BR_VOL   = 520.6e-3 # [m³]

# Our PB coupon
PB_RO       = 830    # [cm] plasma facing outer wall
PB_RI       = 410    # [cm] plasma facing inner wall






"""
MATPLOTLIB SETTINGS
"""
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# Fonts
try:
    font_path = './Python/fonts/DIN-Regular.ttf' # DIN-Regular.ttf' # './Python/arial.ttf'
    fm.fontManager.addfont(font_path)
    prop = fm.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = prop.get_name()
except:
    font_path = './Python/fonts/arial.ttf' # './arial.ttf'
    fm.fontManager.addfont(font_path)
    prop = fm.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = prop.get_name()

# Ticks
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 0.5   # major x‐tick line width
plt.rcParams['ytick.major.width'] = 0.5   # major y‐tick line width
plt.rcParams['xtick.major.size'] = 6      # major x‐tick line length
plt.rcParams['ytick.major.size'] = 6      # major y‐tick line length
plt.rcParams['xtick.minor.size'] = 3      # minor x‐tick line length
plt.rcParams['ytick.minor.size'] = 3      # minor y‐tick line length
plt.rcParams['axes.linewidth']    = 0.5   # axes spines linewidth (use this instead of messing with spine.set_linewidth(1) )

# Font sizes
plt.rcParams['font.size']          = 14   # default text size for labels, legends, etc.
plt.rcParams['axes.titlesize']     = 14   # axes title
plt.rcParams['axes.labelsize']     = 14   # x- and y-axis labels
plt.rcParams['xtick.labelsize']    = 14   # x-tick labels
plt.rcParams['ytick.labelsize']    = 14   # y-tick labels
plt.rcParams['legend.fontsize']    = 12   # legend text
plt.rcParams['figure.titlesize']   = 14   # figure title
