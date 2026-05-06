""" Universal parameters for ALL models """

# Run settings
FERTILE_KGM3 = [0.0, 0.10, 0.50, 1, 10, 25, 50, 75, 100, 150, 250, 500, 750, 999.99]
BLANKETS = ['FLiBe','DCLL','HCPB']
ISOTOPES = ['U238', 'Th232']
BLANKET_COVERAGE = 1.00 # Assume blanket covers 88% of plasma surface, rest for divertor -- removed as irrelevant for our conclusions --ppark 2026-02-11
TEMP_K = 900 # [K] ENDF data generally has: [250, 294, 600, 900, 1200, 2500 K]

N_PARTICLES, N_CYCLES = int(4e6), 25

# Geometry is split into 1/10 sector ... so multiply all tallies by 10x !! -- ppark 2026-05-01
SECTOR            = 1/10        
SECTOR_DEG        = SECTOR * 360
SECTOR_TALLY_MULT = 1/SECTOR

# Fusion power parameters - from JL Ball 24
EV_TO_MJ = 1.602e-25  # 1.602e-19 J/eV × 1e-6 MJ/J
N_PER_MJ = 1 / (17.6 * 1.602e-19)  # 1 MJ / (17.6 MeV × 1.602e-19 J/eV) = 3.546e17 neutrons/MJ 
P_FUS_MW = 1000
NPS_FUS  = P_FUS_MW * N_PER_MJ # n/s


""" Material parameters of breeders """

# FLiBe (2(LiF)-BeF2)
DENSITY_FLIBE   =  1.9505  # [g/cm³] from EOS at 900 K and 101 kPa, calculated in Jupyter/FLiBe/flibe_mats.py, ref. Humrickhouse 17 (INL-44148) --ppark 2026-02-13
ENRICH_FLIBE    =  7.50    # [at%] enrich of Li-6
LIMIT_UF4_FLIBE  = 1000    # 768.0   # [kg/m³]
LIMIT_THF4_FLIBE = 1000    # [kg/m³]

# Lead-lithium (DCLL, 83 at% Pb - 17 at% Li)
DENSITY_DCLL    =  9.40  # [g/cm³]
ENRICH_DCLL     = 90.00  # [at%] enrich of Li-6

# Pebble bed (HCPB, Li4SiO4-Be)
DENSITY_LI4SIO4 =  2.42  # [g/cm³] # should we normalize to ~2.17 for ceramic porosity and 900 K? (pure, room temp g/cm3 = 2.42) --ppark
DENSITY_BE      =  1.85  # [g/cm³] # could normalize to 1.80 for 900 K from 1.85 at STP
ENRICH_HCPB     = 60.00  # [at%] enrich of Li-6


""" Material parameters of fertile material """

# Densities 6.88 g/cm³ (UF4) and 6.61 g/cm³ (ThF4) are computed from molar volume when added to 2(LiF)-BeF2 at 900 K 
# respectively per Cantor (1965) ORNL-TM-3913, p.29 and Cantor (1973) ORNL-TM-4308, p.13 --ppark 2026-02-13
DENSITY_UF4  =  6.88 # [g/cm³] nominally 6.7 for powder mass at room temp 
DENSITY_ThF4 =  6.61 # [g/cm³] nominally 6.3 for powder mass at room temp 

# BISO particle parameters
ENRICH_U     =  0.71 # [wt%]
DENSITY_UO2  = 10.50 # [g/cm³] 
DENSITY_ThO2 = 10.00 # [g/cm³]
DENSITY_SiC  =  3.20 # [g/cm³]

BISO_KERNEL_RADIUS   = 0.04 # [cm]  (400 μm inner radius / 800 μm wide)
BISO_RADIUS          = 0.05 # [cm]  (500 μm inner radius / 100 μm thick)
BISO_VOLUME          = (4/3)*3.14159265359*(BISO_RADIUS)**3         # volume of single BISO particle - trying to avoid np import in this file
KERNEL_VOLUME        = (4/3)*3.14159265359*(BISO_KERNEL_RADIUS)**3  # volume of UO2/ThO2 kernel in single BISO particle
BISO_KERNEL_VOL_FRAC = KERNEL_VOLUME / BISO_VOLUME           # (= 0.512) volume fraction of the kernel in BISO
BISO_COAT_VOL_FRAC   = 1.0 - BISO_KERNEL_VOL_FRAC            # (= 0.488) volume fraction of the SiC coating


""" Geometric parameters for tokamaks """

# --------------------------------------------
# ARC-CLASS TOKAMAK WITH FLIBE BLANKET
# --------------------------------------------

# ARC tokamak specifications from Sorbom et al. (2015)
ARC_R0    = 400     # [cm] major radius
ARC_A     = 100     # [cm] minor radius
ARC_KAPPA = 1.6     #      elongation
ARC_DELTA = 0.5     #      triangularity

# ARC-class tokamak in J. L. Ball et al. (2025)
ARCB_R0     = 400   # [cm] major radius
ARCB_A      = 120   # [cm] minor radius
ARCB_KAPPA  = 1.5   #      elongation
ARCB_DELTA  = 0.5   #      triangularity
ARCB_BL_VOL = 320   # [m³] from ./OpenMC/volume_ARCBall_900K_Li7.5_U000.00kgm3/volume_1.csv

ARC_FW_CM  =   0.3  # [cm] first wall
ARC_ST1_CM =   1.0  # [cm] structural region 1
ARC_BR1_CM =   2.0  # [cm]   breeding region 1
ARC_ST2_CM =   3.0  # [cm] structural region 2
ARC_BR2_CM = 100.0  # [cm]   breeding region 2
ARC_ST3_CM =   3.0  # [cm] structural region 3

ARC_BL_VOL = 291.1  # [m³] from ./OpenMC/volume_ARC_900K_Li7.5_U000.00kgm3/volume_1.csv


# --------------------------------------------
# FLiBe 
# --------------------------------------------

FLIBE_R0     = 600    # [cm] major radius
FLIBE_A      = 200    # [cm] minor radius
FLIBE_KAPPA  = 1.72   #      elongation
FLIBE_DELTA  = 0.4    #      triangularity

FLIBE_FW_CM  =   0.2  # [cm] first wall
FLIBE_BE_CM  =   1.0  # [cm] beryllium neutron multiplier layer
FLIBE_ST1_CM =   1.0  # [cm] structural region 1
FLIBE_BR1_CM =   2.0  # [cm]   breeding region 1
FLIBE_ST2_CM =   3.0  # [cm] structural region 2
FLIBE_BR2_CM = 100.0  # [cm]   breeding region 2
FLIBE_ST3_CM =   3.0  # [cm] structural region 3

FLIBE_BL_VOL = 773.3 # [m³] from ./OpenMC/volume_DCLL_900K_Li90.0_U000.00kgm3/volume_1.csv 
# value IS updated run after adding Be layer, full val: 12.7878381613389 + 760.486252951442 = 773.274091112781


# --------------------------------------------
# DUAL-COOLANT LEAD-LITHIUM (Pb-17Li)
# --------------------------------------------

DCLL_R0       = 620     # [cm] major radius
DCLL_A        = 207     # [cm] minor radius
DCLL_KAPPA    = 1.72    #      elongation
DCLL_DELTA    = 0.4     #      triangularity

DCLL_FW_CM    =   0.2  # [cm]          first wall
DCLL_FWF_CM   =   0.4  # [cm]    first wall front
DCLL_FWC_CM   =   2.0  # [cm]  first wall cooling
DCLL_FWB_CM   =   0.4  # [cm]     first wall back
DCLL_BR1_CM   =  22.5  # [cm]   breeding region 1
DCLL_D1_CM    =   3.2  # [cm]           divider 1
DCLL_BR2_CM   =  21.0  # [cm]   breeding region 2
DCLL_D2_CM    =   3.2  # [cm]           divider 2 (only on outboard)
DCLL_BR3_CM   =  21.0  # [cm]   breeding region 3 (only on outboard)
DCLL_IM_CM    =   8.0  # [cm]      inner manifold 
DCLL_BP_CM    =   1.5  # [cm]          back plate 
DCLL_SS_CM    =  20.0  # [cm]        steel shield

# Volume fractions of material in breeding region of DCLL blanket (Glaser & Goldston 2012, Tb.1, Breeding Ch.1)
DCLL_VF_FS_NOM = 0.019 
DCLL_VF_LL_NOM = 0.808  
DCLL_VF_SI_NOM = 0.076
DCLL_VF_HE_NOM = 0.097

DCLL_BL_VOL = 427.3   # [m³] from ./OpenMC/volume_DCLL_900K_Li90.0_U000.00kgm3/volume_1.csv
DCLL_BR_VOL = DCLL_BL_VOL * DCLL_VF_LL_NOM  


# --------------------------------------------
# HELIUM-COOLED PEBBLE BED BLANKET
# --------------------------------------------

HCPB_R0       = 620    # [cm] major radius
HCPB_A        = 207    # [cm] minor radius
HCPB_KAPPA    = 1.72   #      elongation
HCPB_DELTA    = 0.4    #      triangularity
             
HCPB_FW_CM    =   0.2  # [cm] first wall
HCPB_ST1_CM   =   2.5  # [cm] structural region 1
HCPB_BR1_O_CM =  82.0  # [cm]   breeding region outboard
HCPB_BR1_I_CM =  45.0  # [cm]   breeding region inboard
HCPB_ST2_CM   =   3.0  # [cm] structural region 2

# Volume fractions of material in breeding regions of HCPB blanket (Lu et al. 2017, Tb.2)
HCPB_VF_LI_NOM = 0.1304  # Li4SiO4 (ceramic)
HCPB_VF_BE_NOM = 0.3790  # Be (metal)
HCPB_VF_EU_NOM = 0.1176  # Eurofer
HCPB_VF_HE_NOM = 0.3730  # He-4 (gas)

HCPB_BL_VOL = 520.6 # [m³] from ./OpenMC/volume_HCPB_900K_Li60.0_U000.00kgm3/volume_1.csv
HCPB_BR_VOL = HCPB_BL_VOL * (HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)


HCPB_CONV_U_FPR = [ (   0.1, 0.6092),   # see wedge_data_2026-04-102.xlsx for details  
                    (   0.5, 0.6203),   # -- ppark 2026-04-14
                    (     1, 0.6175),   
                    (    10, 0.6376),
                    (    25, 0.6620),
                    (    50, 0.7051),
                    (    75, 0.7294),
                    (   100, 0.7588),
                    (   150, 0.7955),
                    (   250, 0.8438),
                    (   500, 0.9047),
                    (   750, 0.9343),
                    (999.99, 0.9446), ]

HCPB_CONV_TH_FPR = [(   0.1, 0.7950),
                    (   0.5, 0.7987),
                    (     1, 0.7952),
                    (    10, 0.8035),
                    (    25, 0.8109),
                    (    50, 0.8290),
                    (    75, 0.8369),
                    (   100, 0.8543),
                    (   150, 0.8711),
                    (   250, 0.8976),
                    (   500, 0.9356),
                    (   750, 0.9494),
                    (999.99, 0.9553), ]



"""
MATPLOTLIB SETTINGS
"""
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# Fonts
# try:
#     font_path = './Python/fonts/DIN-Regular.ttf' # DIN-Regular.ttf' # './Python/arial.ttf'
#     fm.fontManager.addfont(font_path)
#     prop = fm.FontProperties(fname=font_path)
#     plt.rcParams['font.family'] = prop.get_name()
# except:
font_path = './Python/fonts/arial.ttf' # './arial.ttf'
fm.fontManager.addfont(font_path)
prop = fm.FontProperties(fname=font_path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['pdf.fonttype'] = 42  # By default, matplotlib converts text to outlines (paths) in PDFs; type 42 embeds text as actual TrueType characters
LONG_DASH = (0, (10, 2))

# Ticks
plt.rcParams['xtick.direction']   = 'in'
plt.rcParams['ytick.direction']   = 'in'
plt.rcParams['xtick.major.width'] = 0.6  # major x‐tick line width
plt.rcParams['ytick.major.width'] = 0.6  # major y‐tick line width
plt.rcParams['xtick.major.size']  = 3.5  # major x‐tick line length
plt.rcParams['ytick.major.size']  = 3.5  # major y‐tick line length
plt.rcParams['xtick.minor.size']  = 2.0  # minor x‐tick line length
plt.rcParams['ytick.minor.size']  = 2.0  # minor y‐tick line length
plt.rcParams['axes.linewidth']    = 0.6  # axes spines linewidth (use this instead of messing with spine.set_linewidth(1) )
plt.rcParams['grid.color'] = '#DBDBDB'

# Font sizes
plt.rcParams['font.size']        =  8   # default text size for labels, legends, etc.
plt.rcParams['axes.titlesize']   =  9   # axes title
plt.rcParams['axes.labelsize']   =  8   # x- and y-axis labels
plt.rcParams['xtick.labelsize']  =  7   # x-tick labels
plt.rcParams['ytick.labelsize']  =  7   # y-tick labels
plt.rcParams['legend.fontsize']  =  7   # legend text
plt.rcParams['figure.titlesize'] = 10   # figure title
