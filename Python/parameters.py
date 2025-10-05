""" USER SPECIFICATIONS
Nominal values used by Emma and Patrick + reasons
  DENSITY_FLIBE : float : 1.94 [g/cc] : used by Ball 25
  DENSITY_UF4   : float : 6.70 [g/cc] : used by Ball 25
  ENRICH_LI : float : 7.5 [wt%] : used by Ball 25 / 20 wt% used by Alex
  ENRICH_U : float : 0.7204 [wt%] : natural uranium from PNNL-15870 
  TEMP_K : float : 900 [K]
  VOL_CC : float : 342,000,000 [cc = 342 m3]
  MASS_U_LIST : list of floats : 0, 0.1, 1, 2.5, 5, 10, 20, 30, 40, 50 [metric tons U]
"""



""" Parameters for ALL models """


FERTILE_BULK_DENSITY_KGM3 = [0, 0.03, 0.3, 0.6, 1.5, 3, 7.5, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150]

DENSITY_FLIBE = 1.94 # g/cm3
ENRICH_FLIBE  = 7.5

DENSITY_UF4   = 6.7 # g/cm3
ENRICH_U  = 0.71  # wt% 
TEMP_K = 900 # K

DENSITY_BISO = 6.93759 #g/cc


""" Parameters for ARC-class tokamak with FLiBe Breeder """
ARC_R0    = 400  # cm - major radius
ARC_A     = 100  # cm - minor radius
ARC_KAPPA = 1.6  #    - elongation
ARC_DELTA = 0.5  #    - triangularity

ARC_FW_CM  =   0.3  # cm - first wall
ARC_ST1_CM =   1.0  # cm - structural region 1
ARC_BR1_CM =   2.0  # cm -   breeding region 1
ARC_ST2_CM =   3.0  # cm - structural region 2
ARC_BR2_CM = 100.0  # cm -   breeding region 2
ARC_ST3_CM =   3.0  # cm - structural region 3

ARC_BR_VOL = 291.1 # m^3 - analytic vol from plot_miller_models.py



""" Parameters for tokamak with FLiBe Breeder """
FLIBE_R0    = 600  # cm - major radius
FLIBE_A     = 200  # cm - minor radius
FLIBE_KAPPA = 1.72 #    - elongation
FLIBE_DELTA = 0.4  #    - triangularity

FLIBE_FW_CM  =   0.2  # cm - first wall
FLIBE_ST1_CM =   1.0  # cm - structural region 1
FLIBE_BR1_CM =   2.0  # cm -   breeding region 1
FLIBE_ST2_CM =   3.0  # cm - structural region 2
FLIBE_BR2_CM = 100.0  # cm -   breeding region 2
FLIBE_ST3_CM =   3.0  # cm - structural region 3

FLIBE_BR_VOL = 771.8 # m^3 - analytic vol from plot_miller_models.py


""" Parameters for tokamak with LL Breeder """


DENSITY_LL  =  9.4 # g/cm3
ENRICH_LL   = 90.0 

LL_R0    = 620  # cm - major radius
LL_A     = 207  # cm - minor radius
LL_KAPPA = 1.72 #    - elongation
LL_DELTA = 0.4  #    - triangularity
#----Outboard Layers----
LL_FW_O_CM  =   0.2  # cm -          first wall
LL_FWF_O_CM =   0.4  # cm -    first wall front
LL_FWC_O_CM =   2.0  # cm -  first wall cooling
LL_FWB_O_CM =   0.4  # cm -     first wall back
LL_BR1_O_CM =  22.5  # cm -   breeding region 1
LL_D1_O_CM  =   3.2  # cm -           divider 1
LL_BR2_O_CM =  21.0  # cm -   breeding region 2
LL_D2_O_CM =   3.2  # cm -           divider 2
LL_BR3_O_CM =  21.0  # cm -   breeding region 3
LL_IM_O_CM =   8.0  # cm -      inner manifold 
#addititional F2H8 back plate

#----Inboard Layers----
# LL_FW_I_CM  =   0.2  # cm -          first wall
# LL_FWF_I_CM =   0.4  # cm -    first wall front
# LL_FWC_I_CM =   2.0  # cm -  first wall cooling
# LL_FWB_I_CM =   0.4  # cm -     first wall back
# LL_BR1_I_CM =  22.5  # cm -   breeding region 1
# LL_D1_I_CM  =   3.2  # cm -           divider 1
# LL_BR2_I_CM =  21.0  # cm -   breeding region 2
LL_IM_I_CM =   8.0  # cm -      inner manifold 

LL_BR_VOL =  440.3757# m^3 - analytic vol from plot_miller_models.py
#Check this value why is it so low?? -ezoccoli






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