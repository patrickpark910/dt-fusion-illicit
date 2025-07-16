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

DENSITY_FLIBE = 1.94 # g/cm3
DENSITY_UF4 = 6.7 # g/cm3

DENSITY_PBLI  = 9.4 # g/cm3
DENSITY_TRISO = 7 # g/cm3

ENRICH_U  = 0.7204 # wt% 
TEMP_K = 900 # K
VOL_CC = 342 * 1e6 # cm3
MASS_U_LIST = [0, 0.0096, 0.1, 0.2, 0.5, 1, 2.5, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50] # [0, 0.1, 1, 5, 10, 20, 30, 40, 50] # metric tons
MASS_U_LIST_FLIBE = [0, 0.0096, 0.1, 0.2, 0.5, 1, 2.5, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50] # [0, 0.1, 1, 5, 10, 20, 30, 40, 50] # metric tons uranium (MTU) 0.01209475,
TEMP_LIST = [294,900]

MASS_U_LIST_PBLI = [0, 0.1, 1, 5, 10, 20, 30, 40, 50]
ENRICH_LI_LIST = [7.5,10.0,12.5,15.0,17.5,20.0] # weight percent

