"""
Code to generate FLiBe number densities
  Alex 2025-06-25
"""

import numpy as np

# ----------------------------------------------------------------------------
# FLIBE
# ----------------------------------------------------------------------------

# Mixture of LiF-BeF2, 66-34 mol%; density = 1.8-2.0 g/cc

avo = 6.02214076E+23

li6_enrichment = 0.900 # CHOOSE YOUR OWN ENRICHMENT IN WEIGHT PERCENT !!!
li6_enrichment = 0.076
li6_enrichment = 0.200

li6_enrichment = 0.200 # CHOOSE YOUR OWN ENRICHMENT IN WEIGHT PERCENT !!!


mass_li6 = 6.0150
mass_li7 = 7.0160

abun_li6 = li6_enrichment * mass_li7 / ((1-li6_enrichment) * mass_li6 + li6_enrichment * mass_li7)
abun_li7 = 1 - abun_li6

mass_be9 = 9.0120
abun_be9 = 1.0000 # Not used

mass_f19 = 18.9984
abun_f19 = 1.0000 # Not used

target_density = 1.80 # CHOOSE YOUR OWN DENSITY !!!
lif_fraction   = 0.66 # CHOOSE YOUR OWN MOLAR FRACTION !!!

bef2_fraction = 1 - lif_fraction

mass_lif  = abun_li6 * mass_li6 + abun_li7 * mass_li7 + mass_f19
mass_bef2 =  mass_be9 + 2 * mass_f19

flibe_per_cc = (avo * target_density * 10**(-24) / 
                
                  (lif_fraction * mass_lif + bef2_fraction * mass_bef2))

nden_lif  =  lif_fraction * flibe_per_cc
nden_bef2 = bef2_fraction * flibe_per_cc

nden_f1 = nden_lif
nden_li6 = nden_lif * abun_li6
nden_li7 = nden_lif * abun_li7

li6_density = nden_li6 * 10**24 / avo * mass_li6 # FOR MANUAL CHECK !!!

# target_value_manual = 0.0041426 / (1.5099223/0.8654072) * 10**24 / avo * mass_li6

nden_be = nden_bef2
nden_f2 = 2 * nden_bef2

nden_ft = nden_f1 + nden_f2

nden_total = nden_ft + nden_be + nden_li6 + nden_li7 

print('Reference numbers for full FLIBE:', nden_li6, nden_li7, nden_be, nden_ft)

# -------------
# DENSITY CHECK
# -------------

actual_density = (nden_li6 * mass_li6 + nden_li7 * mass_li7
                  
                   + nden_ft * mass_f19 + nden_be * mass_be9) * 10**24 / avo

#%%
