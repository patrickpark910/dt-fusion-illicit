"""
Fertile loading limits for FLiBe, DCLL, and HCPB blankets.
Plots fertile kg/m³ vs. vol% of nominal breeder volume for U238 and Th232.
"""
import os, sys

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
os.chdir(project_root)
sys.path.insert(0, project_root)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from Python.parameters import *
from Python.utilities import *


fertile_kgm3 = np.linspace(0.01, 6000, 10000)
x28 = (100 - ENRICH_U) / 100

# --- FLiBe: fertile dissolved as XF4 (cf. calc_xf4_vol_fracs) ---
vf_flibe_u  = fertile_kgm3 * (AMU_UF4  / (x28 * AMU_U238)) / (DENSITY_UF4  * 1000) * 100
vf_flibe_th = fertile_kgm3 * (AMU_ThF4 / AMU_Th232)         / (DENSITY_ThF4 * 1000) * 100

# --- DCLL & HCPB: fertile in BISO particles (cf. calc_biso_vol_fracs) ---
biso_cc_u  = fertile_kgm3 * (AMU_UO2  / (x28 * AMU_U238)) / KERNEL_VOLUME / DENSITY_UO2  / 100**3 * 1000
biso_cc_th = fertile_kgm3 * (AMU_ThO2 / AMU_Th232)         / KERNEL_VOLUME / DENSITY_ThO2 / 100**3 * 1000
vf_biso_u  = biso_cc_u  * BISO_VOLUME * 100
vf_biso_th = biso_cc_th * BISO_VOLUME * 100

# Mask beyond physical loading limit (vol% > 100)
vf_flibe_u  = np.where(vf_flibe_u  <= 100, vf_flibe_u,  np.nan)
vf_flibe_th = np.where(vf_flibe_th <= 100, vf_flibe_th, np.nan)
vf_biso_u   = np.where(vf_biso_u   <= 100, vf_biso_u,   np.nan)
vf_biso_th  = np.where(vf_biso_th  <= 100, vf_biso_th,  np.nan)

# Loading limits (kg/m³ where vol% hits 100)
limit_flibe_u  = (DENSITY_UF4  * 1000) / (AMU_UF4  / (x28 * AMU_U238))
limit_flibe_th = (DENSITY_ThF4 * 1000) / (AMU_ThF4 / AMU_Th232)
limit_biso_u   = 1 / BISO_VOLUME / ((AMU_UO2  / (x28 * AMU_U238)) / KERNEL_VOLUME / DENSITY_UO2  / 100**3 * 1000)
limit_biso_th  = 1 / BISO_VOLUME / ((AMU_ThO2 / AMU_Th232)         / KERNEL_VOLUME / DENSITY_ThO2 / 100**3 * 1000)

print(f"Loading limits [kg/m³ of nominal breeder]:")
print(f"  FLiBe + U238  (UF4):   {limit_flibe_u:.0f}")
print(f"  FLiBe + Th232 (ThF4):  {limit_flibe_th:.0f}")
print(f"  DCLL/HCPB + U238  (BISO):  {limit_biso_u:.0f}")
print(f"  DCLL/HCPB + Th232 (BISO):  {limit_biso_th:.0f}")


# --- Plot ---
fig, ax = plt.subplots(figsize=(3.5, 3.0))

ax.plot(fertile_kgm3, vf_biso_th,  label=r'HCPB/DCLL-ThO$_2$ BISO',  color="#000000", linestyle=LONG_DASH, linewidth=0.75)
ax.plot(fertile_kgm3, vf_biso_u,   label=r'HCPB/DCLL-UO$_2$ BISO',   color="#000000", linestyle='-',       linewidth=0.75)
ax.plot(fertile_kgm3, vf_flibe_th, label=r'FLiBe-ThF$_4$', color='#66b420', linestyle=LONG_DASH, linewidth=0.75)
ax.plot(fertile_kgm3, vf_flibe_u,  label=r'FLiBe-UF$_4$',  color='#66b420', linestyle='-',       linewidth=0.75)
# ax.plot(fertile_kgm3, vf_biso_u,   label=r'DCLL-UO$_2$',   color='#0047ba', linestyle='-',       linewidth=0.75)
# ax.plot(fertile_kgm3, vf_biso_th,  label=r'DCLL-ThO$_2$',  color='#0047ba', linestyle=LONG_DASH, linewidth=0.75)

ax.set_xlabel(r'Fertile isotope density [kg$/$m³]')
ax.set_ylabel('Fertile material fraction [vol% breeding volume]')

x_buf = 0.03 * 5500
y_buf = 0.03 * 100
ax.set_xlim(-x_buf, 5500 + x_buf)
ax.set_ylim(-y_buf, 100 + y_buf)

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(250))

ax.yaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(5))

ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)

plt.tight_layout()
leg = plt.legend(loc='lower right', fancybox=False, edgecolor='none',
                 frameon=True, framealpha=0.0, ncol=2, fontsize=6.5)
leg.get_frame().set_linewidth(0.5)

output_path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
output_name = 'loadings'
plt.savefig(f'{output_path}/{output_name}.png', dpi=200, bbox_inches='tight', pad_inches=0.01)
plt.savefig(f'{output_path}/{output_name}.pdf', dpi=200, bbox_inches='tight', pad_inches=0.01)
print(f"Saved {output_name}.png/pdf")
plt.show()
