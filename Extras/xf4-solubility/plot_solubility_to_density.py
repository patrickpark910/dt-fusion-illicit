#!/usr/bin/env python3
"""
Fertile-isotope density in an XF4 / FLiBe (2LiF-BeF2) molten salt vs. XF4 mol%.

Two fertile fuels are compared:
  - UF4  : fertile isotope U-238  (U naturally enriched -> 0.71 wt% U-235, balance U-238)
  - ThF4 : fertile isotope Th-232 (Th is essentially 100% Th-232 in nature)

FLiBe lithium is enriched to 7.5 at% Li-6.
Volumes are taken as additive (ideal mixing) from the pure-component densities.

Left  subplot : kg of fertile isotope per m^3 of the 2(LiF)-BeF2 carrier ALONE.
Right subplot : kg of fertile isotope per m^3 of the TOTAL fuel salt (XF4 + FLiBe).
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import MultipleLocator

FONT_PATH = '../Wedge/fonts/arial.ttf'
fm.fontManager.addfont(FONT_PATH)
plt.rcParams['font.family'] = fm.FontProperties(fname=FONT_PATH).get_name()
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['pdf.fonttype'] = 42

# ---------------- atomic / isotopic masses (g/mol) ----------------
M_LI6, M_LI7   = 6.015123, 7.016003
M_BE, M_F      = 9.012183, 18.998403
M_U235, M_U238 = 235.043930, 238.050788
M_TH232        = 232.0381

# ---------------- FLiBe carrier: Li2BeF4, Li at 7.5 at% Li-6 ----------------
M_LI      = 0.075 * M_LI6 + 0.925 * M_LI7
M_FLIBE   = 2 * M_LI + M_BE + 4 * M_F
RHO_FLIBE = 1.95  # g/cm^3

# ---------------- uranium: 0.71 wt% U-235, balance U-238 ----------------
W_U235 = 0.0071
W_U238 = 1.0 - W_U235
M_U    = 1.0 / (W_U235 / M_U235 + W_U238 / M_U238)   # mass-avg molar mass
M_UF4  = M_U + 4 * M_F

# ---------------- thorium: natural ~ pure Th-232 ----------------
M_TH   = M_TH232
M_THF4 = M_TH + 4 * M_F

# Per-fuel data: grams of fertile isotope per mole of XF4, salt molar mass, salt density
FUELS = {
    "UF$_4$  (U-238)": dict(
        m_fertile_per_mol=M_U * W_U238,   # g of U-238 per mol UF4
        M_XF4=M_UF4,
        rho_XF4=6.88,
        color="tab:blue",
        ref_molpct=23.0,
    ),
    "ThF$_4$  (Th-232)": dict(
        m_fertile_per_mol=M_TH * 1.0,     # g of Th-232 per mol ThF4
        M_XF4=M_THF4,
        rho_XF4=6.60,
        color="tab:green",
        ref_molpct=17.0,
    ),
}


def carrier_density(x, p):
    """kg of fertile isotope per m^3 of FLiBe carrier.  x = mole fraction of XF4."""
    mass_fertile = x * p["m_fertile_per_mol"]            # g per 1 mol of mixture
    v_flibe      = (1.0 - x) * M_FLIBE / RHO_FLIBE        # cm^3 of carrier
    return mass_fertile / v_flibe * 1000.0               # g/cm^3 -> kg/m^3


def total_density(x, p):
    """kg of fertile isotope per m^3 of TOTAL salt.  x = mole fraction of XF4."""
    mass_fertile = x * p["m_fertile_per_mol"]
    v_flibe      = (1.0 - x) * M_FLIBE / RHO_FLIBE
    v_xf4        = x * p["M_XF4"] / p["rho_XF4"]
    return mass_fertile / (v_flibe + v_xf4) * 1000.0


# ---------------- x-axis: 0 -> 40 mol% XF4 ----------------
molpct = np.linspace(0.0, 40.0, 400)
x = molpct / 100.0

plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8,
                     'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'axes.linewidth': 0.25,
                     'xtick.direction': 'in', 'ytick.direction': 'in',
                     'xtick.major.width': 0.25, 'ytick.major.width': 0.25,
                     'xtick.minor.width': 0.25, 'ytick.minor.width': 0.25,
                     'grid.color': '#DBDBDB'})

fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.5, 3.0), sharex=True)

for name, p in FUELS.items():
    axL.plot(molpct, carrier_density(x, p), color=p["color"], lw=0.75, label=name)
    axR.plot(molpct, total_density(x, p),   color=p["color"], lw=0.75, label=name)

    xr = p["ref_molpct"] / 100.0
    axL.scatter([p["ref_molpct"]], [carrier_density(xr, p)],
                s=15, color=p["color"], edgecolor="k", linewidths=0.25, zorder=5)
    axR.scatter([p["ref_molpct"]], [total_density(xr, p)],
                s=15, color=p["color"], edgecolor="k", linewidths=0.25, zorder=5)

axL.text(0.97, 0.97, r'Per m$^3$ carrier', transform=axL.transAxes,
         fontweight='bold', fontsize=10, ha='right', va='top')
axR.text(0.97, 0.97, r'Per m$^3$ total salt', transform=axR.transAxes,
         fontweight='bold', fontsize=10, ha='right', va='top')

for ax in (axL, axR):
    ax.set_xlabel(r'XF$_4$ content [mol%]')
    ax.set_ylabel(r'Fertile-isotope density [kg/m$^3$]')
    ax.set_xlim(0, 40)
    ax.set_ylim(bottom=0)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(50))
    ax.tick_params(which='both', top=True, right=True)
    ax.grid(axis='x', which='major', linestyle='-', linewidth=0.5)
    ax.grid(axis='y', which='major', linestyle='-', linewidth=0.5)
    ax.legend(loc='lower right', frameon=False, fontsize=8)

plt.tight_layout(pad=0.3)
fig.subplots_adjust(wspace=0.225)
plt.savefig("fertile_density.png", dpi=150, bbox_inches="tight")
plt.savefig("fertile_density.pdf", bbox_inches='tight', pad_inches=0.02)

# ---------------- quick double-check against the worked answers ----------------
print("Reference points (should match the hand calculations):")
for name, p in FUELS.items():
    xr = p["ref_molpct"] / 100.0
    print(f"  {name:18s} at {p['ref_molpct']:4.0f} mol%:  "
          f"carrier = {carrier_density(xr, p):7.1f} kg/m^3   "
          f"total = {total_density(xr, p):7.1f} kg/m^3")

plt.show()