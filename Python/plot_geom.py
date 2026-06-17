"""
Tokamak vertical cross-sections for FLiBe, DCLL, and HCPB blankets.
Run from project root:  python Python/plot_geom.py
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPoly, Patch
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from Python.parameters import *


# ── Colors ────────────────────────────────────────────────────────

C_PLASMA    = '#FFF8E1'
C_FW        = '#616161'
C_BE        = '#A5D6A7'
C_STRUCTURE = '#B0BEC5'
C_BREEDING  = '#FFCC80'
C_COOLANT   = '#90CAF9'
C_SHIELD    = '#78909C'
EC          = '#333333'
LW          = 0.4


# ── Geometry helpers ─────────────────────────────────────────────

def miller(R0, a, kappa, delta, extrude=0, n=500):
    t = np.linspace(0, 2 * np.pi, n)
    R = R0 + a * np.cos(t + delta * np.sin(t))
    Z = kappa * a * np.sin(t)
    if extrude:
        dR = -a * np.sin(t + delta * np.sin(t)) * (1 + delta * np.cos(t))
        dZ = kappa * a * np.cos(t)
        mag = np.hypot(dR, dZ)
        R += extrude * dZ / mag
        Z -= extrude * dR / mag
    return R, Z


def split_at_R0(R, Z, R0):
    """Split closed contour into inboard (R < R0) and outboard (R > R0) halves."""
    crossings = []
    n = len(R)
    for i in range(n - 1):
        if (R[i] - R0) * (R[i + 1] - R0) < 0:
            f = (R0 - R[i]) / (R[i + 1] - R[i])
            crossings.append((i, Z[i] + f * (Z[i + 1] - Z[i])))
    (i1, Z1), (i2, Z2) = crossings
    R_in  = np.concatenate([[R0], R[i1 + 1 : i2 + 1], [R0]])
    Z_in  = np.concatenate([[Z1], Z[i1 + 1 : i2 + 1], [Z2]])
    R_out = np.concatenate([[R0], R[i2 + 1 : n - 1], R[: i1 + 1], [R0]])
    Z_out = np.concatenate([[Z2], Z[i2 + 1 : n - 1], Z[: i1 + 1], [Z1]])
    return (R_in, Z_in), (R_out, Z_out)


def add_ring(ax, Ro, Zo, Ri, Zi, color):
    v = np.column_stack([np.concatenate([Ro, Ri[::-1]]),
                         np.concatenate([Zo, Zi[::-1]])])
    ax.add_patch(MplPoly(v, closed=True, fc=color, ec='none'))


def add_fill(ax, R, Z, color):
    ax.add_patch(MplPoly(np.column_stack([R, Z]), closed=True,
                         fc=color, ec='none'))


def add_line(ax, R, Z):
    ax.plot(R, Z, color=EC, lw=LW)


# ── FLiBe ─────────────────────────────────────────────────────────

def draw_flibe(ax):
    R0, a, kp, dl = FLIBE_R0, FLIBE_A, FLIBE_KAPPA, FLIBE_DELTA
    thk = [FLIBE_FW_CM, FLIBE_BE_CM, FLIBE_ST1_CM, FLIBE_BR1_CM,
           FLIBE_ST2_CM, FLIBE_BR2_CM, FLIBE_ST3_CM]
    d = np.concatenate([[0], np.cumsum(thk)])
    c = [miller(R0, a, kp, dl, e) for e in d]
    col = [C_FW, C_BE, C_STRUCTURE, C_BREEDING,
           C_STRUCTURE, C_BREEDING, C_STRUCTURE]

    add_fill(ax, *c[0], C_PLASMA)
    for i, clr in enumerate(col):
        add_ring(ax, *c[i + 1], *c[i], clr)
    for ci in c:
        add_line(ax, *ci)

    ax.set_title('FLiBe', fontweight='bold')
    return c[-1]


# ── DCLL ──────────────────────────────────────────────────────────

def draw_dcll(ax):
    R0, a, kp, dl = DCLL_R0, DCLL_A, DCLL_KAPPA, DCLL_DELTA
    M = lambda d: miller(R0, a, kp, dl, d)

    # Cumulative extrusion – shared layers
    d_fw  = DCLL_FW_CM
    d_fwf = d_fw  + DCLL_FWF_CM
    d_fwc = d_fwf + DCLL_FWC_CM
    d_fwb = d_fwc + DCLL_FWB_CM
    d_br1 = d_fwb + DCLL_BR1_CM
    d_d1  = d_br1 + DCLL_D1_CM
    d_br2 = d_d1  + DCLL_BR2_CM

    # Inboard only
    d_im_i = d_br2 + DCLL_IM_CM
    d_bp_i = d_im_i + DCLL_BP_CM
    d_ss_i = d_bp_i + DCLL_SS_CM

    # Outboard only
    d_d2   = d_br2  + DCLL_D2_CM
    d_br3  = d_d2   + DCLL_BR3_CM
    d_im_o = d_br3  + DCLL_IM_CM
    d_bp_o = d_im_o + DCLL_BP_CM
    d_ss_o = d_bp_o + DCLL_SS_CM

    c_vc  = M(0);      c_fw  = M(d_fw);   c_fwf = M(d_fwf)
    c_fwc = M(d_fwc);  c_fwb = M(d_fwb);  c_br1 = M(d_br1)
    c_d1  = M(d_d1);   c_br2 = M(d_br2)

    c_im_i = M(d_im_i); c_bp_i = M(d_bp_i); c_ss_i = M(d_ss_i)
    c_d2   = M(d_d2);   c_br3  = M(d_br3)
    c_im_o = M(d_im_o); c_bp_o = M(d_bp_o); c_ss_o = M(d_ss_o)

    # Symmetric layers
    add_fill(ax, *c_vc, C_PLASMA)
    add_ring(ax, *c_fw,  *c_vc,  C_FW)
    add_ring(ax, *c_fwf, *c_fw,  C_STRUCTURE)
    add_ring(ax, *c_fwc, *c_fwf, C_COOLANT)
    add_ring(ax, *c_fwb, *c_fwc, C_STRUCTURE)
    add_ring(ax, *c_br1, *c_fwb, C_BREEDING)
    add_ring(ax, *c_d1,  *c_br1, C_STRUCTURE)
    add_ring(ax, *c_br2, *c_d1,  C_BREEDING)

    # Inboard (R < R0)
    I = lambda c: split_at_R0(*c, R0)[0]
    ib = [I(c_br2), I(c_im_i), I(c_bp_i), I(c_ss_i)]
    ib_col = [C_STRUCTURE, C_STRUCTURE, C_SHIELD]
    for j, clr in enumerate(ib_col):
        add_ring(ax, *ib[j + 1], *ib[j], clr)

    # Outboard (R > R0)
    O = lambda c: split_at_R0(*c, R0)[1]
    ob = [O(c_br2), O(c_d2), O(c_br3), O(c_im_o), O(c_bp_o), O(c_ss_o)]
    ob_col = [C_STRUCTURE, C_BREEDING, C_STRUCTURE, C_STRUCTURE, C_SHIELD]
    for j, clr in enumerate(ob_col):
        add_ring(ax, *ob[j + 1], *ob[j], clr)

    # Boundary lines – shared
    for ci in [c_vc, c_fw, c_fwf, c_fwc, c_fwb, c_br1, c_d1, c_br2]:
        add_line(ax, *ci)
    # Boundary lines – inboard / outboard halves
    for half in ib:
        add_line(ax, *half)
    for half in ob:
        add_line(ax, *half)

    ax.set_title('DCLL', fontweight='bold')
    return c_ss_o


# ── HCPB ──────────────────────────────────────────────────────────

def draw_hcpb(ax):
    R0, a, kp, dl = HCPB_R0, HCPB_A, HCPB_KAPPA, HCPB_DELTA
    M = lambda d: miller(R0, a, kp, dl, d)

    d_fw    = HCPB_FW_CM
    d_st1   = d_fw + HCPB_ST1_CM
    d_br_i  = d_st1 + HCPB_BR1_I_CM
    d_st2_i = d_br_i + HCPB_ST2_CM
    d_br_o  = d_st1 + HCPB_BR1_O_CM
    d_st2_o = d_br_o + HCPB_ST2_CM

    c_vc = M(0);  c_fw = M(d_fw);  c_st1 = M(d_st1)
    c_br_i = M(d_br_i);  c_st2_i = M(d_st2_i)
    c_br_o = M(d_br_o);  c_st2_o = M(d_st2_o)

    # Symmetric layers
    add_fill(ax, *c_vc, C_PLASMA)
    add_ring(ax, *c_fw,  *c_vc,  C_FW)
    add_ring(ax, *c_st1, *c_fw,  C_STRUCTURE)

    # Inboard (R < R0)
    I = lambda c: split_at_R0(*c, R0)[0]
    ist1, ibr, ist2 = I(c_st1), I(c_br_i), I(c_st2_i)
    add_ring(ax, *ibr,  *ist1, C_BREEDING)
    add_ring(ax, *ist2, *ibr,  C_STRUCTURE)

    # Outboard (R > R0)
    O = lambda c: split_at_R0(*c, R0)[1]
    ost1, obr, ost2 = O(c_st1), O(c_br_o), O(c_st2_o)
    add_ring(ax, *obr,  *ost1, C_BREEDING)
    add_ring(ax, *ost2, *obr,  C_STRUCTURE)

    # Boundary lines – shared
    for ci in [c_vc, c_fw, c_st1]:
        add_line(ax, *ci)
    # Boundary lines – inboard / outboard halves
    for half in [ist1, ibr, ist2]:
        add_line(ax, *half)
    for half in [ost1, obr, ost2]:
        add_line(ax, *half)

    ax.set_title('HCPB', fontweight='bold')
    return c_st2_o


# ── Main ──────────────────────────────────────────────────────────

if __name__ == '__main__':

    fig, axes = plt.subplots(1, 3, figsize=(12, 7))

    bounds = []
    for ax, fn in zip(axes, [draw_flibe, draw_dcll, draw_hcpb]):
        outer = fn(ax)
        bounds.append((outer[0].min(), outer[0].max(),
                       outer[1].min(), outer[1].max()))
        ax.set_aspect('equal')

    pad = 30
    r_lo = min(b[0] for b in bounds) - pad
    r_hi = max(b[1] for b in bounds) + pad
    z_lo = min(b[2] for b in bounds) - pad
    z_hi = max(b[3] for b in bounds) + pad

    for ax in axes:
        ax.set_xlim(r_lo, r_hi)
        ax.set_ylim(z_lo, z_hi)
        ax.set_xlabel('R (cm)')
    axes[0].set_ylabel('Z (cm)')

    legend_items = [
        Patch(fc=C_PLASMA,    ec=EC, lw=LW, label='Plasma'),
        Patch(fc=C_FW,        ec=EC, lw=LW, label='First wall (W)'),
        Patch(fc=C_BE,        ec=EC, lw=LW, label='Beryllium'),
        Patch(fc=C_STRUCTURE, ec=EC, lw=LW, label='Structure'),
        Patch(fc=C_BREEDING,  ec=EC, lw=LW, label='Breeding'),
        Patch(fc=C_COOLANT,   ec=EC, lw=LW, label='Coolant'),
        Patch(fc=C_SHIELD,    ec=EC, lw=LW, label='Shield'),
    ]
    fig.legend(handles=legend_items, loc='lower center', ncol=7,
               frameon=False, fontsize=7)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig('./Figures/fig_geom.pdf', bbox_inches='tight')
    plt.show()
