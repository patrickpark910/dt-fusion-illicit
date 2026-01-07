import numpy as np
import pandas as pd
from jinja2 import Template
import os

from wedge_B_pack import *

# GLOBAL PARAMETERS
TEMP_K = 900
BISO_VOLUME = 4/3*np.pi*(0.05)**3
KERNEL_VOLUME = 4/3*np.pi*(0.04)**3

DENSITY_UO2  = 10.50 # [g/cm³]
DENSITY_ThO2 = 10.00 # [g/cm³]

AMU_U235 = 235.0439299
AMU_U238 = 238.05078826
AMU_U = 238.02891 # for natural enrichment
AMU_O = 15.999
AMU_UO2 = AMU_U + 2 * AMU_O
AMU_Th232 = 232.0381

# Volume fractions in HCPB blanket
VF_LI_NOM = 0.1304
VF_BE_NOM = 0.3790
VF_EU_NOM = 0.1176
VF_HE_NOM = 0.3730


def fertile_kgm3_to_biso_per_cc(fertile_kgm3, fertile_isotope='U238'):
    # vol_kernel = KERNEL_VOLUME # 4/3 * np.pi * 0.04**3
    if fertile_isotope == 'U238':
        biso_per_cc = fertile_kgm3 * AMU_UO2 / AMU_U238 / KERNEL_VOLUME / DENSITY_UO2 / 100**3 * 1000
    elif fertile_isotope == 'Th232':
        biso_per_cc = fertile_kgm3 * AMU_ThO2 / AMU_Th232 / KERNEL_VOLUME / DENSITY_ThO2 / 100**3 * 1000
    return biso_per_cc


def main(fertile_kgm3, go_pack=False, go_print=False):
    """
    Performs geometric scaling, particle packing, and MCNP input generation 
    for a "unit prism" representing a rectangular borehole through a HCPB blanket.

    Args:
        fertile_kgm3 (float): The mass of U-238 isotope (kg) per cubic meter of 
            breeding material (Li4SiO4 + Be).
        go_pack (bool, optional): If True, runs the packing algorithm to generate 
            new (x,y,z) coordinates for BISO particles and saves them to CSV. 
            Defaults to False.
        go_print (bool, optional): If True, renders the Jinja2 template and writes 
            the final MCNP .inp file to the ./inputs/ directory. Defaults to False.

    Returns:
        None: The function prints a summary of derived quantities to the console 
            and writes files to disk if the boolean flags are enabled.

    Notes:
        - Function scales the inboard (b1), outboard (b2) face lengths such that 
          their face areas A1, A2 maintain a specific fluence ratio (f1=0.39, f2=0.61) 
          observed in the "full" D-shape HCPB OpenMC model.
        - "Fertile kg/m3" refers to kg of U-238 isotope per m3 of breeding material
        - "Breeding material" refers to Li4SiO4 and Be but NOT Eurofer or He-4(g)
        - "Breeding region" refers to the physical volume containing the breeding
          material AND Eurofer and He-4(g)
        - BISO counts (N1=512, N2=1500) are currently hardcoded to facilitate 3-dimensional 
          prime factorization for lattice packing for Prism C.
        - Expects a directory structure containing `./csv_coords/`, 
          `./inputs/`, and a template file named `wedge_{case}.template`.

    Example:
        >>> main(500.0, go_pack=True, go_print=True)
        # Calculates dimensions for 500 kg/m^3, packs 2012 spheres, and prints .inp
    """

    """ BASIC PARAMETERS
    """
    case = 'B'
    density_str = f"{fertile_kgm3:06.2f}"
    
    # Number of particles for MCNP run
    nps = '5e6' 

    target_total_biso = 2012
    N1, N2 =  512, 1500

    # HCPB tokamak dimensions (along z in model)
    h1, h2 = 45.0, 82.0     # inboard, outboard thicknesses
    d1 = 2.5 + 0.2          # Eurofer + first wall 
    d0 = 400.0              # plasma chamber
    oo = h1 + d0 + 2 * d1   # outboard offset (where outboard starts in z)
    f1, f2 = 0.39, 0.61     # fraction of fluence into inboard vs. outboard in full OpenMC HCPB model


    ''' 
    Nominally, the volume fracs are: Li4SiO4 = 13.04%, Be = 37.90%, Eurofer = 11.76%, He = 37.30%
    We want to deduct the BISO volume we add from the Li4SiO4.

    Here I use: "breeder" = Li4SiO4, Be
                "breeding region" (br) = the physical blanket layer that contains the breeder
                                         this layer usually contains breeder + structure + coolant
                "breeding volume" (bv) = amount of volume nominally occupied in blanket by breeder

    That is, if the breeding region has volume 100 m³ of which 60% is (Li4SiO4 + Be) and the other 
    40% is (Eurofer + Be), then 60 m³ is the breeding volume.
        
    Given fertile_kgm3 [kg U-238 / m³ breeder] we convert this to m³ of BISO.

    **Be careful NOT to renormalize everything. DEDUCT the volume fraction of BISO from 1, and then
    set THAT as the volume fraction of breeder, which then subdivides into Li4SiO4 and Be in a 
    13.04:37.90 ratio.

    We use the fact that 13.04 vol% Li4SiO4 = 0.1304 m³ of Li4SiO4 per m³ of blanket to do our work 
    below in terms of volume fractions.
    '''

    # Nominal fraction the breeding volume is in the breeding region
    VF_BV_BR_NOM = VF_LI_NOM + VF_BE_NOM

    # Nominal volume fractions of Li4SiO4 vs. Be in the breeding volume (bv)
    vf_li_bv_nom = VF_LI_NOM / (VF_LI_NOM + VF_BE_NOM)
    vf_be_bv_nom = VF_BE_NOM / (VF_LI_NOM + VF_BE_NOM)

    # Number of BISO spheres per cc of breeding volume
    biso_per_cc_bv = fertile_kgm3_to_biso_per_cc(fertile_kgm3)
    vf_biso_bv     = biso_per_cc_bv * BISO_VOLUME

    if vf_biso_bv > 1.0:
        print(f"Fatal. Your fertile kg/m³ exceeds what can physically fit in the breeding volume!")
        print(f"Fatal. That is, your volume of BISO per cm³ of breeding volume exceeds 1.")

    # New volume ratios of everything w.r.t. breeders in the breeder volume
    vf_libe_bv   = 1 - vf_biso_bv
    vf_li_bv     = vf_li_bv_nom * vf_libe_bv
    vf_be_bv     = vf_be_bv_nom * vf_libe_bv
    vf_checksum1 = vf_biso_bv + vf_li_bv + vf_be_bv  # should equal 1

    # New volume ratios of everyting w.r.t. everything else in the breeding region
    vf_biso_br = vf_biso_bv * (VF_LI_NOM + VF_BE_NOM)
    vf_li_br   = vf_li_bv * (VF_LI_NOM + VF_BE_NOM)
    vf_be_br   = vf_be_bv * (VF_LI_NOM + VF_BE_NOM)
    # vol frac of Eurofer and He-4 doesn't change
    vf_checksum2 = vf_biso_br + vf_li_br + vf_be_br + VF_EU_NOM + VF_HE_NOM  # should equal 1

    # Number of BISO spheres per cc of breeding region
    biso_per_cc_br = vf_biso_br / BISO_VOLUME

    '''
    NB. Maximum possible value of vf_biso_br is (VF_LI_NOM + VF_BE_NOM) = 50.94%. 
    This means that at most 50.94% of the breeding region can be occupied by BISO
    Per 1 cm³ of breeding region, at most 0.5094 cm³ can be BISO, or 973 spheres/cm³.
    Practically, 400 spheres/cm³ = 1000 kg/m³ so the max possible fertile kg/m³ is 2432.50
    '''

    A_ref = target_total_biso / biso_per_cc_br / (h1+h2) 
    A2 = A_ref * (h1+h2) / ( f1/f2 * h1 + h2)                    # outboard XY face area 
    A1 = A2 * f1/f2                                              # inboard XY face area
    b1, b2 = np.round(np.sqrt(A1), 6), np.round(np.sqrt(A2), 6)  # face lengths
    c1, c2 = b1/2, b2/2                                # face HALF-lengths

    # Volumes of the breeding regions (i.e. area*height)
    V1, V2 = A1*h1, A2*h2


    # Check
    summary = (
            f"c {78*'='}"
            f"\nc Case: {density_str} kg/cm³"
            f"\nc "
            f"\nc Dimensions (inboard, outboard):"
            f"\nc  XY face half-lengths c1, c2: {c1:.6f}, {c2:.6f} [cm]"
            f"\nc  XY face lengths      b1, b2: {b1:.6f}, {b2:.6f} [cm]"
            f"\nc  XY face areas        A1, A2: {A1:.6f}, {A2:.6f} [cm²]"
            f"\nc  Breeding region vols V1, V2: {V1:.6f}, {V2:.6f} [cm³]"
            f"\nc  Fraction of fluence  f1, f2: {f1:.6f}, {f2:.6f} "
            f"\nc "
            f"\nc 'breeder' = Li4SiO4, Be "
            f"\nc 'breeding region' (br) = the physical blanket layer that contains the breeder "
            f"\nc                          this layer usually contains breeder + structure + coolant "
            f"\nc 'breeding volume' (bv) = volume nominally occupied in breeding region by breeder "
            f"\nc "
            f"\nc With respect to the BREEDER material, i.e., per 1 m³ of breeder volume, we have these new volume fractions:"
            f"\nc   vf_biso_breeder_new =  {(vf_biso_bv*100):.6f} vol%"
            f"\nc   vf_li_breeder_new   =  {(vf_li_bv*100):.6f} vol%"
            f"\nc   vf_be_breeder_new   =  {(vf_be_bv*100):.6f} vol%"
            f"\nc   check they add up   = {(vf_checksum1*100):.6f} vol%"
            f"\nc "
            f"\nc With respect to the whole BREEDING REGION, we have these new volume fractions:"
            f"\nc   vf_biso_br =  {(vf_biso_br*100):.6f} vol%"
            f"\nc   vf_li_br   =  {(vf_li_br*100):.6f} vol%"
            f"\nc   vf_be_br   =  {(vf_be_br*100):.6f} vol%"
            f"\nc   VF_EU_NOM  =  {(VF_EU_NOM*100):.6f} vol%"
            f"\nc   VF_HE_NOM  =  {(VF_HE_NOM*100):.6f} vol%"
            f"\nc   check they add up   = {(vf_checksum2*100):.6f} vol%"
            f"\nc   check that BISO + Li4SiO4 + Be adds up to the nominal Li4SiO4 + Be fraction"
            f"\nc   vf_biso_br + vf_li_br + vf_be_br = {((vf_biso_br+vf_li_br+vf_be_br))*100:.6f} : VF_LI_NOM + VF_BE_NOM = {((VF_LI_NOM+VF_BE_NOM))*100:.6f}"
            f"\nc "
            f"\nc Check BISO volume fraction wrt breeding region is correct:"
            f"\nc   Inboard : N1 * BISO_VOLUME / (b1**2 * h1) = {(N1 * BISO_VOLUME / (b1**2 * h1) * 100):.6f} vol%"
            f"\nc   Outboard: N2 * BISO_VOLUME / (b2**2 * h2) = {(N2 * BISO_VOLUME / (b2**2 * h2) * 100):.6f} vol%"
            f"\nc   Total   : (N1+N2) * BISO_VOLUME / (b1**2 * h1 + b2**2 * h2) = {((N1+N2) * BISO_VOLUME / (b1**2 * h1 + b2**2 * h2) * 100):.6f} vol%"
            f"\nc   ...should match                                  vf_biso_br = {(vf_biso_br*100):.6f} vol%"
            f"\nc "
            f"\nc BISO spheres per cc of breeding volume: {biso_per_cc_bv:.6f} spheres/cm³"
            f"\nc              per cc of breeding region: {biso_per_cc_br:.6f} spheres/cm³"
            f"\nc "
            )

    print(summary)

    # Volumes of the breeding region MINUS biso volumes (for MCNP cell vol card)
    V1_br = A1*h1 - N1 * BISO_VOLUME
    V2_br = A2*h2 - N2 * BISO_VOLUME
    V_br  = V1_br + V2_br


    """ PERFORM RANDOM SEQUENTIAL PACKING
    """
    inboard_csv = f"./csv_coords/case{case}_inboard_{density_str}.csv"
    outboard_csv = f"./csv_coords/case{case}_outboard_{density_str}.csv"

    if go_pack:
        biso_inboard  = pack(b1, h1, N1, offset=0,  r=0.05, name='biso')
        biso_inboard.to_csv(f'./csv_coords/case{case}_inboard_{density_str}.csv', index=False)

        biso_outboard = pack(b2, h2, N2, offset=oo, r=0.05, name='biso')
        biso_outboard.to_csv(f'./csv_coords/case{case}_outboard_{density_str}.csv', index=False)


    """ OTHER MCNP SYNTAX GENERATION
    """
    # Calculate surfaces for trapezoidal plasma chamber
    # NB. MCNP plane: 1.0*x + 0.0*y + plane_C*z - plane_D = 0
    z_in_end    = 47.7
    z_out_start = 447.7
    dz = z_out_start - z_in_end
    dx = c2 - c1
    slope = dx / dz
    plane_C = -slope
    plane_D = c1 - (slope * z_in_end)

    # Process CSV coordinates
    df_inboard  = pd.read_csv(inboard_csv)
    df_outboard = pd.read_csv(outboard_csv)

    # Create universe cards and complement strings
    in_univ, out_univ = [], []
    in_comp, out_comp = [], []
    in_tally, out_tally = [], []

    n_in = len(df_inboard)
    n_out = len(df_outboard)
    total_count = n_in + n_out
    total_biso_vol = total_count * BISO_VOLUME
    total_uo2_vol  = total_count * KERNEL_VOLUME

    # Inboard
    for idx, row in df_inboard.iterrows():
        cid = 10000 + idx + 1
        in_univ.append(f"{cid:5d}  0  -151  trcl=({row['x']:.6e} {row['y']:.6e} {row['z']:.6e})  imp:n=1  tmp=7.76e-8  vol={BISO_VOLUME:.3e}  fill=50")
        in_comp.append(f"#{cid}")
        in_tally.append(f"{cid}")

    # Outboard
    for idx, row in df_outboard.iterrows():
        cid = 20000 + idx + 1
        out_univ.append(f"{cid:5d}  0  -151  trcl=({row['x']:.6e} {row['y']:.6e} {row['z']:.6e})  imp:n=1  tmp=7.76e-8  vol={BISO_VOLUME:.3e}  fill=50")
        out_comp.append(f"#{cid}")
        out_tally.append(f"{cid}")

    # Format complement strings (chunks of N for MCNP character limit)
    def chunk(comp_list, N=15):
        return "\n        ".join([" ".join(comp_list[i:i+N]) for i in range(0, len(comp_list), N)])


    """
    Jinja2 template insert
    """
    context = {
        "case": case,
        "fertile_kgm3": fertile_kgm3,
        "density_str" : density_str,
        "summary": summary,
        "b1": b1,
        "b2": b2,
        "c1": c1,
        "c2": c2,
        "V1_fw": f"{(b1**2 * 0.2):.6e}",
        "V1_ei": f"{(b1**2 * 2.5):.6e}",
        "V1_eo": f"{(b1**2 * 2.5):.6e}",
        "V1_br": f"{V1_br:.6e}",
        "V2_fw": f"{(b2**2 * 0.2):.6e}",
        "V2_ei": f"{(b2**2 * 2.5):.6e}",
        "V2_eo": f"{(b2**2 * 2.5):.6e}",
        "V2_br": f"{V2_br:.6e}",
        "V_uo2": f"{KERNEL_VOLUME:.6e}",
        "V_sic": f"{(BISO_VOLUME-KERNEL_VOLUME):.6e}",
        "cells_inboard_complement": chunk(in_comp, N=15),
        "cells_outboard_complement": chunk(out_comp, N=15),
        "inboard_universes": "\n".join(in_univ),
        "outboard_universes": "\n".join(out_univ),
        "plane_C": f"{plane_C:.6e}",
        "plane_D": f"{plane_D:.6e}",
        "cells_inboard_tally": chunk(in_tally, N=15),
        "cells_outboard_tally": chunk(out_tally, N=15),
        "V_tally_u238": f"{total_uo2_vol:.6e}",
        "nps": nps
    }


    """ PRINT TO TEMPLATE
    """
    if go_print:
        with open(f'wedge_{case}.template', 'r', encoding='utf-8') as f:
            template_content = f.read()
        
        template = Template(template_content)
        rendered_mcnp = template.render(context)
        filepath = f'./inputs/wedge_{case}_{density_str}kgm3.inp'

        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(rendered_mcnp)
        
        print(f"Successfully generated: {filepath}\n")
    else:
        print(f"Template was NOT printed\n")

    

if __name__ == "__main__":
    for fertile_kgm3 in [0.10, 0.50, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]:
        main(fertile_kgm3, go_pack=False, go_print=True)