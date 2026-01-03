import numpy as np
import pandas as pd
from jinja2 import Template
import os

from prism_B_pack import *


def fertile_to_biso_density(fertile_kgm3, isotope='U238'):
    vol_kernel = 4/3 * np.pi * 0.04**3
    if isotope == 'U238':
        biso_per_cc = fertile_kgm3 * (238.02+32)/238.05 / vol_kernel / 10.5 / 100**3 * 1000
    return biso_per_cc


def u238_at(r=0.04):
    AMU_U238 = 238.05078826
    AMU_U = 238.02891 # for natural enrichment
    AMU_O = 15.999
    AMU_UO2 = AMU_U + 2 * AMU_O
    # UO2 density 10.5 g/cm3
    mass = 10.50 * (4/3 * np.pi * r**3)
    # assume 0.71 wt% enriched
    u238_mass = mass * AMU_U/AMU_UO2 * (1-0.71)
    # number of atoms = (mass/molar_mass) * Avogadro
    u238_at = u238_mass / AMU_U238 * 6.022e23
    return u238_at


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
          `./inputs/`, and a template file named `prism_{case}.template`.

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

    target_total_biso = 2012  # >> DO NOT CHANGE WITHOUT CHANGING IT IN prism_{A,B,C}.py!! AND CHANGE N1, N2 (WHICH ARE HARDCODED FOR NOW) <<
    
    # Real number of BISO particles in inboard vs. outboard
    N1 =  512  # round(V1 * biso_per_cc_region) # =  523 #  512  # 
    N2 = 1500  # round(V2 * biso_per_cc_region) # = 1489 # 1500  # 
    N_tot = N1 + N2
    # NB. I temporarily hardcoded 512, 1500 to match & make lattice packing easier in Prism C (3-D prime factorization)

    # HCPB tokamak dimensions (along z in model)
    h1, h2 = 45.0, 82.0     # inboard, outboard thicknesses
    d1 = 2.5 + 0.2          # Eurofer + first wall 
    d0 = 400.0              # plasma chamber
    oo = h1 + d0 + 2 * d1   # outboard offset (where outboard starts in z)
    f1, f2 = 0.39, 0.61     # fraction of fluence into inboard vs. outboard in full OpenMC HCPB model

    # Volume of biso sphere and UO2 kernel
    V_BISO   = (4/3) * np.pi * (0.05)**3
    V_KERNEL = (4/3) * np.pi * (0.04)**3

    # Volume fractions of material in breeding regions of HCPB blanket (Lu et al. 2017)
    VF_LI_NOM = 0.1304
    VF_BE_NOM = 0.3790
    VF_EU_NOM = 0.1176
    VF_HE_NOM = 0.3730


    """ 
    DERIVED QUANTITIES
    """

    ''' 
    Nominally, the volume fracs are: Li4SiO4 = 13.04%, Be = 37.90%, Eurofer = 11.76%, He = 37.30%
    We want to deduct the BISO volume we add from the Li4SiO4.

    Here I use: "blanket" = Li4SiO4, Be, Eurofer, He
                "breeder" = Li4SiO4, Be
        
    Given fertile_kgm3 [kg U-238 / m³ breeder] we convert this to m³ of BISO and renormalize everything.
    We use the fact that 13.04 vol% Li4SiO4 = 0.1304 m³ of Li4SiO4 per m³ of blanket to do our work 
    below in terms of volume fractions.
    '''

    biso_per_cc_breeder = fertile_to_biso_density(fertile_kgm3)

    vf_biso_breeder = biso_per_cc_breeder * V_BISO
    vf_li_breeder_nom = VF_LI_NOM / (VF_LI_NOM + VF_BE_NOM)
    vf_be_breeder_nom = VF_BE_NOM / (VF_LI_NOM + VF_BE_NOM)

    # Renormalize everything
    vf_li_breeder_new   = vf_li_breeder_nom / (vf_biso_breeder + vf_li_breeder_nom + vf_be_breeder_nom)
    vf_be_breeder_new   = vf_be_breeder_nom / (vf_biso_breeder + vf_li_breeder_nom + vf_be_breeder_nom)
    vf_biso_breeder_new = vf_biso_breeder   / (vf_biso_breeder + vf_li_breeder_nom + vf_be_breeder_nom)
    vf_checksum1 = vf_biso_breeder_new+vf_li_breeder_new+vf_be_breeder_new # should equal 1

    vf_biso = vf_biso_breeder_new * (VF_LI_NOM + VF_BE_NOM)
    vf_li   = vf_li_breeder_new   * (VF_LI_NOM + VF_BE_NOM)
    vf_be   = vf_be_breeder_new   * (VF_LI_NOM + VF_BE_NOM)
    vf_checksum2 = vf_biso+vf_li+vf_be+VF_EU_NOM+VF_HE_NOM # should equal 1

    biso_per_cc_blanket = vf_biso / V_BISO


    '''
    We derive the nominal face area 'A_ref' above:
        A_ref [cm²] = target_total_biso [spheres] * biso_per_cc_region⁻¹ [cm³/spheres] / total breeder length [cm]
    We now need to scale A_ref so that:
        1. ratio f = A1/A2 = f1/f2 = nominally 0.39/0.61
        2. volume A1*h1 + A2*h2 = A_ref*(h1+h2)
    We solve:
        A2 = A_ref (h1+h2) / ( f1/f2 * h1 + h2 ), A1 = A2 * f1/f2, 
    NB. Technically, the BISO should deduct from the Li4SiO4 and Be vols but NOT Eurofer or He-4(g) vols.
        BUT I do not have the time to script that kind of detail. I expect a very small error from this.
        SO I use the auto-generated 0 kg/m³ breeding region mat from OpenMC that I parse into MCNP.
    '''
    A_ref = target_total_biso / biso_per_cc_blanket / (h1+h2) 
    A2 = A_ref * (h1+h2) / ( f1/f2 * h1 + h2)                    # outboard XY face area 
    A1 = A2 * f1/f2                                              # inboard XY face area
    b1, b2 = np.round(np.sqrt(A1), 6), np.round(np.sqrt(A2), 6)  # face lengths
    c1, c2 = b1/2, b2/2                                          # face HALF-lengths

    # Volumes of the breeding regions (i.e. area*height)
    V1, V2 = A1*h1, A2*h2



    # check
    summary = (
    f"c {78*'='}"
    f"\nc Prism {case}: {fertile_kgm3} kg/cm³"
    f"\nc "
    f"\nc Face lengths b1 = {b1:.4e} cm  : b2 = {b2:.4e} cm"
    f"\nc Face areas   A1 = {A1:.4e} cm² : A2 = {A2:.4e} cm²"
    f"\nc "
    f"\nc Checks:"
    f"\nc   A1/A2 ratio : {(A1/A2):.6f} : should equal {f1}/{f2} = {(f1/f2):.6f}"
    f"\nc   A1*h1 + A2*h2 = {(A1*h1 + A2*h2):.6f} cm³ : should equal A_ref * (h1+h2) = {(A1*h1 + A2*h2):.6f} cm³"
    f"\nc "
    f"\nc Target BISO count : {target_total_biso} spheres"
    f"\nc Actual BISO count : {N_tot} spheres (in: {N1}, out: {N2})"
    f"\nc "
    f"\nc For 1 m³ of breeder, we have:"
    f"\nc   vf_biso_breeder   = {vf_biso_breeder  :.6f} m³ BISO"
    f"\nc   vf_li_breeder_nom = {vf_li_breeder_nom:.6f} m³ Li4SiO4"
    f"\nc   vf_be_breeder_nom = {vf_be_breeder_nom:.6f} m³ Be"
    f"\nc "
    f"\nc With respect to the BREEDER material, we have these new volume fractions:"
    f"\nc   vf_biso_breeder_new =  {(vf_biso_breeder_new*100):.6f} vol%"
    f"\nc   vf_li_breeder_new   =  {(vf_li_breeder_new*100):.6f} vol%"
    f"\nc   vf_be_breeder_new   =  {(vf_be_breeder_new*100):.6f} vol%"
    f"\nc   check they add up   = {(vf_checksum1*100):.6f} vol%"
    f"\nc "
    f"\nc With respect to the BLANKET material, we have these new volume fractions:"
    f"\nc   vf_biso =   {(vf_biso*100):.6f} vol%"
    f"\nc   vf_li   =  {(vf_li*100):.6f} vol%"
    f"\nc   vf_be   =  {(vf_be*100):.6f} vol%"
    f"\nc   vf_eu   =  {(VF_EU_NOM*100):.6f} vol%"
    f"\nc   vf_he   =  {(VF_HE_NOM*100):.6f} vol%"
    f"\nc   check they add up   = {(vf_checksum2*100):.6f} vol%"
    f"\nc   check that BISO + Li4SiO4 + Be adds up to the nominal Li4SiO4 + Be fraction"
    f"\nc   vf_biso + vf_li + vf_be = {((vf_biso+vf_li+vf_be))*100:.6f} : VF_LI_NOM + VF_BE_NOM = {((VF_LI_NOM+VF_BE_NOM))*100:.6f}"
    f"\nc "
    f"\nc Check BISO volume fraction wrt blanket is correct:"
    f"\nc   Inboard : N1 * V_BISO / (b1**2 * h1) = {(N1 * V_BISO / (b1**2 * h1) * 100):.6f} vol%"
    f"\nc   Outboard: N2 * V_BISO / (b2**2 * h2) = {(N2 * V_BISO / (b2**2 * h2) * 100):.6f} vol%"
    f"\nc   Total   : (N1+N2) * V_BISO / (b1**2 * h1 + b2**2 * h2) = {((N1+N2) * V_BISO / (b1**2 * h1 + b2**2 * h2) * 100):.6f} vol%"
    f"\nc   ...should match                                vf_biso = {(vf_biso*100):.6f} vol%"
    f"\nc "
    f"\nc BISO spheres per cc in breeder: {biso_per_cc_breeder:.6f} spheres/cm³ : in blanket : {biso_per_cc_blanket:.6f} spheres/cm³"
    f"\nc "
    f"\nc {78*'='}"
    )

    print(summary)

    # Volumes of the breeding region MINUS biso volumes (for MCNP cell vol card)
    V1_br = A1*h1 - target_total_biso * V_BISO
    V2_br = A2*h2 - target_total_biso * V_BISO
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
    total_biso_vol = total_count * V_BISO
    total_uo2_vol  = total_count * V_KERNEL

    # Inboard
    for idx, row in df_inboard.iterrows():
        cid = 10000 + idx + 1
        in_univ.append(f"{cid:5d}  0  -151  trcl=({row['x']:.6e} {row['y']:.6e} {row['z']:.6e})  imp:n=1  tmp=7.76e-8  vol={V_BISO:.3e}  fill=50")
        in_comp.append(f"#{cid}")
        in_tally.append(f"{cid}")

    # Outboard
    for idx, row in df_outboard.iterrows():
        cid = 20000 + idx + 1
        out_univ.append(f"{cid:5d}  0  -151  trcl=({row['x']:.6e} {row['y']:.6e} {row['z']:.6e})  imp:n=1  tmp=7.76e-8  vol={V_BISO:.3e}  fill=50")
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
        "V_br" : f"{V_br:.6e}",
        "V_uo2": f"{V_KERNEL:.6e}",
        "V_sic": f"{(V_BISO-V_KERNEL):.6e}",
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
        with open(f'prism_{case}.template', 'r', encoding='utf-8') as f:
            template_content = f.read()
        
        template = Template(template_content)
        rendered_mcnp = template.render(context)
        filepath = f'./inputs/prism_{case}_{density_str}kgm3.inp'

        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(rendered_mcnp)
        
        print(f"Successfully generated: {filepath}\n")
    else:
        print(f"Template was NOT printed\n")

    

if __name__ == "__main__":
    for fertile_kgm3 in [0.10, 0.50, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]:
        main(fertile_kgm3, go_pack=False, go_print=True)