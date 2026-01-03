import numpy as np
import pandas as pd
from jinja2 import Template
import os

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


def main(fertile_kgm3, go_print=False):
    """
    Generation an MCNP input file for a HCPB unit prism.
    
    The workflow is:
    1. Defines the HCPB blanket geometry and target particle counts.
    2. Solve for the reference area (A_ref) required to fit 
       exactly 10,000 BISO particles across the inboard and outboard regions.
    3. Calculates the slope and constants for the tilted plasma chamber 
       surfaces (Surfaces 141-144) to ensure geometric consistency.
    4. Reads particle coordinates from CSV files and formatting them into 
       MCNP universe (fill) and complement (exclusion) strings.
    5. Maps all calculated parameters into a Jinja2 context and rendering 
       the 'prism_B.template' into a final '.inp' file.
    """

    """ BASIC PARAMETERS
    """
    case = 'C'
    density_str = f"{fertile_kgm3:07.3f}"
    target_total_biso = 2012  # >> DO NOT CHANGE WITHOUT CHANGING IT IN prism_{A,B,C}.py!! AND CHANGE N1, N2 (WHICH ARE HARDCODED FOR NOW) <<
    
    # HCPB tokamak dimensions (along z in model)
    h1, h2 = 45.0, 82.0     # inboard, outboard thicknesses
    d1 = 2.5 + 0.2          # Eurofer + first wall 
    d0 = 400.0              # plasma chamber
    oo = h1 + d0 + 2 * d1   # outboard offset (where outboard starts in z)
    f1, f2 = 0.39, 0.61     # fraction of fluence into inboard vs. outboard in full OpenMC HCPB model

    # Volume fractions of material in breeding regions of HCPB blanket (Lu et al. 2017)
    vf_li4sio4 = 0.1304
    vf_be      = 0.3790
    vf_eurofer = 0.1176
    vf_he      = 0.3730

    # Number of particles for MCNP run
    nps = '5e6' 

    # Volume of biso sphere and UO2 kernel
    V_biso = (4/3) * np.pi * (0.05)**3
    V_uo2  = (4/3) * np.pi * (0.04)**3
    

    """ DERIVED QUANTITIES
    """

    # NB. Fertile kg/m³ is "kg of U-238 isotope" per "m³ of BREEDING MATERIAL" where 
    #     "breeding material" = Li4SiO4 + Be and EXCLUDING Eurofer, He-4(g) also in breeding REGION.
    #     So, "biso/cm³" is "number of biso spheres" per "cm³ of Li4SiO4 + Be"
    biso_per_cc_breeder = fertile_to_biso_density(fertile_kgm3)  # biso spheres/cm³ of Li4SiO4 + Be
    biso_per_cc_region  = biso_per_cc_breeder * (vf_li4sio4 + vf_be)

    biso_vf_breeder = biso_per_cc_breeder * V_biso
    biso_vf_region  = biso_per_cc_region * V_biso

    # scaling_factor = (f1 * h1) + (f2 * h2)
    A_ref = target_total_biso / biso_per_cc_region / (h1+h2) # ( * scaling_factor)

    '''
    We derive the nominal face area 'A_ref' above:
        A_ref [cm²] = target_total_biso [spheres] * biso_per_cc_region⁻¹ [cm³/spheres] / total breeder length [cm]
    We now need to scale A_ref so that:
        1. ratio f = A1/A2 = f1/f2 = nominally 0.39/0.61
        2. volume A1*h1 + A2*h2 = A_ref*(h1+h2)
    We solve:
        A2 = A_ref (h1+h2) / ( f1/f2 * h1 + h2 ), A1 = A2 * f1/f2, 
    NB. Technically, the BISO should deduct from the Li4SiO4 and Be vols but NOT Eurofer or He-4(g) vols.
        BUT I do not have the time to script that kind of 
        SO I have the auto-generated mats from OpenMC (in ./mat_cards/) that I parse into MCNP.
    '''
    A2 = A_ref * (h1+h2) / ( f1/f2 * h1 + h2)                     # face areas
    A1 = A2 * f1/f2
    b1, b2 = np.round(np.sqrt(A1), 6), np.round(np.sqrt(A2), 6)   # face lengths
    c1, c2 = b1/2, b2/2               # (for MCNP surface cards)  # face HALF-lengths


    # Volumes of the breeding regions (i.e. area*height)
    V1, V2 = A1*h1, A2*h2

    # Real number of BISO particles in inboard vs. outboard
    N1 = round(V1 * biso_per_cc_region) # =  523 #  512  # 
    N2 = round(V2 * biso_per_cc_region) # = 1489 # 1500  # 
    N_tot = N1 + N2
    # NB. I temporarily hardcoded 512, 1500 to match & make lattice packing easier in Prism C (3-D prime factorization)

    # check
    summary = (
    f"c {78*'='}"
    f"\nc Case: {fertile_kgm3} kg/cm³"
    f"\nc "
    f"\nc Per volume of breeder (Li4SiO4 + Be):"
    f"\nc   {round(biso_per_cc_breeder, 6)} spheres/cm³ : {(biso_vf_breeder*100):.6f} vol%"
    f"\nc Per volume of blanket (Li4SiO4 + Be + Eurofer + He):"
    f"\nc   {round(biso_per_cc_region, 6)} spheres/cm³ : {(biso_vf_region*100):.6f} vol%"
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
    f"\nc {78*'='}"
    )

    print(summary)

    # Volumes of the breeding region MINUS biso volumes (for MCNP cell vol card)
    V1_br = A1*h1 - target_total_biso * V_biso
    V2_br = A2*h2 - target_total_biso * V_biso
    V_br  = V1_br + V2_br


    """
    Calculate dimensions of each lattice cell
    """
    lattice_in  = (2,2,128)
    lattice_out = (2,2,375)

    # Lattice cell length
    uli_x = b1 / lattice_in[0] 
    uli_y = b1 / lattice_in[1] 
    uli_z = h1 / lattice_in[2] 
    ulo_x = b2 / lattice_out[0] 
    ulo_y = b2 / lattice_out[1] 
    ulo_z = h2 / lattice_out[2] 

    # Lattice cell half-lengths
    ui_x = uli_x / 2
    ui_y = uli_y / 2
    ui_z = uli_z / 2
    uo_x = ulo_x / 2
    uo_y = ulo_y / 2
    uo_z = ulo_z / 2

    # Lattice offset
    si_x = 3*ui_x
    si_y = 3*ui_y
    si_z = -ui_z + h1
    so_x = 4*uo_x # -b2 + 
    so_y = 4*uo_y
    so_z = -uo_z + h1 + oo + h2


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
        "lattice_in": lattice_in,
        "lattice_out": lattice_out,
        "V1_fw": f"{(b1**2 * 0.2):.6e}",
        "V1_ei": f"{(b1**2 * 2.5):.6e}",
        "V1_eo": f"{(b1**2 * 2.5):.6e}",
        "V1_br": f"{V1_br:.6e}",
        "V2_fw": f"{(b2**2 * 0.2):.6e}",
        "V2_ei": f"{(b2**2 * 2.5):.6e}",
        "V2_eo": f"{(b2**2 * 2.5):.6e}",
        "V2_br": f"{V2_br:.6e}",
        "V_br" : f"{V_br:.6e}",
        "ui_x": f"{ui_x:.6e}",
        "ui_y": f"{ui_y:.6e}",
        "ui_z": f"{ui_z:.6e}",
        "uo_x": f"{uo_x:.6e}",
        "uo_y": f"{uo_y:.6e}",
        "uo_z": f"{uo_z:.6e}",
        "uli_x": f"{(2*ui_x):.8f}",
        "uli_y": f"{(2*ui_y):.8f}",
        "uli_z": f"{(2*ui_z):.8f}",
        "ulo_x": f"{(2*uo_x):.8f}",
        "ulo_y": f"{(2*uo_y):.8f}",
        "ulo_z": f"{(2*uo_z):.8f}",  
        "si_x" : f"{si_x:.6e}",
        "si_y" : f"{si_y:.6e}",
        "si_z" : f"{si_z:.6e}",
        "so_x" : f"{so_x:.6e}",
        "so_y" : f"{so_y:.6e}",
        "so_z" : f"{so_z:.6e}",
        "plane_C": f"{plane_C:.6e}",
        "plane_D": f"{plane_D:.6e}",
        "nps": nps
    }


    """
    Print template
    """
    with open(f'prism_{case}.template', 'r', encoding='utf-8') as f:
        template_content = f.read()
    
    template = Template(template_content)
    rendered_mcnp = template.render(context)
    filepath = f'./inputs/prism_{case}_{density_str}kgm3.inp'

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(rendered_mcnp)
    
    print(f"Successfully generated: {filepath}")

if __name__ == "__main__":
    for fertile_kgm3 in [0.10, 0.50, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.999]: #
        main(fertile_kgm3, go_print=True)