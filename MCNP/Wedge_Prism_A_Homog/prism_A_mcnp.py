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
    case = 'A'
    density_str = f"{fertile_kgm3:06.2f}"
    target_total_biso = 2012  # >> DO NOT CHANGE WITHOUT CHANGING IT IN prism_{A,B,C}.py!! AND CHANGE N1, N2 (WHICH ARE HARDCODED FOR NOW) <<
    
    # HCPB tokamak dimensions (along z in model)
    h1, h2 = 45.0, 82.0     # inboard, outboard thicknesses
    d1 = 2.5 + 0.2          # Eurofer + first wall 
    d0 = 400.0              # plasma chamber
    oo = h1 + d0 + 2 * d1   # outboard offset (where outboard starts in z)
    f1, f2 = 0.39, 0.61     # fraction of fluence into inboard vs. outboard in full OpenMC HCPB model

    # Volume fractions of material in breeding regions of HCPB blanket (Lu et al. 2017)
    VF_LI_NOM = 0.1304
    VF_BE_NOM = 0.3790
    VF_EU_NOM = 0.1176
    VF_HE_NOM = 0.3730

    # Number of particles for MCNP run
    nps = '5e6' 

    # Volume of biso sphere and UO2 kernel
    V_BISO   = (4/3) * np.pi * (0.05)**3
    V_KERNEL = (4/3) * np.pi * (0.04)**3
    
    # Dictionary of homogenized breeding region densities from OpenMC
    rho_dict = {  0.1: 1.8840986616930262,
                  0.5: 1.8843170404651346,
                  1.5: 1.8848629873954048,
                 15.0: 1.8922332709540537,
                 30.0: 1.9004224749081071, 
                 60.0: 1.9168008828162153, 
                 90.0: 1.9331792907243230,
                120.0: 1.9495576986324308,
                150.0: 1.9659361065405383,
                250.0: 2.0205307995675650,
                500.0: 2.1570175321351300,
                750.0: 2.2935042647026950,
                999.999: 2.429990997270261,}
    rho_breeder = rho_dict[fertile_kgm3]

    """
    Geometric parameters for tilted plasma chamber
    """
    z_in_end    = 47.7
    z_out_start = 447.7
    dz = z_out_start - z_in_end
    dx = b2 - b1
    slope = dx / dz
    
    # Plane: 1.0*x + 0.0*y + plane_C*z - plane_D = 0
    plane_C = -slope
    plane_D = b1 - (slope * z_in_end)
    mat_breeder = open(f"./mat_cards/{density_str}.txt").read()

    """
    Jinja2 template insert
    """
    context = {
        "case": case,
        "fertile_kgm3": fertile_kgm3,
        "density_str": density_str,
        "summary": (f"c face lengths scaled for {target_total_biso} BISO in case B\n"
                    f"c inboard (b1): {b1:.6f} cm | outboard: {b2:.6f} cm"),
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
        "rho_breeder": f"{rho_breeder:.6f}",
        "mat_breeder": mat_breeder,
        "plane_C": f"{plane_C:.6e}",
        "plane_D": f"{plane_D:.6e}",
        "nps": nps
    }


    """ PRINT TO TEMPLATE
    """
    if go_print:
        with open(f'prism_{case}.template', 'r') as f:
            template_content = f.read()
        
        template = Template(template_content)
        rendered_mcnp = template.render(context)
        filepath = f'./inputs/prism_{case}_{density_str}kgm3.inp'

        with open(filepath, 'w') as f:
            f.write(rendered_mcnp)
        
        print(f"Successfully generated: {filepath}\n")
    else:
        print(f"Template was NOT printed\n")

    

if __name__ == "__main__":
    for fertile_kgm3 in [0.10, 0.50, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.999]:
        main(fertile_kgm3, go_print=True)