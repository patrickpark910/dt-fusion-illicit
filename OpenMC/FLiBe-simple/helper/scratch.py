"""
Testing out helper functions
  Patrick 2025-06-25
"""

import numpy as np

# Nuclear constants
#   from atom.kaeri.re.kr/nuchart/
AVO = 6.022141076e+23
AMU_LI6, AMU_LI7 = 6.0150, 7.0160 # 6.01512288742, 7.01600343426 # amu = g/mol
AMU_F19 = 18.9984 # 18.99840316207
AMU_BE9 =  9.0120 # 9.012183062

def li_atfracs(w6):
    """
    Convert a lithium-6 enrichment by weight fraction into atomic fractions of Li-6 and Li-7.

    Args: 
        e_li6 (float): Weight fraction of Li-6 in Li (ex. 0.20 for 20 wt% Li-6)
        
    Returns: 
        x6, x7 (floats): Atomic fractions of Li-6, Li-7

    Examples:
        >>> li_atfracs(0.2)
        (0.22576563661046428, 0.7742343633895358)
    """
    if not 0.0 <= w6 <= 1.0: 
        raise ValueError("li6_enrich must be between 0 and 1")

    w7 = 1.0 - w6 # wt frac Li-7
    M6, M7 = w6 / AMU_LI6, w7 / AMU_LI7 # moles
    a6, a7 = M6 / (M6 + M7), M7 / (M6 + M7) # at fracs

    return a6, a7


def flibe_ndens(mol_lif, mol_bef2, dens_flibe, enrich_li):
    """
    Given mols of LiF to BeF2, FLiBe bulk mass density, and Li enrichment (by wt Li-6),
    compute number densities (at/cc) of Li-6, Li-7, Be-9, F-19.

    Args: 
        mol_lif, mol_bef2 (floats): moles of LiF, BeF2 in mixture
        dens_flibe (float): bulk density of FLiBe mixture [g/cc]
        enrich_li (float): wt frac of Li-6 in Li 
        
    Returns: 
        ndens_flibe (dict): number densities {'Li6','Li7','Be9','F19'} in mixture
    """
    
    # Normalize moles [mol frac]
    mol_tot = mol_lif + mol_bef2
    x_lif, x_bef2 = mol_lif / mol_tot, mol_bef2 / mol_tot

    # 
    a6, a7 = li_atfracs(enrich_li)

    # Mols [g/mol]
    M_li = a6 * AMU_LI6 + a7 * AMU_LI7
    M_lif = M_li + AMU_F19 
    M_bef2 = AMU_BE9 + 2*AMU_F19 
    M_flibe = x_lif * M_lif + x_bef2 * M_bef2

    # Number densities [1/cm3]
    N_flibe = dens_flibe / M_flibe
    N_lif = N_flibe * x_lif * AVO
    N_bef2 = N_flibe * x_bef2 * AVO

    ndens_flibe = {'Li6' : 1e-24 * N_lif * a6,       # 1 Li per LiF, a6 Li-6 per Li
                   'Li7' : 1e-24 * N_lif * a7,       # 1 Li per LiF, a7 Li-6 per Li
                   'Be9' : 1e-24 * N_bef2,           # 1 Be per BeF2
                   'F19' : 1e-24 * (N_lif + 2*N_bef2)  # F per LiF + 2 F per BeF2
                    }

    return ndens_flibe



if __name__ == '__main__':
    print(li_atfracs(0.2))
    print(flibe_ndens(mol_lif=0.66, mol_bef2=0.34,dens_flibe=1.80,enrich_li=0.20))


