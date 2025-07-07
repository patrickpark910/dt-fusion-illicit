<<<<<<< Updated upstream
"""
Testing out helper functions
  Patrick 2025-06-25
"""

import numpy as np

# Nuclear constants
#   from atom.kaeri.re.kr/nuchart/
=======

import numpy as np


>>>>>>> Stashed changes
AVO = 6.022141076e+23
AMU_LI6, AMU_LI7 = 6.0150, 7.0160 # 6.01512288742, 7.01600343426 # amu = g/mol
AMU_F19 = 18.9984 # 18.99840316207
AMU_BE9 =  9.0120 # 9.012183062
<<<<<<< Updated upstream

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


=======
AMU_U = 238.02891 # for natural enrichment
AMU_UF4 = AMU_U + 4 * AMU_F19
AMU_FLIBE = 98.89 # g/mol
DENSITY_UF4 = 6.7 # g/cm3


def calc_mix_vol_fracs(mtu, volume=342e6, density_flibe=1.94, displace=False):
    """
    Calculate the volume fractions of FLiBe and UF4. 
    Assumes UF4 concentration small enough that its addition does NOT change FLiBe volume
      after conversaion with J.L. Ball --ppark 2025-07-03

    Args:
        mtu : float : metric tons uranium
        volume : float : cc of system
        density_flibe : float : g/cm3 of FLiBe
        displace : bool : default False, whether to deduct UF4 volume from FLiBe volume (if UF4 is assumed to/not to dissolve in FLiBe)

    Returns:
        (vf_flibe, vf_uf4) : 2-ple of floats : volume fractions
    """
    # Convert inputs to SI units
    mass_u = mtu * 1e3
    density_flibe = density_flibe * 1e3 # kg/m3
    density_uf4   = DENSITY_UF4 * 1e3   # kg/m3  
    volume = float(volume) / 1e6

    # Compute volumes
    mass_uf4 = mass_u * (AMU_UF4 / AMU_U)
    vol_uf4 = mass_uf4 / density_uf4
    
    if displace:
        vol_flibe = volume - vol_uf4
    if not displace:
        vol_flibe = volume

    # Compute volume fractions
    vf_flibe = vol_flibe/(vol_flibe+vol_uf4)
    vf_uf4   = vol_uf4/(vol_flibe+vol_uf4)

    # print(f"volumes: total {volume_m3} m3, flibe {vol_flibe_m3} m3, uf4 {vol_uf4_m3:.4f} m3 | deduct UF4 from FLiBe volume: {subtract}") # For debug

    return vf_flibe, vf_uf4

if __name__ == '__main__':
    print(calc_mix_vol_fracs(50,displace=True))
    print(calc_mix_vol_fracs(50,displace=False))

    rho   = 6.7                     # g/cm³
    M_UF4 = 238.0289 + 4*18.9984    # ≃314.0225 g/mol
    N_A   = 6.02214076e23           # mol⁻¹
    f_U   = 238.0289/M_UF4          # mass frac. of U in UF4
    f238  = 1 - 0.7204/100          # weight frac. U238 (natural U235 = 0.7204 wt%)
    n     = rho/M_UF4 * N_A * f_U * f238   # atoms U238 per cm³
    n_bcm = n/1e24                          # atoms per barn-cm
    print(n_bcm)  # → ≃0.0033651432
>>>>>>> Stashed changes
