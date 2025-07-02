import openmc
import numpy as np

# User specifications
DENSITY_FLIBE = 1.8 # usu 1.8-2.0 g/cc
ENRICH_LI = 20
MOL_LIF, MOL_BEF2 = 0.66, 0.34
TEMP = 900 # K
VOL = 342 * 1e6 # cm3

AMU_F19 = 18.9984 # 18.99840316207
AMU_U = 238.02891 # for natural enrichment
AMU_UF4 = AMU_U + 4* AMU_F19

def calc_uf4_flibe_mass_fracs(mtu, volume=342e6, density_flibe=1.80):
    """
    Calculate the mass fractions of FLiBe and UF4.

    Args:
        mtu : float : metric tons uranium
        volume : float : cubic meters of system (342 m3 in Ball 25)
        density_flibe : float : g/cm3 of FLiBe (1.8 g/cm3)

    Returns:
        (mass_frac_flibe, mass_frac_uf4) : 2-ple of floats : mass fractions
    """
    # convert inputs to SI units
    mass_U_kg = mtu * 1e3
    density_flibe_kg_m3 = density_flibe * 1e3
    volume_m3 = float(volume) / 1e6

    # compute UF4 mass from U mass
    mass_uf4_kg = mass_U_kg * (AMU_UF4 / AMU_U)

    # compute FLiBe mass from density and volume
    mass_flibe_kg = density_flibe_kg_m3 * volume_m3

    # total mass and fractions
    mass_total = mass_flibe_kg + mass_uf4_kg
    frac_flibe = mass_flibe_kg / mass_total
    frac_uf4   = mass_uf4_kg   / mass_total

    return frac_flibe, frac_uf4




if __name__ == "__main__":

    """
    Calculate 

    Args:
        mass_U_list : list of metric tons uranium (MTU) we want to test
        
    """
    MASS_U_LIST = [0, 0.1, 1, 2.5, 5, 10, 20 ,30 ,40, 50]

    """
    MATERIALS
      OpenMC automatically normalizes the fraction of each element in material (like MCNP)
      but the fractions in 'mix_materials' MUST add up to 1
    """
    flibe = openmc.Material(name="FLiBe", temperature=TEMP)
    flibe.set_density('g/cm3', DENSITY_FLIBE) 
    flibe.add_element('Li',MOL_LIF,'ao',ENRICH_LI,'Li6','wo')
    flibe.add_element('Be',MOL_BEF2,'ao')
    flibe.add_element('F',(MOL_LIF+2*MOL_BEF2),'ao')

    uf4 = openmc.Material(name="U",temperature=TEMP)
    uf4.add_elements_from_formula('UF4','wo',0.7204)\

    # Calculate mol ratios of UF4 and FLiBe, ensure they add up to 1
    mass_frac_uf4_list, mass_frac_flibe_list = [], []

    mix_list = []
    for mtu in MASS_U_LIST:
        mass_frac_uf4, mass_frac_flibe = calc_uf4_flibe_mass_fracs(mtu, volume=VOL, density_flibe=DENSITY_FLIBE)
        print(mass_frac_uf4, mass_frac_flibe)
        mix = openmc.Material.mix_materials([flibe, uf4], [mass_frac_flibe, mass_frac_uf4], 'wo')
        mix.name = f"mat-{mtu:.1f}mtu"
        mix.temperature = TEMP
        mix_list.append(mix)

    materials = openmc.Materials(mix_list)



    """
    mix, ind = [], []

    for i,val in enumerate(M_U_val):
        ind.append(i)
        mat_name = f"mat-{val:.1f}tU"
        mix0 = openmc.Material.mix_materials([flibe, uranium], [MR_FLIBE[i], MR_U[i]], 'ao') # 
        mix0.temperature = 900
        mix0.name = mat_name
        mix.append(mix0)
        # print(mix[i].name, mix[i].get_nuclide_atom_densities()) # Returns one or all elements in the material and their atomic densities in units of atom/b-cm

    materials = openmc.Materials(mix)
    """