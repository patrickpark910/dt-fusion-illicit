import openmc

def make_doped_flibe(dopant, dopant_mass, Li6_enrichment=7.5, name='doped_flibe', volume=None, dopant_mass_units="kg"):
    """
    Return openmc material doped with specified fertile material

    Parameters
    ----------
    dopant : str
        "U" for U-238 -> Pu-239, "Th" for Th-232 -> U233
    dopant_mass : float
        mass of fertile material in kilograms
    Li6_enrichment : float
        The percent of the lithium which is Li-6 instead of Li-7.
    name : str
        the name of the material returned
    volume : the volume of the material returned in cubic centimeters
    
    Returns
    -------
    openmc.Material, FLiBe material doped with specified amount of fertile material
    """
    dopant_mass = dopant_mass * 1000 #Dopant mass is provided in units of kg, so here convert to grams

    flibe = openmc.Material()
    flibe.add_elements_from_formula('F4Li2Be', enrichment_target='Li6', enrichment_type='ao', enrichment=Li6_enrichment)
    flibe.set_density('g/cm3', 1.94)
    flibe.depletable = True

    """ Uranium tetrafluroide """
    UF4_molar_mass = 314.02 #g/mol
    uf4 = openmc.Material()
    uf4.add_elements_from_formula('UF4')
    uf4.set_density('g/cm3', 6.7)

    if dopant == 'U':
        tetrafluoride = uf4
    elif dopant == 'Th':
        tetrafluoride = thf4
    else:
        raise ValueError("Invalid dopant passed into blanket liquid function")

    if dopant_mass_units == "kg":
        if volume == None:
            raise ValueError("Volume of blanket specified as None")
        else:
            flibe_mass = flibe.density * volume
            tetrafluoride_mass = get_tetrafluoride_mass(dopant_mass, dopant)
            tetrafluoride_weight_percent = tetrafluoride_mass / (flibe_mass + tetrafluoride_mass)
            print(f"for {dopant_mass} kg U, UF4 mass: {tetrafluoride_mass/1e6:.2f} MT, FLiBe mass: {flibe_mass/1e6:.2f} MT, UF4 wt% in FLiBe + UF4: {tetrafluoride_weight_percent}")
            
    elif dopant_mass_units == "wppm":
        tetrafluoride_weight_percent = dopant_mass/1e6
    
    else:
        raise ValueError("Invalid units given for dopant mass argument")

    doped_mat = openmc.Material.mix_materials([tetrafluoride, flibe], [tetrafluoride_weight_percent, 1 - tetrafluoride_weight_percent], 'wo', name=name)
    doped_mat.volume = volume
    doped_mat.depletable = True
    return doped_mat


def get_tetrafluoride_mass(mass, dopant):
    """
    Computes mass of tetrafluroide from a given mass of pure dopant

    Parameters
    ----------
    mass : float
        mass of fertile material in grams
    dopant : str
        "U" for U-238 -> Pu-239, "Th" for Th-232 -> U233

    Returns
    -------
    float, mass of actinide tetrafluoride containing 'mass' grams of fertile material
    """
    UF4_molar_mass = 314.02 #g/mol

    if dopant == 'U':
        moles = mass / openmc.data.atomic_mass('U238')
        tetrafluoride_mass = moles * UF4_molar_mass

    elif dopant == 'Th':
        moles = mass / openmc.data.atomic_mass('Th232')
        tetrafluoride_mass = moles * ThF4_molar_mass
    else:
        raise ValueError("Not a valid dopant type")

    return tetrafluoride_mass


if __name__ == '__main__':
    make_doped_flibe('U',50e3,volume=342e6)
    