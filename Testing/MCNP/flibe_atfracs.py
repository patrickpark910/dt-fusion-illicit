import numpy as np

def compute_isotopic_number_densities(
    volume_m3: float = 342.0,
    density_salt_gcc: float = 1.94,
    u_to_be_ratio: float = 200e-6,
    u235_wt_frac: float = 0.007204,
    li6_at_frac: float = 0.0759
) -> dict:
    """
    Compute isotopic number densities (atoms/cm³) for U-235, U-238, Li-6, Li-7, Be-9, and F-19
    in a Li₂BeF₄ salt mixture with added metallic uranium.

    Parameters
    ----------
    volume_m3
        Total mixture volume [m³].
    density_salt_gcc
        Density of pure Li₂BeF₄ salt [g/cm³].
    u_to_be_ratio
        Mass ratio of U metal to Be in the mixture (m_U / m_Be).
    u235_wt_frac
        Weight fraction of U-235 in the total uranium.
    li6_at_frac
        Atomic fraction of Li-6 in total lithium.

    Returns
    -------
    densities : dict
        Number densities [atoms/cm³] for each isotope.
    """
    # Atomic masses [g/mol]
    m_u235 = 235.0439299
    m_u238 = 238.05078826
    m_li6 = 6.0151223
    m_li7 = 7.0160030
    m_be9 = 9.0121831
    m_f19 = 18.9984032

    # Avogadro's number [atoms/mol]
    N_A = 6.02214076e23

    # Convert total volume to cm³
    volume_cm3 = volume_m3 * 1e6

    # Total mass of Li₂BeF₄ salt [g]
    mass_salt = density_salt_gcc * volume_cm3

    # Average atomic weight of enriched Li
    A_li = li6_at_frac * m_li6 + (1 - li6_at_frac) * m_li7

    # Formula weight of Li₂BeF₄ [g/mol]
    M_salt = 2 * A_li + m_be9 + 4 * m_f19

    # Mass fractions of elements in the salt
    w_li = 2 * A_li / M_salt
    w_be = m_be9        / M_salt
    w_f  = 4 * m_f19    / M_salt

    # Masses of Li, Be, F in the salt [g]
    mass_li = mass_salt * w_li
    mass_be = mass_salt * w_be
    mass_f  = mass_salt * w_f

    # Mass of U metal added [g]
    mass_u = mass_be * u_to_be_ratio

    # Split uranium mass into isotopes
    mass_u235 = mass_u * u235_wt_frac
    mass_u238 = mass_u * (1 - u235_wt_frac)

    # Compute number densities
    densities = {
        'U235': 1e-24 * (mass_u235 / m_u235) * N_A / volume_cm3,
        'U238': 1e-24 * (mass_u238 / m_u238) * N_A / volume_cm3,
        'Li6' : 1e-24 * ( (mass_li / A_li) * li6_at_frac     * N_A ) / volume_cm3,
        'Li7' : 1e-24 * ( (mass_li / A_li) * (1 - li6_at_frac) * N_A ) / volume_cm3,
        'Be9' : 1e-24 * (mass_be  / m_be9) * N_A / volume_cm3,
        'F19' : 1e-24 * (mass_f   / m_f19) * N_A / volume_cm3
    }

    
    total = 1e-24 * ( (mass_u235 / m_u235) * N_A / volume_cm3 \
              + (mass_u238 / m_u238) * N_A / volume_cm3 \
              + ( (mass_li / A_li) * li6_at_frac     * N_A ) / volume_cm3 \
              + ( (mass_li / A_li) * (1 - li6_at_frac) * N_A ) / volume_cm3 \
              + (mass_be  / m_be9) * N_A / volume_cm3 \
              + (mass_f   / m_f19) * N_A / volume_cm3 )

    densities = {
        'U235': 1e-24/total * (mass_u235 / m_u235) * N_A / volume_cm3,
        'U238': 1e-24/total * (mass_u238 / m_u238) * N_A / volume_cm3,
        'Li6' : 1e-24/total * ( (mass_li / A_li) * li6_at_frac     * N_A ) / volume_cm3,
        'Li7' : 1e-24/total * ( (mass_li / A_li) * (1 - li6_at_frac) * N_A ) / volume_cm3,
        'Be9' : 1e-24/total * (mass_be  / m_be9) * N_A / volume_cm3,
        'F19' : 1e-24/total * (mass_f   / m_f19) * N_A / volume_cm3
    }
    # """

    return densities


def compute_mass_density(
    volume_m3: float = 342.0,
    density_salt_gcc: float = 1.94,
    u_to_be_ratio: float = 200e-6,
    li6_at_frac: float = 0.0759
) -> float:
    """
    Compute the bulk mass density [g/cm³] of a Li₂BeF₄ salt mixture
    with added metallic uranium.

    Parameters
    ----------
    volume_m3
        Total mixture volume [m³].
    density_salt_gcc
        Density of pure Li₂BeF₄ salt [g/cm³].
    u_to_be_ratio
        Mass ratio of U metal to Be in the mixture (m_U / m_Be).
    li6_at_frac
        Atomic fraction of Li-6 in total lithium.

    Returns
    -------
    density_mix
        Bulk mass density of the mixture [g/cm³].
    """
    # Atomic masses [g/mol]
    m_li6 = 6.0151223
    m_li7 = 7.0160030
    m_be9 = 9.0121831
    m_f19 = 18.9984032

    # Convert total volume to cm³
    volume_cm3 = volume_m3 * 1e6

    # Total mass of Li₂BeF₄ salt [g]
    mass_salt = density_salt_gcc * volume_cm3

    # Compute formula weight of enriched Li₂BeF₄
    A_li = li6_at_frac * m_li6 + (1 - li6_at_frac) * m_li7
    M_salt = 2 * A_li + m_be9 + 4 * m_f19

    # Mass fraction of Be in the salt
    w_be = m_be9 / M_salt

    # Mass of Be in the salt [g]
    mass_be = mass_salt * w_be

    # Mass of U metal added [g]
    mass_u = mass_be * u_to_be_ratio

    # Bulk density assuming fixed total volume
    density_mix = (mass_salt + mass_u) / volume_cm3

    return density_mix


if __name__ == "__main__":
    rho = compute_mass_density()
    print(f"Mixture density: {rho:.6f} g/cm³")

    dens = compute_isotopic_number_densities()
    for iso, nd in dens.items():
        print(f"{iso}: {nd:.8e} atoms/b-cm")
