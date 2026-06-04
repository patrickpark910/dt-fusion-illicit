import numpy as np
from scipy.optimize import fsolve

def flibe_density(temp_K):
    """
    Calculates density (rho) for a given T and P using the 
    soft-sphere Helmholtz EOS derived from Humrickhouse 17 (INL-44148)
    """
    # Constants from the document
    P_target  = 101325  # FLiBe blanket is pressurized to ~1 atm (101325 Pa) cf. ORNL/TM-2024/3508 --ppark 2026-02-13
    Rs        = 84.078  # Specific gas constant [J/kg-K]
    Tm        = 732.3   # Melting point temperature [K]
    rhom      = 1991.6  # Melting point density [kg/m^3]
    
    # Parameters from Table 1 (Required for calculation)
    ni = np.array([0.9307, 36.62, 45.17, 32.29, -104.1])
    ti = np.array([0.9599, 0.9663, 1.014, 0.4468, 0.9629])
    di = np.array([0, 0, 3.151, 1.924, 1.998])

    def pressure_eos(rho):
        # Derived from P = - (da/dv)_T = rho^2 * (da/drho)_T
        # Based on Eqs.(1-2) in INL-44148
        # P = rho * Rs * T * (1 + sum(ni * di * (Tm/T)^ti * (rho/rhom)^di))
       
        term_sum = np.sum(ni * di * (Tm / temp_K)**ti * (rho / rhom)**di)
        P = (rho * Rs * temp_K * (1 + term_sum)) - P_target
    
        return P

    # Initial guess: start near melting point density
    rho_guess = 1900.0
    rho_solution, info, ier, msg = fsolve(pressure_eos, rho_guess, full_output=True)

    # print(f"rho_solution = {rho_solution}")
    
    if ier == 1:
        return rho_solution[0]
    else:
        return f"Convergence failed: {msg}"


def uf4_density(temp_K):
    """
    Calculate the density of UF4 as a function of temperature.

    This function uses the correlation for liquid UF4 density provided by the 
    ORNL MSTDB (2024), originally sourced from Kirshenbaum et al. (1961). 
    While the experimental range is 1393--1651 K (UF4 mp = 1309 K), we will use
    this function to extrapolate to the supercooled liquid state (e.g., 900 K).

    The extrapolation is physically justified for molten salt mixtures because 
    dissolved UF4 dissociates into ions (U^4+ and F-) that integrate into 
    the liquid solution structure. As per Cantor (ORNL-4308), molar volume 
    additivity applies, meaning the relevant density is that of the liquid-phase 
    molar volume at the mixture temperature, rather than the room-temperature 
    powder density (approx. 6.7 g/cm³).

    Args:
        temperature (float): The temperature of the salt in Kelvin [K].

    Returns:
        float: The density of UF4 in g/cm³.

    Example:
        >>> uf4_density(900)
        6.4269
    """
    rho = 7.11 - 0.000759 * temp_K
    return rho


def thf4_density(temp_K):
    """
    Density of ThF4 as a function of temperature [K] between (1309, 1614)
    From ORNL MSTDB (2025) which cites A. D. Kirshenbaum et al. (1961) DOI: 10.1016/0022-1902(61)80047-X

    See docstring in uf4_density for more information.

    thf4_density(900) = 6.8872 g/cm^3

    THE PROBLEM IS THAT THIS GIVES YOU A DIFFERENT DENSITY THAN CANTOR (1973, ORNL-4308) DATA FOR
    THF4 MOLAR VOLUME (cm3/mol), which comes out to 6.61 g/cm^3 at 900 K. So, in our 
    OpenMC/parameters.py, we defer to Cantor's experimental data and use DENSITY_ThF4 = 6.61 g/cm^3.
    """
    rho = 7.78 - 0.000992 * temp_K
    return rho


if __name__ == "__main__":
    rho_flibe = flibe_density(900)
    print(f"Density at 900K: {rho_flibe/1000} g/cm^3")

    rho_uf4 = uf4_density(900)
    print(f"Density at 900K: {rho_uf4} g/cm^3")

    rho_thf4 = thf4_density(900)
    print(f"Density at 900K: {rho_thf4} g/cm^3")