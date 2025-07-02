"""
Set of helper functions
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

"""
Nuclear constants -- from atom.kaeri.re.kr/nuchart/
"""
AVO = 6.022141076e+23
AMU_LI6, AMU_LI7 = 6.0150, 7.0160 # 6.01512288742, 7.01600343426 # amu = g/mol
AMU_F19 = 18.9984 # 18.99840316207
AMU_BE9 =  9.0120 # 9.012183062
AMU_U = 238.02891 # for natural enrichment
AMU_UF4 = AMU_U + 4 * AMU_F19
AMU_FLIBE = 98.89 # g/mol
DENSITY_UF4 = 6.7 # g/cm3



"""
Font settings for matplotlib
"""
try:
    font_path = './helper/arial.ttf'  # Your font path goes here
    fm.fontManager.addfont(font_path)
    prop = fm.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = prop.get_name()
except:
    font_path = './arial.ttf'  # Your font path goes here
    fm.fontManager.addfont(font_path)
    prop = fm.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = prop.get_name()


def logspace_per_decade(start, stop, pts_per_decade):
    """
    Returns values from 'start' to 'stop' so that each factor-of-10
    interval contains 'pts_per_decade' points (including its first endpoint).
    Might be a little off if 'stop' isn't precisely at a decade, ex. 20e6 eV

    example: 10 points per decade from 1e-5 → 2e7
    grid = logspace_per_decade(1e-5, 20e6, pts_per_decade=10)
    for i in grid:
        print(np.log10(i))
    # print(np.log10(i) for i in grid)
    """
    log_start = np.log10(start)
    log_stop  = np.log10(stop)
    total_decades = log_stop - log_start
    # how many intervals of size 1 decade we need, as a float
    total_steps = total_decades * pts_per_decade
    # +1 so that we include the very last point at `stop`
    npts = int(np.ceil(total_steps)) + 1
    return np.logspace(log_start, log_stop, num=npts)


def log_midpoints(points):
    """
    Compute logarithmic midpoints between adjacent points in a list.
    The logarithmic midpoint m between two positive values a and b is
    defined as the geometric mean: m = sqrt(a * b).

    Args:
        points: List of floats (must be positive).

    Returns:
        List of floats of length len(points) - 1, where each entry
        is the logarithmic midpoint between points[i] and points[i+1].
    """
    if any(p <= 0 for p in points):
        raise ValueError("All points must be positive to compute log midpoints.")
    return [float(np.sqrt(points[i] * points[i + 1])) for i in range(len(points) - 1)]


def calc_uf4_flibe_mass_fracs(mtu, volume=342e6, density_flibe=1.94):
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
    density_uf4_kg_m3 = DENSITY_UF4 * 1e3
    volume_m3 = float(volume) / 1e6

    # compute UF4 mass from U mass
    mass_uf4_kg = mass_U_kg * (AMU_UF4 / AMU_U)
    vol_uf4_m3 = mass_uf4_kg / density_uf4_kg_m3
    vol_flibe_m3 = volume_m3 - vol_uf4_m3

    # print(f"volumes: total {volume_m3} m3, flibe {vol_flibe_m3} m3, uf4 {vol_uf4_m3}")

    # compute FLiBe mass from density and volume
    mass_flibe_kg = density_flibe_kg_m3 * vol_flibe_m3

    # total mass and fractions
    mass_total = mass_flibe_kg + mass_uf4_kg
    frac_flibe = mass_flibe_kg / mass_total
    frac_uf4   = mass_uf4_kg   / mass_total

    print(mass_uf4_kg)

    return frac_flibe, frac_uf4


def mass_to_molar_fracs(mass_fracs, molar_masses):
    """
    Convert mass fractions to molar fractions.

    Args:
        mass_fracs (dict): component mass fractions, ex. {'flibe': 0.75, 'uf4': 0.25}.
        molar_masses (dict): component molar masses in g/mol, ex. {'flibe': 98.89, 'uf4': 318.02}.

    Returns:
        dict: component molar fractions.
    """
    # Compute "moles per unit mass" for each component
    moles_per_g = {
        comp: mass_fracs[comp] / molar_masses[comp]
        for comp in mass_fracs
    }
    total = sum(moles_per_g.values())

    # Normalize to get molar fractions
    return {
        comp: moles_per_g[comp] / total
        for comp in moles_per_g
    }


if __name__ == '__main__':
    for mtu in [50]: # 0, 0.1, 1, 2.5, 5, 10, 20 ,30 ,40, 50]:
        mf_flibe, mf_uf4 = calc_uf4_flibe_mass_fracs(mtu, volume=342e6, density_flibe=1.94)
        print(mf_flibe, mf_uf4, mf_uf4/(mf_uf4+mf_flibe))
        print(mass_to_molar_fracs({'flibe':mf_flibe,'uf4':mf_uf4}, {'flibe':AMU_FLIBE,'uf4':AMU_UF4}))


