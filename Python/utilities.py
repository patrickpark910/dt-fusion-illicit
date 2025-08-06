"""
Set of helper functions
"""
import os, re, sys
import numpy as np


"""
Nuclear constants -- from atom.kaeri.re.kr/nuchart/
"""
AVO = 6.022141076e+23
AMU_LI6, AMU_LI7 = 6.0150, 7.0160 # 6.01512288742, 7.01600343426 # amu = g/mol
AMU_F19 = 18.9984 # 18.99840316207
AMU_BE9 =  9.0120 # 9.012183062
AMU_U = 238.02891 # for natural enrichment
AMU_U235 = 235.0439299
AMU_U238 = 238.05078826
AMU_O = 15.999
AMU_Th = 232.0381
AMU_UF4 = AMU_U + 4 * AMU_F19
AMU_ThF4 = AMU_Th + 4 * AMU_F19
AMU_FLIBE = 98.89 # g/mol
DENSITY_UF4 = 6.7 # g/cm3
DENSITY_ThF4 = 6.3
AMU_PU239 = 239.0521634
AMU_U233 = 233.039635207
SEC_PER_YR = 3600 * 24 * 365

"""
Fusion parameters - from JL Ball 24
"""
N_PER_MJ = 3.546e17 # 17.6 MeV
P_FUS_MW = 500
NPS_FUS  = P_FUS_MW * N_PER_MJ # n/s





def set_xs_path():
    """
    Temporary solution for finding xs files between WSL and Ubuntu on Computing Cluster without editing PATH --ppark 2025-06-28
    """
    xs_path_ubuntu = '/opt/openmc_data/endfb-viii.0-hdf5/cross_sections.xml'
    xs_path_wsl   = '/mnt/c/openmc/data/endfb-viii.0-hdf5/cross_sections.xml'
    if os.path.isfile(xs_path_ubuntu):
        return xs_path_ubuntu # use this on Zotacs --ppark
    elif os.path.isfile(xs_path_wsl):
        return xs_path_wsl
    else:
        sys.exit(f"Error finding cross section XML!")


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


def calc_mix_mass_fracs(mtu, volume=342e6, density_flibe=1.94, subtract=False):
    """
    Calculate the mass fractions of FLiBe and UF4. 
    Assumes UF4 concentration small enough that its addition does NOT change FLiBe volume
      after conversaion with J.L. Ball --ppark 2025-07-03

    Args:
        mtu : float : metric tons uranium
        volume : float : cc of system
        density_flibe : float : g/cm3 of FLiBe
        subtract: bool : default False, whether to deduct UF4 volume from FLiBe volume (if UF4 is assumed to/not to dissolve in FLiBe)

    Returns:
        (mass_frac_flibe, mass_frac_uf4) : 2-ple of floats : mass fractions
    """
    # Convert inputs to SI units
    mass_U_kg = mtu * 1e3
    density_flibe_kg_m3 = density_flibe * 1e3
    density_uf4_kg_m3 = DENSITY_UF4 * 1e3
    volume_m3 = float(volume) / 1e6

    # Compute UF4 mass from U mass
    mass_uf4_kg = mass_U_kg * (AMU_UF4 / AMU_U)
    vol_uf4_m3 = mass_uf4_kg / density_uf4_kg_m3
    
    if subtract:
        vol_flibe_m3 = volume_m3 - vol_uf4_m3
    else:
        vol_flibe_m3 = volume_m3

    # Compute FLiBe mass from density and volume
    mass_flibe_kg = density_flibe_kg_m3 * vol_flibe_m3

    # Total mass and fractions
    mass_total = mass_flibe_kg + mass_uf4_kg
    frac_flibe = mass_flibe_kg / mass_total
    frac_uf4   = mass_uf4_kg   / mass_total

    # print(f"volumes: total {volume_m3} m3, flibe {vol_flibe_m3} m3, uf4 {vol_uf4_m3:.4f} m3 | deduct UF4 from FLiBe volume: {subtract}") # For debug

    return frac_flibe, frac_uf4


def calc_mix_vol_fracs(mtu, volume=342e6, density_flibe=1.94, displace=False):
    """
    Calculate the volume fractions of FLiBe and UF4. 
    Assumes UF4 concentration small enough that its addition does NOT change FLiBe volume
    after conversaion with J. L. Ball --ppark 2025-07-03

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

    return vf_flibe, vf_uf4
def calc_mix_vol_fracs_th(mtu, volume=342e6, density_flibe=1.94, displace=False):
    """
    Calculate the volume fractions of FLiBe and ThF4. 
    Assumes ThF4 concentration small enough that its addition does NOT change FLiBe volume
    after conversaion with J. L. Ball --ppark 2025-07-03

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
    density_thf4   = DENSITY_ThF4 * 1e3   # kg/m3  
    volume = float(volume) / 1e6

    # Compute volumes
    mass_thf4 = mass_th * (AMU_ThF4 / AMU_Th)
    vol_thf4 = mass_thf4 / density_thf4
    
    if displace:
        vol_flibe = volume - vol_thf4
    if not displace:
        vol_flibe = volume

    # Compute volume fractions
    vf_flibe = vol_flibe/(vol_flibe+vol_thf4)
    vf_thf4   = vol_thf4/(vol_flibe+vol_thf4)

    return vf_flibe, vf_thf4


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


def extract_lie(path: str) -> float:
    """
    Extract the float value following 'Li' in the given path string.

    Args:
        path (str): A filename or identifier containing a substring like 'Li<value>'.

    Returns:
        float: The numeric value that immediately follows 'Li'.

    Raises:
        ValueError: If no valid 'Li<value>' pattern is found.
    """
    pattern = re.compile(r'Li([+-]?\d+(?:\.\d+)?)')
    match = pattern.search(path)
    if not match:
        raise ValueError(f"No 'Li<value>' pattern found in: {path!r}")
    return float(match.group(1))


def set_xs_path():
    """
    Temporary solution for finding xs files between WSL and Ubuntu on Computing Cluster without editing PATH --ppark 2025-06-28
    """
    xs_path_ubuntu = '/opt/openmc_data/endfb-viii.0-hdf5/cross_sections.xml'
    xs_path_wsl   = '/mnt/c/openmc/data/endfb-viii.0-hdf5/cross_sections.xml'

    xs_path = os.environ.get("OPENMC_CROSS_SECTIONS")
    print(f"Checking if 'OPENMC_CROSS_SECTIONS' is set: {xs_path}")

    if xs_path is None:
        if os.path.isfile(xs_path_ubuntu):
            return xs_path_ubuntu # use this on Zotacs --ppark
        elif os.path.isfile(xs_path_wsl):
            return xs_path_wsl
        else:
            sys.exit(f"Error finding cross section XML!")
    else:
        print(f"'OPENMC_CROSS_SECTIONS' found!")
        return xs_path




if __name__ == '__main__':
    for mtu in [0, 0.1, 1, 2.5, 5, 10, 20 ,30 ,40, 50]:
        mf_flibe, mf_uf4 = calc_uf4_flibe_mass_fracs(mtu, volume=342e6, density_flibe=1.94)
        # print(mf_flibe, mf_uf4, mf_uf4/(mf_uf4+mf_flibe))
        # print(mass_to_molar_fracs({'flibe':mf_flibe,'uf4':mf_uf4}, {'flibe':AMU_FLIBE,'uf4':AMU_UF4}))


