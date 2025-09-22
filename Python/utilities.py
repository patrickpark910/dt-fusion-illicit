"""
Set of helper functions
"""
import os, re, sys, time, functools
import numpy as np

from Python.parameters import *

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
DENSITY_ThF4 = 6.3 # g/cm3
AMU_PU239 = 239.0521634
AMU_U233 = 233.039635207
SEC_PER_YR = 3600 * 24 * 365

"""
Fusion parameters - from JL Ball 24
"""
N_PER_MJ = 3.546e17 # 17.6 MeV
P_FUS_MW = 500
NPS_FUS  = P_FUS_MW * N_PER_MJ # n/s



class Colors:
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'  # Reset to default
    


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


def calc_blanket_mass_fracs(fertile_bulk_density_kgm3, breeder_volume_m3, fertile_element='U', fertile_enrich=0.71, breeder_density_kgm3=1.94e3):
    """
    Calculate the mass fractions of FLiBe and (UF4 or ThF4). 
    Assumes UF4/ThF4 dissolution does NOT change FLiBe volume
      after conversaion with J.L. Ball --ppark 2025-07-03

    Args:
        fertile_bulk_density_kgm3 : float : desired bulk density kg/m^3 of U-238 or Th-232 in blanket
        breeder_volume_m3    : float : total volume m^3 of breeding regions in tokamak
        fertile_element      : 'U' or 'Th' : identifies U-238 or Th-232 as the fertile isotope
        fertile_enrich       : float : wt% enrichment of the fertile material
        breeder_density_kgm3 : float : density of FLiBe in kg/m^3

    Returns:
        (breeder_mass_frac, uf4_mass_frac) : 2-ple of floats : mass fractions
    """

    fertile_enrich = fertile_enrich*0.01 # OpenMC uses this in % but we want fraction

    if fertile_element == 'U':
        uf4_bulk_density_kgm3 = fertile_bulk_density_kgm3 / (1-fertile_enrich) * AMU_UF4 / AMU_U
        uf4_mass_kg           = uf4_bulk_density_kgm3 * breeder_volume_m3
        breeder_mass_kg       =  breeder_density_kgm3 * breeder_volume_m3
        blanket_mass_kg       = uf4_mass_kg + breeder_mass_kg
        uf4_mass_frac         = uf4_mass_kg / blanket_mass_kg
        breeder_mass_frac     = breeder_mass_kg /blanket_mass_kg
        return breeder_mass_frac, uf4_mass_frac

    elif fertile_element == 'Th':
        thf4_bulk_density_kgm3 = fertile_bulk_density_kgm3 / (1-fertile_enrich) * AMU_ThF4 / AMU_Th
        thf4_mass_kg           = thf4_bulk_density_kgm3 * breeder_volume_m3
        breeder_mass_kg        =  breeder_density_kgm3 * breeder_volume_m3
        blanket_mass_kg        = thf4_mass_kg + breeder_mass_kg
        thf4_mass_frac         = thf4_mass_kg / blanket_mass_kg
        breeder_mass_frac      = breeder_mass_kg /blanket_mass_kg
        return breeder_mass_frac, thf4_mass_frac


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
            sys.exit(f"{Colors.RED}Error finding cross section XML!{Colors.END}")
    else:
        print(f"{Colors.GREEN}'OPENMC_CROSS_SECTIONS' found!{Colors.END}")
        return xs_path


def miller_model(R0, a, kappa, delta, n=400):
    """
    Calculate parametric coordinates from the Miller local equilibrium model.
    cf. R. L. Miller et al., doi.org/10.1063/1.872666
    cf. (Justin) Ball et al., arxiv.org/pdf/1510.08923
    
    Args:
        t  : array : parameter from 0 to 2pi
        R0 : float : major radius (cm)
        a  : float : minor radius (cm)
        kappa : float : elongation
        delta : float : triangularity
    
    Returns:
        R, Z : arrays of R and Z coordinates
    """
    t = np.linspace(0, 2*np.pi, n)
    R = R0 + a * np.cos(t + delta * np.sin(t))
    Z = kappa * a * np.sin(t)

    # Ensure the contours are closed (first point = last point) for OpenMC
    if not (np.isclose(R[0], R[-1]) and np.isclose(Z[0], Z[-1])):
        R = np.append(R, R[0])
        Z = np.append(Z, Z[0])

    return list(zip(R, Z)) 




def miller_offset(R0, a, kappa, delta, d, n=400):
    """
    Calculate a new contour that is precisely [d] meters extruded from a shape generated by miller_model().
    
    Args:
        t  : array : parameter from 0 to 2pi
        R0 : float : major radius (m)
        a  : float : minor radius (m)
        kappa : float : elongation
        delta : float : triangularity
        d  : float : offset distance (m)
    
    Returns:
        R_offset, Z_offset : arrays of offset R and Z coordinates
    """
    # Original shape
    t = np.linspace(0, 2*np.pi, n)
    R = R0 + a * np.cos(t + delta * np.sin(t))
    Z = kappa * a * np.sin(t)
    
    # Derivatives
    dR_dt = -a * np.sin(t + delta * np.sin(t)) * (1 + delta * np.cos(t))
    dZ_dt = kappa * a * np.cos(t)
    
    # Normal vector components (not normalized)
    N_R = dZ_dt   # kappa * a * cos(t)
    N_Z = -dR_dt  # a * sin(t + delta*sin(t)) * (1 + delta*cos(t))
    
    # Magnitude of normal vector
    N_mag = np.sqrt(N_R**2 + N_Z**2)
    
    # Unit normal components
    N_R_unit = N_R / N_mag
    N_Z_unit = N_Z / N_mag
    
    # Offset coordinates
    R_offset = R + d * N_R_unit
    Z_offset = Z + d * N_Z_unit
    
    # Ensure the contours are closed (first point = last point) for OpenMC
    if not (np.isclose(R_offset[0], R_offset[-1]) and np.isclose(Z_offset[0], Z_offset[-1])):
        R_offset = np.append(R_offset, R_offset[0])
        Z_offset = np.append(Z_offset, Z_offset[0])

    return list(zip(R_offset, Z_offset))


def timer(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        start_cpu = time.process_time()
        
        result = func(*args, **kwargs)
        
        wall_time = time.perf_counter() - start_time
        cpu_time = time.process_time() - start_cpu
        
        print(f"┌─ {func.__name__}")
        print(f"├─ Wall time: {wall_time:.4f}s")
        print(f"└─ CPU time:  {cpu_time:.4f}s")
        # print(f"└─ Args: {len(args)} args, {len(kwargs)} kwargs")
        
        return result
    return wrapper


def has_statepoint_file(directory_path):
    """
    Check if any file starting with 'statepoint' exists in the given directory.
    
    Args:
        directory_path (str): Path to the directory to search
    
    Returns:
        bool: True if a file starting with 'statepoint' is found, False otherwise
    """
    try:
        for filename in os.listdir(directory_path):
            if filename.startswith("statepoint"):
                return True
        return False
    except FileNotFoundError:
        print(f"Directory '{directory_path}' not found.")
        return False
    except PermissionError:
        print(f"Permission denied to access directory '{directory_path}'.")
        return False


if __name__ == "__main__":
    """ Use this to test any of the functions """
    print(calc_blanket_mass_fracs(30, 291.1, fertile_element='U', fertile_enrich=0.71, breeder_density_kgm3=1.94e3))


