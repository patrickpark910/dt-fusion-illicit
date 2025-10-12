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
AMU_UO2 = AMU_U + 2 * AMU_O
AMU_ThO2 = AMU_Th + 2 * AMU_O
AMU_FLIBE = 98.89 # g/mol
DENSITY_UF4 = 6.7 # g/cm3
DENSITY_ThF4 = 6.3 # g/cm3
DENSITY_UO2 = 10.5 # g/cm3
DENSITY_ThO2 = 10.0 # g/cm3
DENSITY_SIC = 3.2 # g/cm3
AMU_PU239 = 239.0521634
AMU_U233 = 233.039635207
SEC_PER_YR = 3600 * 24 * 365

"""
Fusion parameters - from JL Ball 24
"""
N_PER_MJ = 3.546e17 # 17.6 MeV
P_FUS_MW = 1000
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
    UL = '\033[4m'
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
        thf4_bulk_density_kgm3 = fertile_bulk_density_kgm3 * AMU_ThF4 / AMU_Th
        thf4_mass_kg           = thf4_bulk_density_kgm3 * breeder_volume_m3
        breeder_mass_kg        =  breeder_density_kgm3 * breeder_volume_m3
        blanket_mass_kg        = thf4_mass_kg + breeder_mass_kg
        thf4_mass_frac         = thf4_mass_kg / blanket_mass_kg
        breeder_mass_frac      = breeder_mass_kg /blanket_mass_kg
        return breeder_mass_frac, thf4_mass_frac

def calc_biso_blanket_mass_fracs(fertile_bulk_density_kgm3, breeder_volume_m3, fertile_element='U', fertile_enrich=0.71, breeder_density_kgm3=9.4e3):
    """
    Calculate mass fractions of PbLi and BISO particles.
    Following Glaser & Goldston (2012), we assume the BISO/TRISO particles are homogenized 
    throughout the blanket, rather than resolving discrete particles. 
    The kernel/coating geometry is used only to set the relative mass or 
    volume fractions of fertile vs. coating material.
    
    Args:
        fertile_bulk_density_kgm3 : desired bulk density of fertile isotope in blanket [kg/m3]
        breeder_volume_m3         : breeding region volume [m3]
        fertile_element           : 'U' or 'Th'
        fertile_enrich            : enrichment fraction of fertile isotope [wt%]
        breeder_density_kgm3      : density of PbLi [kg/m3]
    """
    fertile_enrich = fertile_enrich*0.01 # OpenMC uses this in % but we want fraction

    if fertile_element == 'U':
        uo2_bulk_density_kgm3 = fertile_bulk_density_kgm3 / (1-fertile_enrich) * AMU_UO2 / AMU_U

        r_uo2 = 400e-4  # r = 400 μm = 0.0400 cm // "800 μm kernel"
        r_sic = 500e-4  # 500 μm = 0.0500 cm // "100 μm thickness"
        V_biso_particle = (4 / 3) * np.pi * (r_sic)**3     # volume of single BISO particle
        V_uo2_in_biso   = (4 / 3) * np.pi * (r_uo2)**3     # volume of UO2 in single BISO particle
        Vf_uo2_in_biso  = V_uo2_in_biso / V_biso_particle  # vol frac UO2 in single BISO
        Vf_sic_in_biso  = 1.0 - Vf_uo2_in_biso  
        # SiC mass per volume of UO2 (ratio inside a particle)
        density_sic     = DENSITY_SIC * 10**3 #kg/m^3
        density_uo2     = DENSITY_UO2 * 10**3#kg/m^3
        m_sic_per_m_uo2 = (density_sic * Vf_sic_in_biso) / (density_uo2 * Vf_uo2_in_biso)

        sic_bulk_density_kgm3 = uo2_bulk_density_kgm3 * m_sic_per_m_uo2

        uo2_mass_kg       = uo2_bulk_density_kgm3 * breeder_volume_m3
        sic_mass_kg       = sic_bulk_density_kgm3 * breeder_volume_m3
        biso_mass_kg      = uo2_mass_kg + sic_mass_kg   # total fertile particle mass
        breeder_mass_kg   = breeder_density_kgm3 * breeder_volume_m3  # PbLi mass
        blanket_mass_kg   = breeder_mass_kg + biso_mass_kg
        breeder_mass_frac = breeder_mass_kg / blanket_mass_kg
        biso_mass_frac    = biso_mass_kg / blanket_mass_kg

        return breeder_mass_frac, biso_mass_frac

    elif fertile_element == 'Th':
        tho2_bulk_density_kgm3 = fertile_bulk_density_kgm3 / (1-fertile_enrich) * AMU_ThO2 / AMU_Th

        r_tho2 = 400e-4  # r = 400 μm = 0.0400 cm // "800 μm kernel"
        r_sic = 500e-4  # 500 μm = 0.0500 cm // "100 μm thickness"
        V_biso_particle = (4 / 3) * np.pi * (r_sic)**3     # volume of single BISO particle
        V_tho2_in_biso   = (4 / 3) * np.pi * (r_tho2)**3     # volume of UO2 in single BISO particle
        Vf_tho2_in_biso  = V_tho2_in_biso / V_biso_particle  # vol frac UO2 in single BISO
        Vf_sic_in_biso  = 1.0 - Vf_tho2_in_biso  
        # SiC mass per volume of UO2 (ratio inside a particle)
        density_sic     = DENSITY_SIC * 10**3 #kg/m^3
        density_tho2     = DENSITY_ThO2 * 10**3#kg/m^3
        m_sic_per_m_tho2 = (density_sic * Vf_sic_in_biso) / (density_tho2 * Vf_tho2_in_biso)

        sic_bulk_density_kgm3 = tho2_bulk_density_kgm3 * m_sic_per_m_tho2

        tho2_mass_kg       = tho2_bulk_density_kgm3 * breeder_volume_m3
        sic_mass_kg       = sic_bulk_density_kgm3 * breeder_volume_m3
        biso_mass_kg      = tho2_mass_kg + sic_mass_kg   # total fertile particle mass
        breeder_mass_kg   = breeder_density_kgm3 * breeder_volume_m3  # PbLi mass
        blanket_mass_kg   = breeder_mass_kg + biso_mass_kg
        breeder_mass_frac = breeder_mass_kg / blanket_mass_kg
        biso_mass_frac    = biso_mass_kg / blanket_mass_kg

        return breeder_mass_frac, biso_mass_frac


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
    xs_path_zotac = '/opt/openmc_data/endfb-viii.0-hdf5/cross_sections.xml'
    xs_path_wsl   = '/mnt/c/OpenMC/data/endfb-viii.0-hdf5/cross_sections.xml'

    xs_path = os.environ.get("OPENMC_CROSS_SECTIONS")
    

    if xs_path is None:
        print(f"{Colors.YELLOW}Warning.{Colors.END} 'OPENMC_CROSS_SECTIONS' is NOT set in PATH.")
        
        if os.path.isfile(xs_path_zotac):
            print(f"{Colors.YELLOW}Warning.{Colors.END} 'OPENMC_CROSS_SECTIONS' is now set to {xs_path_ubuntu}.")
            return xs_path_zotac # use this on Zotacs --ppark

        elif os.path.isfile(xs_path_wsl):
            print(f"{Colors.YELLOW}Warning.{Colors.END} 'OPENMC_CROSS_SECTIONS' is now set to {xs_path_wsl}.")
            return xs_path_wsl
        
        else:
            sys.exit(f"{Colors.RED}Fatal. Cannot find cross section XML!{Colors.END}")
    
    else:
        print(f"{Colors.GREEN}'OPENMC_CROSS_SECTIONS' found!{Colors.END}")
        return xs_path


def miller_model(R0, a, kappa, delta, extrude=0, calc_vol=False, n=100):
    """
    Calculate parametric coordinates from the Miller local equilibrium model.
    cf. R. L. Miller et al., doi.org/10.1063/1.872666
    cf. (Justin) Ball et al., arxiv.org/pdf/1510.08923
    
    Args:
        R0 : float : major radius (cm)
        a  : float : minor radius (cm)
        kappa : float : elongation
        delta : float : triangularity
        extrude : float : thickness by which to extrude the boundary (for blanket layers)
        calc_vol : bool : return volume?
        n  : int : number of (r, z) pairs to generate
    
    Returns:
        R, Z : arrays of R and Z coordinates 
    """
    # Original Miller contour
    t = np.linspace(0, 2*np.pi, n)
    R = R0 + a * np.cos(t + delta * np.sin(t))
    Z = kappa * a * np.sin(t)

    # Derivatives
    dR_dt = -a * np.sin(t + delta * np.sin(t)) * (1 + delta * np.cos(t))
    dZ_dt = kappa * a * np.cos(t)

    # Outward normal
    N_R = dZ_dt
    N_Z = -dR_dt
    N_mag = np.sqrt(N_R**2 + N_Z**2)
    N_R_unit, N_Z_unit = N_R/N_mag, N_Z/N_mag

    # Offset contour
    R_offset = R + extrude * N_R_unit
    Z_offset = Z + extrude * N_Z_unit

    # Close polygon
    if not (np.isclose(R_offset[0], R_offset[-1]) and np.isclose(Z_offset[0], Z_offset[-1])):
        R_offset = np.append(R_offset, R_offset[0])
        Z_offset = np.append(Z_offset, Z_offset[0])

    # if calc_vol:
    #     # Shoelace area
    #     A = 0.5 * np.sum(R_offset[:-1]*Z_offset[1:] - R_offset[1:]*Z_offset[:-1])

    #     # Centroid in R
    #     Cx = (1/(6*A)) * np.sum((R_offset[:-1] + R_offset[1:]) *
    #                                (R_offset[:-1]*Z_offset[1:] - R_offset[1:]*Z_offset[:-1]))
    #     # Torus volume
    #     volume = 2*np.pi*Cx*abs(A)

    #     return R_offset, Z_offset, volume

    # return R_offset, Z_offset
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


def has_statepoint(directory_path):
    """
    Check if any file starting with 'statepoint' exists in the given directory.
    
    Args:
        directory_path (str): Path to the directory to search
    
    Returns:
        bool: True if a file starting with 'statepoint' is found, False otherwise
    """
    found = False
    for filename in os.listdir(directory_path):
        if filename.startswith("statepoint"):
            found = True
    return found
            


if __name__ == "__main__":
    """ Use this to test any of the functions """
    print(calc_blanket_mass_fracs(30, 291.1, fertile_element='U', fertile_enrich=0.71, breeder_density_kgm3=1.94e3))


