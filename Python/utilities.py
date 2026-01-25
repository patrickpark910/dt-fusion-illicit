"""
Set of helper functions
"""
import os, re, sys, time, functools
import numpy as np

from Python.parameters import *
# from parameters import *


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
AMU_Th232 = 232.0381
AMU_UF4 = AMU_U + 4 * AMU_F19
AMU_ThF4 = AMU_Th232 + 4 * AMU_F19
AMU_UO2 = AMU_U + 2 * AMU_O
AMU_ThO2 = AMU_Th232 + 2 * AMU_O
AMU_FLIBE = 98.89 # g/mol
AMU_PU239 = 239.0521634
AMU_U233 = 233.039635207
SEC_PER_YR = 3600 * 24 * 365

"""
Fusion parameters - from JL Ball 24
"""
N_PER_MJ = 3.546e17 # 17.6 MeV
P_FUS_MW = 1000
NPS_FUS  = P_FUS_MW * N_PER_MJ # n/s



class C:
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


def fertile_kgm3_to_biso_per_cc(fertile_kgm3, fertile_isotope='U238'):
    """
    Converts given kg/m³ of fertile material to the number of BISO particles per cm³

    Args:
        fertile_kgm3  (float): kg of fertile isotope per m³ of breeder
        fertile_isotope (str): one of ['U238', 'Th232']

    Returns:
        biso_per_cc (float): number of biso particles per m³ of breeder

    """
    if fertile_isotope == 'U238':
        biso_per_cc = fertile_kgm3 * AMU_UO2 / AMU_U238 / KERNEL_VOLUME / DENSITY_UO2 / 100**3 * 1000
    elif fertile_isotope == 'Th232':
        biso_per_cc = fertile_kgm3 * AMU_ThO2 / AMU_Th232 / KERNEL_VOLUME / DENSITY_ThO2 / 100**3 * 1000
    return biso_per_cc


def calc_biso_breeder_vol_fracs(fertile_kgm3, fertile_isotope='U238'):
    """
    Calculate volume fractions of BISO particles and breeder material.
    Per Glaser & Goldston (2012), we assume the BISO/TRISO particles are homogenized 
    throughout the blanket. The kernel/coating geometry is used only to set the 
    relative mass or volume fractions of fertile vs. coating material.
    
    Args:
        fertile_kgm3  (float): kg of fertile isotope per m³ of breeder
        fertile_isotope (str): one of ['U238', 'Th232']

    Returns:
        vf_biso_br (float): vol frac of BISO relative to nominal breeder volume
        vf_breeder_br (float): new vol frac of breeder relative to nominal breeder volume
        biso_per_cc_br (float): number of BISO spheres per cm³ of breeder
    """
    # Number of BISO spheres per cc of breeder
    biso_per_cc_br = fertile_kgm3_to_biso_per_cc(fertile_kgm3, fertile_isotope=fertile_isotope)
    vf_biso_br     = biso_per_cc_br * BISO_VOLUME

    if vf_biso_br > 1.0:
        print(f"Fatal. Your fertile kg/m³ exceeds what can physically fit in the breeder volume!")
        print(f"Fatal. That is, your volume of BISO per cm³ of breeder volume exceeds 1.")
        sys.exit()

    # New volume ratio of breeder relative to nominal breeder volume
    vf_breeder_br = 1 - vf_biso_br

    return vf_biso_br, vf_breeder_br, biso_per_cc_br


def set_xs_path():
    """
    Temporary solution for finding xs files between WSL and Ubuntu on Computing Cluster without editing PATH --ppark 2025-06-28
    """
    xs_path_zotac = '/opt/openmc_data/endfb-viii.0-hdf5/cross_sections.xml'
    xs_path_wsl   = '/mnt/c/OpenMC/data/endfb-viii.0-hdf5/cross_sections.xml'

    xs_path = os.environ.get("OPENMC_CROSS_SECTIONS")
    

    if xs_path is None:
        print(f"{C.YELLOW}Warning.{C.END} 'OPENMC_CROSS_SECTIONS' is NOT set in PATH.")
        
        if os.path.isfile(xs_path_zotac):
            print(f"{C.YELLOW}Warning.{C.END} 'OPENMC_CROSS_SECTIONS' is now set to {xs_path_ubuntu}.")
            return xs_path_zotac # use this on Zotacs --ppark

        elif os.path.isfile(xs_path_wsl):
            print(f"{C.YELLOW}Warning.{C.END} 'OPENMC_CROSS_SECTIONS' is now set to {xs_path_wsl}.")
            return xs_path_wsl
        
        else:
            sys.exit(f"{C.RED}Fatal. Cannot find cross section XML!{C.END}")
    
    else:
        print(f"{C.GREEN}'OPENMC_CROSS_SECTIONS' found!{C.END}")
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
    # calc_biso_blanket_vol_fracs(fertile_kgm3, breeder_volume_m3, fertile_isotope='U238', fertile_enrich=0.71)
    print(calc_biso_blanket_vol_fracs(1000, LL_BR_VOL, fertile_isotope='Th232'))
    # for 1000, LL_BR_VOL, 'U', 0.71: (0.7874738114760764, 0.21252618852392352)


