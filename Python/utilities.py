"""
Set of helper functions
"""
import os, re, sys
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


def calc_blanket_mass_fracs(fertile_density_kgm3, fertile_element='U', fertile_enrich=0.71, breeder_density_kgm3=1.94e3):
    """
    Calculate the mass fractions of FLiBe and (UF4 or ThF4). 
    Assumes UF4/ThF4 dissolution does NOT change FLiBe volume
      after conversaion with J.L. Ball --ppark 2025-07-03

    Args:
        fertile_density_kgm3 : float : desired density of U-238 or Th-232 in kg/m^3
        fertile_element      : 'U' or 'Th' : identifies U-238 or Th-232 as the fertile isotope
        fertile_enrich       : float : wt% enrichment of the fertile material
        breeder_density_kgm3 : float : density of FLiBe in kg/m^3

    Returns:
        (mass_frac_flibe, mass_frac_uf4) : 2-ple of floats : mass fractions
    """

    if fertile_element == 'U':
        uf4_density_kgm3   = fertile_density_kgm3 / (1-fertile_enrich) * AMU_UF4 / AMU_U
        blanket_mass_in_m3 = breeder_density_kgm3*1 + uf4_density_kgm3*1 # adding *1 for clarity that density = mass for 1 m^3
        # ^ this heavily assumes uf4 "solute" doesn't bigly change volume of molten salt "solvent"
        breeder_mass_frac  = breeder_density_kgm3/blanket_mass_in_m3 
        uf4_mass_frac      = uf4_density_kgm3/blanket_mass_in_m3
        return breeder_mass_frac, uf4_mass_frac

    if fertile_element == 'Th':
        thf4_density_kgm3  = fertile_density_kgm3 / (1-fertile_enrich) * AMU_ThF4 / AMU_Th
        blanket_mass_in_m3 = breeder_density_kgm3*1 + thf4_density_kgm3*1 # adding *1 for clarity that density = mass for 1 m^3
        # ^ this heavily assumes thf4 "solute" doesn't bigly change volume of molten salt "solvent"
        breeder_mass_frac  = breeder_density_kgm3/blanket_mass_in_m3 
        thf4_mass_frac     = thf4_density_kgm3/blanket_mass_in_m3
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
            sys.exit(f"Error finding cross section XML!")
    else:
        print(f"'OPENMC_CROSS_SECTIONS' found!")
        return xs_path


def miller_points(R0_m: float, a_m: float, delta: float, kappa: float, n: int = 400) -> np.ndarray:
    """Return Nx2 array of (r,z) points [cm] for a Miller D-shape.
    R0_m : major radius [m]
    a_m  : minor radius [m]
    delta: triangularity delta
    kappa: elongation kappa
    n    : number of points around the boundary
    """
    theta = np.linspace(0, 2*np.pi, n, endpoint=False)
    R = R0_m + a_m*np.cos(theta + delta*np.sin(theta))
    Z = kappa*a_m*np.sin(theta)
    pts_m = np.column_stack([R, Z])
    pts_cm = 100.0*pts_m
    # Close the loop (Polygon expects a closed path)
    return np.vstack([pts_cm, pts_cm[0]])


def polygon_area_centroid_cm(pts_cm: np.ndarray):
    """Compute signed area [cm^2] and centroid (Cx,Cy) [cm] of a polygon.
    Accepts Nx2 with optional repeated first point at the end.
    """
    if pts_cm.shape[1] != 2:
        raise ValueError("pts_cm must be of shape (N, 2) for (R,Z).")
    # Drop duplicated closing point if present
    if np.allclose(pts_cm[0], pts_cm[-1]):
        pts_cm = pts_cm[:-1]
    x = pts_cm[:, 0]
    y = pts_cm[:, 1]
    x2 = np.roll(x, -1)
    y2 = np.roll(y, -1)
    cross = x*y2 - x2*y
    A = 0.5*np.sum(cross)                    # signed area [cm^2]
    if abs(A) < 1e-20:
        raise ValueError("Degenerate polygon (area ~ 0).")
    Cx = (1.0/(6.0*A)) * np.sum((x + x2) * cross)  # centroid x [cm]
    Cy = (1.0/(6.0*A)) * np.sum((y + y2) * cross)  # centroid y [cm]
    return A, Cx, Cy

def torus_volume(R0_m: float, a_m: float, delta: float, kappa: float, n: int = 2048) -> float:
    """Return volume [m^3] of a Miller D-shaped torus via Pappus–Guldinus: V = 2π R̄ A.
    Uses `miller_points` (cm) → polygon area & centroid (cm) → convert cm^3 → m^3.
    Raises if the cross-section touches the axis (min R ≤ 0 cm).
    """
    pts_cm = miller_points(R0_m, a_m, delta, kappa, n=n)
    R_cm = pts_cm[:, 0]
    if np.min(R_cm) <= 0.0:
        raise ValueError("Cross-section intersects the rotation axis (min R ≤ 0).")
    A_cm2, Rbar_cm, _ = polygon_area_centroid_cm(pts_cm)
    A_cm2 = abs(A_cm2)
    V_cm3 = 2.0 * np.pi * Rbar_cm * A_cm2
    V_m3  = V_cm3 * 1e-6  # (cm^3 → m^3)
    return float(V_m3)

# ---------- Example ----------
if __name__ == "__main__":
    # Example: R0=6.2 m, a=2.0 m, δ=0.4, κ=1.6
    V = torus_volume(6.2, 2.0, 0.4, 1.6, n=10000)
    print(f"Volume ≈ {V:.3f} m^3")



