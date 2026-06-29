"""
Set of helper functions
"""
import os, sys, time, functools
import numpy as np
import pandas as pd
import openmc
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from Python.parameters import *
# from parameters import *


"""
Nuclear constants -- from atom.kaeri.re.kr/nuchart/
"""
AVO = 6.022141076e+23
AMU_LI6, AMU_LI7 = 6.0150, 7.0160 # 6.01512288742, 7.01600343426 # amu = g/mol
AMU_F19 = 18.9984 # 18.99840316207
AMU_BE9 =  9.0120 # 9.012183062

ENRICH_U = 0.71 # [wt%] needs to be in wt% for OpenMC argument
AMU_U235 = 235.0439299
AMU_U238 = 238.05078826
AMU_U    = (ENRICH_U * AMU_U235 + (100-ENRICH_U) * AMU_U238)/100
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


def log_buffer(lo, hi, buf=0.03):
    f = (hi / lo) ** buf
    return (lo / f, hi * f)


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


def read_and_sort(path:str):
    """
    Reads a CSV file and sorts it by 'fertile_kg/m3' if the column exists.
    """
    df = pd.read_csv(path)
    if 'fertile_kg/m3' in df.columns:
        df = df.sort_values(by='fertile_kg/m3').reset_index(drop=True)
    return df
    

def readtxtFile(path): 
    energy, microxs = [], []

    with open(path, 'r') as file:
        file.readline()
        file.readline()

        for line in file:
            values = line.split()
            energy.append(float(values[0]))
            microxs.append(float(values[1]))

    return np.array(energy), np.log(np.array(microxs))


def fitquartic(x,y):
    '''
    Fits [a, b, c, d] for ax^4 + bx^3 + cx^2 + dx = y (force zero production at zero enrichment)
    Return approximation expression, derivative of approximation
    
    x: array, enrichment amounts
    y: array, Pu-239 production per year 
    '''
    def quartic(x,a,b,c,d):
        return a*x**4 + b*x**3 + c*x**2 + d*x
    
    params, _ = curve_fit(quartic, x, y)
    a_fit,b_fit,c_fit,d_fit = params

    f = lambda x: a_fit*x**4 + b_fit*x**3 + c_fit*x**2 + d_fit*x
    fprime = lambda x: 4*a_fit*x**3 + 3*b_fit*x**2 + 2*c_fit*x + d_fit

    return f, fprime


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
    biso_per_cc     = fertile_kgm3_to_biso_per_cc(fertile_kgm3, fertile_isotope=fertile_isotope)
    vf_biso_breeder = biso_per_cc * BISO_VOLUME

    if vf_biso_breeder > 1.0:
        print(f"Fatal. Your fertile kg/m³ exceeds what can physically fit in the breeder volume!")
        print(f"Fatal. That is, your volume of BISO per cm³ of breeder volume exceeds 1.")
        sys.exit()

    """
    vf_biso_breeder tells us the X cm³ of BISO per 1 cm³ of breeder material.
    We want the volume fractions of BISO and breeder relative to the NOMINAL breeder volume, (1 + X) cm³.
    """
    vf_biso    = vf_biso_breeder / (vf_biso_breeder + 1)
    vf_breeder = 1 / (vf_biso_breeder + 1)

    return vf_biso, vf_breeder, biso_per_cc


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


def miller_model(R0, a, kappa, delta, extrude=0, calc=None, n=100):
    """Parametric coordinates from the Miller local equilibrium model,
    with optional toroidal volume or surface area via Pappus's theorem.

    Refs: R.L. Miller et al., doi.org/10.1063/1.872666
          (Justin) Ball et al., arxiv.org/pdf/1510.08923

    Args:
        R0 (float): Major radius [cm].
        a (float): Minor radius [cm].
        kappa (float): Elongation.
        delta (float): Triangularity.
        extrude (float): Outward extrusion thickness along surface normals [cm].
        calc (str | None): 'vol' for toroidal volume [cm³], 'area' for surface area [cm²], None for coordinates only.
        n (int): Number of (R, Z) pairs along the contour.

    Returns:
        list[tuple[float, float]]: Closed (R, Z) contour [cm].
        float: Toroidal volume [cm³], if calc='vol'.
        float: Toroidal surface area [cm²], if calc='area'.
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

    if calc == 'vol':
        # Shoelace area
        A = 0.5 * np.sum(R_offset[:-1]*Z_offset[1:] - R_offset[1:]*Z_offset[:-1])

        # Centroid in R
        Cx = (1/(6*A)) * np.sum((R_offset[:-1] + R_offset[1:]) *
                                (R_offset[:-1]*Z_offset[1:] - R_offset[1:]*Z_offset[:-1]))
        # Torus volume
        volume = 2*np.pi*Cx*abs(A)

        return R_offset, Z_offset, volume

    if calc == 'area':
        dR = np.diff(R_offset)
        dZ = np.diff(Z_offset)
        dl = np.sqrt(dR**2 + dZ**2)
        R_mid = 0.5 * (R_offset[:-1] + R_offset[1:])
        area = 2 * np.pi * np.sum(R_mid * dl)
        
        return R_offset, Z_offset, area

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


def has_statepoint(directory_path,cycle=None):
    """
    Check if any file starting with 'statepoint' exists in the given directory.
    
    Args:
        directory_path (str): Path to the directory to search
    
    Returns:
        bool: True if a file starting with 'statepoint' is found, False otherwise
    """
    found = False
    for filename in os.listdir(directory_path):
        if cycle:
            cycle_str = str(cycle).zfill(2)  # Pad cycle number with zeros to get 5 -> '05'
            if filename.startswith(f"statepoint.{cycle_str}"):
                found = True
        else:
            if filename.startswith("statepoint"):
                found = True
    return found


def sum_over_cells(df, nuclides, score):
    """Sum reaction rates over all cells per energy bin for given nuclide(s) and score.

    Args:
        df (pd.DataFrame): Spectrum tally DataFrame with 'nuclide', 'score', energy bin, 'mean', and 'std. dev.' columns.
        nuclides (str or list[str]): Nuclide name(s) (e.g. 'U238' or ['U235', 'U238']).
        score (str): Reaction score to filter on (e.g. 'fission').

    Returns:
        pd.DataFrame: One row per energy bin with columns 'mean' and 'std_dev'.

    Example:
        >>> fission_spectrum = sum_over_cells(fertile_spec, ['U235', 'U238'], 'fission')
    """
    if isinstance(nuclides, str):
        nuclides = [nuclides]
    mask = df['nuclide'].isin(nuclides) & (df['score'] == score)
    grouped = (df[mask]
                .groupby(['energy low [eV]', 'energy high [eV]', 'energy mid [eV]'])
                .agg(mean=('mean', 'sum'),
                    std_dev=('std. dev.', lambda x: np.sqrt((x**2).sum())))  # this names the column std_dev instead of std. dev. which gets fucky with OpenMC's outputs
                .rename(columns={'std_dev': 'std. dev.'})                    # so we change it back here --ppark 2026-05-20
                .reset_index()
                .sort_values('energy mid [eV]'))
    return grouped


def sum_over_nuclides(df, score):
    """Sum reaction rates over all nuclides per cell for a given score.

    Args:
        df (pd.DataFrame): Tally DataFrame with 'score', 'cell', 'mean', and 'std. dev.' columns.
        score (str): Reaction score to filter on (e.g. 'fission').

    Returns:
        pd.DataFrame: One row per cell with columns 'mean' and 'std_dev'.

    Example:
        >>> total_fission = sum_over_nuclides(fertile, 'fission')
    """
    return (df[df['score'] == score]
            .groupby('cell')
            .agg(mean=('mean', 'sum'),
                 std_dev=('std. dev.', lambda x: np.sqrt((x**2).sum()))) # this names the column std_dev instead of std. dev. which gets fucky with OpenMC's outputs
            .rename(columns={'std_dev': 'std. dev.'})                    # so we change it back here --ppark 2026-05-20
            .reset_index()
            .sort_values('cell'))


def find_sp(sp_dir, n_cycles):
    """
    Load the OpenMC statepoint from `path`. Tries the exact cycle count
    first, then falls back to the highest available statepoint.
 
    Args:
        sp_dir (str): Directory containing the statepoint HDF5 file(s).
        n_cycles (int): Expected number of batches (used to build the default filename).
 
    Returns:
        openmc.StatePoint
    """
        
    # First try finding statepoint that exactly matches the programmed n_cycles

    sp_path = os.path.join(sp_dir, f"statepoint.{str(n_cycles).zfill(2)}.h5")

    try:
        print(f"{C.YELLOW}Comment.{C.END} Loading statepoint: {sp_path}")
        return openmc.StatePoint(sp_path) 

    except Exception as e:
        print(f"{C.YELLOW}Warning.{C.END} {e}. statepoint.{n_cycles}.h5 missing or could not be read at: {sp_path}")
    
        # Second try finding highest statepoint
        try:
            print(f"{C.YELLOW}Comment.{C.END} Looking for the highest statepoint in: {sp_dir}")

            sp_files = [f for f in os.listdir(sp_dir) if f.startswith('statepoint.') and f.endswith('.h5')]

            if not sp_files:
                print(f"{C.RED}Warning.{C.END} <utilities.py/find_sp> No statepoints found at all in: {sp_dir}")
                print(f"But, reactor.py/openmc() should have run OpenMC if: (1) there are no statepoints (2) you have self.run_openmc=True (current state: {self.run_openmc})."
                        f"If you see this error, self.run_openmc=False or you might have broken something in our logic.")
                sys.exit(2)
                
            # Get the latest statepoint by batch number
            sp_path = os.path.join(sp_dir, max(sp_files, key=lambda x: int(x.split('.')[1])))
            return openmc.StatePoint(sp_path) 

        except Exception as e:
            print(f"{C.RED}  Fatal.{C.END} <reactor.py/extract_tallies> {e}. No valid statepoints found at all in: {sp_dir}")
            sys.exit(2)



if __name__ == "__main__":
    """ Use this to test any of the functions """
    # calc_biso_blanket_vol_fracs(fertile_kgm3, breeder_volume_m3, fertile_isotope='U238', fertile_enrich=0.71)
    # print(calc_biso_blanket_vol_fracs(1000, LL_BR_VOL, fertile_isotope='Th232'))
    # for 1000, LL_BR_VOL, 'U', 0.71: (0.7874738114760764, 0.21252618852392352)

    A_FW_FLiBe  = miller_model(FLIBE_R0, FLIBE_A, FLIBE_KAPPA, FLIBE_DELTA, extrude=(FLIBE_FW_CM), calc='area', n=10000)
    A_FW_DCLL   = miller_model(DCLL_R0, DCLL_A, DCLL_KAPPA, DCLL_DELTA,     extrude=(DCLL_FW_CM),  calc='area', n=10000)
    A_FW_HCPB   = miller_model(HCPB_R0, HCPB_A, HCPB_KAPPA, HCPB_DELTA,     extrude=(HCPB_FW_CM),  calc='area', n=10000)

    print(f"FW surface area for FLiBe: {(A_FW_FLiBe[2]/1e4):.2f} m²")   
    print(f"FW surface area for DCLL:  {(A_FW_DCLL[2]/1e4):.2f}  m²")
    print(f"FW surface area for HCPB:  {(A_FW_HCPB[2]/1e4):.2f}  m²")

    P_NEUTRON_MW = P_FUS_MW * 14.06 / 17.6
    print(P_NEUTRON_MW)

    print(f"Neutron wall loading for FLiBe: {(P_NEUTRON_MW / (A_FW_FLiBe[2]/1e4))} MW/m²")   
    print(f"Neutron wall loading for DCLL:  {(P_NEUTRON_MW / (A_FW_DCLL[2]/1e4))} MW/m²")
    print(f"Neutron wall loading for HCPB:  {(P_NEUTRON_MW / (A_FW_HCPB[2]/1e4))} MW/m²")


def spectrum_stats(E, y, mode_weight=None, quantiles=(0.10, 0.25, 0.75, 0.90)):
    """
    Characteristic energies of a binned spectrum, in the units of E.

    E           : bin-midpoint energies
    y           : per-bin weighting for mean, median & quantiles
                  (flux, leakage current, reaction contribution, ...)
    mode_weight : optional separate weighting whose peak defines the mode
                  (e.g. per-lethargy flux). Defaults to y.
    quantiles   : cumulative-y fractions whose energies are returned in `q`
                  (default 10/25/75/90 %).

      mean   = sum(E*y)/sum(y)            (weighted average)
      median = E at 50% of cumulative y   (interpolated)
      mode   = E of the largest mode_weight bin
      q[p]   = E at fraction p of cumulative y (interpolated), for each p

    Returns (mean_E, median_E, mode_E, q), where q is a dict {p: E_p}
    keyed by the requested fractions. All values are NaN if there is no
    positive y. (median_E equals the 0.50 point.)
    """
    E = np.asarray(E, dtype=float)
    y = np.asarray(y, dtype=float)
    w = y if mode_weight is None else np.asarray(mode_weight, dtype=float)

    good = np.isfinite(E) & np.isfinite(y) & (y > 0)
    if not good.any() or y[good].sum() == 0:
        return np.nan, np.nan, np.nan, {p: np.nan for p in quantiles}
    E_g, y_g, w_g = E[good], y[good], w[good]

    mean_E = np.sum(E_g * y_g) / np.sum(y_g)

    order = np.argsort(E_g)
    E_sorted = E_g[order]
    cum = np.cumsum(y_g[order])
    cum /= cum[-1]
    cum = np.concatenate(([0.0], cum))
    E_for_interp = np.concatenate(([E_sorted[0]], E_sorted))
    median_E = float(np.interp(0.5, cum, E_for_interp))
    q = {p: float(np.interp(p, cum, E_for_interp)) for p in quantiles}

    mode_E = float(E_g[np.argmax(w_g)])
    return mean_E, median_E, mode_E, q


def second_mode(E, y, source_eV=14.07e6, guard=1.3, min_prominence=0.1):
    """
    Second mode of a spectrum: the tallest *genuine* local maximum of y(E)
    lying below the D-T source peak near `source_eV`.

    Built for fusion-blanket flux spectra whose global mode is the
    14.07 MeV source peak; this returns the next peak down in energy
    (the slowing-down / scattered-neutron feature inside the plotted range).

    Peaks are detected on log10(y), so a peak's prominence is measured in
    decades to match a log-scaled axis. Any peak at E >= source_eV/guard is
    taken to be the source peak (or its shoulder) and discarded.
    `min_prominence` is the minimum topographic prominence, in decades, for
    a bump to count (this rejects Monte-Carlo noise). Among the surviving
    sub-source peaks the tallest (largest y) is returned.

    Note: maxima sitting exactly at the array edges (e.g. a thermal upturn
    at the lowest bin) are not detected by design.

    Returns (E2, y2) of the chosen peak, or (nan, nan) if none qualifies.
    """
    E = np.asarray(E, dtype=float)
    y = np.asarray(y, dtype=float)

    order = np.argsort(E)
    E, y = E[order], y[order]
    good = np.isfinite(E) & np.isfinite(y) & (y > 0)
    E, y = E[good], y[good]
    if E.size < 3:
        return np.nan, np.nan

    peaks, _ = find_peaks(np.log10(y), prominence=min_prominence)
    if peaks.size == 0:
        return np.nan, np.nan

    peaks = peaks[E[peaks] < source_eV / guard]      # keep only sub-source peaks
    if peaks.size == 0:
        return np.nan, np.nan

    best = peaks[np.argmax(y[peaks])]                # tallest secondary peak
    return float(E[best]), float(y[best])


def strip_latex(s):
    """Render simple matplotlib LaTeX markup as plain/unicode text for the console."""
    return (s.replace('$_4$', '4').replace('$_2$', '2')
             .replace('$^6$', '⁶').replace('$^7$', '⁷')
             .replace(r'$\gamma$', 'γ')
             .replace('$', ''))          # extend as more markup appears