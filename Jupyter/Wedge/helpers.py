import os
import numpy as np


# GLOBAL PARAMETERS
TEMP_K = 900
BISO_VOLUME = 4/3*np.pi*(0.05)**3
KERNEL_VOLUME = 4/3*np.pi*(0.04)**3

DENSITY_UO2  = 10.50 # [g/cm³]
DENSITY_ThO2 = 10.00 # [g/cm³]

AMU_U235 = 235.0439299
AMU_U238 = 238.05078826
AMU_U = 238.02891 # for natural enrichment
AMU_O = 15.999
AMU_UO2 = AMU_U + 2 * AMU_O
AMU_Th232 = 232.0381

# Volume fractions in HCPB blanket
VF_LI_NOM = 0.1304
VF_BE_NOM = 0.3790
VF_EU_NOM = 0.1176
VF_HE_NOM = 0.3730


def fertile_kgm3_to_biso_per_cc(fertile_kgm3, fertile_isotope='U238'):
    # vol_kernel = V_KERNEL # 4/3 * np.pi * 0.04**3
    if fertile_isotope == 'U238':
        biso_per_cc = fertile_kgm3 * AMU_UO2 / AMU_U238 / KERNEL_VOLUME / DENSITY_UO2 / 100**3 * 1000
    elif fertile_isotope == 'Th232':
        biso_per_cc = fertile_kgm3 * AMU_ThO2 / AMU_Th232 / KERNEL_VOLUME / DENSITY_ThO2 / 100**3 * 1000
    return biso_per_cc


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


def lattice_coords(lower_left, shape, pitch):
    """
    Calculates the center coordinates of all cells in a rectangular lattice
    and returns them as a list of tuples.
    
    Args:
        lower_left (iterable of float): (x, y, z) coordinates of the lower-left corner.
        shape (iterable of int): number of cells (Nx, Ny, Nz).
        pitch (iterable of float): width, height, and depth of each cell (Dx, Dy, Dz).
        
    Returns: 
        coords (list of tuple): flat list containing (x, y, z) center coordinates.
    """
    nx, ny, nz = shape
    dx, dy, dz = pitch
    lx, ly, lz = lower_left

    coords = []
    
    # Iterate through each dimension to calculate the center of each voxel
    for i in range(nx):
        x = lx + (i + 0.5) * dx
        for j in range(ny):
            y = ly + (j + 0.5) * dy
            for k in range(nz):
                z = lz + (k + 0.5) * dz
                coords.append((x, y, z))
                
    return coords


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