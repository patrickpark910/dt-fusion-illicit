import os
import sys
import numpy as np
import openmc

# export OPENMC_CROSS_SECTIONS=/Users/gretali/Desktop/research2025/FissileDependence/endfb-viii.0-hdf5/cross_sections.xml

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

# Lead-lithium (DCLL, 83 at% Pb - 17 at% Li)
DENSITY_DCLL    =  9.40  # [g/cm³]
ENRICH_DCLL     = 90.00  # [at%] enrich of Li-6

# Volume fractions in DCLL blanket (Glaser & Goldston 2012, Tb.1)
DCLL_VF_FS_NOM = 0.019 
DCLL_VF_LL_NOM = 0.808  
DCLL_VF_SI_NOM = 0.076
DCLL_VF_HE_NOM = 0.097


def fertile_kgm3_to_biso_per_cc(fertile_kgm3, fertile_isotope='U238'):
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
    
    for i in range(nx):
        x = lx + (i + 0.5) * dx
        for j in range(ny):
            y = ly + (j + 0.5) * dy
            for k in range(nz):
                z = lz + (k + 0.5) * dz
                coords.append((x, y, z))
                
    return coords



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


class Wedge():

    def __init__(self, case, fertile_kgm3, isotope='U238', breeder_enrich=90.0, write_openmc=True, run_openmc=False,):

        self.case = case
        self.fertile_kgm3 = fertile_kgm3
        self.fertile_str = f"{fertile_kgm3:06.2f}"
        self.breeder_enr_str = f"{breeder_enrich:04.1f}"
        self.isotope = isotope 
        self.name = f"wedge{self.case}_hcpb_{self.isotope}_{self.fertile_str}kgm3"         
        self.path = f"./OpenMC_{self.isotope}/{self.name}"

        os.makedirs(self.path, exist_ok=True)


    def materials(self, debug=False):

        # Tungsten | First wall
        self.tungsten = openmc.Material(name='firstwall', temperature=TEMP_K)
        self.tungsten.set_density('g/cm3',19.0)
        self.tungsten.add_element('W',1)
        self.tungsten.depletable = False

        # Structure | Eurofer
        self.eurofer = openmc.Material(name='Eurofer', temperature=TEMP_K)
        self.eurofer.depletable = False
        self.eurofer.set_density('g/cm3', 7.80)
        self.eurofer.add_element('Fe', 89.36, percent_type='wo')
        self.eurofer.add_element('C' ,  0.11, percent_type='wo')
        self.eurofer.add_element('Cr',  9.00, percent_type='wo')
        self.eurofer.add_element('W' ,  1.10, percent_type='wo')
        self.eurofer.add_element('Mn',  0.40, percent_type='wo')
        self.eurofer.add_element('N' ,  0.03, percent_type='wo')
        # self.eurofer.add_s_alpha_beta('c_Fe56', 0.89006) # 89.36 wt% Fe = 89.006 at% Fe
        # OpenMC v.0.15.3. NotImplementedError: Currently we do not support mixing materials containing S(a,b) tables

        # He-4 (gas) | coolant
        he = openmc.Material(name='Helium') 
        he.set_density('g/cm3', 0.004279) # Helium density at 900 K ~80 bar 
        he.add_element('He', 1) 


        ''' Breeder material '''
        # Li4SiO4 | Tritium breeder
        li4sio4 = openmc.Material(name='Li4SiO4', temperature=TEMP_K) 
        li4sio4.set_density('g/cm3', 2.42)
        li4sio4.add_elements_from_formula('Li4SiO4', enrichment_target='Li6', enrichment_type='ao', enrichment=60.0)

        # Beryllium | Neutron multiplier
        be = openmc.Material(name='Beryllium', temperature=TEMP_K) 
        be.set_density('g/cm3', 1.85) 
        be.add_element('Be', 1, percent_type='wo')     
        # be.add_s_alpha_beta('c_Be')  # NotImplementedError   


        ''' BISO particle '''
        # Fertile material
        if self.isotope == 'U238':
            self.kernel = openmc.Material(name='UO2', temperature=TEMP_K)
            self.kernel.add_elements_from_formula('UO2')
            self.kernel.set_density('g/cm3', DENSITY_UO2)  

        elif self.isotope == 'Th232': 
            self.kernel = openmc.Material(name='ThO2', temperature=TEMP_K) 
            self.kernel.add_elements_from_formula('ThO2') 
            self.kernel.set_density('g/cm3', DENSITY_ThO2)   

        # SiC coating
        self.sic = openmc.Material(name='SiC', temperature=TEMP_K)
        self.sic.add_elements_from_formula('SiC')
        self.sic.set_density('g/cm3', 3.2)

        if case == 'A': # homogeneous 
            biso = openmc.Material.mix_materials([self.kernel, self.sic], [KERNEL_VOLUME/BISO_VOLUME, (1-(KERNEL_VOLUME/BISO_VOLUME))], 'vo') 
            self.blanket = openmc.Material.mix_materials([biso, pbli, self.f82h, self.sic, he], [self.vf_biso_br, self.vf_pbli_br, DCLL_VF_FS_NOM, DCLL_VF_SI_NOM, DCLL_VF_HE_NOM], 'vo') 
            self.blanket.temperature = TEMP_K
            self.blanket.name = (f"{self.fertile_str} kg/m3"
                                 f" | {self.biso_per_cc_bv:.4f} spheres/cc = {(self.vf_biso_bv*100):.4f} vol% in breeder volume"
                                 f" | {self.biso_per_cc_br:.4f} spheres/cc = {(self.vf_biso_br*100):.4f} vol% in breeding region")
            self.blanket.temperature = TEMP_K

        ''' Adjust materials according to case '''

        if case == 'A':

            biso = openmc.Material.mix_materials([self.kernel, self.sic], [KERNEL_VOLUME/BISO_VOLUME, (1-(KERNEL_VOLUME/BISO_VOLUME))], 'vo') 

            breeder = openmc.Material.mix_materials([li4sio4, be, biso], [self.vf_li_bv, self.vf_be_bv, self.vf_biso_bv], 'vo') 
            self.blanket = openmc.Material.mix_materials([breeder, self.eurofer, he], [(VF_LI_NOM+VF_BE_NOM), VF_EU_NOM, VF_HE_NOM], 'vo') 
            self.blanket.name = (f"{self.fertile_str} kg/m3"
                                 f" | {self.biso_per_cc_bv:.4f} spheres/cc = {(self.vf_biso_bv*100):.4f} vol% in breeder volume"
                                 f" | {self.biso_per_cc_br:.4f} spheres/cc = {(self.vf_biso_br*100):.4f} vol% in breeding region")
            self.blanket.temperature = TEMP_K

            self.materials = openmc.Materials([self.tungsten, self.eurofer, self.blanket]) 

            if debug:
                print(f"BISO materials mix")
                print(f"  Vol frac kernel : {(KERNEL_VOLUME/BISO_VOLUME):.6f}")
                print(f"  Vol frac coat   : {(1-(KERNEL_VOLUME/BISO_VOLUME)):.6f}")
                print(f"")
                print(f"Breeding material mix")
                print(f"  Vol frac Li4SiO4: {self.vf_li_bv:.6f}")
                print(f"  Vol frac Be(met): {self.vf_be_bv:.6f}")
                print(f"  Vol frac BISO   : {self.vf_biso_bv:.6f}")
                print(f"")
                print(f"Blanket mix")
                print(f"  Vol frac breeder: {(VF_LI_NOM+VF_BE_NOM):.6f}")
                print(f"  Vol frac Eurofer: {VF_EU_NOM:.6f}")
                print(f"  Vol frac He-4(g): {VF_HE_NOM:.6f}")
                print(f"")


        elif case in ['B','C']:

            self.blanket = openmc.Material.mix_materials([li4sio4, be, self.eurofer, he], [VF_LI_NOM, VF_BE_NOM, VF_EU_NOM, VF_HE_NOM], 'vo') 
            self.blanket.temperature = TEMP_K

            self.materials = openmc.Materials([self.tungsten, self.eurofer, self.blanket, self.kernel, self.sic]) 



    def geometry(self, debug=False):

        # -----------------------------------
        # SURFACES
        # -----------------------------------

        # X-planes
        s111 = openmc.XPlane(surface_id=111, x0= -self.geom['c_in'], boundary_type="reflective")
        s112 = openmc.XPlane(surface_id=112, x0= +self.geom['c_in'], boundary_type="reflective")
        s113 = openmc.XPlane(surface_id=113, x0= -self.geom['c_out'], boundary_type="reflective")
        s114 = openmc.XPlane(surface_id=114, x0= +self.geom['c_out'], boundary_type="reflective")

        # Y-planes
        s121 = openmc.YPlane(surface_id=121, y0= -self.geom['c_in'], boundary_type="reflective")
        s122 = openmc.YPlane(surface_id=122, y0= +self.geom['c_in'], boundary_type="reflective")
        s123 = openmc.YPlane(surface_id=123, y0= -self.geom['c_out'], boundary_type="reflective")
        s124 = openmc.YPlane(surface_id=124, y0= +self.geom['c_out'], boundary_type="reflective")

        # Z-planes
        s130 = openmc.ZPlane(surface_id=130, z0=  -2.50, boundary_type="vacuum")  # 
        s131 = openmc.ZPlane(surface_id=131, z0=   0.00)  #   2.5 cm - inboard outer eurofer 
        s132 = openmc.ZPlane(surface_id=132, z0=  45.00)  #  45.0 cm - inboard breeding region
        s133 = openmc.ZPlane(surface_id=133, z0=  47.50)  #   2.5 cm - inboard inner eurofer
        s134 = openmc.ZPlane(surface_id=134, z0=  47.70)  #   0.2 cm - inboard W first wall
        s135 = openmc.ZPlane(surface_id=135, z0= 447.70)  # 400.0 cm - plasma chamber
        s136 = openmc.ZPlane(surface_id=136, z0= 447.90)  #   0.2 cm - outboard W first wall 
        s137 = openmc.ZPlane(surface_id=137, z0= 450.40)  #   2.5 cm - outboard inner eurofer
        s138 = openmc.ZPlane(surface_id=138, z0= 532.40)  #  82.0 cm - outboard breeding region 
        s139 = openmc.ZPlane(surface_id=139, z0= 534.90, boundary_type="vacuum")  #   2.5 cm - outboard outer eurofer

        # Plasma chamber boundaries
        # Plane normal card entries: Ax + By + Cz - D = 0
        s171 = openmc.Plane(surface_id=171, a= +1.0, b=  0.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s172 = openmc.Plane(surface_id=172, a= -1.0, b=  0.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s173 = openmc.Plane(surface_id=173, a=  0.0, b= +1.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s174 = openmc.Plane(surface_id=174, a=  0.0, b= -1.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")

        # BISO universe
        s150 = openmc.Sphere(surface_id=150, r=0.040) # Kernel
        s151 = openmc.Sphere(surface_id=151, r=0.050) # Outer Coating
        # s152 = openmc.Sphere(surface_id=152, r=0.0501) # Outer universe boundary


        # -----------------------------------
        # CELLS
        # -----------------------------------

        # inboard cells 
        c20 = openmc.Cell(cell_id=20, region= +s111 & -s112 & +s121 & -s122 & +s130 & -s131, fill=self.shield) # shield 
        c21 = openmc.Cell(cell_id=21, region= +s111 & -s112 & +s121 & -s122 & +s131 & -s132, fill=self.f82h) # plate 
        c22 = openmc.Cell(cell_id=22, region= +s111 & -s112 & +s121 & -s122 & +s132 & -s133, fill=self.manifold) # manifold
        # c23 breeding region 1 
        c24 = openmc.Cell(cell_id=24, region= +s111 & -s112 & +s121 & -s122 & +s134 & -s135, fill=self.divider) # divider
        # c25 breeding region 2
        c26 = openmc.Cell(cell_id=26, region= +s111 & -s112 & +s121 & -s122 & +s136 & -s137, fill=self.f82h) # f82h
        c27 = openmc.Cell(cell_id=27, region= +s111 & -s112 & +s121 & -s122 & +s137 & -s138, fill=self.coolant) # coolant 
        c28 = openmc.Cell(cell_id=28, region= +s111 & -s112 & +s121 & -s122 & +s138 & -s139, fill=self.f82h) # f82h
        c29 = openmc.Cell(cell_id=29, region= +s111 & -s112 & +s121 & -s122 & +s139 & -s140, fill=self.firstwall) # tungsten

        # plasma chamber void 
        c10 = openmc.Cell(cell_id=10, region= -s171 & -s172 & -s173 & -s174 & +s140 & -s141)
        c10.importance = {'neutron':1}

        # outboard cells
        c30 = openmc.Cell(cell_id=30, region= +s113 & -s114 & +s123 & -s124 & +s141 & -s150, fill=self.firstwall) # tungsten
        c31 = openmc.Cell(cell_id=31, region= +s113 & -s114 & +s123 & -s124 & +s150 & -s151, fill=self.f82h) # f82h
        c32 = openmc.Cell(cell_id=32, region= +s113 & -s114 & +s123 & -s124 & +s151 & -s152, fill=self.coolant) # coolant
        c33 = openmc.Cell(cell_id=33, region= +s113 & -s114 & +s123 & -s124 & +s152 & -s153, fill=self.f82h) # f82h
        # c34 breeding region 1
        c35 = openmc.Cell(cell_id=35, region= +s113 & -s114 & +s123 & -s124 & +s154 & -s155, fill=self.divider) # divider 1
        # c36 breeding region 2
        c37 = openmc.Cell(cell_id=37, region= +s113 & -s114 & +s123 & -s124 & +s156 & -s157, fill=self.divider) # divider 2
        # c38 breeding region 3
        c39 = openmc.Cell(cell_id=39, region= +s113 & -s114 & +s123 & -s124 & +s158 & -s159, fill=self.manifold) # manifold
        c40 = openmc.Cell(cell_id=40, region= +s113 & -s114 & +s123 & -s124 & +s159 & -s160, fill=self.f82h) # plate
        c41 = openmc.Cell(cell_id=41, region= +s113 & -s114 & +s123 & -s124 & +s160 & -s161, fill=self.shield) # shield

        if case in ['A']:
            # inboard 
            c23 = openmc.Cell(cell_id=23, region= +s111 & -s112 & +s121 & -s122 & +s133 & -s134, fill=self.blanket)
            c25 = openmc.Cell(cell_id=25, region= +s111 & -s112 & +s121 & -s122 & +s135 & -s136, fill=self.blanket)

        if case in ['A']:

            c13 = openmc.Cell(cell_id=13, region= +s111 & -s112 & +s121 & -s122 & +s131 & -s132, fill=self.blanket)
            c23 = openmc.Cell(cell_id=23, region= +s113 & -s114 & +s123 & -s124 & +s137 & -s138, fill=self.blanket)
            

        elif case in ['B', 'C']:
            # Common setup for B and C: Define the BISO Universe
            c13 = openmc.Cell(cell_id=13, region= +s111 & -s112 & +s121 & -s122 & +s131 & -s132)
            c23 = openmc.Cell(cell_id=23, region= +s113 & -s114 & +s123 & -s124 & +s137 & -s138)

            # outboard 
            c34 = openmc.Cell(cell_id=34, region= +s113 & -s114 & +s123 & -s124 & +s153 & -s154)
            c36 = openmc.Cell(cell_id=36, region= +s113 & -s114 & +s123 & -s124 & +s155 & -s156)
            c38 = openmc.Cell(cell_id=38, region= +s113 & -s114 & +s123 & -s124 & +s157 & -s158)

            # single BISO
            c50 = openmc.Cell(cell_id=50, fill=self.kernel,  region= -s180)
            c51 = openmc.Cell(cell_id=51, fill=self.sic,     region= +s180 & -s181) 
            u50 = openmc.Universe(name='BISO Universe', cells=[c50, c51]) 

            # Lower left (ll), upper right (ur) of inboard, outboard bounding boxes
            ll_in_23, ur_in_23   = c23.region.bounding_box
            ll_in_25, ur_in_25   = c25.region.bounding_box
            ll_out_34, ur_out_34 = c34.region.bounding_box
            ll_out_36, ur_out_36 = c36.region.bounding_box
            ll_out_38, ur_out_38 = c38.region.bounding_box 

            pitch_in_23 = (ur_in_23 - ll_in_23) / self.geom['shape_in_per']
            pitch_in_25 = (ur_in_25 - ll_in_25) / self.geom['shape_in_per']
            pitch_out_34 = (ur_out_34 - ll_out_34) / self.geom['shape_out_per']
            pitch_out_36 = (ur_out_36 - ll_out_36) / self.geom['shape_out_per']
            pitch_out_38 = (ur_out_38 - ll_out_38) / self.geom['shape_out_per']

            shape_in  = (4, 4, 256)  # (Nx, Ny, Nz)
            shape_out = (4, 4, 750)

            pitch_in  = (ur_in - ll_in) / shape_in
            pitch_out = (ur_out - ll_out) / shape_out


            if case == 'B':

                # Random packing -- NB. 'pack_spheres' returns a list of coordinates, not cell instances
                centers_in  = openmc.model.pack_spheres(radius=0.05, region= +s111 & -s112 & +s121 & -s122 & +s131 & -s132, num_spheres=self.geom['N1'])
                centers_out = openmc.model.pack_spheres(radius=0.05, region= +s113 & -s114 & +s123 & -s124 & +s137 & -s138, num_spheres=self.geom['N2'])
                
                # openmc.model.TRISO(outer_radius, fill, center=(0.0, 0.0, 0.0)) -- each model.TRISO instance IS a cell
                biso_in  = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_in]
                biso_out = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_out]

                c13.fill = openmc.model.create_triso_lattice(biso_in, ll_in, pitch_in, shape_in, self.blanket)
                c23.fill = openmc.model.create_triso_lattice(biso_out, ll_out, pitch_out, shape_out, self.blanket)


            elif case == 'C':
                
                # Rectangular lattice
                centers_in  = lattice_coords(ll_in, shape_in, pitch_in)
                centers_out = lattice_coords(ll_out, shape_out, pitch_out)
                
                # openmc.model.TRISO(outer_radius, fill, center=(0.0, 0.0, 0.0)) -- each model.TRISO instance IS a cell
                biso_in_23  = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_23]
                biso_in_25  = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_25]
                biso_out_34 = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_34]
                biso_out_36 = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_36]
                biso_out_38 = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_38]

                # fill this mega cells with blanket material 
                c23.fill = openmc.model.create_triso_lattice(biso_in_23, ll_in_23, pitch_in_23, self.geom['shape_in_per'], self.blanket)
                c25.fill = openmc.model.create_triso_lattice(biso_in_25, ll_in_25, pitch_in_25, self.geom['shape_in_per'], self.blanket)
                c34.fill = openmc.model.create_triso_lattice(biso_out_34, ll_out_34, pitch_out_34, self.geom['shape_out_per'], self.blanket)
                c36.fill = openmc.model.create_triso_lattice(biso_out_36, ll_out_36, pitch_out_36, self.geom['shape_out_per'], self.blanket)
                c38.fill = openmc.model.create_triso_lattice(biso_out_38, ll_out_38, pitch_out_38, self.geom['shape_out_per'], self.blanket)

                if debug:
                    print(f"Case C: Fixed lattice BISO placement")
                    print(f"  Inboard lattice:")
                    print(f"    Shape: {self.geom['shape_in_per']}")
                    print(f"    Pitch: {pitch_in_23} cm")
                    print(f"    Lower left: {ll_in_23}")
                    print(f"    Upper right: {ur_in_23}")
                    print(f"    Subtotal one inboard breeding region BISO particles: {np.prod(self.geom['shape_in_per'])}, x2 = {np.prod(self.geom['shape_in_per']) * 2}")
                    print(f"  Outboard lattice:")
                    print(f"    Shape: {self.geom['shape_out_per']}")
                    print(f"    Pitch: {pitch_out_34} cm")
                    print(f"    Lower left: {ll_out_34}")
                    print(f"    Upper right: {ur_out_34}")
                    print(f"    Subtotal one outboard breeding region BISO particles: {np.prod(self.geom['shape_out_per'])}, x3 = {np.prod(self.geom['shape_out_per']) * 3}")
                    print(f"  Total BISO particles: {np.prod(self.geom['shape_in_per'])*2 + np.prod(self.geom['shape_out_per'])*3}")

        self.cells = [c10, c20, c21, c22, c23, c24, c25, c26, c27, c28, c29, c30, c31, c32, c33, c34, c35, c36, c37, c38, c39, c40, c41]
        self.geometry = openmc.Geometry(openmc.Universe(cells=self.cells))


    def tallies(self, debug=False):

        self.tallies = openmc.Tallies() # initialize
        cell_filter = openmc.CellFilter([cell for cell in self.cells if cell.fill is not None])
        energy_filter = openmc.EnergyFilter(logspace_per_decade(1e-5, 20e6, 100)) # './helpers/utilities.py'

        # Flux tally 
        flux_tally        = self.make_tally('flux', ['flux'], filters=[cell_filter])
        flux_energy_tally = self.make_tally('flux spectrum', ['flux'], filters=[cell_filter, energy_filter])

        # Fertile element reaction rates
        fertile_tally_tot    = self.make_tally(f'Total fertile rxn rate',     ['(n,gamma)', 'fission', '(n,2n)', 'elastic'], nuclides=['U238', 'U235', 'Th232'])
        fertile_tally        = self.make_tally(f'Fertile rxn rates',          ['(n,gamma)', 'fission', '(n,2n)', 'elastic'], filters=[cell_filter], nuclides=['U238', 'U235', 'Th232'])
        fertile_energy_tally = self.make_tally(f'Fertile rxn rates spectrum', ['(n,gamma)', 'fission', '(n,2n)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['U238', 'U235', 'Th232'])

        # Lithium reaction rates
        Li_tally_tot    = self.make_tally('Total (n,t) rxn rate',  ['(n,Xt)'])
        Li_tally        = self.make_tally('Li rxn rates',          ['(n,Xt)', '(n,gamma)', 'elastic'], filters=[cell_filter], nuclides=['Li6', 'Li7'])
        Li_energy_tally = self.make_tally('Li rxn rates spectrum', ['(n,Xt)', '(n,gamma)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['Li6', 'Li7'])

        # Put in order you want it to print in... recommend fertile_tally's first
        self.tallies.extend([fertile_tally_tot, Li_tally_tot, fertile_tally, Li_tally, flux_tally, fertile_energy_tally, flux_energy_tally, Li_energy_tally]) 


    def make_tally(self, name, scores, filters:list=None, nuclides:list=None):
        tally = openmc.Tally(name=name)
        tally.scores = scores
        if filters is not None:
            tally.filters = filters
        if nuclides is not None:
            tally.nuclides = nuclides
        return tally


    def settings(self, debug=False):

        # 14-MeV point source at center of plasma chamber
        source = openmc.IndependentSource()
        # source.space = openmc.stats.Point((0,0,247.7))                                                           
        source.space = openmc.stats.CartesianIndependent(x=openmc.stats.Uniform(a= -self.geom['c_in'], b= +self.geom['c_in']),
                                                         y=openmc.stats.Uniform(a= -self.geom['c_in'], b= +self.geom['c_in']),
                                                         z=openmc.stats.Discrete([277.7], [1.0]) )
        source.particle = 'neutron'
        source.energy   = openmc.stats.Discrete([14.0e6], [1.0])  # 14 MeV
        
        # Define the angular distribution -- mu=1 is straight up (+z), mu=-1 is straight down (-z)
        # We give each a 50% (0.5) probability.
        source.angle = openmc.stats.PolarAzimuthal(
                           mu=openmc.stats.Discrete([1.0, -1.0], [0.5, 0.5]),
                           phi=openmc.stats.Uniform(0, 2*np.pi) # azimuthal angle doesn't matter for mu=+/-1
                           )
        source.angle = openmc.stats.Isotropic()

        # Create settings and assign source
        self.settings = openmc.Settings()
        self.settings.source = source

        # Run type
        # if self.fertile_kgm3 < 0.09:
        #     nps, batches = int(1e5), int(10)
        # elif 0.09 < self.fertile_kgm3 < 10.0:
        #     nps, batches = int(1e5), int(10)
        # else:
        #     nps, batches = int(1e5), int(10)

        self.settings.run_mode  = 'fixed source'
        self.settings.particles = int(4e4) # nps
        self.settings.batches   = int(25)  # batches

        # self.settings.trace = (1,1,1)
        self.settings.max_tracks = 4


    def prism_helpers(self, debug=True):

        fertile_kgm3 = self.fertile_kgm3

        target_total_biso = 16096 # 20120 
        N1, N2 = 4096, 12000      # 5120, 15000

        # N_in/N_out = breeder_length_in/breeder_length_out * f1/f2
        N_in, N_out = 6240, 9856     # inboard is 39% of total blanket volume, just like the fluence fraction
        N_in_per    = int(N_in)  / 2      # per single breeding inboard region
        N_out_per   = int(N_out) / 3      # per single breeding outboard region

        shape_in_per  = (4, 4, int(N_in_per/16))  # per inboard breeding region
        shape_out_per = (4, 4, int(N_out_per/16)) # per outboard breeding region

        # DCLL tokamak dimensions (along z in model)
        f1, f2        = 0.39, 0.61  # fraction of fluence into inboard vs. outboard in full OpenMC DCLL model
        thickness_per = 21          # thickness of each breeding region in cm
        thickness_in  = 2 * thickness_per
        thickness_out = 3 * thickness_per


        biso_per_cc_bv = fertile_kgm3_to_biso_per_cc(fertile_kgm3) # Number of BISO spheres per cc of breeding volume
        vf_biso_bv = biso_per_cc_bv * BISO_VOLUME # volume of BISO spheres per cc of breeding volume 
        if vf_biso_bv > 1.0:
            print(f"Fatal. Your fertile kg/m³ exceeds what can physically fit in the breeding volume!")
            print(f"Fatal. That is, your volume of BISO per cm³ of breeding volume exceeds 1.")

        vf_pbli_bv = 1 - vf_biso_bv # new volume ratio wrt breeders in breeder volume 
        vf_checksum1 = vf_biso_bv + vf_pbli_bv # should equal 1
        vf_biso_br, vf_pbli_br = vf_biso_bv * DCLL_VF_LL_NOM, vf_pbli_bv * DCLL_VF_LL_NOM
        vf_checksum2 = vf_biso_br + vf_pbli_br + DCLL_VF_FS_NOM + DCLL_VF_SI_NOM + DCLL_VF_HE_NOM  # should equal 1 

        biso_per_cc_br = vf_biso_br / BISO_VOLUME # Number of BISO spheres per cc of breeding region

        A_ref = target_total_biso / biso_per_cc_br / (thickness_in + thickness_out) 
        A_out = A_ref * (thickness_in+thickness_out) / ( f1/f2 * thickness_in + thickness_out) # outboard XY face area 
        A_in = A_out * f1/f2    # inboard XY face area                                            
        b_in, b_out = np.round(np.sqrt(A_in), 6), np.round(np.sqrt(A_out), 6)  # face lengths
        c_in, c_out = b_in/2, b_out/2  # face HALF-lengths

        V_in, V_in_per = A_in * thickness_in, A_in * thickness_per
        V_out, V_out_per = A_out * thickness_out, A_out * thickness_per # volumes of breeding regions

        if debug:
            print(40*f"=")
            print(f"Case: {fertile_kgm3} kg/cm³")
            print(f"")
            print(f"Dimensions (inboard, outboard):")
            print(f" XY face half-lengths c_in, c_out: {c_in:.6f}, {c_out:.6f} [cm]")
            print(f" XY face lengths      b_in, b_out: {b_in:.6f}, {b_out:.6f} [cm]")
            print(f" XY face areas        A_in, A_out: {A_in:.6f}, {A_out:.6f} [cm²]")
            print(f" Breeding region vols V_in, V_out: {V_in:.6f}, {V_out:.6f} [cm³]")
            print(f" Fraction of fluence  f1, f2: {f1:.6f}, {f2:.6f} ")
            print(f"")
            print(f"'breeder' = Pb-Li ")
            print(f"'breeding region' (br) = the physical blanket layer that contains the breeder ")
            print(f"                         this layer usually contains breeder + structure + coolant ")
            print(f"'breeding volume' (bv) = volume nominally occupied in breeding region by breeder ")
            print(f"")
            print(f"With respect to the BREEDER material, i.e., per 1 m³ of breeder volume, we have these new volume fractions:")
            print(f"  vf_biso_breeder_new =  {(vf_biso_bv*100):.6f} vol%")
            print(f"  vf_pbli_breeder_new   =  {(vf_pbli_bv*100):.6f} vol%")
            print(f"  check they add up   = {(vf_checksum1*100):.6f} vol%")
            print(f"")
            print(f"With respect to the whole BREEDING REGION, we have these new volume fractions:")
            print(f"  vf_biso_br =  {(vf_biso_br*100):.6f} vol%")
            print(f"  vf_pbli_br =  {(vf_pbli_br*100):.6f} vol%")
            print(f"  VF_FS_NOM  =  {(DCLL_VF_FS_NOM*100):.6f} vol%")
            print(f"  VF_SI_NOM  =  {(DCLL_VF_SI_NOM*100):.6f} vol%")
            print(f"  VF_HE_NOM  =  {(DCLL_VF_HE_NOM*100):.6f} vol%")
            print(f"  check they add up   = {(vf_checksum2*100):.6f} vol%")
            print(f"  check that BISO + PbLi adds up to the nominal PbLi fraction")
            print(f"  vf_biso_br + vf_pbli_br = {((vf_biso_br+vf_pbli_br))*100:.6f} : VF_PBLI_NOM = {DCLL_VF_LL_NOM*100:.6f}")
            print(f"")
            print(f"Check BISO volume fraction wrt breeding region is correct:")
            print(f"  Inboard : N_in * BISO_VOLUME / (b_in**2 * thickness_in) = {(N_in * BISO_VOLUME / (b_in**2 * thickness_in) * 100):.6f} vol%")
            print(f"  Outboard: N_out * BISO_VOLUME / (b_out**2 * thickness_out) = {(N_out * BISO_VOLUME / (b_out**2 * thickness_out) * 100):.6f} vol%")
            print(f"  Total   : (N1+N2) * BISO_VOLUME / (b1**2 * h1 + b2**2 * h2) = {((N_in+N_out) * BISO_VOLUME / (b_in**2 * thickness_in + b_out**2 * thickness_out) * 100):.6f} vol%")
            print(f"  ...should match                                  vf_biso_br = {(vf_biso_br*100):.6f} vol%")
            print(f"")
            print(f"BISO spheres per cc of breeding volume: {biso_per_cc_bv:.6f} spheres/cm³")
            print(f"             per cc of breeding region: {biso_per_cc_br:.6f} spheres/cm³")
            print(f"")
        
        # Calculate surfaces for trapezoidal plasma chamber
        # NB. MCNP plane: 1.0*x + 0.0*y + plane_C*z - plane_D = 0
        z_in_end    = 77.7
        z_out_start = 400 + z_in_end
        dz = z_out_start - z_in_end
        dx = c_out - c_in
        slope = dx / dz
        
        plane_C = -slope
        plane_D = c_in - (slope * z_in_end)

        self.biso_per_cc_bv = biso_per_cc_bv
        self.biso_per_cc_br = biso_per_cc_br
        self.vf_biso_bv = vf_biso_bv
        self.vf_biso_br = vf_biso_br
        self.vf_pbli_bv = vf_pbli_bv
        self.vf_pbli_br = vf_pbli_br

        self.geom = {'shape_in_per':shape_in_per, 'shape_out_per':shape_out_per, 'N_in':N_in, 'N_out':N_out, 'b_in':b_in, 'b_out':b_out, 'c_in':c_in, 'c_out':c_out, 'thickness_per':thickness_per, 'pC':plane_C, 'pD':plane_D}


    def openmc(self, debug=False, write=True, run=False):

        if write:
            self.prism_helpers(debug=debug)
            self.materials(debug=debug)
            self.geometry(debug=debug)
            self.settings()
            self.tallies()

            self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
            self.model.export_to_model_xml(self.path)

        if run:
            if has_statepoint(self.path):
                print(f"Warning. File {self.path}/statepoint.h5 already exists, so this OpenMC run will be aborted...")
            else:
                self.model.run(cwd=self.path, tracks=True)


if __name__ == '__main__':

    os.makedirs(f"./OpenMC/", exist_ok=True)

    for case in ['A','C',]:
        for isotope in ['U238',]: # 
            for fertile_kgm3 in reversed([0.10, 0.50, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]): #     
                current_run = Wedge(case, fertile_kgm3, isotope=isotope)
                current_run.openmc(debug=True, write=True, run=True)