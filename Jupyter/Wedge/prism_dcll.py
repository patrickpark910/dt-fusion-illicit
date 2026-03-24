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



class Prism():

    def __init__(self, case, fertile_kgm3, isotope='U238', breeder_enrich=90.0, write_openmc=True, run_openmc=False,):

        self.case = case
        self.fertile_kgm3 = fertile_kgm3
        self.fertile_str = f"{fertile_kgm3:06.2f}"
        self.breeder_enr_str = f"{breeder_enrich:04.1f}"
        self.isotope = isotope 
        self.breeder_enrich = breeder_enrich
        self.name = f"wedge{self.case}_dcll_Li{self.breeder_enr_str}_{self.isotope}_{self.fertile_str}kgm3"         
        self.path = f"./OpenMC/{self.name}"

        os.makedirs(self.path, exist_ok=True)


    def materials(self, debug=False):

        # Tungsten | First wall
        self.firstwall = openmc.Material(name='firstwall', temperature=TEMP_K)
        self.firstwall.set_density('g/cm3',19.3)
        self.firstwall.add_element('W',1)

        # F82H steel (first wall front, back + part of other structural components)
        self.f82h = openmc.Material(name='f82h', temperature=TEMP_K)
        self.f82h.add_element('Fe', 89.3686, percent_type='wo')
        self.f82h.add_element('C' ,  0.1000, percent_type='wo')
        self.f82h.add_element('Si',  0.1000, percent_type='wo')
        self.f82h.add_element('Mn',  0.1300, percent_type='wo')
        self.f82h.add_element('Cr',  8.1600, percent_type='wo')
        self.f82h.add_element('W' ,  1.9400, percent_type='wo')
        self.f82h.add_element('V' ,  0.2000, percent_type='wo')
        self.f82h.add_element('N' ,  0.0014, percent_type='wo')
        self.f82h.set_density('g/cm3', 7.78)

        # He-4 (gas) - COOLANT
        he = openmc.Material(name='helium') 
        he.set_density('g/cm3', 0.004) # Helium density at 900 K ~80 bar 
        he.add_element('He', 1) 

        # Coolant 
        self.coolant = openmc.Material.mix_materials([self.f82h, he], [0.170, 0.830], 'vo')
        self.coolant.set_density('atom/b-cm', 0.01471892350) # from Glaser et al. (2025) MCNP
        self.coolant.temperature = TEMP_K
        self.coolant.name = "coolant (17.0 vol% F82H, 83.0 vol% He-4)" 

        self.divider = openmc.Material.mix_materials([self.f82h, he], [0.512, 0.488], 'vo')
        self.divider.temperature = TEMP_K
        self.divider.name        = "divider (51.2 vol% F82H, 48.8 vol% He-4)"

        self.manifold = openmc.Material.mix_materials([self.f82h, he], [0.453, 0.547], 'vo')
        self.manifold.temperature = TEMP_K
        self.manifold.name        = "inner manifold (45.3 vol% F82H, 54.7 vol% He-4)"

        self.shield = openmc.Material.mix_materials([self.f82h, he], [0.80, 0.20], 'vo')
        self.shield.temperature = TEMP_K
        self.shield.name = "steel shield (80 vol% F82H, 20 vol% He-4)"

        # Breeder material
        pbli = openmc.Material(name='breeder', temperature=TEMP_K)
        pbli.set_density('g/cm3', DENSITY_DCLL)
        pbli.add_element('Pb', 0.83, percent_type='ao') 
        pbli.add_element('Li', 0.17, percent_type='ao', enrichment_target='Li6', enrichment_type='ao', enrichment=self.breeder_enrich) # Li-6 enrichment to 90 at% 

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

            self.materials = openmc.Materials([self.manifold, self.blanket, self.divider, self.f82h, self.coolant, self.firstwall, self.shield])

            if debug:
                print(f"Case A: Homogeneous BISO")
                print(f"BISO materials mix")
                print(f"  Vol frac kernel : {(KERNEL_VOLUME/BISO_VOLUME):.6f}")
                print(f"  Vol frac coat   : {(1-(KERNEL_VOLUME/BISO_VOLUME)):.6f}")
                print(f"")
                print(f"Breeding material mix")
                print(f"  Vol frac PbLi: {self.vf_pbli_bv:.6f}")
                print(f"  Vol frac BISO   : {self.vf_biso_bv:.6f}")
                print(f"")
                print(f"Blanket mix")
                print(f"  Vol frac breeder: {(self.vf_biso_br + self.vf_pbli_br):.6f}")
                print(f"  Vol frac Ferritic Steel: {DCLL_VF_FS_NOM:.6f}")
                print(f"  Vol frac Helium Coolant: {DCLL_VF_HE_NOM:.6f}")
                print(f"  Vol frac Silicon: {DCLL_VF_SI_NOM:.6f}")
                print(f"")

        elif case in ['B','C']:
            # homogenize everything else in nominal fractions without the BISO
            self.blanket = openmc.Material.mix_materials([pbli, self.f82h, self.sic, he], [DCLL_VF_LL_NOM, DCLL_VF_FS_NOM, DCLL_VF_SI_NOM, DCLL_VF_HE_NOM], 'vo')
            self.blanket.temperature = TEMP_K

            self.materials = openmc.Materials([self.manifold, self.blanket, self.kernel, self.sic, self.divider, self.f82h, self.coolant, self.firstwall, self.shield])


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

        # Z-planes - inboard 
        s130 = openmc.ZPlane(surface_id=130, z0=  -0.00, boundary_type="vacuum")  # 0.0 cm - inboard - vacuum-shield boundary
        s131 = openmc.ZPlane(surface_id=131, z0=   20.00)   #   20.0 cm - inboard - shield-back plate boundary
        s132 = openmc.ZPlane(surface_id=132, z0=  21.5)     #   21.5 cm - inboard - back plate-manifold boundary 
        s133 = openmc.ZPlane(surface_id=133, z0=  29.5)     #   29.5 cm - inboard - manifold-breeding region 1 boundary 
        s134 = openmc.ZPlane(surface_id=134, z0=  50.5)     #   50.5 cm - inboard - breeding region 1-divider boundary 
        s135 = openmc.ZPlane(surface_id=135, z0= 53.7)      #   53.7 cm - inboard - divider-breeding region 2 boundary
        s136 = openmc.ZPlane(surface_id=136, z0= 74.7)      #   74.7 cm - inbodard - breeding region 2 - first wall back boundary 
        s137 = openmc.ZPlane(surface_id=137, z0= 75.1)      #   75.1 cm - inboard - first wall back-coolant boundary 
        s138 = openmc.ZPlane(surface_id=138, z0= 77.1)      #   77.1 cm - inboard - coolant-first wall front boundary 
        s139 = openmc.ZPlane(surface_id=139, z0= 77.5)      #   77.5 cm - inboard - first wall front-tungsten boundary 

        # Z-planes - plasma chamber
        s140 = openmc.ZPlane(surface_id=140, z0= 77.7)      #   77.7 cm - inboard - tungsten - plasma chamber boundary 
        s141 = openmc.ZPlane(surface_id=141, z0= 477.7)     #   477.7 cm - outboard - plasma chamber-tungsten boundary 

        # Z-planes - outboard
        s150 = openmc.ZPlane(surface_id=150, z0= 477.9)     #   477.9 cm - outboard - tungsten-first wall front boundary
        s151 = openmc.ZPlane(surface_id=151, z0= 478.3)     #   478.3 cm - outboard - first wall front-coolant boundary
        s152 = openmc.ZPlane(surface_id=152, z0= 480.3)     #   480.3 cm - outboard - coolant-first wall back boundary
        s153 = openmc.ZPlane(surface_id=153, z0= 480.7)     #   480.7 cm - outboard - first wall back-breeding region 1 boundary
        s154 = openmc.ZPlane(surface_id=154, z0= 501.7)     #   501.7 cm - outboard - breeding region 1-divider 1 boundary
        s155 = openmc.ZPlane(surface_id=155, z0= 504.9)     #   504.9 cm - outboard - divider 1-breeding region 2 boundary
        s156 = openmc.ZPlane(surface_id=156, z0= 525.9)     #   525.9 cm - outboard - breeding region 2-divider 2 boundary
        s157 = openmc.ZPlane(surface_id=157, z0= 529.1)     #   529.1 cm - outboard - divider 2-breeding region 3 boundary
        s158 = openmc.ZPlane(surface_id=158, z0= 550.1)     #   550.1 cm - outboard - breeding region 3-manifold boundary
        s159 = openmc.ZPlane(surface_id=159, z0= 558.1)     #   558.1 cm - outboard - manifold-back plate boundary
        s160 = openmc.ZPlane(surface_id=160, z0= 559.6)     #   559.6 cm - outboard - back plate-shield boundary
        s161 = openmc.ZPlane(surface_id=161, z0= 579.6, boundary_type='vacuum')     #   579.6 cm - outboard - shield-vacuum boundary

        # Plasma chamber boundaries
        # Plane normal card entries: Ax + By + Cz - D = 0
        s171 = openmc.Plane(surface_id=171, a= +1.0, b=  0.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s172 = openmc.Plane(surface_id=172, a= -1.0, b=  0.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s173 = openmc.Plane(surface_id=173, a=  0.0, b= +1.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s174 = openmc.Plane(surface_id=174, a=  0.0, b= -1.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")

        # BISO universe
        s180 = openmc.Sphere(surface_id=180, r=0.040) # Kernel
        s181 = openmc.Sphere(surface_id=181, r=0.050) # Outer Coating

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

            # outboard 
            c34 = openmc.Cell(cell_id=34, region= +s113 & -s114 & +s123 & -s124 & +s153 & -s154, fill=self.blanket)
            c36 = openmc.Cell(cell_id=36, region= +s113 & -s114 & +s123 & -s124 & +s155 & -s156, fill=self.blanket)
            c38 = openmc.Cell(cell_id=38, region= +s113 & -s114 & +s123 & -s124 & +s157 & -s158, fill=self.blanket)
        
        elif case in ['B', 'C']:
            # Common setup for B and C: Define the BISO Universe
            # inboard 
            c23 = openmc.Cell(cell_id=23, region= +s111 & -s112 & +s121 & -s122 & +s133 & -s134)
            c25 = openmc.Cell(cell_id=25, region= +s111 & -s112 & +s121 & -s122 & +s135 & -s136)

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

            if case == 'B': # Random packing
                # -- NB. 'pack_spheres' returns a list of coordinates, not cell instances
                r23 = +s111 & -s112 & +s121 & -s122 & +s133 & -s134
                r25 = +s111 & -s112 & +s121 & -s122 & +s135 & -s136
                r34 = +s113 & -s114 & +s123 & -s124 & +s153 & -s154
                r36 = +s113 & -s114 & +s123 & -s124 & +s155 & -s156
                r38 = +s113 & -s114 & +s123 & -s124 & +s157 & -s158  

                N_in_per = int(self.geom['N_in']/2)
                N_out_per = int(self.geom['N_out']/3)

                centers_23 = openmc.model.pack_spheres(0.05, r23, num_spheres=N_in_per)
                centers_25 = openmc.model.pack_spheres(0.05, r25, num_spheres=N_in_per)
                centers_34 = openmc.model.pack_spheres(0.05, r34, num_spheres=N_out_per)
                centers_36 = openmc.model.pack_spheres(0.05, r36, num_spheres=N_out_per)
                centers_38 = openmc.model.pack_spheres(0.05, r38, num_spheres=N_out_per)
                
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

            elif case == 'C':
                # Rectangular lattice
                centers_23  = lattice_coords(ll_in_23, self.geom['shape_in_per'], pitch_in_23)
                centers_25  = lattice_coords(ll_in_25, self.geom['shape_in_per'], pitch_in_25)
                centers_34  = lattice_coords(ll_out_34, self.geom['shape_out_per'], pitch_out_34)
                centers_36  = lattice_coords(ll_out_36, self.geom['shape_out_per'], pitch_out_36)
                centers_38  = lattice_coords(ll_out_38, self.geom['shape_out_per'], pitch_out_38)
                
                # openmc.model.TRISO(outer_radius, fill, center=(0.0, 0.0, 0.0)) -- each model.TRISO instance IS a cell
                biso_in_23  = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_23]
                biso_in_25  = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_25]
                biso_out_34 = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_34]
                biso_out_36 = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_36]
                biso_out_38 = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_38]

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
        # source.angle = openmc.stats.Isotropic()

        # Create settings and assign source
        self.settings = openmc.Settings()
        self.settings.source = source

        # Run type
        self.settings.run_mode  = 'fixed source'
        self.settings.particles = self.nps
        self.settings.batches   = self.batches

        # note to self 2026-01-22: change to 20 or 50 batches

        # self.settings.trace = (1,1,1)
        self.settings.max_tracks = 4


    def prism_helpers(self, debug=False):

        fertile_kgm3 = self.fertile_kgm3

        self.nps, self.batches = int(4e2), int(25)
        # if fertile_kgm3 < 10.0:
        #     self.nps = int(1e7)
        #     self.batches = int(10)
        # else:
        #     self.nps = int(1e5)
        #     self.batches = int(10)

        target_total_biso = 16096 # maybe change, maybe it's not 2012 particles. 

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

    for iso in ['U238','Th232']:
        for case in ['A','C']:
            for fertile_kgm3 in [0.10, 0.50, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]: # [0.10, 0.50, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99] : #  [0.10, 0.50, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]
                openmc.reset_auto_ids()
                current_run = Prism(case, fertile_kgm3, isotope=iso)
                current_run.openmc(debug=True, write=True, run=True)