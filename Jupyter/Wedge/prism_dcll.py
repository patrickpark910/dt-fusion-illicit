import os
import numpy as np
import openmc


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

# Volume fractions in DCLL blanket (Glaser & Goldston 2012, Tb.1)
VF_FS_NOM = 0.019   # F82H (ferritic steel)
VF_LL_NOM = 0.808   # Pb-17Li (lead-lithium)
VF_SI_NOM = 0.076   # SiC (silicon carbide)
VF_HE_NOM = 0.097   # He-4 (gas)


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


def logspace_per_decade(start, stop, pts_per_decade):
    """
    Returns values from 'start' to 'stop' so that each factor-of-10
    interval contains 'pts_per_decade' points (including its first endpoint).
    Might be a little off if 'stop' isn't precisely at a decade, ex. 20e6 eV

    example: 10 points per decade from 1e-5 → 2e7
    grid = logspace_per_decade(1e-5, 20e6, pts_per_decade=10)
    for i in grid:
        print(np.log10(i))
    """
    log_start = np.log10(start)
    log_stop  = np.log10(stop)
    total_decades = log_stop - log_start
    total_steps = total_decades * pts_per_decade
    npts = int(np.ceil(total_steps)) + 1
    return np.logspace(log_start, log_stop, num=npts)


class Wedge():

    def __init__(self, case, fertile_kgm3, isotope='U238', write_openmc=True, run_openmc=False,):

        self.case = case
        self.fertile_kgm3 = fertile_kgm3
        self.fertile_str = f"{fertile_kgm3:06.2f}"
        self.isotope = isotope 

        self.n_particles, self.n_batches = int(4e4), int(25)
        n = f"{self.n_particles:.0e}".replace("+0", "").replace("+", "")
        
        self.name = f"dcll_wedge{self.case}_{self.isotope}_{self.fertile_str}kgm3_{n}x{self.n_batches}"         
        self.path = f"./OpenMC/{self.name}"

        os.makedirs(self.path, exist_ok=True)


    def materials(self, debug=False):

        # Tungsten | First wall
        self.tungsten = openmc.Material(name='firstwall', temperature=TEMP_K)
        self.tungsten.set_density('g/cm3', 19.3)
        self.tungsten.add_element('W', 1)
        self.tungsten.depletable = False

        # Structure | F82H steel (Hirose 2014)
        self.f82h = openmc.Material(name='F82H', temperature=TEMP_K)
        self.f82h.depletable = False
        self.f82h.set_density('g/cm3', 7.78)
        self.f82h.add_element('Fe', 89.3686, percent_type='wo')
        self.f82h.add_element('C' ,  0.1000, percent_type='wo')
        self.f82h.add_element('Si',  0.1000, percent_type='wo')
        self.f82h.add_element('Mn',  0.1300, percent_type='wo')
        self.f82h.add_element('Cr',  8.1600, percent_type='wo')
        self.f82h.add_element('W' ,  1.9400, percent_type='wo')
        self.f82h.add_element('V' ,  0.2000, percent_type='wo')
        self.f82h.add_element('N' ,  0.0014, percent_type='wo')

        # He-4 (gas) | coolant
        he = openmc.Material(name='Helium') 
        he.set_density('g/cm3', 0.004)
        he.add_element('He', 1) 

        # Coolant mix (17 vol% F82H + 83 vol% He-4) | first wall cooling channel
        self.coolant = openmc.Material.mix_materials([self.f82h, he], [0.170, 0.830], 'vo')
        self.coolant.set_density('atom/b-cm', 0.01471892350)  # from Glaser et al. (2025) MCNP
        self.coolant.temperature = TEMP_K
        self.coolant.name = "coolant (17.0 vol% F82H, 83.0 vol% He-4)" 

        # Divider material (51.2 vol% F82H + 48.8 vol% He-4)
        self.divider = openmc.Material.mix_materials([self.f82h, he], [0.512, 0.488], 'vo')
        self.divider.temperature = TEMP_K
        self.divider.name = "divider (51.2 vol% F82H, 48.8 vol% He-4)"

        # Manifold material (45.3 vol% F82H + 54.7 vol% He-4)
        self.manifold = openmc.Material.mix_materials([self.f82h, he], [0.453, 0.547], 'vo')
        self.manifold.temperature = TEMP_K
        self.manifold.name = "inner manifold (45.3 vol% F82H, 54.7 vol% He-4)"

        # SiC (structural insert in blanket)
        sic_struct = openmc.Material(name='SiC_struct', temperature=TEMP_K)
        sic_struct.add_elements_from_formula('SiC')
        sic_struct.set_density('g/cm3', 3.2)


        ''' Breeder material '''
        # Pb-17Li (83 at% Pb, 17 at% Li enriched to 90 at% Li-6)
        pbli = openmc.Material(name='PbLi', temperature=TEMP_K) 
        pbli.set_density('g/cm3', 9.40) 
        pbli.add_element('Pb', 0.83, percent_type='ao') 
        pbli.add_element('Li', 0.17, percent_type='ao', enrichment_target='Li6', enrichment_type='ao', enrichment=90.0)


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


        ''' Adjust materials according to case '''

        if case == 'A':

            biso = openmc.Material.mix_materials([self.kernel, self.sic], [KERNEL_VOLUME/BISO_VOLUME, (1-(KERNEL_VOLUME/BISO_VOLUME))], 'vo') 

            breeder = openmc.Material.mix_materials([pbli, biso], [self.vf_pbli_bv, self.vf_biso_bv], 'vo') 
            self.blanket = openmc.Material.mix_materials([breeder, self.f82h, sic_struct, he], [VF_LL_NOM, VF_FS_NOM, VF_SI_NOM, VF_HE_NOM], 'vo') 
            self.blanket.name = (f"{self.fertile_str} kg/m3"
                                 f" | {self.biso_per_cc_bv:.4f} spheres/cc = {(self.vf_biso_bv*100):.4f} vol% in breeder volume"
                                 f" | {self.biso_per_cc_br:.4f} spheres/cc = {(self.vf_biso_br*100):.4f} vol% in breeding region")
            self.blanket.temperature = TEMP_K

            self.materials = openmc.Materials([self.tungsten, self.f82h, self.coolant, self.blanket, self.divider, self.manifold]) 

            if debug:
                print(f"BISO materials mix")
                print(f"  Vol frac kernel : {(KERNEL_VOLUME/BISO_VOLUME):.6f}")
                print(f"  Vol frac coat   : {(1-(KERNEL_VOLUME/BISO_VOLUME)):.6f}")
                print(f"")
                print(f"Breeding material mix (PbLi + BISO)")
                print(f"  Vol frac PbLi   : {self.vf_pbli_bv:.6f}")
                print(f"  Vol frac BISO   : {self.vf_biso_bv:.6f}")
                print(f"")
                print(f"Blanket mix")
                print(f"  Vol frac breeder: {VF_LL_NOM:.6f}")
                print(f"  Vol frac F82H   : {VF_FS_NOM:.6f}")
                print(f"  Vol frac SiC    : {VF_SI_NOM:.6f}")
                print(f"  Vol frac He-4(g): {VF_HE_NOM:.6f}")
                print(f"")


        elif case in ['B','C']:

            self.blanket = openmc.Material.mix_materials([pbli, self.f82h, sic_struct, he], [VF_LL_NOM, VF_FS_NOM, VF_SI_NOM, VF_HE_NOM], 'vo') 
            self.blanket.temperature = TEMP_K

            self.materials = openmc.Materials([self.tungsten, self.f82h, self.coolant, self.blanket, self.divider, self.manifold, self.kernel, self.sic]) 


    def geometry(self, debug=False):

        # -----------------------------------
        # SURFACES
        # -----------------------------------

        # X-planes
        s111 = openmc.XPlane(surface_id=111, x0= -self.geom['c1'], boundary_type="reflective")
        s112 = openmc.XPlane(surface_id=112, x0= +self.geom['c1'], boundary_type="reflective")
        s113 = openmc.XPlane(surface_id=113, x0= -self.geom['c2'], boundary_type="reflective")
        s114 = openmc.XPlane(surface_id=114, x0= +self.geom['c2'], boundary_type="reflective")

        # Y-planes
        s121 = openmc.YPlane(surface_id=121, y0= -self.geom['c1'], boundary_type="reflective")
        s122 = openmc.YPlane(surface_id=122, y0= +self.geom['c1'], boundary_type="reflective")
        s123 = openmc.YPlane(surface_id=123, y0= -self.geom['c2'], boundary_type="reflective")
        s124 = openmc.YPlane(surface_id=124, y0= +self.geom['c2'], boundary_type="reflective")

        # Z-planes -- DCLL layered structure along z
        # 20 planes (z0..z19) define 19 cells arranged as:
        #
        # INBOARD (back → plasma):
        #   z0  ──[manifold]── z1 ──[BR2]── z2 ──[D1]── z3 ──[BR1]── z4
        #   z4  ──[FW-back]── z5 ──[FW-cool]── z6 ──[FW-front]── z7 ──[W-FW]── z8
        #
        # PLASMA:
        #   z8  ──[plasma chamber]── z9
        #
        # OUTBOARD (plasma → back):
        #   z9  ──[W-FW]── z10 ──[FW-front]── z11 ──[FW-cool]── z12 ──[FW-back]── z13
        #   z13 ──[BR1]── z14 ──[D1]── z15 ──[BR2]── z16 ──[D2]── z17 ──[BR3]── z18
        #   z18 ──[manifold]── z19
        #
        zv = self.geom['zv']  # list of 20 z-plane values

        sz = []  # list of OpenMC ZPlane surfaces
        for i, zval in enumerate(zv):
            btype = "vacuum" if (i == 0 or i == len(zv)-1) else "transmission"
            sz.append(openmc.ZPlane(surface_id=200+i, z0=zval, boundary_type=btype))

        # Plasma chamber boundaries (trapezoidal)
        # Plane normal card entries: Ax + By + Cz - D = 0
        s141 = openmc.Plane(surface_id=141, a= +1.0, b=  0.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s142 = openmc.Plane(surface_id=142, a= -1.0, b=  0.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s143 = openmc.Plane(surface_id=143, a=  0.0, b= +1.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s144 = openmc.Plane(surface_id=144, a=  0.0, b= -1.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")

        # BISO universe
        s150 = openmc.Sphere(surface_id=150, r=0.040) # Kernel
        s151 = openmc.Sphere(surface_id=151, r=0.050) # Outer Coating


        # -----------------------------------
        # CELLS
        # -----------------------------------

        # Shorthand for inboard and outboard XY regions
        xy_in  = +s111 & -s112 & +s121 & -s122
        xy_out = +s113 & -s114 & +s123 & -s124

        # Plasma chamber
        c10 = openmc.Cell(cell_id=10, region= -s142 & -s141 & -s144 & -s143 & +sz[8] & -sz[9])
        c10.importance = {'neutron':1}

        # ---- INBOARD (z0→z8, back to plasma) ----
        c40 = openmc.Cell(cell_id=40, region= xy_in & +sz[0]  & -sz[1],  fill=self.manifold)   #  8.0 cm manifold
        c31 = openmc.Cell(cell_id=31, region= xy_in & +sz[1]  & -sz[2])                        # 21.0 cm BR2 (filled below)
        c20 = openmc.Cell(cell_id=20, region= xy_in & +sz[2]  & -sz[3],  fill=self.divider)    #  3.2 cm D1
        c32 = openmc.Cell(cell_id=32, region= xy_in & +sz[3]  & -sz[4])                        # 22.5 cm BR1 (filled below)
        c21 = openmc.Cell(cell_id=21, region= xy_in & +sz[4]  & -sz[5],  fill=self.f82h)       #  0.4 cm FW back
        c22 = openmc.Cell(cell_id=22, region= xy_in & +sz[5]  & -sz[6],  fill=self.coolant)    #  2.0 cm FW cooling
        c23 = openmc.Cell(cell_id=23, region= xy_in & +sz[6]  & -sz[7],  fill=self.f82h)       #  0.4 cm FW front
        c11 = openmc.Cell(cell_id=11, region= xy_in & +sz[7]  & -sz[8],  fill=self.tungsten)   #  0.2 cm W first wall

        # ---- OUTBOARD (z9→z19, plasma to back) ----
        c12 = openmc.Cell(cell_id=12, region= xy_out & +sz[9]  & -sz[10], fill=self.tungsten)  #  0.2 cm W first wall
        c24 = openmc.Cell(cell_id=24, region= xy_out & +sz[10] & -sz[11], fill=self.f82h)      #  0.4 cm FW front
        c25 = openmc.Cell(cell_id=25, region= xy_out & +sz[11] & -sz[12], fill=self.coolant)   #  2.0 cm FW cooling
        c26 = openmc.Cell(cell_id=26, region= xy_out & +sz[12] & -sz[13], fill=self.f82h)      #  0.4 cm FW back
        c33 = openmc.Cell(cell_id=33, region= xy_out & +sz[13] & -sz[14])                      # 22.5 cm BR1 (filled below)
        c27 = openmc.Cell(cell_id=27, region= xy_out & +sz[14] & -sz[15], fill=self.divider)   #  3.2 cm D1
        c34 = openmc.Cell(cell_id=34, region= xy_out & +sz[15] & -sz[16])                      # 21.0 cm BR2 (filled below)
        c28 = openmc.Cell(cell_id=28, region= xy_out & +sz[16] & -sz[17], fill=self.divider)   #  3.2 cm D2
        c35 = openmc.Cell(cell_id=35, region= xy_out & +sz[17] & -sz[18])                      # 21.0 cm BR3 (filled below)
        c41 = openmc.Cell(cell_id=41, region= xy_out & +sz[18] & -sz[19], fill=self.manifold)  #  8.0 cm manifold

        # Breeding region cells list (for fill logic below)
        br_cells_in  = [c31, c32]           # inboard BR2, BR1
        br_cells_out = [c33, c34, c35]      # outboard BR1, BR2, BR3


        if case in ['A']:

            for c in br_cells_in + br_cells_out:
                c.fill = self.blanket


        elif case in ['B', 'C']:
            # Define the BISO Universe
            c_kern = openmc.Cell(cell_id=51, fill=self.kernel,  region= -s150)
            c_coat = openmc.Cell(cell_id=52, fill=self.sic,     region= +s150 & -s151)
            u50 = openmc.Universe(name='BISO Universe', cells=[c_kern, c_coat])

            # For each breeding region, compute lattice and BISO placement
            all_br_cells = br_cells_in + br_cells_out
            all_br_Nz    = self.geom['Nz_list']   # Nz for each breeding region

            for idx, (cell, nz) in enumerate(zip(all_br_cells, all_br_Nz)):

                ll, ur = cell.region.bounding_box
                shape  = (4, 4, nz)
                pitch  = (ur - ll) / shape

                if case == 'B':
                    n_spheres = int(np.prod(shape))
                    centers = openmc.model.pack_spheres(radius=0.05, region=cell.region, num_spheres=n_spheres)
                    trisos  = [openmc.model.TRISO(0.05, u50, center=c) for c in centers]
                    cell.fill = openmc.model.create_triso_lattice(trisos, ll, pitch, shape, self.blanket)

                elif case == 'C':
                    centers = lattice_coords(ll, shape, pitch)
                    trisos  = [openmc.model.TRISO(0.05, u50, center=c) for c in centers]
                    cell.fill = openmc.model.create_triso_lattice(trisos, ll, pitch, shape, self.blanket)

            if debug:
                print(f"Case {case}: BISO placement in {len(all_br_cells)} breeding regions")
                for idx, (cell, nz) in enumerate(zip(all_br_cells, all_br_Nz)):
                    ll, ur = cell.region.bounding_box
                    shape = (4, 4, nz)
                    pitch = (ur - ll) / shape
                    print(f"  Region {idx} (cell {cell.id}):")
                    print(f"    Shape: {shape}")
                    print(f"    Pitch: {pitch} cm")
                    print(f"    Lower left: {ll}")
                    print(f"    Upper right: {ur}")
                    print(f"    Subtotal BISO particles: {np.prod(shape)}")
                total_biso = sum(4*4*nz for nz in all_br_Nz)
                print(f"  Total BISO particles: {total_biso}")


        self.cells = [c10,
                      c40, c31, c20, c32, c21, c22, c23, c11,                  # inboard
                      c12, c24, c25, c26, c33, c27, c34, c28, c35, c41]        # outboard
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
        source.space = openmc.stats.CartesianIndependent(x=openmc.stats.Uniform(a= -self.geom['c1'], b= +self.geom['c1']),
                                                         y=openmc.stats.Uniform(a= -self.geom['c1'], b= +self.geom['c1']),
                                                         z=openmc.stats.Discrete([self.geom['z_src']], [1.0]) )
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

        self.settings.run_mode  = 'fixed source'
        self.settings.particles = self.n_particles
        self.settings.batches   = self.n_batches


    def prism_helpers(self, debug=True):

        fertile_kgm3 = self.fertile_kgm3

        target_total_biso = 16096
        N1, N2 = 4096, 12000      # inboard, outboard

        # DCLL blanket layer thicknesses [cm]
        t_mf   =  8.0    # manifold
        t_br1  = 22.5    # breeding region 1
        t_br2  = 21.0    # breeding region 2
        t_br3  = 21.0    # breeding region 3 (outboard only)
        t_d    =  3.2    # divider
        t_fwb  =  0.4    # first wall back (F82H)
        t_fwc  =  2.0    # first wall cooling channel
        t_fwf  =  0.4    # first wall front (F82H)
        t_fw   =  0.2    # tungsten first wall
        d0     = 400.0   # plasma chamber

        # Total breeding height (PbLi regions only, excluding dividers)
        h1 = t_br1 + t_br2                    # inboard = 43.5 cm
        h2 = t_br1 + t_br2 + t_br3            # outboard = 64.5 cm

        # Structure thickness between plasma and breeding (FW complex)
        d1 = t_fwb + t_fwc + t_fwf + t_fw   # = 3.0 cm

        f1, f2 = 0.39, 0.61   # fraction of fluence into inboard vs. outboard -- DCLL tokamak has same R0=620, a=207 as HCPB

        ''' 
        Nominally, the volume fracs are: PbLi = 80.80%, F82H = 1.90%, SiC = 7.60%, He = 9.70%
        We want to deduct the BISO volume we add from the PbLi.

        Here I use: "breeder" = PbLi
                    "breeding region" (br) = the physical blanket layer that contains the breeder
                                                this layer usually contains PbLi + F82H + SiC + He
                    "breeding volume" (bv) = amount of volume nominally occupied in blanket by breeder (PbLi)

        That is, if the breeding region has volume 100 m³ of which 80.8% is PbLi and the other 
        19.2% is (F82H + SiC + He), then 80.8 m³ is the breeding volume.
            
        Given fertile_kgm3 [kg U-238 / m³ breeder] we convert this to m³ of BISO.

        **Be careful NOT to renormalize everything. DEDUCT the volume fraction of BISO from 1, and then
        set THAT as the volume fraction of PbLi, keeping F82H/SiC/He constant.

        We use the fact that 80.8 vol% PbLi = 0.808 m³ of PbLi per m³ of blanket to do our work 
        below in terms of volume fractions.  --ppark
        '''

        # Number of BISO spheres per cc of breeding volume
        biso_per_cc_bv = fertile_kgm3_to_biso_per_cc(fertile_kgm3)
        vf_biso_bv     = biso_per_cc_bv * BISO_VOLUME

        if vf_biso_bv > 1.0:
            print(f"Fatal. Your fertile kg/m³ exceeds what can physically fit in the breeding volume!")
            print(f"Fatal. That is, your volume of BISO per cm³ of breeding volume exceeds 1.")

        # New volume ratios of everything w.r.t. PbLi in the breeder volume
        vf_pbli_bv   = 1 - vf_biso_bv
        vf_checksum1 = vf_biso_bv + vf_pbli_bv  # should equal 1

        # New volume ratios of everything w.r.t. everything else in the breeding region
        vf_biso_br = vf_biso_bv * VF_LL_NOM
        vf_pbli_br = vf_pbli_bv * VF_LL_NOM
        # vol frac of F82H, SiC, and He-4 doesn't change
        vf_checksum2 = vf_biso_br + vf_pbli_br + VF_FS_NOM + VF_SI_NOM + VF_HE_NOM  # should equal 1

        # Number of BISO spheres per cc of breeding region
        biso_per_cc_br = vf_biso_br / BISO_VOLUME

        '''
        NB. Maximum possible value of vf_biso_br is VF_LL_NOM = 80.80%. 
        This means that at most 80.80% of the breeding region can be occupied by BISO.
        Per 1 cm³ of breeding region, at most 0.808 cm³ can be BISO, or 1543 spheres/cm³.
        Practically, 400 spheres/cm³ ≈ 1000 kg/m³ so the max possible fertile kg/m³ is ~3858.
        '''

        A_ref = target_total_biso / biso_per_cc_br / (h1+h2) 
        A2 = A_ref * (h1+h2) / ( f1/f2 * h1 + h2)                    # outboard XY face area 
        A1 = A2 * f1/f2                                              # inboard XY face area
        b1, b2 = np.round(np.sqrt(A1), 6), np.round(np.sqrt(A2), 6)  # face lengths
        c1, c2 = b1/2, b2/2                                # face HALF-lengths

        # Volumes of the breeding regions (i.e. area*height, PbLi zones only)
        V1, V2 = A1*h1, A2*h2


        # ---- Compute z-plane positions sequentially ----
        # z0..z19 define 19 cells + 2 vacuum boundaries
        z = 0.0
        zv = []

        # Inboard (back → plasma)
        zv.append(z);  z += t_mf     #  0: inboard manifold bottom (vacuum)
        zv.append(z);  z += t_br2    #  1: top of manifold = bottom of BR2
        zv.append(z);  z += t_d      #  2: top of BR2 = bottom of D1
        zv.append(z);  z += t_br1    #  3: top of D1 = bottom of BR1
        zv.append(z);  z += t_fwb    #  4: top of BR1 = bottom of FW-back
        zv.append(z);  z += t_fwc    #  5: top of FW-back = bottom of FW-cool
        zv.append(z);  z += t_fwf    #  6: top of FW-cool = bottom of FW-front
        zv.append(z);  z += t_fw     #  7: top of FW-front = bottom of W-FW
        zv.append(z);  z += d0       #  8: top of W-FW = plasma bottom

        # Outboard (plasma → back)
        zv.append(z);  z += t_fw     #  9: plasma top = bottom of W-FW
        zv.append(z);  z += t_fwf    # 10: top of W-FW = bottom of FW-front
        zv.append(z);  z += t_fwc    # 11: top of FW-front = bottom of FW-cool
        zv.append(z);  z += t_fwb    # 12: top of FW-cool = bottom of FW-back
        zv.append(z);  z += t_br1    # 13: top of FW-back = bottom of BR1
        zv.append(z);  z += t_d      # 14: top of BR1 = bottom of D1
        zv.append(z);  z += t_br2    # 15: top of D1 = bottom of BR2
        zv.append(z);  z += t_d      # 16: top of BR2 = bottom of D2
        zv.append(z);  z += t_br3    # 17: top of D2 = bottom of BR3
        zv.append(z);  z += t_mf     # 18: top of BR3 = bottom of manifold
        zv.append(z)                  # 19: top of manifold (vacuum)

        # Source position (center of plasma chamber)
        z_src = (zv[8] + zv[9]) / 2.0

        # Calculate surfaces for trapezoidal plasma chamber
        # NB. MCNP plane: 1.0*x + 0.0*y + plane_C*z - plane_D = 0
        z_in_end    = zv[8]    # inboard plasma-side boundary
        z_out_start = zv[9]    # outboard plasma-side boundary
        dz = z_out_start - z_in_end
        dx = c2 - c1
        slope = dx / dz
        
        plane_C = -slope
        plane_D = c1 - (slope * z_in_end)


        # ---- Compute Nz for each breeding region (for cases B/C) ----
        # Distribute BISO proportional to breeding region thickness
        # Order: inboard BR2, inboard BR1, outboard BR1, outboard BR2, outboard BR3

        Nz_in_br2  = int(np.round(N1 / 16 * t_br2 / h1))
        Nz_in_br1  = int(np.round(N1 / 16 * t_br1 / h1))
        Nz_out_br1 = int(np.round(N2 / 16 * t_br1 / h2))
        Nz_out_br2 = int(np.round(N2 / 16 * t_br2 / h2))
        Nz_out_br3 = int(np.round(N2 / 16 * t_br3 / h2))

        Nz_list = [Nz_in_br2, Nz_in_br1, Nz_out_br1, Nz_out_br2, Nz_out_br3]


        if debug:
            print(40*f"=")
            print(f"Case: {fertile_kgm3} kg/m³")
            print(f"")
            print(f"Dimensions (inboard, outboard):")
            print(f" XY face half-lengths c1, c2: {c1:.6f}, {c2:.6f} [cm]")
            print(f" XY face lengths      b1, b2: {b1:.6f}, {b2:.6f} [cm]")
            print(f" XY face areas        A1, A2: {A1:.6f}, {A2:.6f} [cm²]")
            print(f" Breeding heights     h1, h2: {h1:.1f}, {h2:.1f} [cm]")
            print(f" Breeding vols        V1, V2: {V1:.6f}, {V2:.6f} [cm³]")
            print(f" Fraction of fluence  f1, f2: {f1:.6f}, {f2:.6f} ")
            print(f"")
            print(f"'breeder' = Pb-17Li ")
            print(f"'breeding region' (br) = the physical blanket layer that contains the breeder ")
            print(f"                         this layer usually contains PbLi + F82H + SiC + He ")
            print(f"'breeding volume' (bv) = volume nominally occupied in breeding region by PbLi ")
            print(f"")
            print(f"With respect to the BREEDER material, i.e., per 1 m³ of breeder (PbLi) volume, we have these new volume fractions:")
            print(f"  vf_biso_breeder_new =  {(vf_biso_bv*100):.6f} vol%")
            print(f"  vf_pbli_breeder_new =  {(vf_pbli_bv*100):.6f} vol%")
            print(f"  check they add up   = {(vf_checksum1*100):.6f} vol%")
            print(f"")
            print(f"With respect to the whole BREEDING REGION, we have these new volume fractions:")
            print(f"  vf_biso_br =  {(vf_biso_br*100):.6f} vol%")
            print(f"  vf_pbli_br =  {(vf_pbli_br*100):.6f} vol%")
            print(f"  VF_FS_NOM  =  {(VF_FS_NOM*100):.6f} vol%")
            print(f"  VF_SI_NOM  =  {(VF_SI_NOM*100):.6f} vol%")
            print(f"  VF_HE_NOM  =  {(VF_HE_NOM*100):.6f} vol%")
            print(f"  check they add up   = {(vf_checksum2*100):.6f} vol%")
            print(f"  check that BISO + PbLi adds up to the nominal PbLi fraction")
            print(f"  vf_biso_br + vf_pbli_br = {((vf_biso_br+vf_pbli_br))*100:.6f} : VF_LL_NOM = {(VF_LL_NOM)*100:.6f}")
            print(f"")
            in_vol  = b1**2 * h1
            out_vol = b2**2 * h2
            actual_N1 = sum(4*4*nz for nz in Nz_list[:2])
            actual_N2 = sum(4*4*nz for nz in Nz_list[2:])
            print(f"Check BISO volume fraction wrt breeding region is correct:")
            print(f"  Inboard : N1_actual * BISO_VOLUME / (b1² * h1) = {(actual_N1 * BISO_VOLUME / in_vol * 100):.6f} vol%")
            print(f"  Outboard: N2_actual * BISO_VOLUME / (b2² * h2) = {(actual_N2 * BISO_VOLUME / out_vol * 100):.6f} vol%")
            print(f"  Total   : (N1+N2)_actual * BISO_VOLUME / (b1² * h1 + b2² * h2) = {((actual_N1+actual_N2) * BISO_VOLUME / (in_vol + out_vol) * 100):.6f} vol%")
            print(f"  ...should match                                      vf_biso_br = {(vf_biso_br*100):.6f} vol%")
            print(f"")
            print(f"BISO spheres per cc of breeding volume: {biso_per_cc_bv:.6f} spheres/cm³")
            print(f"             per cc of breeding region: {biso_per_cc_br:.6f} spheres/cm³")
            print(f"")
            print(f"Nz per breeding sub-region: {Nz_list}")
            print(f"  Inboard  BR2={Nz_list[0]}, BR1={Nz_list[1]}")
            print(f"  Outboard BR1={Nz_list[2]}, BR2={Nz_list[3]}, BR3={Nz_list[4]}")
            print(f"")


        # Store needed variables in class
        self.biso_per_cc_bv = biso_per_cc_bv
        self.biso_per_cc_br = biso_per_cc_br
        self.vf_biso_bv = vf_biso_bv
        self.vf_biso_br = vf_biso_br
        self.vf_pbli_bv = vf_pbli_bv

        self.geom = {'N1':N1, 'N2':N2, 'b1':b1, 'b2':b2, 'c1':c1, 'c2':c2, 
                      'h1':h1, 'h2':h2, 'pC':plane_C, 'pD':plane_D,
                      'zv':zv, 'z_src':z_src, 'Nz_list':Nz_list}


    def openmc(self, debug=False, write=True, run=False):

        if write:

            self.prism_helpers()
            self.materials()
            self.geometry(debug=debug)
            self.settings()
            self.tallies()

            self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
            self.model.export_to_model_xml(self.path)

        if run:
            if has_statepoint(self.path):
                print(f"Warning. File {self.path}/statepoint.h5 already exists, so this OpenMC run will be aborted...")
            else:
                self.model.run(cwd=self.path, tracks=False)


if __name__ == '__main__':

    os.makedirs(f"./OpenMC/", exist_ok=True)

    for case in ['A','C']:
        for isotope in ['U238',]: # 
            for fertile_kgm3 in [0.10, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]: # [0.10, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]: #     
                current_run = Wedge(case, fertile_kgm3, isotope=isotope)
                current_run.openmc(debug=True, write=True, run=True)