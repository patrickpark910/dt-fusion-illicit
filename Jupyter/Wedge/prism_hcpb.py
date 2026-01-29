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


class Prism():

    def __init__(self, case, fertile_kgm3, isotope='U238', write_openmc=True, run_openmc=False,):

        self.case = case
        self.fertile_kgm3 = fertile_kgm3
        self.fertile_str = f"{fertile_kgm3:06.2f}"
        self.isotope = isotope 
        self.name = f"prism_{self.case}_{self.isotope}_{self.fertile_str}kgm3"         
        self.path = f"./OpenMC/{self.name}"

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
            # self.kernel.add_s_alpha_beta('c_U_in_UO2')  # NotImplementedError
            # self.kernel.add_s_alpha_beta('c_O_in_UO2')  # NotImplementedError

        elif self.isotope == 'Th232': 
            self.kernel = openmc.Material(name='ThO2', temperature=TEMP_K) 
            self.kernel.add_elements_from_formula('ThO2') 
            self.kernel.set_density('g/cm3', DENSITY_ThO2)   

        # SiC coating
        self.sic = openmc.Material(name='SiC', temperature=TEMP_K)
        self.sic.add_elements_from_formula('SiC')
        self.sic.set_density('g/cm3', 3.2)
        # self.sic.add_s_alpha_beta('c_Si_in_SiC')  # NotImplementedError


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


        if case in ['B','C']:

            self.blanket = openmc.Material.mix_materials([li4sio4, be, self.eurofer, he], [VF_LI_NOM, VF_BE_NOM, VF_EU_NOM, VF_HE_NOM], 'vo') 
            self.blanket.temperature = TEMP_K

            self.materials = openmc.Materials([self.tungsten, self.eurofer, self.blanket, self.kernel, self.sic]) 


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
        s141 = openmc.Plane(surface_id=141, a= +1.0, b=  0.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s142 = openmc.Plane(surface_id=142, a= -1.0, b=  0.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s143 = openmc.Plane(surface_id=143, a=  0.0, b= +1.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s144 = openmc.Plane(surface_id=144, a=  0.0, b= -1.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")

        # BISO universe
        s150 = openmc.Sphere(surface_id=150, r=0.040) # Kernel
        s151 = openmc.Sphere(surface_id=151, r=0.050) # Outer Coating
        # s152 = openmc.Sphere(surface_id=152, r=0.0501) # Outer universe boundary


        # -----------------------------------
        # CELLS
        # -----------------------------------

        c10 = openmc.Cell(cell_id=10, region= -s142 & -s141 & -s144 & -s143 & +s134 & -s135)  # plasma chamber void
        c10.importance = {'neutron':1}

        c11 = openmc.Cell(cell_id=11, region= +s111 & -s112 & +s121 & -s122 & +s133 & -s134, fill=self.tungsten) # first wall 
        c12 = openmc.Cell(cell_id=12, region= +s111 & -s112 & +s121 & -s122 & +s132 & -s133, fill=self.eurofer)  # after first wall
        c14 = openmc.Cell(cell_id=14, region= +s111 & -s112 & +s121 & -s122 & +s130 & -s131, fill=self.eurofer)  # manifold

        c21 = openmc.Cell(cell_id=21, region= +s113 & -s114 & +s123 & -s124 & +s135 & -s136, fill=self.tungsten)
        c22 = openmc.Cell(cell_id=22, region= +s113 & -s114 & +s123 & -s124 & +s136 & -s137, fill=self.eurofer)
        c24 = openmc.Cell(cell_id=24, region= +s113 & -s114 & +s123 & -s124 & +s138 & -s139, fill=self.eurofer)

        if case == 'A':

            c13 = openmc.Cell(cell_id=13, region= +s111 & -s112 & +s121 & -s122 & +s131 & -s132, fill=self.blanket)
            c23 = openmc.Cell(cell_id=23, region= +s113 & -s114 & +s123 & -s124 & +s137 & -s138, fill=self.blanket)
            

        elif case in ['B', 'C']:
            # Common setup for B and C: Define the BISO Universe
            c13 = openmc.Cell(cell_id=13, region= +s111 & -s112 & +s121 & -s122 & +s131 & -s132)
            c23 = openmc.Cell(cell_id=23, region= +s113 & -s114 & +s123 & -s124 & +s137 & -s138)

            c31 = openmc.Cell(cell_id=31, fill=self.kernel,  region= -s150)
            c32 = openmc.Cell(cell_id=32, fill=self.sic,     region= +s150 & -s151) # 
            # c33 = openmc.Cell(cell_id=33, fill=self.blanket, region= +s151)
            u50 = openmc.Universe(name='BISO Universe', cells=[c31, c32]) # c33

            # Lower left (ll), upper right (ur) of inboard, outboard bounding boxes
            # Example values:
            #   Lower left : [-12.01134016 -12.01134016   0.01    ]
            #   Upper right: [ 12.01134016  12.01134016  45.45    ]
            
            ll_in, ur_in   = c13.region.bounding_box
            ll_out, ur_out = c23.region.bounding_box

            # ll_in, ur_in   = ll_in*np.array([1.01, 1.01, 0.99]), ur_in*np.array([1.01, 1.01, 1.01])
            # ll_out, ur_out = ll_out*np.array([1.01, 1.01, 0.99]), ur_out*np.array([1.01, 1.01, 1.01])

            shape_in  = (2, 2, 128)  # (Nx, Ny, Nz)
            shape_out = (2, 2, 375)

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
                biso_in  = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_in]
                biso_out = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_out]

                c13.fill = openmc.model.create_triso_lattice(biso_in, ll_in, pitch_in, shape_in, self.blanket)
                c23.fill = openmc.model.create_triso_lattice(biso_out, ll_out, pitch_out, shape_out, self.blanket)


                # # Inboard lattice
                # lattice_in = openmc.RectLattice()
                # lattice_in.lower_left = ll_in
                # lattice_in.pitch = pitch_in
                # lattice_in.universes = np.full(shape_in[::-1], u50) # NB. lattices indices are [z, y, x]
                # lattice_in.outer = openmc.Universe(cells=[openmc.Cell(fill=self.blanket)])

                # # Outboard lattice
                # lattice_out = openmc.RectLattice()
                # lattice_out.lower_left = ll_out
                # lattice_out.pitch = pitch_out
                # lattice_out.universes = np.full(shape_out[::-1], u50)
                # lattice_out.outer = openmc.Universe(cells=[openmc.Cell(fill=self.blanket)])
                
                # c13.fill = lattice_in
                # c23.fill = lattice_out

                if debug:
                    print(f"Case C: Fixed lattice BISO placement")
                    print(f"  Inboard lattice:")
                    print(f"    Shape: {shape_in}")
                    print(f"    Pitch: {pitch_in} cm")
                    print(f"    Lower left: {ll_in}")
                    print(f"    Upper right: {ur_in}")
                    print(f"    Subtotal BISO particles: {np.prod(shape_in)}")
                    print(f"  Outboard lattice:")
                    print(f"    Shape: {shape_out}")
                    print(f"    Pitch: {pitch_out} cm")
                    print(f"    Lower left: {ll_out}")
                    print(f"    Upper right: {ur_out}")
                    print(f"    Subtotal BISO particles: {np.prod(shape_out)}")
                    print(f"  Total BISO particles: {np.prod(shape_in) + np.prod(shape_out)}")


        self.cells = [c10, c11, c12, c13, c14, c21, c22, c23, c24]
        self.geometry = openmc.Geometry(openmc.Universe(cells=self.cells))


    def tallies(self, debug=False):

        self.tallies = openmc.Tallies() # initialize
        cell_filter = openmc.CellFilter([cell for cell in self.cells if cell.fill is not None])
        energy_filter = openmc.EnergyFilter(logspace_per_decade(1e-5, 20e6, 100)) # './helpers/utilities.py'

        # Flux tally 
        flux_tally        = self.make_tally('flux', ['flux'], filters=[cell_filter])
        flux_energy_tally = self.make_tally('flux spectrum', ['flux'], filters=[cell_filter, energy_filter])

        # Fertile element reaction rates
        fertile_tally_tot    = self.make_tally(f'Total fertile rxn rate',     ['(n,gamma)', 'fission', 'elastic'], nuclides=['U238', 'U235', 'Th232'])
        fertile_tally        = self.make_tally(f'Fertile rxn rates',          ['(n,gamma)', 'fission', 'elastic'], filters=[cell_filter], nuclides=['U238', 'U235', 'Th232'])
        fertile_energy_tally = self.make_tally(f'Fertile rxn rates spectrum', ['(n,gamma)', 'fission', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['U238', 'U235', 'Th232'])

        # Lithium reaction rates
        Li_tally_tot    = self.make_tally('Total Li rxn rate',     ['(n,gamma)', '(n,Xt)', 'elastic'], nuclides=['Li6', 'Li7'])
        Li_tally        = self.make_tally('Li rxn rates',          ['(n,gamma)', '(n,Xt)', 'elastic'], filters=[cell_filter], nuclides=['Li6', 'Li7'])
        Li_energy_tally = self.make_tally('Li rxn rates spectrum', ['(n,gamma)', '(n,Xt)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['Li6', 'Li7'])

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
                                                         z=openmc.stats.Discrete([247.7], [1.0]) )
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

        if fertile_kgm3 < 10.0:
            self.nps = int(1e7)
            self.batches = int(10)
        else:
            self.nps = int(1e5)
            self.batches = int(10)

        target_total_biso = 2012 
        N1, N2 = 512, 1500

        # HCPB tokamak dimensions (along z in model)
        h1, h2 = 45.0, 82.0     # inboard, outboard thicknesses
        d1 = 2.5 + 0.2          # Eurofer + first wall 
        d0 = 400.0              # plasma chamber
        oo = h1 + d0 + 2 * d1   # outboard offset (where outboard starts in z)
        f1, f2 = 0.39, 0.61     # fraction of fluence into inboard vs. outboard in full OpenMC HCPB model


        ''' 
        Nominally, the volume fracs are: Li4SiO4 = 13.04%, Be = 37.90%, Eurofer = 11.76%, He = 37.30%
        We want to deduct the BISO volume we add from the Li4SiO4.

        Here I use: "breeder" = Li4SiO4, Be
                    "breeding region" (br) = the physical blanket layer that contains the breeder
                                             this layer usually contains breeder + structure + coolant
                    "breeding volume" (bv) = amount of volume nominally occupied in blanket by breeder

        That is, if the breeding region has volume 100 m³ of which 60% is (Li4SiO4 + Be) and the other 
        40% is (Eurofer + Be), then 60 m³ is the breeding volume.
            
        Given fertile_kgm3 [kg U-238 / m³ breeder] we convert this to m³ of BISO.

        **Be careful NOT to renormalize everything. DEDUCT the volume fraction of BISO from 1, and then
        set THAT as the volume fraction of breeder, which then subdivides into Li4SiO4 and Be in a 
        13.04:37.90 ratio.

        We use the fact that 13.04 vol% Li4SiO4 = 0.1304 m³ of Li4SiO4 per m³ of blanket to do our work 
        below in terms of volume fractions.  --ppark
        '''

        # Nominal fraction the breeding volume is in the breeding region
        VF_BV_BR_NOM = VF_LI_NOM + VF_BE_NOM

        # Nominal volume fractions of Li4SiO4 vs. Be in the breeding volume (bv)
        vf_li_bv_nom = VF_LI_NOM / (VF_LI_NOM + VF_BE_NOM)
        vf_be_bv_nom = VF_BE_NOM / (VF_LI_NOM + VF_BE_NOM)

        # Number of BISO spheres per cc of breeding volume
        biso_per_cc_bv = fertile_kgm3_to_biso_per_cc(fertile_kgm3)
        vf_biso_bv     = biso_per_cc_bv * BISO_VOLUME

        if vf_biso_bv > 1.0:
            print(f"Fatal. Your fertile kg/m³ exceeds what can physically fit in the breeding volume!")
            print(f"Fatal. That is, your volume of BISO per cm³ of breeding volume exceeds 1.")

        # New volume ratios of everything w.r.t. breeders in the breeder volume
        vf_libe_bv   = 1 - vf_biso_bv
        vf_li_bv     = vf_li_bv_nom * vf_libe_bv
        vf_be_bv     = vf_be_bv_nom * vf_libe_bv
        vf_checksum1 = vf_biso_bv + vf_li_bv + vf_be_bv  # should equal 1

        # New volume ratios of everyting w.r.t. everything else in the breeding region
        vf_biso_br = vf_biso_bv * (VF_LI_NOM + VF_BE_NOM)
        vf_li_br   = vf_li_bv * (VF_LI_NOM + VF_BE_NOM)
        vf_be_br   = vf_be_bv * (VF_LI_NOM + VF_BE_NOM)
        # vol frac of Eurofer and He-4 doesn't change
        vf_checksum2 = vf_biso_br + vf_li_br + vf_be_br + VF_EU_NOM + VF_HE_NOM  # should equal 1

        # Number of BISO spheres per cc of breeding region
        biso_per_cc_br = vf_biso_br / BISO_VOLUME

        '''
        NB. Maximum possible value of vf_biso_br is (VF_LI_NOM + VF_BE_NOM) = 50.94%. 
        This means that at most 50.94% of the breeding region can be occupied by BISO
        Per 1 cm³ of breeding region, at most 0.5094 cm³ can be BISO, or 973 spheres/cm³.
        Practically, 400 spheres/cm³ = 1000 kg/m³ so the max possible fertile kg/m³ is 2432.50
        '''

        A_ref = target_total_biso / biso_per_cc_br / (h1+h2) 
        A2 = A_ref * (h1+h2) / ( f1/f2 * h1 + h2)                    # outboard XY face area 
        A1 = A2 * f1/f2                                              # inboard XY face area
        b1, b2 = np.round(np.sqrt(A1), 6), np.round(np.sqrt(A2), 6)  # face lengths
        c1, c2 = b1/2, b2/2                                # face HALF-lengths

        # Volumes of the breeding regions (i.e. area*height)
        V1, V2 = A1*h1, A2*h2


        if debug:
            print(40*f"=")
            print(f"Case: {fertile_kgm3} kg/cm³")
            print(f"")
            print(f"Dimensions (inboard, outboard):")
            print(f" XY face half-lengths c1, c2: {c1:.6f}, {c2:.6f} [cm]")
            print(f" XY face lengths      b1, b2: {b1:.6f}, {b2:.6f} [cm]")
            print(f" XY face areas        A1, A2: {A1:.6f}, {A2:.6f} [cm²]")
            print(f" Breeding region vols V1, V2: {V1:.6f}, {V2:.6f} [cm³]")
            print(f" Fraction of fluence  f1, f2: {f1:.6f}, {f2:.6f} ")
            print(f"")
            print(f"'breeder' = Li4SiO4, Be ")
            print(f"'breeding region' (br) = the physical blanket layer that contains the breeder ")
            print(f"                         this layer usually contains breeder + structure + coolant ")
            print(f"'breeding volume' (bv) = volume nominally occupied in breeding region by breeder ")
            print(f"")
            print(f"With respect to the BREEDER material, i.e., per 1 m³ of breeder volume, we have these new volume fractions:")
            print(f"  vf_biso_breeder_new =  {(vf_biso_bv*100):.6f} vol%")
            print(f"  vf_li_breeder_new   =  {(vf_li_bv*100):.6f} vol%")
            print(f"  vf_be_breeder_new   =  {(vf_be_bv*100):.6f} vol%")
            print(f"  check they add up   = {(vf_checksum1*100):.6f} vol%")
            print(f"")
            print(f"With respect to the whole BREEDING REGION, we have these new volume fractions:")
            print(f"  vf_biso_br =  {(vf_biso_br*100):.6f} vol%")
            print(f"  vf_li_br   =  {(vf_li_br*100):.6f} vol%")
            print(f"  vf_be_br   =  {(vf_be_br*100):.6f} vol%")
            print(f"  VF_EU_NOM  =  {(VF_EU_NOM*100):.6f} vol%")
            print(f"  VF_HE_NOM  =  {(VF_HE_NOM*100):.6f} vol%")
            print(f"  check they add up   = {(vf_checksum2*100):.6f} vol%")
            print(f"  check that BISO + Li4SiO4 + Be adds up to the nominal Li4SiO4 + Be fraction")
            print(f"  vf_biso_br + vf_li_br + vf_be_br = {((vf_biso_br+vf_li_br+vf_be_br))*100:.6f} : VF_LI_NOM + VF_BE_NOM = {((VF_LI_NOM+VF_BE_NOM))*100:.6f}")
            print(f"")
            print(f"Check BISO volume fraction wrt breeding region is correct:")
            print(f"  Inboard : N1 * BISO_VOLUME / (b1**2 * h1) = {(N1 * BISO_VOLUME / (b1**2 * h1) * 100):.6f} vol%")
            print(f"  Outboard: N2 * BISO_VOLUME / (b2**2 * h2) = {(N2 * BISO_VOLUME / (b2**2 * h2) * 100):.6f} vol%")
            print(f"  Total   : (N1+N2) * BISO_VOLUME / (b1**2 * h1 + b2**2 * h2) = {((N1+N2) * BISO_VOLUME / (b1**2 * h1 + b2**2 * h2) * 100):.6f} vol%")
            print(f"  ...should match                                  vf_biso_br = {(vf_biso_br*100):.6f} vol%")
            print(f"")
            print(f"BISO spheres per cc of breeding volume: {biso_per_cc_bv:.6f} spheres/cm³")
            print(f"             per cc of breeding region: {biso_per_cc_br:.6f} spheres/cm³")
            print(f"")


        # Calculate surfaces for trapezoidal plasma chamber
        # NB. MCNP plane: 1.0*x + 0.0*y + plane_C*z - plane_D = 0
        z_in_end    = 47.7
        z_out_start = 447.7
        dz = z_out_start - z_in_end
        dx = c2 - c1
        slope = dx / dz
        
        plane_C = -slope
        plane_D = c1 - (slope * z_in_end)

        # I moved the needed variables to these class dictionaries bc the amount of 'self.' 
        # in the arithmetic above was really screwing readability --ppark 2026-01-02

        self.biso_per_cc_bv = biso_per_cc_bv
        self.biso_per_cc_br = biso_per_cc_br
        self.vf_biso_bv = vf_biso_bv
        self.vf_biso_br = vf_biso_br
        self.vf_li_bv = vf_li_bv
        self.vf_be_bv = vf_be_bv

        self.geom = {'N1':N1, 'N2':N2, 'b1':b1, 'b2':b2, 'c1':c1, 'c2':c2, 'h1':h1, 'h2':h2, 'pC':plane_C, 'pD':plane_D}


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
                self.model.run(cwd=self.path, tracks=True)


if __name__ == '__main__':

    os.makedirs(f"./OpenMC/", exist_ok=True)

    for case in ['C',]:
        for fertile_kgm3 in [0.10, 0.50, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99] : #  [0.10, 0.50, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]   # [0.10,60,150]
            current_run = Prism(case, fertile_kgm3, isotope='Th232')
            current_run.openmc(debug=True, write=True, run=True)