import openmc
import os, sys
import numpy as np
from datetime import date


# Import helper functions
from Python.parameters import *
from Python.utilities import *

class Reactor:

    @timer
    def __init__(self, breeder='flibe', fertile_element='U', fertile_density_kgm3=0.0, run_openmc=False):

        self.fertile_density_kgm3 = fertile_density_kgm3
        self.fertile_element = fertile_element.capitalize()
        self.run = run_openmc # i'm calling this self.run bc self.run_openmc() is the function, which main.py will confuse!

        if breeder.lower() == 'flibe':
            self.temp_k          = TEMP_FLIBE_K
            self.breeder_name         = 'FLiBe'
            self.breeder_name_density = DENSITY_FLIBE
            self.breeder_name_enrich  = ENRICH_FLIBE
            # self.breeder_name_volume  = VOLUME_FLIBE
            # print(self.blanket_volume)

        # elif blanket.lower() == 'll':
        #     self.blanket         = 'LL'
        #     self.blanket_density = DENSITY_LL
        #     self.blanket_enrich  = ENRICH_LL

        # elif blanket.lower() == 'pb':
        #     self.blanket         = 'PB'
        #     self.blanket_density = DENSITY_PB
        #     self.blanket_enrich  = ENRICH_PB


        """
        Name file based on reactor config
        """
        today = date.today().strftime("%Y-%m-%d")
        self.name = f"{self.breeder_name}_{self.fertile_element}{self.fertile_density_kgm3:.2f}kgm3_Li{self.breeder_name_enrich:04.1f}_{self.temp_k}K" # _{today}
        # should come out to smth like: FLiBe_U010kgm3_Li7.5_900K
        self.path = f"./OpenMC/{self.name}/"
        
        # os.makedirs(self.path, exist_ok=True)

        start_msg = f"="*42+f"\nWriting an OpenMC input for {self.name}"
        print(f"{Colors.RED}{start_msg}{Colors.END}")

        self.materials()
        self.geometry()
        self.settings()
        self.tallies()


    @timer
    def materials(self):
        """
        OpenMC Materials
        """

        self.firstwall = openmc.Material(name='firstwall', temperature=self.temp_k)
        self.firstwall.depletable = False

        if self.breeder_name == 'FLiBe':
            self.firstwall.add_element('O',5/1e6,percent_type='wo')
            self.firstwall.add_element('N',5/1e6,percent_type='wo')
            self.firstwall.add_element('C',5/1e6,percent_type='wo')
            self.firstwall.add_element('Na',4/1e6,percent_type='wo')
            self.firstwall.add_element('K',2.5/1e6,percent_type='wo')
            self.firstwall.add_element('Al',3/1e6,percent_type='wo')
            self.firstwall.add_element('Ca',0.5/1e6,percent_type='wo')
            self.firstwall.add_element('Cr',0.5/1e6,percent_type='wo')
            self.firstwall.add_element('Cu',0.5/1e6,percent_type='wo')
            self.firstwall.add_element('Fe',5/1e6,percent_type='wo')
            self.firstwall.add_element('W',1-(5+5+5+4+2.5+3+0.5+0.5+0.5+5)/1e6,percent_type='wo')
            self.firstwall.set_density('g/cm3',19.3)
        



        self.structure = openmc.Material(name='structure', temperature=self.temp_k)
        self.structure.depletable = False

        if self.breeder_name == 'FLiBe':
            """ V-4Cr-4Ti """
            # This code is from jlball but V-4Cr-4Ti specs also can be found
            # from ANL material id BL-47, p.3 (1994) altho these are slightly different
            # cf. osti.gov/servlets/purl/10194461 --ppark 2025-07-22
            self.structure.add_element('Cr',0.04,percent_type='wo')
            self.structure.add_element('Ti',0.04,percent_type='wo')
            self.structure.add_element('C',56/1e6,percent_type='wo')
            self.structure.add_element('O',181/1e6,percent_type='wo')
            self.structure.add_element('N',103/1e6,percent_type='wo')
            self.structure.add_element('B',7/1e6,percent_type='wo')
            self.structure.add_element('Na',17/1e6,percent_type='wo')
            self.structure.add_element('Mg',0.5/1e6,percent_type='wo')
            self.structure.add_element('Al',119/1e6,percent_type='wo')
            self.structure.add_element('Si',280/1e6,percent_type='wo')
            self.structure.add_element('Mn',0.5/1e6,percent_type='wo')
            self.structure.add_element('Fe',80/1e6,percent_type='wo')
            self.structure.add_element('Ni',13/1e6,percent_type='wo')
            self.structure.add_element('Cu',4/1e6,percent_type='wo')
            self.structure.add_element('V',1-0.04-0.04-(56+181+103+7+17+0.5+119+280+0.5+80+13+4)/1e6,percent_type='wo')
            self.structure.set_density('g/cm3',6.05) 
            # This density value is sus and needs a good source --jlball 
            # This value is from Metals Handbook, 9th ed, vol 2: "Properties and Selection: Nonferrous Alloys and Pure Metals" (1979) --ppark 2025-07-22



        # Breeder
        self.breeder = openmc.Material(name='breeder', temperature=self.temp_k)
        self.breeder.set_density('g/cm3', self.breeder_name_density)

        if self.breeder_name == 'FLiBe':
            self.breeder.add_elements_from_formula('F4Li2Be', 'ao', enrichment_target='Li6', enrichment_type='wo', enrichment=self.breeder_name_enrich)

        # elif self.breeder_name.lower() == 'll':
        #     breeder.add_element('Pb', 0.83)
        #     breeder.add_element('Li', 0.17, enrichment=self.breeder_name_enrich, enrichment_target='Li6', enrichment_type='wo')  # Li-6 enrichment to 90%

        # elif self.blanket.lower() == 'pb':
        #     self.blanket_density = DENSITY_PB
        #     self.blanket_enrich  = ENRICH_PB
        

        # Fertile compound
        self.fertile = openmc.Material(name='fertile', temperature=self.temp_k)

        if self.breeder_name == 'FLiBe':
            if self.fertile_element == 'U':
                self.fertile.add_elements_from_formula('UF4','ao') # 'ao' here refers to 1:4 atomic ratio of U:F in UF4
                self.fertile.set_density('g/cm3', DENSITY_UF4) 
            elif self.fertile_element == 'Th':
                self.fertile.add_elements_from_formula('ThF4','ao') # 'ao' here refers to 1:4 atomic ratio of U:F in UF4
                self.fertile.set_density('g/cm3', DENSITY_ThF4) 

        # elif self.breeder_name in ['LL', 'PB']:
        #     if self.fertile_element == 'U':
        #         pass  # UO2 pebbles
        #     elif self.fertile_element == 'Th':
        #         pass  # ThO2 pebbles


        # Mix fertile compound to blanket
        if self.breeder_name == 'FLiBe':
            breeder_mass_frac, fertile_compound_mass_frac = calc_blanket_mass_fracs(self.fertile_density_kgm3, 
                                                                                    fertile_element=self.fertile_element, 
                                                                                    fertile_enrich=0.71, 
                                                                                    breeder_density_kgm3=1.94e3)
            self.blanket = openmc.Material.mix_materials([self.breeder, self.fertile], [breeder_mass_frac, fertile_compound_mass_frac], 'wo') # fractions in 'mix_materials' MUST add up to 1
            self.blanket.name, self.blanket.temperature = self.name, self.temp_k
            

        # elif self.breeder_name in ['LL', 'PB']:
        #     pass

        self.materials = openmc.Materials([self.firstwall, self.structure, self.blanket])
        # self.materials.export_to_xml(self.path)


    @timer
    def geometry(self):

        # Tokamak parameters
        if self.breeder_name == 'FLiBe':
            R0, a, kappa, delta = TOKAMAK_R0, TOKAMAK_A, TOKAMAK_KAPPA, TOKAMAK_DELTA 
            d_fw  = FLIBE_FW_CM 
            d_st1 = d_fw  + FLIBE_ST1_CM
            d_br1 = d_st1 + FLIBE_BR1_CM
            d_st2 = d_br1 + FLIBE_ST2_CM
            d_br2 = d_st2 + FLIBE_BR2_CM
            d_st3 = d_br2 + FLIBE_ST3_CM


        """ Numerically create points for poloidal cross-section contours """
        # Parameter t array
        t = np.linspace(0, 2*np.pi, 100)

        # Arrays of a (minor r) and z points
        points_vc  =  miller_model(t, R0, a, kappa, delta)         # coords around vacuum chamber
        points_fw  = miller_offset(t, R0, a, kappa, delta, d_fw)   # outer coords around first wall
        points_st1 = miller_offset(t, R0, a, kappa, delta, d_st1)  # outer coords around structural region 1
        points_br1 = miller_offset(t, R0, a, kappa, delta, d_br1)  # outer coords around breeding region 1
        points_st2 = miller_offset(t, R0, a, kappa, delta, d_st2)  # outer coords around structural region 2
        points_br2 = miller_offset(t, R0, a, kappa, delta, d_br2)  # outer coords around breeding region 2
        points_st3 = miller_offset(t, R0, a, kappa, delta, d_st3)  # outer coords around structural region 3

        # Create OpenMC surfaces
        surfaces = [openmc.model.Polygon(points_vc , basis='rz'), # ZTorus is created by revolving an RZ contour around the Z-axis
                    openmc.model.Polygon(points_fw , basis='rz'),
                    openmc.model.Polygon(points_st1, basis='rz'),
                    openmc.model.Polygon(points_br1, basis='rz'),
                    openmc.model.Polygon(points_st2, basis='rz'),
                    openmc.model.Polygon(points_br2, basis='rz'),
                    openmc.model.Polygon(points_st3, basis='rz') ]

        # Add boundary surfaces
        r_max = max([point[0] for point in points_st3]) + 1.0  # 1 cm beyond outermost surface
        z_max = max([point[1] for point in points_st3]) + 1.0
        z_min = min([point[1] for point in points_st3]) - 1.0

        outer_cylinder = openmc.ZCylinder(r=r_max, boundary_type='vacuum')
        top_plane      = openmc.ZPlane(z0=z_max, boundary_type='vacuum')
        bottom_plane   = openmc.ZPlane(z0=z_min, boundary_type='vacuum')

        """ Create and assign cells """
        cell_vc   = openmc.Cell(region= -surfaces[0])
        cell_vc.importance = {'neutron':1}
        cell_fw   = openmc.Cell(region= +surfaces[0] & -surfaces[1] , fill=self.firstwall) 
        cell_st1  = openmc.Cell(region= +surfaces[1] & -surfaces[2] , fill=self.structure)
        cell_br1  = openmc.Cell(region= +surfaces[2] & -surfaces[3] , fill=self.blanket)
        cell_st2  = openmc.Cell(region= +surfaces[3] & -surfaces[4] , fill=self.structure)
        cell_br2  = openmc.Cell(region= +surfaces[4] & -surfaces[5] , fill=self.blanket)
        cell_st3  = openmc.Cell(region= +surfaces[5] & -surfaces[6] , fill=self.structure) 

        # Void cell with proper boundaries (otherwise causes error with just Polygons)
        cell_void = openmc.Cell(region= +surfaces[6] & -outer_cylinder & +bottom_plane & -top_plane)
        cell_void.importance = {'neutron':0}

        self.cells = [cell_vc, cell_fw, cell_st1, cell_br1, cell_st2, cell_br2, cell_st3, cell_void]
        self.geometry = openmc.Geometry(openmc.Universe(cells=self.cells))
        # self.geometry.export_to_xml(self.path)


    @timer
    def plot(self):
        # ---------- Plots.xml + generate images ----------
        # RZ poloidal cross-section (best for axisymmetric Miller shapes)
        xy = openmc.Plot()
        xy.filename = "tokamak_xy"
        xy.basis = "xy"
        xy.width  = (1800, 1800)
        xy.pixels = (18000, 18000)
        xy.color_by = "material"
        xy.colors = {
            self.firstwall: (255, 0, 0),    # Red for first wall
            self.structure: (0, 255, 0),    # Green for structure
            self.blanket: (0, 0, 255),    # Blue for breeder
            # Void regions will be white by default
        }

        # XZ Cartesian slice (also useful)
        xz = openmc.Plot()
        xz.filename = "tokamak_xz"
        xz.basis = "xz"
        xz.width  = (1800, 1000)
        xz.pixels = (18000, 10000)
        xz.color_by = "material"
        xz.colors = {
            self.firstwall: (255, 0, 0),    # Red for first wall
            self.structure: (0, 255, 0),    # Green for structure
            self.blanket: (0, 0, 255),    # Blue for breeder
            # Void regions will be white by default
        }

        openmc.Plots([xy, xz]).export_to_xml(self.path)
        openmc.plot_geometry(path_input=self.path)  # writes tokamak_rz.ppm and tokamak_xz.ppm
        print("Done. Files: tokamak_rz.ppm, tokamak_xz.ppm")


    @timer
    def tallies(self):
        """ TALLIES """
        self.tallies = openmc.Tallies() # initialize

        # Filters
        cell_filter = openmc.CellFilter(self.cells)

        E_bin_edges = logspace_per_decade(1e-5, 20e6, 100) # './helpers/utilities.py'
        energy_filter = openmc.EnergyFilter(E_bin_edges)
        filters = [energy_filter, cell_filter]

        # Flux tally 
        flux_tally = openmc.Tally(name='flux')
        flux_tally.scores = ['flux'] # specific names required
        flux_tally.filters = filters

        # Uranium reaction rates
        U_tally = openmc.Tally(name='uranium rxn rates')
        U_tally.scores = ['(n,gamma)','fission', 'elastic'] # specific names required
        U_tally.nuclides = ['U238', 'U235']
        U_tally.filters = filters

        # Lithium reaction rates
        Li_tally = openmc.Tally(name='lithium rxn rates')
        Li_tally.scores = ['(n,gamma)','(n,Xt)', 'elastic'] # specific names required
        Li_tally.nuclides = ['Li6', 'Li7']
        Li_tally.filters = filters

        # Fluorine reaction rates
        F_tally = openmc.Tally(name='fluorine rxn rates')
        F_tally.scores = ['(n,gamma)', 'elastic'] # specific names required
        F_tally.nuclides = ['F19']
        F_tally.filters = filters

        # Beryllium reaction rates
        Be_tally = openmc.Tally(name='beryllium rxn rates')
        Be_tally.scores = ['(n,gamma)','(n,2n)', 'elastic'] # specific names required
        Be_tally.nuclides = ['Be9']
        Be_tally.filters = filters

        self.tallies.extend([flux_tally, U_tally, Li_tally, F_tally, Be_tally])


    @timer
    def settings(self):
        """ SETTINGS """
       
        source = openmc.IndependentSource()
        source.space = openmc.stats.CylindricalIndependent(
                                      r=openmc.stats.Discrete([600.0], [1.0]),  # r   = 600 [cm]
                                    phi=openmc.stats.Uniform(0.0, 2*np.pi)   ,  # phi = 0 to 2pi
                                      z=openmc.stats.Discrete([0.0], [1.0])   ) # z   = 0
                                                          
        source.particle = 'neutron'
        source.energy   = openmc.stats.Discrete([14.0e6], [1.0])  # 14 MeV
        source.angle    = openmc.stats.Isotropic()

        # Create settings and assign source
        self.settings = openmc.Settings()
        self.settings.source = source

        """ Run type """
        self.settings.run_mode = 'fixed source'
        self.settings.particles = int(1e6)
        self.settings.batches = 10

    @timer
    def run_openmc(self):
        self.materials.cross_sections = set_xs_path()
        self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
        self.model.export_to_model_xml(f"./OpenMC/{self.name}/")
        self.model.run(cwd=f"./OpenMC/{self.name}/")