import openmc
import os, sys
import numpy as np
from datetime import date
from shapely.geometry import Polygon as SPoly


# Import helper functions
from Python.parameters import *
from Python.utilities import *

class Reactor:

    def __init__(self, breeder='flibe', fertile_element='U', fertile_density_kgm3=0.0, run_openmc=False):

        self.fertile_density_kgm3 = fertile_density_kgm3
        self.fertile_element = fertile_element.capitalize()
        self.run_openmc = run_openmc

        if breeder.lower() == 'flibe':
            self.temp_k          = TEMP_FLIBE_K
            self.breeder         = 'FLiBe'
            self.breeder_density = DENSITY_FLIBE
            self.breeder_enrich  = ENRICH_FLIBE
            # self.breeder_volume  = VOLUME_FLIBE
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
        self.name = f"{self.breeder}_{self.fertile_element}{self.fertile_density_kgm3:.2f}kgm3_Li{self.breeder_enrich:04.1f}_{self.temp_k}K" # _{today}
        # should come out to smth like: FLiBe_U010kgm3_Li7.5_900K
        self.path = f"./OpenMC/{self.name}/"
        
        os.makedirs(self.path, exist_ok=True)

        start_msg = f"""
        Writing an OpenMC input for {self.name}
        """
        print(start_msg)

        self.materials()
        print('1')
        self.geometry()
        print('2')
        # self.settings()
        print('3')
        # self.tallies()
        print('4')
        # self.run_openmc()
        print('5')


    def materials(self):
        """
        OpenMC Materials
        """

        tungsten = openmc.Material(name='W')
        tungsten.add_element('O',5/1e6,percent_type='wo')
        tungsten.add_element('N',5/1e6,percent_type='wo')
        tungsten.add_element('C',5/1e6,percent_type='wo')
        tungsten.add_element('Na',4/1e6,percent_type='wo')
        tungsten.add_element('K',2.5/1e6,percent_type='wo')
        tungsten.add_element('Al',3/1e6,percent_type='wo')
        tungsten.add_element('Ca',0.5/1e6,percent_type='wo')
        tungsten.add_element('Cr',0.5/1e6,percent_type='wo')
        tungsten.add_element('Cu',0.5/1e6,percent_type='wo')
        tungsten.add_element('Fe',5/1e6,percent_type='wo')

        tungsten.add_element('W',1-(5+5+5+4+2.5+3+0.5+0.5+0.5+5)/1e6,percent_type='wo')
        tungsten.set_density('g/cm3',19.3)
        tungsten.depletable = False

        """ V-4Cr-4Ti """
        # This code is from jlball but V-4Cr-4Ti specs also can be found
        # from ANL material id BL-47, p.3 (1994) altho these are slightly different
        # https://www.osti.gov/servlets/purl/10194461 --ppark 2025-07-22

        vcrti = openmc.Material(name='V-4Cr-4Ti VV')
        vcrti.depletable = False

        vcrti.add_element('Cr',0.04,percent_type='wo')
        vcrti.add_element('Ti',0.04,percent_type='wo')

        vcrti.add_element('C',56/1e6,percent_type='wo')
        vcrti.add_element('O',181/1e6,percent_type='wo')
        vcrti.add_element('N',103/1e6,percent_type='wo')
        vcrti.add_element('B',7/1e6,percent_type='wo')
        vcrti.add_element('Na',17/1e6,percent_type='wo')
        vcrti.add_element('Mg',0.5/1e6,percent_type='wo')
        vcrti.add_element('Al',119/1e6,percent_type='wo')
        vcrti.add_element('Si',280/1e6,percent_type='wo')
        vcrti.add_element('Mn',0.5/1e6,percent_type='wo')
        vcrti.add_element('Fe',80/1e6,percent_type='wo')
        vcrti.add_element('Ni',13/1e6,percent_type='wo')
        vcrti.add_element('Cu',4/1e6,percent_type='wo')
        vcrti.add_element('V',1-0.04-0.04-(56+181+103+7+17+0.5+119+280+0.5+80+13+4)/1e6,percent_type='wo')
        vcrti.set_density('g/cm3',6.05) 
        # This density value is sus and needs a good source --jlball 
        # This value is from Metals Handbook, 9th ed, vol 2: "Properties and Selection: Nonferrous Alloys and Pure Metals" (1979) --ppark 2025-07-22



        # Breeder
        breeder = openmc.Material()
        breeder.set_density('g/cm3', self.breeder_density)

        if self.breeder == 'FLiBe':
            breeder.add_elements_from_formula('F4Li2Be', 'ao', enrichment_target='Li6', enrichment_type='wo', enrichment=self.breeder_enrich)

        # elif self.breeder.lower() == 'll':
        #     breeder.add_element('Pb', 0.83)
        #     breeder.add_element('Li', 0.17, enrichment=self.breeder_enrich, enrichment_target='Li6', enrichment_type='wo')  # Li-6 enrichment to 90%

        # elif self.blanket.lower() == 'pb':
        #     self.blanket_density = DENSITY_PB
        #     self.blanket_enrich  = ENRICH_PB
        

        # Fertile compound
        fertile = openmc.Material()

        if self.breeder == 'FLiBe':
            if self.fertile_element == 'U':
                fertile.add_elements_from_formula('UF4','ao') # 'ao' here refers to 1:4 atomic ratio of U:F in UF4
                fertile.set_density('g/cm3', DENSITY_UF4) 
            elif self.fertile_element == 'Th':
                fertile.add_elements_from_formula('ThF4','ao') # 'ao' here refers to 1:4 atomic ratio of U:F in UF4
                fertile.set_density('g/cm3', DENSITY_ThF4) 

        elif self.breeder in ['LL', 'PB']:
            if self.fertile_element == 'U':
                pass  # UO2 pebbles
            elif self.fertile_element == 'Th':
                pass  # ThO2 pebbles


        # Mix fertile compound to blanket
        if self.breeder == 'FLiBe':
            breeder_mass_frac, fertile_compound_mass_frac = calc_blanket_mass_fracs(self.fertile_density_kgm3, 
                                                                                    fertile_element=self.fertile_element, 
                                                                                    fertile_enrich=0.71, 
                                                                                    breeder_density_kgm3=1.94e3)
            blanket = openmc.Material.mix_materials([breeder, fertile], [breeder_mass_frac, fertile_compound_mass_frac], 'wo') # fractions in 'mix_materials' MUST add up to 1
            blanket.name, blanket.temperature = self.name, self.temp_k
            self.materials = openmc.Materials([tungsten,vcrti,blanket])

        # elif self.breeder in ['LL', 'PB']:
        #     pass

        self.materials.export_to_xml(self.path)


    def geometry(self):

        R0, a, kappa, delta = 4.0, 1.2, 1.84, 0.5

        L = [
            ('FW',  FLIBE_FW_CM,  'steel'),
            ('ST1', FLIBE_ST1_CM, 'steel'),
            ('BR1', FLIBE_BR1_CM, 'flibe'),
            ('ST2', FLIBE_ST2_CM, 'steel'),
            ('BR2', FLIBE_BR2_CM, 'flibe'),
            ('ST3', FLIBE_ST3_CM, 'steel'),
        ]

        # --- helpers ---
        def omc_poly(spoly):  # shapely -> OpenMC polygon (closed loop)
            return openmc.Polygon(np.asarray(spoly.exterior.coords))

        # base plasma boundary (in cm), shapely polygon
        p_in = SPoly(miller_points(R0, a, delta, kappa))

        # Z slab (extrude to simple 3D)
        zmin = openmc.ZPlane(z0=-50.0, boundary_type='vacuum')
        zmax = openmc.ZPlane(z0=+50.0, boundary_type='vacuum')

        # bounding box with vacuum boundaries
        tot = sum(t for _, t, _ in L); aout = a + tot/100.0
        xmin = (R0 - aout)*100 - 50; xmax = (R0 + aout)*100 + 50
        ymin = -kappa*aout*100 - 50; ymax =  kappa*aout*100 + 50
        xlo = openmc.XPlane(x0=xmin, boundary_type='vacuum'); xhi = openmc.XPlane(x0=xmax, boundary_type='vacuum')
        ylo = openmc.YPlane(y0=ymin, boundary_type='vacuum'); yhi = openmc.YPlane(y0=ymax, boundary_type='vacuum')
        box = +xlo & -xhi & +ylo & -yhi

        # build layers by true outward offsets
        cells = []
        for name, t_cm, mname in L:
            p_out = p_in.buffer(t_cm, resolution=48)             # exact offset by t_cm (cm)
            inner, outer = omc_poly(p_in), omc_poly(p_out)
            region = box & -outer & +inner & +zmin & -zmax       # between inner & outer
            cells.append(openmc.Cell(name=name, fill=mat[mname], region=region))
            p_in = p_out                                         # advance outward

        # void region between outermost layer and box
        outermost = omc_poly(p_in)
        cells.append(openmc.Cell(name='void_outside', region=box & +outermost & +zmin & -zmax))

        self.geometry = openmc.Geometry(openmc.Universe(cells=cells))
        self.geometry.export_to_xml(self.path)

        # ---------- Plots.xml + generate images ----------
        # RZ poloidal cross-section (best for axisymmetric Miller shapes)
        rz = openmc.Plot()
        rz.filename = "tokamak_xz"
        rz.basis = "xz"
        rz.width  = (2*R_world_cm, 2*R_world_cm)  # show everything
        rz.pixels = (1200, 1200)
        rz.color_by = "cell"

        # XZ Cartesian slice (also useful)
        xz = openmc.Plot()
        xz.filename = "tokamak_xz"
        xz.basis = "xz"
        xz.width  = (2*R_world_cm, 2*R_world_cm)
        xz.pixels = (1200, 1200)
        xz.color_by = "cell"

        openmc.Plots([rz, xz]).export_to_xml(self.path)
        openmc.plot_geometry(path_input=self.path)  # writes tokamak_rz.ppm and tokamak_xz.ppm
        print("Done. Files: tokamak_rz.ppm, tokamak_xz.ppm")



    def tallies(self):
        """ TALLIES """
        self.tallies = openmc.Tallies() # initialize

        # Filters
        cell_filter = openmc.CellFilter(cells)

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


    def settings(self):
        """ SETTINGS """
        self.settings = openmc.Settings()

        """ Source from First Wall Effects """
        # Determine a flux spectrum from the neutrons leaving the surface of a two layer first wall model './Python/FirstWall.py'
        # The neutrons added to each MTU box in our model reflect the outgoing current spectrum of out vanadium shell
        sp = openmc.StatePoint(f'./OpenMC/FirstWall_V4cm_900K_2025-07-22/statepoint.100.h5')
        out_tally = sp.get_tally(name='outgoing_spectrum')

        energy_bins = np.array(out_tally.filters[1].bins) # shape (N, 2)
        energy_midpoints = 0.5 * (energy_bins[:, 0] + energy_bins[:, 1])

        current_spectrum = out_tally.get_values(scores=['current']).flatten() # Neutron generation at each current in the spectrum is weighted by its probability of occurrence
        total_current = current_spectrum.sum()
        probabilities = current_spectrum / total_current

        energies, weights = energy_midpoints.tolist(), probabilities.tolist()
        source = []

        src = openmc.IndependentSource()
        src.space  = openmc.stats.Point((p,0,0))
        src.angle  = openmc.stats.Isotropic()
        src.energy = openmc.stats.Discrete(energies, weights)
        source.append(src)
        self.settings.source = source


        """ Run type """
        self.settings.run_mode = 'fixed source'
        self.settings.particles = int(1e6)
        self.settings.batches = 100


    def run_openmc(self):
        self.materials.cross_sections = set_xs_path()
        self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
        self.model.export_to_model_xml(f"./OpenMC/{self.name}/")
        self.model.run(cwd=f"./OpenMC/{self.name}/")