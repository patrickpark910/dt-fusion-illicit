import openmc
import os, sys, re
import numpy as np
import pandas as pd
from PIL import Image
from abc import ABC, abstractmethod

# Import helper functions
from Python.parameters import *
from Python.utilities import *
from Python.output import OutputMixin


class Reactor(OutputMixin, ABC):

    def __init__(self, blanket_name, fertile_kgm3=0.0, fertile_isotope='U238', lithium_enrich=None, run_type='tallies', n_particles=N_PARTICLES, n_cycles=N_CYCLES, print_xml=True, run_openmc=True, run_debug=True):

        self.fertile_kgm3 = fertile_kgm3
        self.fertile_isotope = fertile_isotope
        
        if self.fertile_isotope == 'U238': 
            self.fissile_isotope = 'Pu239'
        elif self.fertile_isotope == 'Th232': 
            self.fissile_isotope = 'U233'

        self.run_type    = run_type
        self.print_xml   = print_xml
        self.run_openmc  = run_openmc
        self.run_debug   = run_debug
        self.n_particles = n_particles
        self.n_cycles    = n_cycles

        # All class templates must define these variables:
        self.temp_k          = None           # TEMP_K
        self.blanket_name    = None           # 'ARC'
        self.blanket_volume  = None           # ARC_BL_VOL    # m³
        self.breeder_density = None           # DENSITY_FLIBE # g/cm³
        self.breeder_enrich  = lithium_enrich # ENRICH_FLIBE  # wt%
        self.name = None # f"{self.run_type}_{self.blanket_name}_{self.temp_k}K_Li{self.breeder_enrich:04.1f}_{self.fertile_isotope}{self.fertile_kgm3:06.2f}kgm3"
        self.path = None # f"./OpenMC/{self.name}"


    @abstractmethod
    def materials(self):
        pass


    @abstractmethod
    def geometry(self):
        pass
    

    def tallies(self):
        """ 
        TALLIES 
        """
        self.tallies = openmc.Tallies() # initialize

        # Filters
        void_cell_filter    = openmc.CellFilter([cell for cell in self.cells if cell.fill is None])
        cellfrom_filter     = openmc.CellFromFilter(self.cells)
        cellwithvoid_filter = openmc.CellFilter(self.cells)
        cell_filter         = openmc.CellFilter([cell for cell in self.cells if cell.fill is not None])  # 
        energy_filter       = openmc.EnergyFilter(logspace_per_decade(1e-3, 20e6, 50)) # './helpers/utilities.py'
                

        # ---------------------------------------------------------------
        # Flux and current tallies 
        # - Our Miller-shape Polygon doesn't play well with the surface_filter
        #   (bc it is the union of like 1500 conic sections) so we get currents
        #   by defining cell_filter and cellfrom_filter per docs below
        # 
        # cf. Misc scores > current > description @ docs.openmc.org/en/stable/usersguide/tallies.html 
        # ---------------------------------------------------------------

        # Current tally
        current_tally        = self.make_tally('current', ['current'], filters=[cellwithvoid_filter, cellfrom_filter])

        # Leakage current spectrum - scope CellFilter to void cells only to keep output manageable
        current_energy_tally = self.make_tally('current spectrum', ['current'],
                                               filters=[void_cell_filter, cellfrom_filter, energy_filter])

        # Flux tally 
        flux_tally        = self.make_tally('flux', ['flux'], filters=[cell_filter])
        flux_energy_tally = self.make_tally('flux spectrum', ['flux'], filters=[cell_filter, energy_filter])


        # ---------------------------------------------------------------
        # Power tallies
        # - Tally heating-local throughout the blanket, then tally fission-q-recoverable to isolate fission power contribution
        # ---------------------------------------------------------------

        heating_tally        = self.make_tally('heating', ['heating-local'], filters=[cell_filter])
        heating_energy_tally = self.make_tally('heating spectrum', ['heating-local'], filters=[cell_filter, energy_filter])  # not having cell_filter on purpose

        fisq_tally           = self.make_tally('fission-q', ['fission-q-recoverable'], filters=[cell_filter])
        fisq_energy_tally    = self.make_tally('fission-q spectrum', ['fission-q-recoverable'], filters=[cell_filter, energy_filter])  # not having cell_filter on purpose


        # ---------------------------------------------------------------
        # Reaction rate tallies 
        # ---------------------------------------------------------------

        # Fertile element reaction rates // MT 4 = inelastic scattering
        fertile_tally_tot    = self.make_tally('Total fertile rxn rate',     ['(n,gamma)', 'fission', 'nu-fission','(n,2n)', '(n,3n)', '4', 'elastic'], nuclides=['U238', 'U235', 'Th232'])
        fertile_tally        = self.make_tally('Fertile rxn rates',          ['(n,gamma)', 'fission', 'nu-fission','(n,2n)', '(n,3n)', '4', 'elastic'], nuclides=['U238', 'U235', 'Th232'], filters=[cell_filter])
        fertile_energy_tally = self.make_tally('Fertile rxn rates spectrum', ['(n,gamma)', 'fission', 'nu-fission','(n,2n)', '(n,3n)', '4', 'elastic'], nuclides=['U238', 'U235', 'Th232'], filters=[cell_filter, energy_filter])

        # Lithium reaction rates
        Li_tally_tot    = self.make_tally('Total Li rxn rate',     ['(n,gamma)', '(n,Xt)', 'elastic'], nuclides=['Li6', 'Li7'])
        Li_tally        = self.make_tally('Li rxn rates',          ['(n,gamma)', '(n,Xt)', 'elastic'], nuclides=['Li6', 'Li7'], filters=[cell_filter], )
        Li_energy_tally = self.make_tally('Li rxn rates spectrum', ['(n,gamma)', '(n,Xt)', 'elastic'], nuclides=['Li6', 'Li7'], filters=[cell_filter, energy_filter], )

        # Fluorine reaction rates
        # F_tally        = self.make_tally('F rxn rates', ['(n,gamma)', 'elastic'], filters=[cell_filter], nuclides=['F19'])
        # F_energy_tally = self.make_tally('F rxn rates spectrum', ['(n,gamma)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['F19'])

        # Beryllium reaction rates
        Be_tally_tot    = self.make_tally('Total Be rxn rate',     ['(n,gamma)', '(n,2n)', '(n,a)', 'elastic'], nuclides=['Be9'])
        Be_tally        = self.make_tally('Be rxn rates',          ['(n,gamma)', '(n,2n)', '(n,a)', 'elastic'], nuclides=['Be9'], filters=[cell_filter],                )
        Be_energy_tally = self.make_tally('Be rxn rates spectrum', ['(n,gamma)', '(n,2n)', '(n,a)', 'elastic'], nuclides=['Be9'], filters=[cell_filter, energy_filter], )

        # Lead reaction rates
        Pb_tally_tot    = self.make_tally('Total Pb rxn rate',     ['(n,gamma)', '(n,2n)', '4', 'elastic'], nuclides=['Pb204', 'Pb206', 'Pb207', 'Pb208'])
        Pb_tally        = self.make_tally('Pb rxn rates',          ['(n,gamma)', '(n,2n)', '4', 'elastic'], nuclides=['Pb204', 'Pb206', 'Pb207', 'Pb208'], filters=[cell_filter],               )
        Pb_energy_tally = self.make_tally('Pb rxn rates spectrum', ['(n,gamma)', '(n,2n)', '4', 'elastic'], nuclides=['Pb204', 'Pb206', 'Pb207', 'Pb208'], filters=[cell_filter, energy_filter],)


        self.tallies.extend([fertile_tally_tot, Li_tally_tot, Be_tally_tot, Pb_tally_tot])
        self.tallies.extend([fertile_tally, Li_tally, Be_tally, Pb_tally])
        self.tallies.extend([current_tally, flux_tally])
        self.tallies.extend([heating_tally, fisq_tally])
        self.tallies.extend([fertile_energy_tally, Li_energy_tally, Be_energy_tally, Pb_energy_tally])
        self.tallies.extend([current_energy_tally, flux_energy_tally]) 
        self.tallies.extend([heating_energy_tally, fisq_energy_tally])

    
    @staticmethod
    def make_tally(name, scores, filters:list=None, nuclides:list=None):
        tally = openmc.Tally(name=name)
        tally.scores = scores
        if filters is not None:
            tally.filters = filters
        if nuclides is not None:
            tally.nuclides = nuclides
        return tally


    def settings(self):
        """ SETTINGS """
       
        source = openmc.IndependentSource()
        source.space = openmc.stats.CylindricalIndependent(
                                      r=openmc.stats.Discrete([self.R0], [1.0]),  # r = R0 major radius [cm]
                                    phi=openmc.stats.Uniform(0.0, 2*np.pi)   ,  # phi = 0 to 2pi
                                      z=openmc.stats.Discrete([0.0], [1.0])   ) # z   = 0
                                                          
        source.particle = 'neutron'
        source.energy   = openmc.stats.Discrete([14.0e6], [1.0])  # 14 MeV
        source.angle    = openmc.stats.Isotropic()

        # Create settings and assign source
        self.settings = openmc.Settings()
        self.settings.source = source

        """ Run type """
        self.settings.run_mode  = 'fixed source'
        self.settings.particles = self.n_particles
        self.settings.batches   = self.n_cycles

        """ Dump data to statepoint intermittently to allow restarts later """
        save_sp_every = 5
        self.settings.statepoint = {'batches': range(save_sp_every, self.n_cycles + 1, save_sp_every)} 
        self.settings.output     = {'summary': False, 'tallies': False}


    def compile(self):
        
        openmc.reset_auto_ids()
        self.materials()
        self.geometry()
        self.settings()
        self.tallies()

        self.materials.cross_sections = set_xs_path()
        self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)

        if self.print_xml:
            self.model.export_to_model_xml(self.path)
            print(f"{C.YELLOW}Comment.{C.END} OpenMC model XML files printed to: {self.path}")
        else:
            print(f"{C.RED}Warning.{C.END} OpenMC model XML files NOT printed (print_xml=False) to: {self.path}")


    def openmc(self):

        if self.run_openmc:
            if has_statepoint(self.path,cycle=self.n_cycles):
                print(f"{C.YELLOW}Warning.{C.END} File {self.path}/statepoint.{str(self.n_cycles).zfill(2)}.h5 already exists, so this OpenMC run will be skipped...")
            else:
                self.model.run(cwd=self.path)

                # Clean up intermediate statepoint files -- KEEP UNDER THIS LAST 'ELSE'!! --ppark 2026-05-01
                for filename in os.listdir(self.path):
                    if filename.startswith('statepoint.') and filename.endswith('.h5'):
                        try:
                            # Extract the integer batch number from the filename
                            batch_num = int(filename.split('.')[1])  
                            # Delete if it's less than the final programmed cycle
                            if batch_num < self.n_cycles:
                                print(f"{C.YELLOW}Comment.{C.END} Removing intermediate statepoint.{batch_num}.h5 (< {self.n_cycles} cycles)...")
                                file_path = os.path.join(self.path, filename)
                                os.remove(file_path)
                        except ValueError:
                            # Safely ignore any improperly formatted statepoint files
                            pass
        else:    
            print(f"{C.RED}Warning.{C.END} OpenMC calculation skipped (run_openmc=False) for: {self.path}")


    def volumes(self, samples=int(1e10)):
        """
        Calculate volumes of all cells using OpenMC stochastic volume calculation.
        """

        self.materials()
        self.geometry()
        
        # Get all cells that have materials (exclude voids)
        cells_to_calc = self.cells # [cell for cell in self.cells if cell.fill == self.blanket] #  is not None]
        
        # Create volume calculation settings
        vol_calc = openmc.VolumeCalculation(cells_to_calc, samples,
                                            [-1*self.extent_r, -1*self.extent_r, -1*self.extent_z], # lower left of bounding box
                                            [self.extent_r, self.extent_r, self.extent_z]) # upper right of bounding box
        
        # vol_calc.set_trigger(1e-02, 'std_dev')

        settings_vol_calc = openmc.Settings()
        settings_vol_calc.volume_calculations = [vol_calc]
        
        print(f"{C.GREEN}Running stochastic volume calculation with {samples} samples...{C.END}")

        self.materials.cross_sections = set_xs_path()
        self.model_vol_calc = openmc.model.Model(self.geometry, self.materials, settings_vol_calc)
        self.model_vol_calc.calculate_volumes(cwd=self.path, export_model_xml=False, apply_volumes=False)


        # ---------------------------
        # Process volume calc results
        # ---------------------------
        vol_results = openmc.VolumeCalculation.from_hdf5(f"{self.path}/volume_1.h5")
        
        #print(f"{C.GREEN}Stochastic Volume Calculation Results:{C.END}")

        vol_dict = {}
        for k, v in vol_results.volumes.items():
            vol_dict[k] = (v.nominal_value/1e6, v.std_dev/1e6)
        # vol_dict['sum'] = ( sum(v.nominal_value/1e6 for v in vol_results.volumes.values()),
        #                     sum(v.std_dev/1e6 for v in vol_results.volumes.values())       )

        df = pd.DataFrame.from_dict(vol_dict, orient='index', columns=['volume_m3', 'stdev_m3'])
        df.reset_index(inplace=True)
        df.rename(columns={'index': 'cells'}, inplace=True)
        df.to_csv(f'{self.path}/volume_1.csv',index=False)

        print(f"{C.YELLOW}Comment.{C.END} CSV file 'volume_1.csv' printed to {self.path}")
        # print(df)

        """
        vol_results.volumes is a dictionary that looks like:
          {2: "909560.9858392784+/-140318.8308700002",
           3: "5543990.770829888+/-346055.6196989608"}
        so we will rewrite it to separate the value from the standev, like this:
          {2: (909560.9858392784, 140318.8308700002),
           3: (5543990.770829888, 346055.6196989608)}

        But also whoever designed the output of VolueCalculation.volumes() 
        honestly needs to WRITE in the documentation that the data is stored/written 
        using the 'uncertainties' package. I was going crazy using all sorts of 
        regex match patterns to try to understand why splitting the vol_results.volumes.items()
        into k, v was causing issues vs. what v looked like when I was printing it. 
        It was 'uncertainties' changing the formatting of v when you go print it. 
          --ppark  2025-09-20
        """


    def plot(self):

        self.materials()
        self.geometry()

        Image.MAX_IMAGE_PIXELS = None # suppreses DecompressionBombWarning lmao

        settings_plot = openmc.Settings()
        settings_plot.run_mode = 'plot'

        if self.blanket_name in ['ARC', 'ARCBall', 'FLiBe']:
            C = {self.firstwall: (30, 27, 41),    # very dark gray
                      self.structure: (109, 110, 113), # gray
                      self.blanket: (129, 204, 185),   # teal
                      # Void regions will be white by default
                      }

        elif self.blanket_name in ['DCLL']:
            C = {self.firstwall: (30, 27, 41),           # very dark gray
                      self.structure: (109, 110, 113),        # gray
                      self.coolant: (37, 150, 190), # cobalt blue 
                      self.blanket: (129, 204, 185),          # teal
                      self.divider: (109, 110, 113),          # gray
                      self.inner_manifold:(176, 123, 76),     # wood-color
                      # Void regions will be white by default
                      }

        elif self.blanket_name in ['HCPB']:
            C = {self.firstwall: (30, 27, 41),           # very dark gray
                      self.structure: (109, 110, 113),        # gray
                      self.blanket: (129, 204, 185),          # teal
                      # Void regions will be white by default
                      }

        x_width = round(2*self.extent_r)
        z_width = round(2*self.extent_z)
        
        # XY toroidal slice
        xy = openmc.Plot()
        xy.filename = f"{self.blanket_name}_xy" # {self.path}/
        xy.basis = "xy"
        xy.width  = (x_width, x_width)
        xy.pixels = (10*x_width, 10*x_width)
        xy.color_by = "material"
        xy.C = C

        # XZ poloidal slice
        xz = openmc.Plot()
        xz.filename = f"{self.blanket_name}_xz" # {self.path}/
        xz.basis = "xz"
        xz.width  = (x_width, z_width)
        xz.pixels = (10*x_width, 10*z_width)
        xz.color_by = "material"
        xz.C = C

        plots = openmc.Plots([xy, xz])
        model_plot = openmc.model.Model(self.geometry, self.materials, settings_plot)
        model_plot.plots = plots 
        model_plot.plot_geometry(cwd=f"{self.path}")  # writes tokamak_rz.ppm and tokamak_xz.ppm

        for basename in [f"{self.blanket_name}_xy", f"{self.blanket_name}_xz"]:
            ppm_file = os.path.join(self.path, f"{basename}.ppm")
            png_file = os.path.join(self.path, f"{basename}.png")
            if os.path.exists(ppm_file):
                with Image.open(ppm_file) as im:
                    im.save(png_file)
                print(f"{C.YELLOW}Comment.{C.END} Plots '{self.blanket_name}_xy.png', '{self.blanket_name}_xz.png' printed to {self.path}")
            else:
                print(f"{C.YELLOW}Error.{C.END} OpenMC did not print '{self.blanket_name}_xy.ppm', '{self.blanket_name}_xz.ppm' to {self.path}!")