import openmc
import os, sys, re
import numpy as np
import pandas as pd
from PIL import Image
from abc import ABC, abstractmethod

# Import helper functions
from Python.parameters import *
from Python.utilities import *


class Reactor(ABC):

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

        # Fertile element reaction rates
        fertile_tally_tot    = self.make_tally('Total fertile rxn rate',     ['(n,gamma)', 'fission', '(n,2n)', 'elastic'], nuclides=['U238', 'U235', 'Th232'])
        fertile_tally        = self.make_tally('Fertile rxn rates',          ['(n,gamma)', 'fission', '(n,2n)', 'elastic'], filters=[cell_filter], nuclides=['U238', 'U235', 'Th232'])
        fertile_energy_tally = self.make_tally('Fertile rxn rates spectrum', ['(n,gamma)', 'fission', '(n,2n)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['U238', 'U235', 'Th232'])

        # Lithium reaction rates
        Li_tally_tot    = self.make_tally('Total Li rxn rate',     ['(n,gamma)', '(n,Xt)', 'elastic'], nuclides=['Li6', 'Li7'])
        Li_tally        = self.make_tally('Li rxn rates',          ['(n,gamma)', '(n,Xt)', 'elastic'], filters=[cell_filter], nuclides=['Li6', 'Li7'])
        Li_energy_tally = self.make_tally('Li rxn rates spectrum', ['(n,gamma)', '(n,Xt)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['Li6', 'Li7'])

        # Fluorine reaction rates
        # F_tally        = self.make_tally('F rxn rates', ['(n,gamma)', 'elastic'], filters=[cell_filter], nuclides=['F19'])
        # F_energy_tally = self.make_tally('F rxn rates spectrum', ['(n,gamma)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['F19'])

        # Beryllium reaction rates
        Be_tally        = self.make_tally('Be rxn rates',          ['(n,gamma)', '(n,2n)', 'elastic'], filters=[cell_filter],                nuclides=['Be9'])
        Be_energy_tally = self.make_tally('Be rxn rates spectrum', ['(n,gamma)', '(n,2n)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['Be9'])
        Pb_tally        = self.make_tally('Pb rxn rates',          ['(n,gamma)', '(n,2n)', 'elastic'], filters=[cell_filter],                nuclides=['Pb204', 'Pb206', 'Pb207', 'Pb208'])
        Pb_energy_tally = self.make_tally('Pb rxn rates spectrum', ['(n,gamma)', '(n,2n)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['Pb204', 'Pb206', 'Pb207', 'Pb208'])


        self.tallies.extend([fertile_tally_tot, Li_tally_tot])
        self.tallies.extend([fertile_tally, Li_tally, Be_tally, Pb_tally])
        self.tallies.extend([current_tally, flux_tally])
        self.tallies.extend([heating_tally, fisq_tally])
        self.tallies.extend([fertile_energy_tally, Li_energy_tally, Be_energy_tally, Pb_energy_tally])
        self.tallies.extend([current_energy_tally, flux_energy_tally]) 
        self.tallies.extend([heating_energy_tally, fisq_energy_tally])


    def make_tally(self, name, scores, filters:list=None, nuclides:list=None):
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
            if has_statepoint(self.path):
                print(f"{C.YELLOW}Warning.{C.END} File {self.path}/statepoint.XX.h5 already exists, so this OpenMC run will be skipped...")
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


    def extract_tallies(self):

        """ Load tallies """
        sp_path = f'{self.path}/statepoint.{self.n_cycles}.h5'
        
        # First try finding statepoint that exactly matches the programmed n_cycles
        try:
            print(f"{C.YELLOW}Comment.{C.END} Loading statepoint: {sp_path}")
            sp = openmc.StatePoint(sp_path) 

        except Exception as e:
            print(f"{C.YELLOW}Warning.{C.END} {e}. statepoint.{self.n_cycles}.h5 missing or could not be read at: {sp_path}")
        
            # Second try finding highest statepoint
            try:
                print(f"{C.YELLOW}Comment.{C.END} Looking for the highest statepoint in: {self.path}")

                sp_files = [f for f in os.listdir(self.path) if f.startswith('statepoint.') and f.endswith('.h5')]

                if not sp_files:
                    print(f"{C.RED}Warning.{C.END} <reactor.py/extract_tallies> No statepoints found at all in: {self.path}")
                    print(f"But, reactor.py/openmc() should have run OpenMC if: (1) there are no statepoints (2) you have self.run_openmc=True (current state: {self.run_openmc})."
                          f"If you see this error, self.run_openmc=False or you might have broken something in our logic.")
                    sys.exit(2)
                    
                # Get the latest statepoint by batch number
                sp_path = os.path.join(self.path, max(sp_files, key=lambda x: int(x.split('.')[1])))
                sp = openmc.StatePoint(sp_path) 

            except Exception as e:
                print(f"{C.RED}  Fatal.{C.END} <reactor.py/extract_tallies> {e}. No valid statepoints found at all in: {self.path}")
                sys.exit(2)

        
        print(f"Reading tallies...")

        """
        Convert tallies into usable forms. These dataframes have headers:
        cell  energy low [eV]  energy high [eV] nuclide score     mean  std. dev.  energy mid [eV]
        """

        flux         = sp.get_tally(name='flux').get_pandas_dataframe()
        flux_spec    = sp.get_tally(name='flux spectrum').get_pandas_dataframe()
        Li           = sp.get_tally(name='Li rxn rates').get_pandas_dataframe()
        Li_spec      = sp.get_tally(name='Li rxn rates spectrum').get_pandas_dataframe()
        Be           = sp.get_tally(name='Be rxn rates').get_pandas_dataframe()
        fertile      = sp.get_tally(name='Fertile rxn rates').get_pandas_dataframe()
        fertile_spec = sp.get_tally(name='Fertile rxn rates spectrum').get_pandas_dataframe()
        current_df   = sp.get_tally(name='current').get_pandas_dataframe()
        heating      = sp.get_tally(name='heating').get_pandas_dataframe()
        fisq         = sp.get_tally(name='fission-q').get_pandas_dataframe()

        # Add new column for energy bin midpoint (for plotting)
        for df in [flux_spec, fertile_spec, Li_spec]:
            df['energy mid [eV]'] = (df['energy low [eV]'] + df['energy high [eV]'])/ 2


        # =====================================================================
        # PROCESS LEAKAGE CURRENT INTO VOID CELLS
        # =====================================================================
        # Identify void vs material cells by fill
        void_ids    = [cell.id for cell in self.cells if cell.fill is None and cell.id != 10]  # exclude plasma
        material_ids = [23, 22, 42, 44]  # these are the outermost blanket cells for FLiBe (23), HCPB (22, 23), DCLL (42, 44)

        # Filter: destination is a void cell, origin is a material cell
        void_current_df = current_df[
            (current_df['cell'].isin(void_ids)) &
            (current_df['cellfrom'].isin(material_ids))
        ].copy()

        # Total leakage
        total_leakage_current = void_current_df['mean'].sum()
        total_leakage_err = np.sqrt((void_current_df['std. dev.']**2).sum())

        # Per-void-cell breakdown (inboard vs outboard for DCLL)
        leakage_rows = []
        for vid in void_ids:
            mask = void_current_df['cell'] == vid
            leak_val = void_current_df.loc[mask, 'mean'].sum()
            leak_err = np.sqrt((void_current_df.loc[mask, 'std. dev.']**2).sum())
            # Which material cell(s) actually contribute
            neighbors = void_current_df.loc[mask & (void_current_df['mean'] > 0), 'cellfrom'].unique().tolist()
            leakage_rows.append({'void_cell': vid, 'cellfrom': str(neighbors),
                                 'leakage': leak_val, 'leakage_stdev': leak_err})

        leakage_rows.append({'void_cell': 'total', 'cellfrom': '',
                             'leakage': total_leakage_current, 'leakage_stdev': total_leakage_err})

        leakage_df = pd.DataFrame(leakage_rows)
        leakage_df.insert(0, 'fertile_kgm3', self.fertile_kgm3)


        # =====================================================================
        # PROCESS LEAKAGE CURRENT SPECTRUM INTO VOID (CELL 99)
        # =====================================================================
        current_spec_df = sp.get_tally(name='current spectrum').get_pandas_dataframe()
        current_spec_df['energy mid [eV]'] = (current_spec_df['energy low [eV]'] + current_spec_df['energy high [eV]']) / 2

        # Filter: destination is cell 99, exclude self-crossing (cellfrom != 99)
        leakage_spec_df = current_spec_df[
            (current_spec_df['cell'] == 99) &
            (current_spec_df['cellfrom'] != 99)
        ].copy()

        # Sum partial currents over all cellfrom origins per energy bin
        leakage_Ebin_df = ( leakage_spec_df.groupby(['energy low [eV]', 'energy high [eV]', 'energy mid [eV]'])
                                           .agg(mean=('mean', 'sum'),std_dev=('std. dev.', lambda x: np.sqrt((x**2).sum())))
                                           .reset_index() )
        leakage_Ebin_df.columns = ['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'mean', 'std. dev.']


        # =====================================================================
        # PROCESS STANDARD REACTION RATES
        # =====================================================================
        # Flux spectrum dataframe
        flux_Ebin_df     = flux_spec[flux_spec['cell'].between(30, 39, inclusive='both')]  

        U238_ng_Ebin_df  = fertile_spec[(fertile_spec['nuclide'] == self.fertile_isotope) & (fertile_spec['score'] == '(n,gamma)')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        U238_fis_Ebin_df = fertile_spec[(fertile_spec['nuclide'] == self.fertile_isotope) & (fertile_spec['score'] == 'fission')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        Li6_nt_Ebin_df   = Li_spec[(Li_spec['nuclide'] == 'Li6') & (Li_spec['score'] == '(n,Xt)')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        Li7_nt_Ebin_df   = Li_spec[(Li_spec['nuclide'] == 'Li7') & (Li_spec['score'] == '(n,Xt)')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]


        # Flux
        flux_list = flux['mean'].tolist()

        # Lithium reaction rates 
        Li6_nt = Li[(Li['nuclide']=='Li6') & (Li['score']=='(n,Xt)')][['cell','score','mean','std. dev.']] 
        Li7_nt = Li[(Li['nuclide']=='Li7') & (Li['score']=='(n,Xt)')][['cell','score','mean','std. dev.']] 
        
        Li6_nt_list     = Li6_nt['mean'].tolist()
        Li7_nt_list     = Li7_nt['mean'].tolist() 
        Li6_nt_err_list = Li6_nt['std. dev.'].tolist()
        Li7_nt_err_list = Li7_nt['std. dev.'].tolist()

        # Beryllium (n, 2n) reaction rates
        Be9_n2n = Be[(Be['nuclide'] == 'Be9') & (Be['score'] == '(n,2n)')][['cell', 'score', 'mean', 'std. dev.']]
        Be9_n2n_list     = Be9_n2n['mean'].tolist()
        Be9_n2n_err_list = Be9_n2n['std. dev.'].tolist()
        
        tbr_list = [x + y for x, y in zip(Li6_nt_list, Li7_nt_list)] 
        tbr_err_list = [x + y for x, y in zip(Li6_nt_err_list, Li7_nt_err_list)]

        # U-238 or Th-232 reaction rates 
        U238              = fertile[(fertile['nuclide']==self.fertile_isotope)][['cell','score','mean','std. dev.']]
        U238_fis_list     = U238[U238['score'] == 'fission'][['mean']]['mean'].tolist()
        U238_fis_err_list = U238[U238['score'] == 'fission'][['std. dev.']]['std. dev.'].tolist()
        U238_ng_list      = U238[U238['score'] == '(n,gamma)'][['mean']]['mean'].tolist()
        U238_ng_err_list  = U238[U238['score'] == '(n,gamma)'][['std. dev.']]['std. dev.'].tolist()
        
        # Heating and fission-q-recoverable, converted to MW
        heat_list     = (heating['mean'] * NPS_FUS * EV_TO_MJ).tolist()
        heat_err_list = (heating['std. dev.'] * NPS_FUS * EV_TO_MJ).tolist()
        fisq_list     = (fisq['mean'] * NPS_FUS * EV_TO_MJ).tolist()
        fisq_err_list = (fisq['std. dev.'] * NPS_FUS * EV_TO_MJ).tolist()


        # Determine the atomic mass based on the isotope
        amu = AMU_PU239 if self.fertile_isotope == 'U238' else AMU_U233

        # Scaling factor
        scaling_factor = NPS_FUS * SEC_PER_YR * amu / AVO / 1e3  # kg/yr

        fissile_per_yr_list = [Pu_per_srcn * scaling_factor for Pu_per_srcn  in U238_ng_list]
        fissile_per_yr_err_list = [err * scaling_factor for err in U238_ng_err_list]


        """Create list of cell IDs randomly assigned by OpenMC to match mtu loading to cell ID"""
        cell_ids = [str(x) for x in flux['cell'].unique().tolist()] 

        df = pd.DataFrame({'cell': cell_ids,
                           'flux': flux_list,
                           'tbr'               : tbr_list,
                           'tbr_stdev'         : tbr_err_list,
                           f'{self.fertile_isotope}(n,g)'         : U238_ng_list,
                           f'{self.fertile_isotope}(n,g)_stdev'   : U238_ng_err_list,
                           f'{self.fissile_isotope}_kg/yr'        : fissile_per_yr_list,
                           f'{self.fissile_isotope}_kg/yr_stdev'  : fissile_per_yr_err_list,
                           'Li6(n,t)'          : Li6_nt_list,
                           'Li6(n,t)_stdev'    : Li6_nt_err_list,
                           'Li7(n,t)'          : Li7_nt_list,
                           'Li7(n,t)_stdev'    : Li7_nt_err_list,
                           'Be9(n,2n)'         : Be9_n2n_list,
                           'Be9(n,2n)_stdev'   : Be9_n2n_err_list,
                           f'{self.fertile_isotope}(n,fis)'       : U238_fis_list,
                           f'{self.fertile_isotope}(n,fis)_stdev' : U238_fis_err_list,
                           'heating [MW]'      : heat_list ,
                           'heating_stdev'     : heat_err_list,
                           'fisq [MW]'         : fisq_list,
                           'fisq_stdev'        : fisq_err_list })

        totals = df.sum(numeric_only=True)
        totals['cell'] = 'total'
        df = pd.concat([df, pd.DataFrame([totals])], ignore_index=True)

        print(df.iloc[:, :8])  # prints first 8 columns 
        
        # Export the Dataframes
        df.to_csv(f'./{self.path}/tallies_summary.csv', index=False)
        leakage_df.to_csv(f'./{self.path}/leakage.csv', index=False)
        U238_ng_Ebin_df.to_csv(f'./{self.path}/{self.fertile_isotope}_n-gamma_Ebins.csv', index=False)
        flux_Ebin_df.to_csv(f'./{self.path}/{self.fertile_isotope}_flux_Ebins.csv', index=False)
        leakage_Ebin_df.to_csv(f'./{self.path}/{self.fertile_isotope}_leakage_Ebins.csv', index=False)
