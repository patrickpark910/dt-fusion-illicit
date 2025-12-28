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

    def __init__(self, breeder_name='FLiBe', fertile_element='U', fertile_bulk_density_kgm3=0.0, run_type='tallies', run_openmc=False):

        self.fertile_bulk_density_kgm3 = fertile_bulk_density_kgm3


        self.fertile_element = fertile_element
        if fertile_element == 'U': 
            self.fertile_isotope = 'U238'
            self.fissile_isotope = 'Pu239'
        elif fertile_element == 'Th': 
            self.fertile_isotope = 'Th232'
            self.fissile_isotope = 'U233'

        self.run_type    = run_type
        self.run_openmc  = run_openmc
        self.n_particles = N_PARTICLES # int(1e6)
        self.n_cycles    = N_CYCLES    # 10

        # All class templates must define these variables:
        self.temp_k          = None # TEMP_K
        self.breeder_name    = None # 'ARC'
        self.breeder_density = None # DENSITY_FLIBE # g/cm³
        self.breeder_enrich  = None # ENRICH_FLIBE  # wt%
        self.breeder_volume  = None # ARC_BR_VOL    # m³
        self.name = None # f"{run_type}_{breeder_name}_{self.fertile_element}{fertile_bulk_density_kgm3:3.2f}kgm3_Li{self.breeder_enrich:04.1f}_{self.temp_k}K"         
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

        energy_filter = openmc.EnergyFilter(logspace_per_decade(1e-5, 20e6, 100)) # './helpers/utilities.py'
        
        # For PBHet and PBHom: use material-based tallies instead of cell-based
        if self.breeder_name == 'PBHet' or self.breeder_name == 'PBHom':
            # Material filters
            blanket_material_filter = openmc.MaterialFilter([self.blanket])
            
            # For PBHet: U/Th is in kernel material (BISO particles)
            # For PBHom: U/Th is homogenized in blanket material (no BISO, no separate kernel)
            if self.breeder_name == 'PBHet' and hasattr(self, 'kernel'):
                # PBHet has BISO particles with kernel material
                fertile_material_filter = openmc.MaterialFilter([self.kernel])
            else:
                # PBHom: U/Th is homogenized in blanket material (no separate kernel material)
                fertile_material_filter = blanket_material_filter
            
            # Flux tally - filtered by blanket material (where most reactions occur)
            flux_tally = self.make_tally('flux', ['flux'], filters=[blanket_material_filter])
            flux_energy_tally = self.make_tally('flux spectrum', ['flux'], filters=[blanket_material_filter, energy_filter])
            self.tallies.extend([flux_tally, flux_energy_tally])
            
            # Fertile element (n,gamma) reaction rates
            if self.fertile_element == 'U':
                fertile_ng_tally = self.make_tally('U n-gamma', ['(n,gamma)'], 
                                                   filters=[fertile_material_filter], 
                                                   nuclides=['U238', 'U235'])
                fertile_ng_energy_tally = self.make_tally('U n-gamma spectrum', ['(n,gamma)'], 
                                                          filters=[fertile_material_filter, energy_filter], 
                                                          nuclides=['U238', 'U235'])
                self.tallies.extend([fertile_ng_tally, fertile_ng_energy_tally])
            
            elif self.fertile_element == 'Th':
                fertile_ng_tally = self.make_tally('Th n-gamma', ['(n,gamma)'], 
                                                   filters=[fertile_material_filter], 
                                                   nuclides=['Th232'])
                fertile_ng_energy_tally = self.make_tally('Th n-gamma spectrum', ['(n,gamma)'], 
                                                           filters=[fertile_material_filter, energy_filter], 
                                                           nuclides=['Th232'])
                self.tallies.extend([fertile_ng_tally, fertile_ng_energy_tally])
            
            # Lithium tritium breeding reactions - filtered by blanket material
            # Li6 (n,Xt) for tritium production
            Li6_nt_tally = self.make_tally('Li6 n-Xt', ['(n,Xt)'], 
                                         filters=[blanket_material_filter], 
                                         nuclides=['Li6'])
            Li6_nt_energy_tally = self.make_tally('Li6 n-Xt spectrum', ['(n,Xt)'], 
                                                  filters=[blanket_material_filter, energy_filter], 
                                                  nuclides=['Li6'])
            
            # Li7 (n,Xt) for tritium production
            Li7_nt_tally = self.make_tally('Li7 n-Xt', ['(n,Xt)'], 
                                          filters=[blanket_material_filter], 
                                          nuclides=['Li7'])
            Li7_nt_energy_tally = self.make_tally('Li7 n-Xt spectrum', ['(n,Xt)'], 
                                                  filters=[blanket_material_filter, energy_filter], 
                                                  nuclides=['Li7'])
            
            self.tallies.extend([Li6_nt_tally, Li6_nt_energy_tally, Li7_nt_tally, Li7_nt_energy_tally])
        
        else:
            cell_filter = openmc.CellFilter([cell for cell in self.cells if cell.fill is not None])
        
            heat_tally        = self.make_tally('volumetric heating',          ['heating-local'], filters=[cell_filter])
            heat_energy_tally = self.make_tally('volumetric heating spectrum', ['heating-local'], filters=[cell_filter, energy_filter])
            self.tallies.extend([heat_tally, heat_energy_tally])
            



            # ---------------------------------------------------------------
            # Flux and reaction rate tallies 
            # ---------------------------------------------------------------

            # Flux tally 
            flux_tally        = self.make_tally('flux', ['flux'], filters=[cell_filter])
            flux_energy_tally = self.make_tally('flux spectrum', ['flux'], filters=[cell_filter, energy_filter])

            # Fertile element reaction rates
            if self.fertile_element == 'U':
                fertile_tally        = self.make_tally(f'U rxn rates',          ['(n,gamma)', 'fission', 'elastic'], filters=[cell_filter], nuclides=['U238', 'U235'])
                fertile_energy_tally = self.make_tally(f'U rxn rates spectrum', ['(n,gamma)', 'fission', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['U238', 'U235'])

            elif self.fertile_element == 'Th':
                fertile_tally        = self.make_tally(f'Th rxn rates',          ['(n,gamma)', 'fission', 'elastic'], filters=[cell_filter], nuclides=['Th232'])
                fertile_energy_tally = self.make_tally(f'Th rxn rates spectrum', ['(n,gamma)', 'fission', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['Th232'])

            # Lithium reaction rates
            Li_tally        = self.make_tally('Li rxn rates',          ['(n,gamma)', '(n,Xt)', 'elastic'], filters=[cell_filter], nuclides=['Li6', 'Li7'])
            Li_energy_tally = self.make_tally('Li rxn rates spectrum', ['(n,gamma)', '(n,Xt)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['Li6', 'Li7'])

            # Fluorine reaction rates
            F_tally        = self.make_tally('F rxn rates', ['(n,gamma)', 'elastic'], filters=[cell_filter], nuclides=['F19'])
            F_energy_tally = self.make_tally('F rxn rates spectrum', ['(n,gamma)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['F19'])

            # Beryllium reaction rates
            Be_tally        = self.make_tally('Be rxn rates', ['(n,gamma)', '(n,2n)', 'elastic'], filters=[cell_filter], nuclides=['Be9'])
            Be_energy_tally = self.make_tally('Be rxn rates spectrum', ['(n,gamma)', '(n,2n)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['Be9'])

            self.tallies.extend([flux_tally, fertile_tally, Li_tally, F_tally, Be_tally])
            self.tallies.extend([flux_energy_tally, fertile_energy_tally, Li_energy_tally, F_energy_tally, Be_energy_tally])


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
        if self.breeder_name == 'PBHet' or self.breeder_name == 'PBHom':
            # Create plane sources on first wall plasma-facing surfaces
            # Split neutron count 50/50 between inner and outer sources
            inner_weight = 0.5
            outer_weight = 0.5
            
            # Inner first wall plane source (a1 x a1)
            source_inner = openmc.IndependentSource()
            # Use Box distribution: x fixed at plane, y and z uniform over the area
            a1_half = (self.inner_fw_size - 0.001) / 2.0
            source_inner.space = openmc.stats.Box(
                [self.inner_fw_pf_x - 0.01, -a1_half, -a1_half],  # Lower corner (small x thickness for plane)
                [self.inner_fw_pf_x - 0.001,  a1_half,  a1_half]   # Upper corner
            )
            source_inner.particle = 'neutron'
            source_inner.energy = openmc.stats.Discrete([14.0e6], [1.0])  # 14 MeV
            source_inner.angle = openmc.stats.Isotropic()
            source_inner.strength = inner_weight
            
            # Outer first wall plane source (a2 x a2)
            source_outer = openmc.IndependentSource()
            a2_half = self.outer_fw_size / 2.0
            source_outer.space = openmc.stats.Box(
                [self.outer_fw_pf_x + 0.001, -a2_half, -a2_half],  # Lower corner (small x thickness for plane)
                [self.outer_fw_pf_x + 0.01,  a2_half,  a2_half]   # Upper corner
            )
            source_outer.particle = 'neutron'
            source_outer.energy = openmc.stats.Discrete([14.0e6], [1.0])  # 14 MeV
            source_outer.angle = openmc.stats.Isotropic()
            source_outer.strength = outer_weight
            
            # Use list of sources for multiple sources
            sources = [source_inner, source_outer]
        else:
            source = openmc.IndependentSource()
            source.space = openmc.stats.CylindricalIndependent(
                                          r=openmc.stats.Discrete([self.R0], [1.0]),  # r = R0 major radius [cm]
                                        phi=openmc.stats.Uniform(0.0, 2*np.pi)   ,  # phi = 0 to 2pi
                                          z=openmc.stats.Discrete([0.0], [1.0])   ) # z   = 0
                                                          
            source.particle = 'neutron'
            source.energy   = openmc.stats.Discrete([14.0e6], [1.0])  # 14 MeV
            source.angle    = openmc.stats.Isotropic()
            sources = source

        # Create settings and assign source(s)
        self.settings = openmc.Settings()
        self.settings.source = sources

        """ Run type """
        self.settings.run_mode  = 'fixed source'
        self.settings.particles = self.n_particles
        self.settings.batches   = self.n_cycles


    def openmc(self):
        
        self.materials()
        self.geometry()
        self.settings()
        self.tallies()

        self.materials.cross_sections = set_xs_path()
        self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
        self.model.export_to_model_xml(self.path)
        self.model.run(cwd=self.path)


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
        
        print(f"{Colors.GREEN}Running stochastic volume calculation with {samples} samples...{Colors.END}")

        self.materials.cross_sections = set_xs_path()
        self.model_vol_calc = openmc.model.Model(self.geometry, self.materials, settings_vol_calc)
        self.model_vol_calc.calculate_volumes(cwd=self.path, export_model_xml=False, apply_volumes=False)


        # ---------------------------
        # Process volume calc results
        # ---------------------------
        vol_results = openmc.VolumeCalculation.from_hdf5(f"{self.path}/volume_1.h5")
        
        #print(f"{Colors.GREEN}Stochastic Volume Calculation Results:{Colors.END}")

        vol_dict = {}
        for k, v in vol_results.volumes.items():
            vol_dict[k] = (v.nominal_value/1e6, v.std_dev/1e6)
        # vol_dict['sum'] = ( sum(v.nominal_value/1e6 for v in vol_results.volumes.values()),
        #                     sum(v.std_dev/1e6 for v in vol_results.volumes.values())       )

        df = pd.DataFrame.from_dict(vol_dict, orient='index', columns=['volume_m3', 'stdev_m3'])
        df.reset_index(inplace=True)
        df.rename(columns={'index': 'cells'}, inplace=True)
        df.to_csv(f'{self.path}/volume_1.csv',index=False)

        print(f"{Colors.YELLOW}Comment.{Colors.END} CSV file 'volume_1.csv' printed to {self.path}")
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

        if self.breeder_name in ['ARC', 'ARCBall', 'FLiBe']:
            colors = {self.firstwall: (30, 27, 41),    # very dark gray
                      self.structure: (109, 110, 113), # gray
                      self.blanket: (129, 204, 185),   # teal
                      # Void regions will be white by default
                      }

        elif self.breeder_name in ['LL']:
            colors = {self.firstwall: (30, 27, 41),           # very dark gray
                      self.structure: (109, 110, 113),        # gray
                      self.coolant: (37, 150, 190), # cobalt blue 
                      self.blanket: (129, 204, 185),          # teal
                      self.divider: (109, 110, 113),          # gray
                      self.inner_manifold:(176, 123, 76),     # wood-color
                      # Void regions will be white by default
                      }

        elif self.breeder_name in ['PB']:
            colors = {self.firstwall: (30, 27, 41),           # very dark gray
                      self.structure: (109, 110, 113),        # gray
                      self.blanket: (129, 204, 185),          # teal
                      # Void regions will be white by default
                      }

        x_width = round(2*self.extent_r)
        z_width = round(2*self.extent_z)
        
        # XY toroidal slice
        xy = openmc.Plot()
        xy.filename = f"{self.breeder_name}_xy" # {self.path}/
        xy.basis = "xy"
        xy.width  = (x_width, x_width)
        xy.pixels = (10*x_width, 10*x_width)
        xy.color_by = "material"
        xy.colors = colors

        # XZ poloidal slice
        xz = openmc.Plot()
        xz.filename = f"{self.breeder_name}_xz" # {self.path}/
        xz.basis = "xz"
        xz.width  = (x_width, z_width)
        xz.pixels = (10*x_width, 10*z_width)
        xz.color_by = "material"
        xz.colors = colors

        plots = openmc.Plots([xy, xz])
        model_plot = openmc.model.Model(self.geometry, self.materials, settings_plot)
        model_plot.plots = plots 
        model_plot.plot_geometry(cwd=f"{self.path}")  # writes tokamak_rz.ppm and tokamak_xz.ppm

        for basename in [f"{self.breeder_name}_xy", f"{self.breeder_name}_xz"]:
            ppm_file = os.path.join(self.path, f"{basename}.ppm")
            png_file = os.path.join(self.path, f"{basename}.png")
            if os.path.exists(ppm_file):
                with Image.open(ppm_file) as im:
                    im.save(png_file)
                print(f"{Colors.YELLOW}Comment.{Colors.END} Plots '{self.breeder_name}_xy.png', '{self.breeder_name}_xz.png' printed to {self.path}")
            else:
                print(f"{Colors.YELLOW}Error.{Colors.END} OpenMC did not print '{self.breeder_name}_xy.ppm', '{self.breeder_name}_xz.ppm' to {self.path}!")


    def extract_tallies(self):

        """ Load tallies """
        sp_path = f'{self.path}/statepoint.{self.n_cycles}.h5'
        print(f"Loading statepoint: {sp_path}")
        
        try:
            sp = openmc.StatePoint(sp_path) 
        except Exception as e:
            print(f"\n{e}\n")
            sys.exit(f"can't read the sp")


        """ Convert tallies into usable forms """
        # Read tallies
        print(f"Reading tallies...")

        # For PBHet and PBHom: use material-based tallies
        if self.breeder_name == 'PBHet' or self.breeder_name == 'PBHom':
            # Read tallies - similar structure to cell-based
            flux = sp.get_tally(name='flux').get_pandas_dataframe()
            flux_spec = sp.get_tally(name='flux spectrum').get_pandas_dataframe()
            
            # Read fertile element tallies (U or Th) - using similar naming
            if self.fertile_element == 'U':
                fertile = sp.get_tally(name='U n-gamma').get_pandas_dataframe()
                fertile_spec = sp.get_tally(name='U n-gamma spectrum').get_pandas_dataframe()
            elif self.fertile_element == 'Th':
                fertile = sp.get_tally(name='Th n-gamma').get_pandas_dataframe()
                fertile_spec = sp.get_tally(name='Th n-gamma spectrum').get_pandas_dataframe()
            
            # Read Li tallies - combine into single dataframe like cell-based
            Li6_nt_df = sp.get_tally(name='Li6 n-Xt').get_pandas_dataframe()
            Li7_nt_df = sp.get_tally(name='Li7 n-Xt').get_pandas_dataframe()
            # Combine Li6 and Li7 into single dataframe (like cell-based 'Li rxn rates')
            Li = pd.concat([Li6_nt_df, Li7_nt_df], ignore_index=True)
            
            Li6_nt_spec = sp.get_tally(name='Li6 n-Xt spectrum').get_pandas_dataframe()
            Li7_nt_spec = sp.get_tally(name='Li7 n-Xt spectrum').get_pandas_dataframe()
            # Combine spectra
            Li_spec = pd.concat([Li6_nt_spec, Li7_nt_spec], ignore_index=True)

            # Add new column for energy bin midpoint (for plotting)
            for df in [flux_spec, fertile_spec, Li_spec]:
                df['energy mid [eV]'] = (df['energy low [eV]'] + df['energy high [eV]'])/ 2

            """ Reaction rates for every fertile loading and for every energy bin """
            # Use 'material' column instead of 'cell', but same structure
            fertile_ng_Ebin_df = fertile_spec[(fertile_spec['nuclide'] == self.fertile_isotope)][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'material', 'mean', 'std. dev.']]
            Li6_nt_Ebin_df = Li_spec[(Li_spec['nuclide'] == 'Li6')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'material', 'mean', 'std. dev.']]
            Li7_nt_Ebin_df = Li_spec[(Li_spec['nuclide'] == 'Li7')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'material', 'mean', 'std. dev.']]

            # Lithium reaction rates - filter by nuclide like cell-based
            Li6_nt = Li[(Li['nuclide']=='Li6')][['material','nuclide','mean','std. dev.']]
            Li7_nt = Li[(Li['nuclide']=='Li7')][['material','nuclide','mean','std. dev.']]
            
            Li6_nt_list = Li6_nt['mean'].tolist()
            Li7_nt_list = Li7_nt['mean'].tolist()
            Li6_nt_err_list = Li6_nt['std. dev.'].tolist()
            Li7_nt_err_list = Li7_nt['std. dev.'].tolist()
            
            tbr_list = [x + y for x, y in zip(Li6_nt_list, Li7_nt_list)]
            tbr_err_list = [np.sqrt(x**2 + y**2) for x, y in zip(Li6_nt_err_list, Li7_nt_err_list)]

            # Fertile element reaction rates - same structure as cell-based
            fertile_ng = fertile[(fertile['nuclide']==self.fertile_isotope)][['material','nuclide','mean','std. dev.']]
            fertile_ng_list = fertile_ng['mean'].tolist()
            fertile_ng_err_list = fertile_ng['std. dev.'].tolist()

            # Fissile production
            fissile_per_yr_list = []
            for ng_per_srcn in fertile_ng_list:
                if self.fertile_element == 'U':
                    fissile_per_yr_list.append(ng_per_srcn * NPS_FUS * SEC_PER_YR * AMU_PU239 / AVO / 1e3) # kg/yr
                elif self.fertile_element == 'Th':
                    fissile_per_yr_list.append(ng_per_srcn * NPS_FUS * SEC_PER_YR * AMU_U233 / AVO / 1e3) # kg/yr

            """Create list of material IDs - similar to cell IDs"""
            all_material_ids = set(flux['material'].unique().tolist())
            all_material_ids.update(Li['material'].unique().tolist())
            all_material_ids.update(fertile['material'].unique().tolist())
            material_ids = [str(x) for x in sorted(all_material_ids)]
            
            print("unique materials in flux   :", len(flux['material'].unique()),   " rows:", len(flux))
            print("unique materials in Li     :", len(Li['material'].unique()),   " rows:", len(Li))
            print("unique materials in fertile:", len(fertile['material'].unique()), " rows:", len(fertile))
            print("total unique materials     :", len(material_ids))

            # Create dictionaries keyed by material ID - same structure as cell-based
            flux_dict = dict(zip(flux['material'], flux['mean']))
            Li6_nt_dict = dict(zip(Li6_nt['material'], Li6_nt['mean']))
            Li7_nt_dict = dict(zip(Li7_nt['material'], Li7_nt['mean']))
            Li6_nt_err_dict = dict(zip(Li6_nt['material'], Li6_nt['std. dev.']))
            Li7_nt_err_dict = dict(zip(Li7_nt['material'], Li7_nt['std. dev.']))
            fertile_ng_dict = dict(zip(fertile_ng['material'], fertile_ng['mean']))
            fertile_ng_err_dict = dict(zip(fertile_ng['material'], fertile_ng['std. dev.']))

            # Build lists aligned by material_id, using 0.0 as default for missing values
            flux_list = [flux_dict.get(mid, 0.0) for mid in material_ids]
            Li6_nt_list = [Li6_nt_dict.get(mid, 0.0) for mid in material_ids]
            Li7_nt_list = [Li7_nt_dict.get(mid, 0.0) for mid in material_ids]
            Li6_nt_err_list = [Li6_nt_err_dict.get(mid, 0.0) for mid in material_ids]
            Li7_nt_err_list = [Li7_nt_err_dict.get(mid, 0.0) for mid in material_ids]
            tbr_list = [x + y for x, y in zip(Li6_nt_list, Li7_nt_list)]
            tbr_err_list = [np.sqrt(x**2 + y**2) for x, y in zip(Li6_nt_err_list, Li7_nt_err_list)]
            fertile_ng_list = [fertile_ng_dict.get(mid, 0.0) for mid in material_ids]
            fertile_ng_err_list = [fertile_ng_err_dict.get(mid, 0.0) for mid in material_ids]

            # Fissile production
            fissile_per_yr_list = []
            for ng_per_srcn in fertile_ng_list:
                if self.fertile_element == 'U':
                    fissile_per_yr_list.append(ng_per_srcn * NPS_FUS * SEC_PER_YR * AMU_PU239 / AVO / 1e3) # kg/yr
                elif self.fertile_element == 'Th':
                    fissile_per_yr_list.append(ng_per_srcn * NPS_FUS * SEC_PER_YR * AMU_U233 / AVO / 1e3) # kg/yr

            print("\n-- list lengths (should all match) --")
            print("material_ids       :", len(material_ids))
            print("flux_list          :", len(flux_list))
            print("Li6_nt_list        :", len(Li6_nt_list))
            print("Li7_nt_list        :", len(Li7_nt_list))
            print("tbr_list           :", len(tbr_list))
            print("fertile_ng_list    :", len(fertile_ng_list))
            print("fissile_per_yr_list:", len(fissile_per_yr_list))
            print("=== END DEBUG ===\n")

            # Create DataFrame - use 'cell' column name for compatibility, but contains material IDs
            df = pd.DataFrame({'cell': material_ids,
                           'flux': flux_list,
                           'Li6(n,t)': Li6_nt_list,
                           'Li6(n,t)_stdev': Li6_nt_err_list,
                           'Li7(n,t)': Li7_nt_list,
                           'Li7(n,t)_stdev': Li7_nt_err_list,
                           'tbr': tbr_list,
                           'tbr_stdev': tbr_err_list,
                           f'{self.fertile_isotope}(n,g)': fertile_ng_list,
                           f'{self.fertile_isotope}(n,g)_stdev': fertile_ng_err_list,
                           f'{self.fissile_isotope}_kg/yr': fissile_per_yr_list})

            # Sum to get totals - same as cell-based
            totals = df.sum(numeric_only=True)
            totals['cell'] = 'total'
            df = pd.concat([df, pd.DataFrame([totals])], ignore_index=True)
            
            df.to_csv(f'./{self.path}/tallies_summary.csv', index=False)
            
            # Energy-binned data - group by energy and sum over materials
            fertile_ng_Ebin_df = fertile_ng_Ebin_df.groupby('energy mid [eV]', as_index=False).agg({
                'energy low [eV]': 'first',
                'energy high [eV]': 'first',
                'mean': 'sum',
                'std. dev.': lambda x: np.sqrt((x**2).sum())
            })
            fertile_ng_Ebin_df.to_csv(f'./{self.path}/{self.fertile_isotope}_n-gamma_Ebins.csv', index=False)
            
        else:
            # Original cell-based tallies for other reactor types
            flux      = sp.get_tally(name='flux').get_pandas_dataframe()
            flux_spec = sp.get_tally(name='flux spectrum').get_pandas_dataframe()
            Li        = sp.get_tally(name='Li rxn rates').get_pandas_dataframe()
            Li_spec   = sp.get_tally(name='Li rxn rates spectrum').get_pandas_dataframe()
            fertile      = sp.get_tally(name=f'{self.fertile_element} rxn rates').get_pandas_dataframe()
            fertile_spec = sp.get_tally(name=f'{self.fertile_element} rxn rates spectrum').get_pandas_dataframe()

            # Add new column for energy bin midpoint (for plotting)
            for df in [flux_spec, fertile_spec, Li_spec]:
                df['energy mid [eV]'] = (df['energy low [eV]'] + df['energy high [eV]'])/ 2

            """ Reaction rates for every fertile loading and for every energy bin """
            U238_ng_Ebin_df  = fertile_spec[(fertile_spec['nuclide'] == self.fertile_isotope) & (fertile_spec['score'] == '(n,gamma)')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]
            U238_fis_Ebin_df = fertile_spec[(fertile_spec['nuclide'] == self.fertile_isotope) & (fertile_spec['score'] == 'fission')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]
            Li6_nt_Ebin_df   = Li_spec[(Li_spec['nuclide'] == 'Li6') & (Li_spec['score'] == '(n,Xt)')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]
            Li7_nt_Ebin_df   = Li_spec[(Li_spec['nuclide'] == 'Li7') & (Li_spec['score'] == '(n,Xt)')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]

            # Lithium reaction rates 
            Li6_nt = Li[(Li['nuclide']=='Li6') & (Li['score']=='(n,Xt)')][['cell','score','mean','std. dev.']]
            Li7_nt = Li[(Li['nuclide']=='Li7') & (Li['score']=='(n,Xt)')][['cell','score','mean','std. dev.']] 
            
            Li6_nt_list     = Li6_nt['mean'].tolist()
            Li7_nt_list     = Li7_nt['mean'].tolist() 
            Li6_nt_err_list = Li6_nt['std. dev.'].tolist()
            Li7_nt_err_list = Li7_nt['std. dev.'].tolist()
            
            tbr_list = [x + y for x, y in zip(Li6_nt_list, Li7_nt_list)]
            tbr_err_list = [x + y for x, y in zip(Li6_nt_err_list, Li7_nt_err_list)]

            # U-238 or Th-232 reaction rates for each mtu loading summed over all energies
            U238              = fertile[(fertile['nuclide']==self.fertile_isotope)][['cell','score','mean','std. dev.']]
            U238_fis_list     = U238[U238['score'] == 'fission'][['mean']]['mean'].tolist()
            U238_fis_err_list = U238[U238['score'] == 'fission'][['std. dev.']]['std. dev.'].tolist()
            U238_ng_list      = U238[U238['score'] == '(n,gamma)'][['mean']]['mean'].tolist()
            U238_ng_err_list  = U238[U238['score'] == '(n,gamma)'][['std. dev.']]['std. dev.'].tolist()

            # Plutonium kg per year
            fissile_per_yr_list = []
            for Pu_per_srcn in U238_ng_list:
                if self.fertile_element == 'U':
                    fissile_per_yr_list.append( Pu_per_srcn * NPS_FUS * SEC_PER_YR * AMU_PU239 / AVO / 1e3) # kg/yr
                elif self.fertile_element == 'Th':
                    fissile_per_yr_list.append( Pu_per_srcn * NPS_FUS * SEC_PER_YR * AMU_U233 / AVO / 1e3) # kg/yr

            """Create list of cell IDs - get union of all cells from all tallies"""
            all_cell_ids = set(flux['cell'].unique().tolist())
            all_cell_ids.update(Li['cell'].unique().tolist())
            all_cell_ids.update(fertile['cell'].unique().tolist())
            cell_ids = [str(x) for x in sorted(all_cell_ids)]
            
            print("unique cells in flux   :", len(flux['cell'].unique()),   " rows:", len(flux))
            print("unique cells in Li    :", len(Li['cell'].unique()),   " rows:", len(Li))
            print("unique cells in fertile:", len(fertile['cell'].unique()), " rows:", len(fertile))
            print("total unique cells     :", len(cell_ids))

            # Create dictionaries keyed by cell ID to ensure proper alignment
            flux_dict = dict(zip(flux['cell'], flux['mean']))
            Li6_nt_dict = dict(zip(Li6_nt['cell'], Li6_nt['mean']))
            Li7_nt_dict = dict(zip(Li7_nt['cell'], Li7_nt['mean']))
            Li6_nt_err_dict = dict(zip(Li6_nt['cell'], Li6_nt['std. dev.']))
            Li7_nt_err_dict = dict(zip(Li7_nt['cell'], Li7_nt['std. dev.']))
            U238_fis_dict = dict(zip(U238[U238['score'] == 'fission']['cell'], 
                                      U238[U238['score'] == 'fission']['mean']))
            U238_fis_err_dict = dict(zip(U238[U238['score'] == 'fission']['cell'], 
                                          U238[U238['score'] == 'fission']['std. dev.']))
            U238_ng_dict = dict(zip(U238[U238['score'] == '(n,gamma)']['cell'], 
                                    U238[U238['score'] == '(n,gamma)']['mean']))
            U238_ng_err_dict = dict(zip(U238[U238['score'] == '(n,gamma)']['cell'], 
                                         U238[U238['score'] == '(n,gamma)']['std. dev.']))

            # Build lists aligned by cell_id, using 0.0 as default for missing values
            flux_list = [flux_dict.get(cid, 0.0) for cid in cell_ids]
            Li6_nt_list = [Li6_nt_dict.get(cid, 0.0) for cid in cell_ids]
            Li7_nt_list = [Li7_nt_dict.get(cid, 0.0) for cid in cell_ids]
            Li6_nt_err_list = [Li6_nt_err_dict.get(cid, 0.0) for cid in cell_ids]
            Li7_nt_err_list = [Li7_nt_err_dict.get(cid, 0.0) for cid in cell_ids]
            tbr_list = [x + y for x, y in zip(Li6_nt_list, Li7_nt_list)]
            tbr_err_list = [x + y for x, y in zip(Li6_nt_err_list, Li7_nt_err_list)]
            U238_fis_list = [U238_fis_dict.get(cid, 0.0) for cid in cell_ids]
            U238_fis_err_list = [U238_fis_err_dict.get(cid, 0.0) for cid in cell_ids]
            U238_ng_list = [U238_ng_dict.get(cid, 0.0) for cid in cell_ids]
            U238_ng_err_list = [U238_ng_err_dict.get(cid, 0.0) for cid in cell_ids]

            # Plutonium kg per year
            fissile_per_yr_list = []
            for Pu_per_srcn in U238_ng_list:
                if self.fertile_element == 'U':
                    fissile_per_yr_list.append( Pu_per_srcn * NPS_FUS * SEC_PER_YR * AMU_PU239 / AVO / 1e3) # kg/yr
                elif self.fertile_element == 'Th':
                    fissile_per_yr_list.append( Pu_per_srcn * NPS_FUS * SEC_PER_YR * AMU_U233 / AVO / 1e3) # kg/yr

            print("\n-- list lengths (should all match) --")
            print("cell_ids           :", len(cell_ids))
            print("Li6_nt_list        :", len(Li6_nt_list))
            print("Li7_nt_list        :", len(Li7_nt_list))
            print("tbr_list           :", len(tbr_list))
            print("U238_ng_list       :", len(U238_ng_list))
            print("fissile_per_yr_list:", len(fissile_per_yr_list))
            print("=== END DEBUG ===\n")

            df = pd.DataFrame({'cell': cell_ids,
                           'flux': flux_list,
                           'Li6(n,t)'          : Li6_nt_list,
                           'Li6(n,t)_stdev'    : Li6_nt_err_list,
                           'Li7(n,t)'          : Li7_nt_list,
                           'Li7(n,t)_stdev'    : Li7_nt_err_list,
                           'tbr'               : tbr_list,
                           'tbr_stdev'         : tbr_err_list,
                           f'{self.fertile_isotope}(n,fis)'       : U238_fis_list,
                           f'{self.fertile_isotope}(n,fis)_stdev' : U238_fis_err_list,
                           f'{self.fertile_isotope}(n,g)'         : U238_ng_list,
                           f'{self.fertile_isotope}(n,g)_stdev'   : U238_ng_err_list,
                           f'{self.fissile_isotope}_kg/yr'        : fissile_per_yr_list,   })

            totals = df.sum(numeric_only=True)
            totals['cell'] = 'total'
            df = pd.concat([df, pd.DataFrame([totals])], ignore_index=True)
            
            df.to_csv(f'./{self.path}/tallies_summary.csv', index=False)
            
            U238_ng_Ebin_df.to_csv(f'./{self.path}/{self.fertile_isotope}_n-gamma_Ebins.csv', index=False)



            
