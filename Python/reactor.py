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
        self.n_particles = int(1e5)
        self.n_cycles    = 10

        # All class templates must define these variables:
        self.temp_k          = None # TEMP_K
        self.breeder_name    = None # 'ARC'
        self.breeder_density = None # DENSITY_FLIBE # g/cm^3
        self.breeder_enrich  = None # ENRICH_FLIBE  # wt%
        self.breeder_volume  = None # ARC_BR_VOL    # m^3
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

        # Filters
        cell_filter = openmc.CellFilter([cell for cell in self.cells if cell.fill is not None])
        energy_filter = openmc.EnergyFilter(logspace_per_decade(1e-5, 20e6, 100)) # './helpers/utilities.py'
        

        # ---------------------------------------------------------------
        # Neutron wall loading and heating tallies 
        # ---------------------------------------------------------------
        """
        nwl_surface_filter = openmc.SurfaceFilter([self.surface_vc.id]) # Surface filter for the plasma-facing surface (inner surface of first wall)
        nwl_tally        = self.make_tally('neutron wall loading',          ['current'], filters=[nwl_surface_filter])
        nwl_energy_tally = self.make_tally('neutron wall loading spectrum', ['current'], filters=[nwl_surface_filter, energy_filter])
        # NWL = neutron current through surface, Using 'current' gives net current, 'partial-current' gives directional
        self.tallies.extend([nwl_tally, nwl_energy_tally])
        # I omit these bc our surfaces are openmc.model.Polygon's which SurfaceFilter can NOT read --ppark 2025-09-21
        """

        # -----------------------------------------------------------------------
        # Heating tallies (heating-local estimates n+gamma heating using only n)
        # -----------------------------------------------------------------------

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


    def openmc(self):
        
        self.materials()
        self.geometry()
        self.settings()
        self.tallies()

        self.materials.cross_sections = set_xs_path()
        self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
        self.model.export_to_model_xml(self.path)
        self.model.run(cwd=self.path)


    def volumes(self, samples=int(1e8)):
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

        flux      = sp.get_tally(name='flux').get_pandas_dataframe()
        flux_spec = sp.get_tally(name='flux spectrum').get_pandas_dataframe()
        Li        = sp.get_tally(name='Li rxn rates').get_pandas_dataframe()
        Li_spec   = sp.get_tally(name='Li rxn rates spectrum').get_pandas_dataframe()
        fertile      = sp.get_tally(name=f'{self.fertile_element} rxn rates').get_pandas_dataframe()
        fertile_spec = sp.get_tally(name=f'{self.fertile_element} rxn rates spectrum').get_pandas_dataframe()

        # Add new column for energy bin midpoint (for plotting)
        for df in [flux_spec, fertile_spec, Li_spec]:
            df['energy mid [eV]'] = (df['energy low [eV]'] + df['energy high [eV]'])/ 2

        """
        flux, U, and Li dataframes should have column headers: 
            (index), cell, nuclide, score, mean, std. dev.
        flux_spec, U_spec, and Li_spec should have column headers: 
            (index), cell, energy low [eV], energy high [eV], nuclide, score, mean, std. dev., energy mid [eV]
        """

        """ Reaction rates for every fertile loading and for every energy bin """
        U238_ng_Ebin_df  = fertile_spec[(fertile_spec['nuclide'] == self.fertile_isotope) & (fertile_spec['score'] == '(n,gamma)')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        U238_fis_Ebin_df = fertile_spec[(fertile_spec['nuclide'] == self.fertile_isotope) & (fertile_spec['score'] == 'fission')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        Li6_nt_Ebin_df   = Li_spec[(Li_spec['nuclide'] == 'Li6') & (Li_spec['score'] == '(n,Xt)')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]
        Li7_nt_Ebin_df   = Li_spec[(Li_spec['nuclide'] == 'Li7') & (Li_spec['score'] == '(n,Xt)')][['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'cell', 'mean', 'std. dev.']]


        # Flux
        flux_list = flux['mean'].tolist()

        # Lithium reaction rates 
        Li6_nt = Li[(Li['nuclide']=='Li6') & (Li['score']=='(n,Xt)')][['cell','score','mean','std. dev.']] # Li6_rr.get_pandas_dataframe()
        Li7_nt = Li[(Li['nuclide']=='Li7') & (Li['score']=='(n,Xt)')][['cell','score','mean','std. dev.']] 
        
        Li6_nt_list     = Li6_nt['mean'].tolist()
        Li7_nt_list     = Li7_nt['mean'].tolist() 
        Li6_nt_err_list = Li6_nt['std. dev.'].tolist()
        Li7_nt_err_list = Li7_nt['std. dev.'].tolist()
        
        tbr_list = [x + y for x, y in zip(Li6_nt_list, Li7_nt_list)] # will look like: [0.0, 0.0, 0.094, 0.0, 1.074, 0.0]
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


        """Create list of cell IDs randomly assigned by OpenMC to match mtu loading to cell ID"""
        cell_ids = [str(x) for x in flux['cell'].unique().tolist()] 
        # turn into str to add 'total' at end -- otherwise gets pandas error for type mismatch -- ppark 2025-10-07

        df = pd.DataFrame({'cell': cell_ids,
                           'flux': flux_list,
                           'Li6(n,t)'          : Li6_nt_list,
                           'Li6(n,t)_stdev'    : Li6_nt_err_list,
                           'Li7(n,t)'          : Li7_nt_list,
                           'Li7(n,t)_stdev'    : Li7_nt_err_list,
                           'tbr'               : tbr_list,
                           'tbr_stdev'         : tbr_err_list,
                           f'{self.fertile_isotope}(n,fis)'       : Li7_nt_list,
                           f'{self.fertile_isotope}(n,fis)_stdev' : U238_fis_list,
                           f'{self.fertile_isotope}(n,g)'         : U238_ng_list,
                           f'{self.fertile_isotope}(n,g)_stdev'   : U238_ng_err_list,
                           f'{self.fissile_isotope}_kg/yr'        : fissile_per_yr_list,   })

        totals = df.sum(numeric_only=True)
        totals['cell'] = 'total'
        df = pd.concat([df, pd.DataFrame([totals])], ignore_index=True)
        
        df.to_csv(f'./{self.path}/tallies_summary.csv', index=False)
        
        U238_ng_Ebin_df.to_csv(f'./{self.path}/{self.fertile_isotope}_n-gamma_Ebins.csv', index=False)



            
