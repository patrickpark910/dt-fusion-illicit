import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from scipy.interpolate import Akima1DInterpolator
from scipy.interpolate import PchipInterpolator

# Import helper functions
from Python.reactor import *
from Python.parameters import *
from Python.utilities import *


class Plot(Reactor):

    def __init__(self, reactor: Reactor, *, save=False, show=False):
        # Adopt the already-built Reactor's state without re-running Reactor.__init__
        self.__dict__.update(reactor.__dict__)
        # Plot-specific flags (your plotting methods already expect these)
        self.save = save
        self.show = show


    @timer
    def volumes(self, samples=int(1e8)):
        """
        Calculate volumes of all cells using OpenMC stochastic volume calculation.
        """
        
        # Get all cells that have materials (exclude voids)
        cells_to_calc = self.cells # [cell for cell in self.cells if cell.fill == self.blanket] #  is not None]
        
        # Create volume calculation settings
        vol_calc = openmc.VolumeCalculation(cells_to_calc, samples,
                                            [-1*self.extent_r, -1*self.extent_r, -1*self.extent_z], # lower left of bounding box
                                            [self.extent_r, self.extent_r, self.extent_z]) # upper right of bounding box
        
        # vol_calc.set_trigger(1e-02, 'std_dev')

        settings_vol_calc = openmc.Settings()
        settings_vol_calc.volume_calculations = [vol_calc]
        
        print(f"{Colors.BLUE}Running stochastic volume calculation with {samples} samples...{Colors.END}")

        self.materials.cross_sections = set_xs_path()
        self.model_vol_calc = openmc.model.Model(self.geometry, self.materials, settings_vol_calc)
        self.model_vol_calc.calculate_volumes(cwd=self.path, export_model_xml=False, apply_volumes=False)
        print('did openmc run')

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


    @timer
    def plot(self):

        Image.MAX_IMAGE_PIXELS = None # suppreses DecompressionBombWarning lmao

        settings_plot = openmc.Settings()
        settings_plot.run_mode = 'plot'

        if self.breeder_name in ['ARC', 'FLiBe']:
            colors = {self.firstwall: (30, 27, 41),    # very dark gray
                      self.structure: (109, 110, 113), # gray
                      self.blanket: (129, 204, 185),   # teal
                      # Void regions will be white by default
                      }

        elif self.breeder_name in ['LL']:
            colors = {self.firstwall: (30, 27, 41),           # very dark gray
                      self.structure: (109, 110, 113),        # gray
                      self.firstwall_cooling: (37, 150, 190), # cobalt blue 
                      self.blanket: (129, 204, 185),          # teal
                      self.divider: (109, 110, 113),          # gray
                      self.inner_manifold:(176, 123, 76),     # wood-color
                      # Void regions will be white by default
                      }
        elif self.breeder_name in ['PB']:
            colors = {self.firstwall: (30, 27, 41),           # very dark gray
                      self.eurofer:   (109, 110, 113),        # gray
                      self.blanket: (129, 204, 185),          # teal
                      }

        x_width = round(2*self.extent_r)
        z_width = round(2*self.extent_z)
        
        # XY toroidal slice
        xy = openmc.Plot()
        xy.filename = f"{self.breeder_name}_xy" # {self.path}/
        xy.basis = "xy"
        xy.width  = (x_width, x_width)
        xy.pixels = (8*x_width, 8*x_width)
        xy.color_by = "material"
        xy.colors = colors

        # XZ poloidal slice
        xz = openmc.Plot()
        xz.filename = f"{self.breeder_name}_xz" # {self.path}/
        xz.basis = "xz"
        xz.width  = (x_width, z_width)
        xz.pixels = (8*x_width, 8*z_width)
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