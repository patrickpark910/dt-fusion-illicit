import openmc
import os, sys
import numpy as np

# Import helper functions
sys.path.insert(0, f"{os.getcwd()}/Python/")
from parameters import *
from utilities import *


def main():

    print("\n\n")
    print("="*42)
    print(f"Running First Wall")
        
    current_run = FirstWall()

    if current_run.run:
        if os.path.isdir(f"./OpenMC/{current_run.name}/"):
            print(f"Warning. Directory {f"./OpenMC/{current_run.name}/"} already exists, so running OpenMC will fail. Skipping...")

        else:
            current_run.run_openmc()


class FirstWall:

    def __init__(self, run_openmc=False):
        self.name = 'FirstWall_50'
        self.run  = run_openmc

        """ MATERIALS """
        """ TUNGSTEN from ball"""
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
   
        self.materials = openmc.Materials([tungsten,vcrti])
        # self.materials.cross_sections = set_xs_path()

        """GEOMETRY"""
        '''first wall slab thickness approximations from balls simplified model
        but we are ignoring the effects of imbeded structural components within the flibe mixture
        such as the channel structures, RF heating, and vacuum systems.'''
        plasma_thickness = 100
        tung_thickness = 0.3 #cm
        vcrti_thickness = 4

        r0 = openmc.Sphere(r=plasma_thickness)
        r1 = openmc.Sphere(r=tung_thickness + plasma_thickness)  # inner radius 10 cm
        r2 = openmc.Sphere(r=vcrti_thickness+tung_thickness+plasma_thickness, boundary_type='vacuum')

        plasma_region = -r0
        tung_region = +r0 & -r1
        vcrti_region = +r1 & -r2

        plasma_cell = openmc.Cell(region=plasma_region)
        tung_cell = openmc.Cell(fill=tungsten, region=tung_region)
        vcrti_cell = openmc.Cell(fill=vcrti, region=vcrti_region)

        self.geometry = openmc.Geometry([plasma_cell, tung_cell, vcrti_cell])


        """ TALLIES """
        self.tallies = openmc.Tallies() # initialize

        # Filters
        E_bin_edges = logspace_per_decade(1e-5, 20e6, 100) # './helpers/utilities.py'
        energy_filter = openmc.EnergyFilter(E_bin_edges)

        surface_filter = openmc.SurfaceFilter(r2)

        # Current tally to measure outgoing neutron spectrum
        out_tally = openmc.Tally(name='outgoing_spectrum')
        out_tally.filters = [surface_filter, energy_filter]
        out_tally.scores = ['current']

        self.tallies = openmc.Tallies([out_tally])
        

        """ SETTINGS """
        self.settings = openmc.Settings()

        """ Source
        Isotropic 14.07 MeV 
        """
        """ Run type """
        self.settings = openmc.Settings()
        self.settings.run_mode = 'fixed source'
        self.settings.particles = int(1e7)
        self.settings.batches = 100

        source = openmc.IndependentSource()
        source.space = openmc.stats.Point((0,0,0)) 
        source.angle = openmc.stats.Isotropic()     
        source.energy = openmc.stats.Discrete([14.07e6], [1.0])  

        self.settings.source = source

    def run_openmc(self):
        self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
        self.model.export_to_model_xml(f"./OpenMC/{self.name}/")
        self.model.run(cwd=f"./OpenMC/{self.name}/")





if __name__ == "__main__":
    main()
