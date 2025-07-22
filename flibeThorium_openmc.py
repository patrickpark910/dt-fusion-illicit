import openmc
import os, sys
import numpy as np

# Import helper functions
from Python.parameters import *
from Python.utilities import *


def main():
    print("\n\n")
    print("="*42)
    print(f"Running FLiBe Thorium OpenMC model for Li-6 enrichment: {7.5} wt%")

    current_run = FLIBE_Th(enrich_li=7.5, temp_k=900)

    print(f"Check if '{current_run.path}' exists: {os.path.isdir(current_run.path)}")

    if os.path.isdir(current_run.path):
        print(f"Warning. Directory {current_run.path} already exists, so running OpenMC will fail. Skipping...")
        return
    else:
        current_run.set_xs_path()
        current_run.run_openmc()


class FLIBE_Th:

    def __init__(self, u_list=MASS_U_LIST_FLIBE, enrich_li=7.5, temp_k=900):

        self.lie = enrich_li
        self.temp = temp_k
        self.u_list = u_list
        self.name = f"FLiBe_Th_Li{self.lie:04.1f}_7_22"
        self.path = f"./OpenMC/{self.name}/"

        flibe = openmc.Material()
        flibe.add_elements_from_formula('F4Li2Be', 'ao', enrichment_target='Li6', enrichment_type='wo', enrichment=self.lie)
        flibe.set_density('g/cm3', DENSITY_FLIBE) 

        # Uranium tetrafluoride -- assumed to dissolve in FLiBe
        # Modeled Thorium due to concerns about Pu florination for reprocessing, it has not been proven that this method can work experimentally
        """ Thorium tetrafluoride """
        thf4 = openmc.Material()
        thf4.add_elements_from_formula('ThF4')
        thf4.set_density('g/cm3', 6.3)

        # Calculate volume ratios of UF4 and FLiBe, ensure they add up to 1
        mix_list = []
        for mtu in self.u_list:
            vf_flibe, vf_thf4 = calc_mix_vol_fracs_th(mtu, volume=VOL_CC, density_flibe=DENSITY_FLIBE, displace=True) # './helper/utilities.py'
            mix = openmc.Material.mix_materials([flibe, thf4], [vf_flibe, vf_thf4], 'vo') # fractions in 'mix_materials' MUST add up to 1
            mix.name, mix.volume, mix.temperature = f"FLiBe {mtu:.1f} MTU at {self.temp} K", VOL_CC, self.temp
            mix_list.append(mix)

        self.materials = openmc.Materials(mix_list)


        """ GEOMETRY """
        cells, pitch, half_box = [], 120, 50 # +/- 50 cm bounds
        box_centers = [pitch * i for i in range(len(self.u_list))] # used in Sources

        for i, material in enumerate(mix_list):
            x_min = openmc.XPlane(x0= -half_box + box_centers[i], boundary_type='reflective')
            x_max = openmc.XPlane(x0=  half_box + box_centers[i], boundary_type='reflective')
            y_min, y_max = openmc.YPlane(-50, boundary_type='reflective'), openmc.YPlane( 50, boundary_type='reflective')
            z_min, z_max = openmc.ZPlane(-50, boundary_type='reflective'), openmc.ZPlane( 50, boundary_type='reflective')
            region = +x_min & -x_max & +y_min & -y_max & +z_min & -z_max
            
            cell = openmc.Cell(fill=material, region=region)
            cell.name = f"mix-{i+1}"
            cells.append(cell)

        root_univ = openmc.Universe(cells=cells) # Create root universe with all material cells
        root_univ.plot(width=(pitch * len(self.u_list), 110), origin=(pitch * (len(self.u_list) - 1) / 2, 0.0, 0.0)) # Visualize
        self.geometry = openmc.Geometry(root_univ) # Set geometry


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
        Th_tally = openmc.Tally(name='thorium rxn rates')
        Th_tally.scores = ['(n,gamma)','fission', 'elastic'] # specific names required
        Th_tally.nuclides = ['Th232']
        Th_tally.filters = filters

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

        self.tallies.extend([flux_tally, Th_tally, Li_tally, F_tally, Be_tally])
        # tallies.export_to_xml("./xml/tallies.xml") --don't need bc 'model.export_to_model_xml' below --ppark 2025-06-28


        
        """First Wall Effects"""
        sp = openmc.StatePoint(f'./OpenMC/FirstWall/statepoint.100.h5')  
        out_tally = sp.get_tally(name='outgoing_spectrum')

        energy_bins = out_tally.filters[1].bins
        energy_bins = np.array(energy_bins)  # shape (N, 2)

        energy_midpoints = 0.5 * (energy_bins[:, 0] + energy_bins[:, 1])

        current_spectrum = out_tally.get_values(scores=['current']).flatten()

        total_current = current_spectrum.sum()
        probabilities = current_spectrum / total_current

        """ SETTINGS """
        self.settings = openmc.Settings()

        """ Source
        Isotropic 14.07 MeV point source at center of each cube
        """
        energies = energy_midpoints.tolist()
        weights = probabilities.tolist()
        source = []
        for p in box_centers:
            src = openmc.IndependentSource()
            src.space  = openmc.stats.Point((p,0,0))
            src.angle  = openmc.stats.Isotropic()
            src.energy = openmc.stats.Discrete(energies, weights)
            source.append(src)
        self.settings.source = source

        """ Run type """
        self.settings.run_mode = 'fixed source'
        self.settings.particles = len(MASS_U_LIST) * int(1e6)  #  
        self.settings.batches = 100

    def run_openmc(self):
        self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
        self.model.export_to_model_xml(f"./OpenMC/{self.name}/")  # _Li6-20wt
        self.model.run(cwd=f"./OpenMC/{self.name}/")  # _Li6-20wt


    def set_xs_path(self):
        """ 
        Temporary solution for finding xs files between WSL and Ubuntu on Computing Cluster without editing PATH --ppark 2025-06-28 
        """
        xs_path_ubuntu = '/opt/openmc_data/endfb-viii.0-hdf5/cross_sections.xml'
        xs_path_wsl   = '/mnt/c/openmc/data/endfb-viii.0-hdf5/cross_sections.xml'
        if os.path.isfile(xs_path_ubuntu):
            self.materials.cross_sections = xs_path_ubuntu # use this on Zotacs --ppark
        elif os.path.isfile(xs_path_wsl):
            self.materials.cross_sections = xs_path_wsl
            print(self.materials.cross_sections )
        else:
            sys.exit(f"Error finding cross section XML!")


if __name__ == "__main__":
    main()



