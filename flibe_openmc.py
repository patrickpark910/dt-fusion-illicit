import openmc
import sys
import numpy as np

# Import helper functions
from Python.parameters import *
from Python.utilities import *


def main():
    for e in [7.5]:
        print("\n\n")
        print("="*42)
        print(f"Running FLiBe OpenMC model for Li-6 enrichment: {e} wt%")

        current_run = FLIBE(enrich_li=e, temp_k=900, run_openmc=True)

        print(f"Check if '{current_run.path}' exists: {os.path.isdir(current_run.path)}")

        if current_run.run:
            if os.path.isdir(current_run.path):
                print(f"Warning. Directory {current_run.path} already exists, so running OpenMC will fail. Skipping...")
                continue
            else:
                current_run.run_openmc()


class FLIBE:

    def __init__(self, u_list=MASS_U_LIST_FLIBE, enrich_li=7.5, temp_k=900, run_openmc=False):

        self.lie = enrich_li
        self.temp = temp_k
        self.u_list = u_list
        self.name = f"FLiBe_FW_Li{self.lie:04.1f}_{self.temp}K_2025-07-22"
        self.path = f"./OpenMC/{self.name}/"
        self.run = run_openmc

        flibe = openmc.Material()
        flibe.add_elements_from_formula('F4Li2Be', 'ao', enrichment_target='Li6', enrichment_type='wo', enrichment=self.lie)
        flibe.set_density('g/cm3', DENSITY_FLIBE) 

        # Uranium tetrafluoride -- assumed to dissolve in FLiBe
        uf4 = openmc.Material()
        uf4.add_elements_from_formula('UF4','ao',ENRICH_U) # 'ao' here refers to 1:4 atomic ratio of U:F in UF4
        uf4.set_density('g/cm3', DENSITY_UF4) 

        # Calculate volume ratios of UF4 and FLiBe, ensure they add up to 1
        mix_list = []
        for mtu in self.u_list:
            vf_flibe, vf_uf4 = calc_mix_vol_fracs(mtu, volume=VOL_CC, density_flibe=DENSITY_FLIBE, displace=True) # './helper/utilities.py'
            mix = openmc.Material.mix_materials([flibe, uf4], [vf_flibe, vf_uf4], 'vo') # fractions in 'mix_materials' MUST add up to 1
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
        # energy_filter = openmc.EnergyFilter([0., 0.625, 20.0e6])
        # --Default thermal, intermediate, fast energy cutoffs in MCNP
        # energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')
        # --These have extra bins in key energy ranges. A full list of energy structures is available here: --ppark 2025-06-27
        #   https://github.com/openmc-dev/openmc/blob/6254be37582e09acff038f5656332b89e53e4eae/openmc/mgxs/__init__.py#L50-L420
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


        """ SETTINGS """
        self.settings = openmc.Settings()

        """Source from First Wall Effects"""
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

        for p in box_centers:
            src = openmc.IndependentSource()
            src.space  = openmc.stats.Point((p,0,0))
            src.angle  = openmc.stats.Isotropic()
            src.energy = openmc.stats.Discrete(energies, weights)
            source.append(src)
        self.settings.source = source


        """ Run type """
        self.settings.run_mode = 'fixed source'
        self.settings.particles = len(MASS_U_LIST) * int(1e6)
        self.settings.batches = 100


    def run_openmc(self):
        self.materials.cross_sections = set_xs_path()
        self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
        self.model.export_to_model_xml(f"./OpenMC/{self.name}/")
        self.model.run(cwd=f"./OpenMC/{self.name}/")





if __name__ == "__main__":
    main()



