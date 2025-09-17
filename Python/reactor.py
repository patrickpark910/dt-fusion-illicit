import openmc
import os, sys
import numpy as np
from datetime import date

# Import helper functions
from Python.parameters import *
from Python.utilities import *

class Reactor:

    def __init__(self, blanket='flibe', fertile_element='U', fertile_density_kgm3=0.0, temp_k=900, run_openmc=False):

        self.fertile_density_gcm3 = fertile_density_kgm3*0.001  #  g/cm3
        self.fertile_element = fertile_element.capitalize()
        self.temp_k = int(temp_k)
        self.run_openmc = run_openmc

        if blanket.lower() == 'flibe':
            self.blanket         = 'FLiBe'
            self.blanket_density = DENSITY_FLIBE
            self.blanket_enrich  = ENRICH_FLIBE
            self.blanket_volume  = VOLUME_FLIBE
            print(self.blanket_volume)

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
        self.name = f"{self.blanket}_{self.fertile_element}{fertile_density_kgm3:.2f}kgm3_Li{self.blanket_enrich:04.1f}_{self.temp_k}K" # _{today}
        self.path = f"./OpenMC/{self.name}/"
        self.run = run_openmc

        self.materials()
        print('1')
        self.geometry()
        print('2')
        self.settings()
        print('3')
        self.tallies()
        print('4')
        self.run_openmc()
        print('5')

    def materials(self):
        """
        OpenMC Materials
        """
        # Blanket
        blanket = openmc.Material()
        blanket.set_density('g/cm3', self.blanket_density)

        if self.blanket == 'FLiBe':
            blanket.add_elements_from_formula('F4Li2Be', 'ao', enrichment_target='Li6', enrichment_type='wo', enrichment=self.blanket_enrich)

        elif self.blanket.lower() == 'll':
            blanket.add_element('Pb', 0.83)
            blanket.add_element('Li', 0.17, enrichment=self.blanket_enrich, enrichment_target='Li6', enrichment_type='wo')  # Li-6 enrichment to 90%

        # elif self.blanket.lower() == 'pb':
        #     self.blanket_density = DENSITY_PB
        #     self.blanket_enrich  = ENRICH_PB
        

        # Fertile compound
        fertile = openmc.Material()

        if self.blanket == 'FLiBe':
            if self.fertile_element == 'U':
                fertile.add_elements_from_formula('UF4','ao') # 'ao' here refers to 1:4 atomic ratio of U:F in UF4
                fertile.set_density('g/cm3', DENSITY_UF4) 
            elif self.fertile_element == 'Th':
                fertile.add_elements_from_formula('ThF4','ao') # 'ao' here refers to 1:4 atomic ratio of U:F in UF4
                fertile.set_density('g/cm3', DENSITY_ThF4) 

        elif self.blanket in ['LL', 'PB']:
            if self.fertile_element == 'U':
                pass  # UO2 pebbles
            elif self.fertile_element == 'Th':
                pass  # ThO2 pebbles


        # Mix fertile compound to blanket
        if self.blanket == 'FLiBe':
            pass 
            # mass_uf4 = 


            vf_flibe, vf_uf4 = calc_mix_vol_fracs(mtu, volume=VOL_CC, density_flibe=DENSITY_FLIBE, displace=True) # './helper/utilities.py'
            mix = openmc.Material.mix_materials([blanket, fertile], [vf_flibe, vf_uf4], 'vo') # fractions in 'mix_materials' MUST add up to 1
            mix.name, mix.volume, mix.temperature = f"FLiBe {mtu:.1f} MTU at {self.temp} K", VOL_CC, self.temp
            mix_list.append(mix)
            self.materials = openmc.Materials(mix_list)

        elif self.blanket in ['LL', 'PB']:
            pass


    def geometry(self):

        root_univ = openmc.Universe(cells=cells) # Create root universe with all material cells
        root_univ.plot(width=(pitch * len(self.u_list), 110), origin=(pitch * (len(self.u_list) - 1) / 2, 0.0, 0.0)) # Visualize
        self.geometry = openmc.Geometry(root_univ) # Set geometry


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