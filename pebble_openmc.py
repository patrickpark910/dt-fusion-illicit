import openmc
import os, sys
import numpy as np

# Import helper functions
from Python.parameters import *
from Python.utilities import *


def main():

    print("\n\n")
    print("="*42)
    print(f"Running Helium-Cooled Pebble Bed OpenMC model for Li-6 enrichment: 60 wt%")
        
    current_run = HCPB(enrich_li=60, u_list=MASS_U_LIST_HCPB, run_openmc=False)

    if current_run.run:
        if os.path.isdir(current_run.path):
            print(f"Warning. Directory {current_run.path} already exists, so running OpenMC will fail. Skipping...")
        else:
            current_run.run_openmc()



class HCPB:

    def __init__(self, enrich_li=60, u_list=MASS_U_LIST_HCPB, temp_k=900, run_openmc=False):

        self.e = enrich_li
        self.u_list = u_list
        self.temp = temp_k
        self.run  = run_openmc
        self.name = f'HCPB_FW4.3cm_Li{self.e}_NUO2_900K_2025-07-22'
        self.path = f"./OpenMC/{self.name}/"

        """ MATERIALS
        Primary reference: P. Lu et al., "Activation analysis for the HCPB blanket module in the European DEMO", Fus. Engr. Des. (2017)
        """

        # Eurofer -- structural
        eurofer = openmc.Material(name='Eurofer')
        eurofer.set_density('g/cm3', 7.8)
        eurofer.add_element('Fe', 89.0026, percent_type='wo')
        eurofer.add_element('B', 0.001, percent_type='wo')
        eurofer.add_element('C', 0.1049, percent_type='wo')
        eurofer.add_element('N', 0.04, percent_type='wo')
        eurofer.add_element('O', 0.001, percent_type='wo')
        eurofer.add_element('Al', 0.004, percent_type='wo')
        eurofer.add_element('Si', 0.026, percent_type='wo')
        eurofer.add_element('P', 0.002, percent_type='wo')
        eurofer.add_element('S', 0.003, percent_type='wo')
        eurofer.add_element('Ti', 0.001, percent_type='wo')
        eurofer.add_element('V', 0.01963, percent_type='wo')
        eurofer.add_element('Cr', 9.00, percent_type='wo')
        eurofer.add_element('Mn', 0.55, percent_type='wo')
        eurofer.add_element('Co', 0.005, percent_type='wo')
        eurofer.add_element('Ni', 0.01, percent_type='wo')
        eurofer.add_element('Cu', 0.003, percent_type='wo')
        eurofer.add_element('Nb', 0.005, percent_type='wo')
        eurofer.add_element('Mo', 0.003, percent_type='wo')
        eurofer.add_element('Ta', 0.12, percent_type='wo')
        eurofer.add_element('W', 1.0987, percent_type='wo')

        # Li4SiO4 -- lithium ceramic breeder
        li4sio4 = openmc.Material(name='Li4SiO4')
        '''density from Activation calculations for multiple recycling of breeder ceramics by
        melt processing'''
        li4sio4.set_density('g/cm3', 2.42) 
        li4sio4.add_element('Li', 22.415, percent_type='wo', enrichment_target='Li6', enrichment_type='wo', enrichment=enrich_li)
        li4sio4.add_element('Si', 24.077, percent_type='wo')
        li4sio4.add_element('O', 53.39, percent_type='wo')
        li4sio4.add_element('Al', 0.003, percent_type='wo')
        li4sio4.add_element('C', 0.1, percent_type='wo')
        li4sio4.add_element('Ca', 0.003, percent_type='wo')
        li4sio4.add_element('Co', 0.0002, percent_type='wo')
        li4sio4.add_element('Cr', 0.0001, percent_type='wo')
        li4sio4.add_element('Cu', 0.0001, percent_type='wo')
        li4sio4.add_element('Fe', 0.0005, percent_type='wo')
        li4sio4.add_element('K', 0.0001, percent_type='wo')
        li4sio4.add_element('Mg', 0.0005, percent_type='wo')
        li4sio4.add_element('Mn', 0.0001, percent_type='wo')
        li4sio4.add_element('Pt', 0.009, percent_type='wo')
        li4sio4.add_element('Na', 0.002, percent_type='wo')
        li4sio4.add_element('Ni', 0.0002, percent_type='wo')
        li4sio4.add_element('Ti', 0.0005, percent_type='wo')
        li4sio4.add_element('Zn', 0.0002, percent_type='wo')
        li4sio4.add_element('Zr', 0.0001, percent_type='wo')

        # Beryllium -- neutron multiplier
        be = openmc.Material(name='Beryllium')
        be.set_density('g/cm3', 1.85)
        be.add_element('Be', 98.749, percent_type='wo')
        be.add_element('O', 0.9, percent_type='wo')
        be.add_element('Al', 0.09, percent_type='wo')
        be.add_element('Fe', 0.1, percent_type='wo')
        be.add_element('Mg', 0.08, percent_type='wo')
        be.add_element('Si', 0.06, percent_type='wo')
        be.add_element('Mn', 0.01, percent_type='wo')
        #Beryllium has a natural amount of Uranium that is already an activation risk
        #when I add the number of triso's i need to achieve my desired U238 mass
        # The whole sustem has a U238 mass > than expected due to the Uranium in Beryllium --ezoccoli1
        be.add_element('U', 0.01, percent_type='wo')
        be.add_element('Co', 0.001, percent_type='wo')
        be.add_element('Cu', 0.001, percent_type='wo')
        be.add_element('Fe', 0.003, percent_type='wo')
        be.add_element('K', 0.001, percent_type='wo')
        be.add_element('Mg', 0.0005, percent_type='wo')
        be.add_element('Mn', 0.0005, percent_type='wo')
        be.add_element('Na', 0.001, percent_type='wo')
        be.add_element('Nb', 0.001, percent_type='wo')
        be.add_element('Ni', 0.0005, percent_type='wo')
        be.add_element('Pb', 0.0005, percent_type='wo')
        be.add_element('Ta', 0.002, percent_type='wo')

        # Helium coolant
        he = openmc.Material(name='Helium')
        he.set_density('g/cm3', 0.0001785)
        he.add_element('He', 1)

        # UO2 kernel
        uo2 = openmc.Material(name='UO2')
        uo2.set_density('g/cm3', 10.5)
        uo2.add_elements_from_formula('UO2',enrichment=ENRICH_U) # (nominal ENRICH_U = 0.7204) --ppark 2025-07-23
        mf_uo2 = ((100-ENRICH_U)*AMU_U238+ENRICH_U*AMU_U235)/((100-ENRICH_U)*AMU_U238+ENRICH_U*AMU_U235+2*AMU_O)
        # 'mf_uo2' = mass frac of U in UO2 = 0.88 at nat U // ENRICH_U in % so subtract from 100 not 1! --ppark 2025-07-23

        # SiC
        sic = openmc.Material(name='SiC')
        sic.set_density('g/cm3', 3.2)
        sic.add_elements_from_formula('SiC')

        # BISO -- specs and geometry from Glaser & Goldston (2012)
        r_uo2 = 400e-4                                     # r = 400 μm = 0.0400 cm // "800 μm kernel"
        r_sic = 500e-4                                     # 500 μm = 0.0500 cm // "100 μm thickness"
        V_biso_particle = (4 / 3) * np.pi * (r_sic)**3     # volume of single BISO particle
        V_uo2_in_biso   = (4 / 3) * np.pi * (r_uo2)**3     # volume of UO2 in single BISO particle
        Vf_uo2_in_biso  = V_uo2_in_biso / V_biso_particle  # vol frac UO2 in single BISO
        Vf_sic_in_biso  = 1.0 - Vf_uo2_in_biso             # vol frac SiC in single BISO
        biso = openmc.Material.mix_materials([uo2, sic], [Vf_uo2_in_biso, Vf_sic_in_biso], 'vo')
        biso.set_density('g/cm3', 6.94)                    # 6.93759 g/cc from 10.5 g/cc UO2, 3.2 g/cc SiC
                                                           # -- ~7 g/cc from Glaser & Goldston --ppark 2025-07-22




        # Calculate volume ratios BISO, Eurofer, Li4SiO4, Be, He -- ensure they add up to 1
        mix_list = []
        for mtu in self.u_list:
            # Calculate volume of BISO to put in
            m_uo2 = mtu*1e6 * (1/mf_uo2)
            V_uo2 = m_uo2 / 10.5            # total vol of uo2 to be put in model
            V_biso = V_uo2 / Vf_uo2_in_biso # total vol of biso to be put in model
            Vf_biso = V_biso / VOL_CC
            print(f"V_uo2 {V_uo2:.4f} // V_uo2_in_biso {V_uo2_in_biso:.4f}")

            # Nominal volfracs minus BISO volfrac
            Vf_eurofer = 0.1176 # Lu et al. (2017)
            Vf_li4sio4 = 0.1304 - Vf_biso
            Vf_be      = 0.3790
            Vf_he      = 0.3730

            mix = openmc.Material.mix_materials([eurofer, li4sio4, be, he, biso],[Vf_eurofer, Vf_li4sio4, Vf_be, Vf_he, Vf_biso],'vo')
            mix.name = 'Full Lithium Ceramic + TRISO Blanket Homogenized'
            mix.temperature = self.temp
            mix_list.append(mix)

        self.materials = openmc.Materials(mix_list)
        self.materials.cross_sections = set_xs_path()


        """ GEOMETRY """
        cells, pitch, half_box = [], 120, 50  # +/- 50 cm bounds
        box_centers = [pitch * i for i in range(len(self.u_list))]  # used in Sources

        for i, material in enumerate(mix_list):
            x_min = openmc.XPlane(x0=-half_box + box_centers[i], boundary_type='reflective')
            x_max = openmc.XPlane(x0=half_box + box_centers[i], boundary_type='reflective')
            y_min, y_max = openmc.YPlane(-50, boundary_type='reflective'), openmc.YPlane(50, boundary_type='reflective')
            z_min, z_max = openmc.ZPlane(-50, boundary_type='reflective'), openmc.ZPlane(50, boundary_type='reflective')
            region = +x_min & -x_max & +y_min & -y_max & +z_min & -z_max

            cell = openmc.Cell(fill=material, region=region)
            cell.name = f"mix-{i + 1}"
            cells.append(cell)

        root_univ = openmc.Universe(cells=cells)  # Create root universe with all material cells
        root_univ.plot(width=(pitch * len(self.u_list), 110), origin=(pitch * (len(self.u_list) - 1) / 2, 0.0, 0.0))  # Visualize
        self.geometry = openmc.Geometry(root_univ)  # Set geometry


        """ TALLIES """
        self.tallies = openmc.Tallies() # initialize

        # Filters
        cell_filter = openmc.CellFilter(cells)

        E_bin_edges = logspace_per_decade(1e-5, 20e6, 100) # './Python/utilities.py'
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

        # Beryllium reaction rates
        Be_tally = openmc.Tally(name='beryllium rxn rates')
        Be_tally.scores = ['(n,gamma)','(n,2n)', 'elastic'] # specific names required
        Be_tally.nuclides = ['Be9']
        Be_tally.filters = filters

        self.tallies.extend([flux_tally, U_tally, Li_tally, Be_tally])

        
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
        self.settings.particles = len(self.u_list) * int(1e6)
        self.settings.batches = 100



    def run_openmc(self):
        self.materials.cross_sections = set_xs_path()
        self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
        self.model.export_to_model_xml(f"./OpenMC/{self.name}/")
        self.model.run(cwd=f"./OpenMC/{self.name}/")





if __name__ == "__main__":
    main()



