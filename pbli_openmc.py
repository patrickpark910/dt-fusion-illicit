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
    print(f"Running PbLi OpenMC model for Li-6 enrichment: 90 wt%")
        
    current_run = DCLL(enrich_li=90, u_list=MASS_U_LIST_PBLI)

    if os.path.isdir(current_run.path):
        print(f"Warning. Directory {current_run.path} already exists, so running OpenMC will fail. Skipping...")
        return
    else:
        current_run.set_xs_path()
        current_run.run_openmc()



class DCLL:

    def __init__(self, enrich_li=90, u_list=MASS_U_LIST_PBLI):

        self.e = enrich_li
        self.name = f'PbLi_Li{self.e}_7_22'
        self.u_list = u_list
        self.path = f"./OpenMC/{self.name}/"

        # PbLi
        pbli = openmc.Material()  # 83% Pb, 17% Li by atomic fraction
        pbli.set_density('g/cm3', DENSITY_PBLI)
        pbli.add_element('Pb', 0.83)
        pbli.add_element('Li', 0.17, enrichment=self.e, enrichment_target='Li6', enrichment_type='wo')  # Li-6 enrichment to 90%

        # ceramic TRISO -- neglect buoyancy effects --ezoccoli 2025-07-07
        # TRISO model from OpenMC examples
        fuel = openmc.Material(name='Fuel')
        fuel.set_density('g/cm3', 10.5)
        fuel.add_nuclide('U235', 4.6716e-02)
        fuel.add_nuclide('U238', 2.8697e-01)
        fuel.add_nuclide('O16', 5.0000e-01)
        fuel.add_element('C', 1.6667e-01)
        Mu238 = fuel.get_mass_density('U238')
        Mfuel = fuel.density  # used to compute mass fraction of U238 in fuel
        # calculate total fuel mass given a mass of U238 (self.u_list)

        # SiC FCI or structural inserts AND coating for TRISO
        SiC = openmc.Material(name='SiC')
        SiC.set_density('g/cm3', 3.2)
        SiC.add_element('C', 0.5)
        SiC.add_element('Si', 0.5)

        radius_fuel = 400e-4  # 400 μm = 0.0400 cm
        radius_sic = 500e-4  # 500 μm = 0.0500 cm
        vol_tot = (4 / 3) * np.pi * (radius_sic) ** 3
        vol_fuel = (4 / 3) * np.pi * (radius_fuel) ** 3
        vf_fuel = vol_fuel / vol_tot
        vf_sic = 1.0 - vf_fuel

        triso = openmc.Material.mix_materials([fuel, SiC], [vf_fuel, vf_sic], 'vo')
        triso.set_density('g/cm3', DENSITY_TRISO)

        f82h = openmc.Material(name='F82H Steel')  # low-activation ferritic steel FS
        f82h.set_density('g/cm3', 7.87)  # typical density
        # https://www.tandfonline.com/doi/abs/10.1179/mst.1998.14.6.573
        # Mass fractions from fusion blanket design literature
        f82h.add_element('Fe', 0.88675)
        f82h.add_element('Cr', 0.0780)
        f82h.add_element('W', 0.020)
        f82h.add_element('V', 0.0016)
        f82h.add_element('Ti', 0.0002)
        f82h.add_element('Mn', 0.0018)
        f82h.add_element('C', 0.0009)
        f82h.add_element('Si', 0.0013)
        f82h.add_element('N', 0.00006)
        f82h.add_element('P', 0.00004)
        f82h.add_element('S', 0.00003)
        f82h.add_element('Ni', 0.0004)
        f82h.add_element('Mo', 0.0001)
        f82h.add_element('Al', 0.0001)
        f82h.add_element('Nb', 0.0001)
        # f82h.add_element('Ta', 0.0002) # no Ta180_m1 in ENDF8

        he = openmc.Material(name='Helium')
        he.set_density('g/cm3', 0.0001785)
        he.add_element('He', 1)

        # Calculate volume ratios of TRISO vs PbLi+structure, ensure they add up to 1
        mix_list = []
        for mtu in self.u_list:
            mass_fuel = mtu * 10 ** 6 * (Mfuel / Mu238)
            V_fuel = mass_fuel / 10.5  # mass over fuel density
            V_triso = V_fuel / vf_fuel
            vf_triso = V_triso / (VOL_CC * 0.805)
            vf_pbli = 1 - vf_triso
            pbliTRISO = openmc.Material.mix_materials([pbli, triso], [vf_pbli, vf_triso], 'vo')  # fractions in 'mix_materials' MUST add up to 1
            pbliTRISO.volume = 0.805 * VOL_CC
            vf_pbliTRISO = 0.805
            vf_he = 0.097
            vf_f82h = 0.019
            vf_sic = 0.079  # volume fractions from Breeding Ch. 2 Table 1. Glaser and Goldston
            mix = openmc.Material.mix_materials([pbliTRISO, he, f82h, SiC], [vf_pbliTRISO, vf_he, vf_f82h, vf_sic], 'vo')
            mix.name = 'Full PbLi + TRISO Blanket Homogenized'
            mix.temperature = TEMP_K
            mix.volume = VOL_CC
            mix_list.append(mix)

            #print(f"\nProperties of PbLi-TRISO with {mtu} MTU:")
            #print(mix)

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

        # Lead reaction rates
        Pb_tally = openmc.Tally(name='lead rxn rates')
        Pb_tally.scores = ['(n,gamma)', '(n,2n)', 'elastic']  # specific names required
        Pb_tally.nuclides = ['Pb204', 'Pb206', 'Pb207', 'Pb208']
        Pb_tally.filters = filters

        self.tallies.extend([flux_tally, U_tally, Li_tally, Pb_tally])

 
        
        """First Wall Effects"""
        sp = openmc.StatePoint(f'./OpenMC/FirstWall_50/statepoint.100.h5')  
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
        self.settings.particles = len(MASS_U_LIST) * int(1e1)  #  
        self.settings.batches = 100


    def run_openmc(self):
        self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
        self.model.export_to_model_xml(f"./OpenMC/{self.name}/")
        self.model.run(cwd=f"./OpenMC/{self.name}/")

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
            #print(self.materials.cross_sections )
        else:
            sys.exit(f"Error finding cross section XML!")





if __name__ == "__main__":
    main()



