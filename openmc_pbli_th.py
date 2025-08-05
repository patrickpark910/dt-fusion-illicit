import openmc
import os, sys
import numpy as np

# Import helper functions
from Python.parameters import *
from Python.utilities import *


def main():

    print("\n\n")
    print("="*42)
    print(f"Running PbLi OpenMC model for Li-6 enrichment: 90 wt%")
        
    current_run = DCLL(enrich_li=90, u_list=MASS_U_LIST_PBLI, run_openmc=True)

    if current_run.run:
        if os.path.isdir(current_run.path):
            print(f"Warning. Directory {current_run.path} already exists, so running OpenMC will fail. Skipping...")
        else:
            current_run.run_openmc()


class DCLL:

    def __init__(self, enrich_li=90, u_list=MASS_U_LIST_PBLI, run_openmc=False):

        self.e = enrich_li
        self.name = f'DCLL_Li{self.e:04.1f}_NUO2_900K_2025-07-22'
        self.u_list = u_list
        self.path = f"./OpenMC/{self.name}/"
        self.run = run_openmc

        # Pb-17Li
        pbli = openmc.Material()  # 83% Pb, 17% Li by atomic fraction
        pbli.set_density('g/cm3', DENSITY_PBLI)
        pbli.add_element('Pb', 0.83)
        pbli.add_element('Li', 0.17, enrichment=self.e, enrichment_target='Li6', enrichment_type='wo')  # Li-6 enrichment to 90%

        # 13.8 wt%-enriched U2O3C (half UO2 + half UCO) -- from OpenMC demo --ppark 2025-07-23
        # uco = openmc.Material(name='kernel')
        # uco.set_density('g/cm3', 10.5)
        # uco.add_nuclide('U235', 4.6716e-02)
        # uco.add_nuclide('U238', 2.8697e-01)
        # uco.add_nuclide('O16', 5.0000e-01)
        # uco.add_element('C', 1.6667e-01)
        # Mu238 = uco.get_mass_density('U238')
        # Mfuel = uco.density  # used to compute mass fraction of U238 in fuel
        # calculate total fuel mass given a mass of U238 as MTU (self.u_list)

        # UO2 kernel (nominal ENRICH_U = 0.7204) --ppark 2025-07-23
        uo2 = openmc.Material(name='UO2')
        uo2.set_density('g/cm3', 10.5)
        uo2.add_elements_from_formula('UO2',enrichment=ENRICH_U)
        mf_uo2 = ((100-ENRICH_U)*AMU_U238+ENRICH_U*AMU_U235)/((100-ENRICH_U)*AMU_U238+ENRICH_U*AMU_U235+2*AMU_O)
        # 'mf_uo2' = mass frac of U in UO2 = 0.88 at nat U // ENRICH_U in % so subtract from 100 not 1! --ppark 2025-07-23

        # SiC FCI or structural inserts AND coating for TRISO
        sic = openmc.Material(name='SiC')
        sic.set_density('g/cm3', 3.2)
        sic.add_elements_from_formula('SiC')

        # Glaser & Goldston BISO particles used for simplicity
        # and to validate our results with their model
        # BISO particles include our fuel pellet coated in one layer of SiC

        r_uo2 = 400e-4  # r = 400 μm = 0.0400 cm // "800 μm kernel"
        r_sic = 500e-4  # 500 μm = 0.0500 cm // "100 μm thickness"
        V_biso_particle = (4 / 3) * np.pi * (r_sic)**3     # volume of single BISO particle
        V_uo2_in_biso   = (4 / 3) * np.pi * (r_uo2)**3     # volume of UO2 in single BISO particle
        Vf_uo2_in_biso  = V_uo2_in_biso / V_biso_particle  # vol frac UO2 in single BISO
        Vf_sic_in_biso  = 1.0 - Vf_uo2_in_biso             # vol frac SiC in single BISO

        biso = openmc.Material.mix_materials([uo2, sic], [Vf_uo2_in_biso, Vf_sic_in_biso], 'vo')
        biso.set_density('g/cm3', 6.94) # ~7 g/cc from Glaser & Goldston // 6.93759 g/cc from 10.5 g/cc UO2, 3.2 g/cc SiC --ppark 2025-07-22

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

        # helium coolant
        he = openmc.Material(name='Helium')
        he.set_density('g/cm3', 0.0001785)
        he.add_element('He', 1)

        # Calculate volume ratios of TRISO vs PbLi+structure, ensure they add up to 1
        mix_list = []
        for mtu in self.u_list:
            # Calculate volume of BISO to put in
            m_uo2 = mtu*1e6 * (1/mf_uo2)
            V_uo2 = m_uo2 / 10.5  # total vol of uo2 to be put in model
            V_biso = V_uo2 / Vf_uo2_in_biso # total vol of biso to be put in model
            print(f"V_uo2 {V_uo2:.4f} // V_uo2_in_biso {V_uo2_in_biso:.4f}")
            # Deduct BISO volume from Lead-Lithium's volume
            Vf_fs, Vf_ll, Vf_sic, Vf_he = 0.019, 0.805, 0.079, 0.097  # volfracs Glaser & Goldston 12, Tb.1, breeding channel 2
            Vf_biso = V_biso / VOL_CC
            Vf_ll = 0.805 - Vf_biso
            print(f"volume BISO: {V_biso:.4f}")
            print(f"volume fractions FS/LL/SiC/He/BISO: {Vf_fs:.4f}, {Vf_ll:.4f}, {Vf_sic:.4f}, {Vf_he:.4f}, {Vf_biso:.4f}")

            mix = openmc.Material.mix_materials([f82h,pbli,sic,he,biso], [Vf_fs, Vf_ll, Vf_sic, Vf_he, Vf_biso], 'vo')
            mix.name = 'FS+LL+SiC+He+BISO Homogenized'
            mix.temperature = TEMP_K
            mix_list.append(mix)
            # print(f"\nProperties of PbLi-BISO with {mtu} MTU:")
            # print(mix)

        self.materials = openmc.Materials(mix_list)


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
        flux_tally.scores = ['flux'] 
        flux_tally.filters = filters

        # Uranium reaction rates
        U_tally = openmc.Tally(name='uranium rxn rates')
        U_tally.scores = ['(n,gamma)','fission', 'elastic'] 
        U_tally.nuclides = ['U238', 'U235']
        U_tally.filters = filters

        # Lithium reaction rates
        Li_tally = openmc.Tally(name='lithium rxn rates')
        Li_tally.scores = ['(n,gamma)','(n,Xt)', 'elastic'] 
        Li_tally.nuclides = ['Li6', 'Li7']
        Li_tally.filters = filters

        # Lead reaction rates
        Pb_tally = openmc.Tally(name='lead rxn rates')
        Pb_tally.scores = ['(n,gamma)', '(n,2n)', 'elastic']  
        Pb_tally.nuclides = ['Pb204', 'Pb206', 'Pb207', 'Pb208']
        Pb_tally.filters = filters

        self.tallies.extend([flux_tally, U_tally, Li_tally, Pb_tally])


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
        self.settings.particles = len(MASS_U_LIST_PBLI) * int(1e6)
        self.settings.batches = 100


    def run_openmc(self):
        self.materials.cross_sections = set_xs_path()
        self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
        self.model.export_to_model_xml(f"./OpenMC/{self.name}/")
        self.model.run(cwd=f"./OpenMC/{self.name}/")


if __name__ == "__main__":
    main()



