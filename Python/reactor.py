import openmc
import os, sys
import numpy as np
import pandas as pd
from datetime import date


# Import helper functions
from Python.parameters import *
from Python.utilities import *

class Reactor:

    @timer
    def __init__(self, breeder='flibe', fertile_element='U', fertile_bulk_density_kgm3=0.0, calc_volumes=False, run_openmc=False):

        self.fertile_bulk_density_kgm3 = fertile_bulk_density_kgm3
        self.fertile_element = fertile_element.capitalize()
        self.run_openmc   = run_openmc    # i'm calling this self.run_openmc and self.openmc() is the function
        self.calc_volumes = calc_volumes

        if breeder.lower() == 'flibe':
            self.temp_k          = TEMP_K
            self.breeder_name    = 'FLiBe'
            self.breeder_density = DENSITY_FLIBE # g/cm^3
            self.breeder_enrich  = ENRICH_FLIBE  # wt% 
            self.breeder_volume  = FLIBE_BR_VOL  # m^3

        elif breeder.lower() == 'll':
            self.temp_k          = TEMP_K
            self.breeder_name    = 'LL'
            self.breeder_density = DENSITY_LL # g/cm^3
            self.breeder_enrich  = ENRICH_LL  # wt% 
            self.breeder_volume  = LL_BR_VOL  # m^3 must check could be fishy -ezoccoli

        # elif blanket.lower() == 'pb':
        #     self.blanket         = 'PB'
        #     self.blanket_density = DENSITY_PB
        #     self.blanket_enrich  = ENRICH_PB

        elif breeder.lower() == 'arc':
            self.temp_k          = TEMP_K
            self.breeder_name    = 'ARC'
            self.breeder_density = DENSITY_FLIBE # g/cm^3
            self.breeder_enrich  = ENRICH_FLIBE  # wt%
            self.breeder_volume  = ARC_BR_VOL    # m^3


        """
        Name file based on reactor config
        """
        today     = date.today().strftime("%Y-%m-%d")
        self.name = f"{self.breeder_name}_{self.fertile_element}{self.fertile_bulk_density_kgm3:3.2f}kgm3_Li{self.breeder_enrich:04.1f}_{self.temp_k}K" # _{today}
        # should come out to smth like: FLiBe_U010kgm3_Li7.5_900K
        self.path = f"./OpenMC/{self.name}"
        
        # os.makedirs(self.path, exist_ok=True)

        start_msg = f"="*42+f"\nWriting an OpenMC input for {self.name}"
        print(f"{Colors.RED}{start_msg}{Colors.END}")

        self.materials()
        self.geometry()
        self.settings()
        self.tallies()

    @timer
    def materials(self):
        """
        OpenMC Materials
        """

        # ------------------------------------------------------------------
        # Air
        # ------------------------------------------------------------------
        self.air = openmc.Material(name='air')
        self.air.set_density('g/cm3', 0.001225)

        # Atom fractions for dry air
        self.air.add_element('N', 0.78084, percent_type='ao')   # Nitrogen
        self.air.add_element('O', 0.20946, percent_type='ao')   # Oxygen
        self.air.add_element('Ar', 0.00934, percent_type='ao')  # Argon
        self.air.add_element('C', 0.00036, percent_type='ao')   # Carbon from CO2

        # ------------------------------------------------------------------
        # First wall 
        # ------------------------------------------------------------------

        self.firstwall = openmc.Material(name='firstwall', temperature=self.temp_k)
        self.firstwall.depletable = False

        if self.breeder_name in ['FLiBe','ARC']:
            self.firstwall.add_element('O',5/1e6,percent_type='wo')
            self.firstwall.add_element('N',5/1e6,percent_type='wo')
            self.firstwall.add_element('C',5/1e6,percent_type='wo')
            self.firstwall.add_element('Na',4/1e6,percent_type='wo')
            self.firstwall.add_element('K',2.5/1e6,percent_type='wo')
            self.firstwall.add_element('Al',3/1e6,percent_type='wo')
            self.firstwall.add_element('Ca',0.5/1e6,percent_type='wo')
            self.firstwall.add_element('Cr',0.5/1e6,percent_type='wo')
            self.firstwall.add_element('Cu',0.5/1e6,percent_type='wo')
            self.firstwall.add_element('Fe',5/1e6,percent_type='wo')
            self.firstwall.add_element('W',1-(5+5+5+4+2.5+3+0.5+0.5+0.5+5)/1e6,percent_type='wo')
            self.firstwall.set_density('g/cm3',19.3)

        elif self.breeder_name in ['LL']:
            self.firstwall.add_nuclide("W180", 0.00007586640, "ao")
            self.firstwall.add_nuclide("W182", 0.01675383000, "ao")
            self.firstwall.add_nuclide("W183", 0.00904706820, "ao")
            self.firstwall.add_nuclide("W184", 0.01937122080, "ao")
            self.firstwall.add_nuclide("W186", 0.01797401460, "ao")
            self.firstwall.set_density("g/cm3", 19.3)  # tungsten density from mcnp 0.06322200000 atomic density in atoms/barn-cm
        

        # ------------------------------------------------------------------
        # Structure 
        # ------------------------------------------------------------------


        if self.breeder_name in ['FLiBe','ARC']:
            """ V-4Cr-4Ti """
            self.structure = openmc.Material(name='structure', temperature=self.temp_k)
            self.structure.depletable = False
            # This code is from jlball but V-4Cr-4Ti specs also can be found
            # from ANL material id BL-47, p.3 (1994) altho these are slightly different
            # cf. osti.gov/servlets/purl/10194461 --ppark 2025-07-22
            self.structure.add_element('Cr',0.04,percent_type='wo')
            self.structure.add_element('Ti',0.04,percent_type='wo')
            self.structure.add_element('C',56/1e6,percent_type='wo')
            self.structure.add_element('O',181/1e6,percent_type='wo')
            self.structure.add_element('N',103/1e6,percent_type='wo')
            self.structure.add_element('B',7/1e6,percent_type='wo')
            self.structure.add_element('Na',17/1e6,percent_type='wo')
            self.structure.add_element('Mg',0.5/1e6,percent_type='wo')
            self.structure.add_element('Al',119/1e6,percent_type='wo')
            self.structure.add_element('Si',280/1e6,percent_type='wo')
            self.structure.add_element('Mn',0.5/1e6,percent_type='wo')
            self.structure.add_element('Fe',80/1e6,percent_type='wo')
            self.structure.add_element('Ni',13/1e6,percent_type='wo')
            self.structure.add_element('Cu',4/1e6,percent_type='wo')
            self.structure.add_element('V',1-0.04-0.04-(56+181+103+7+17+0.5+119+280+0.5+80+13+4)/1e6,percent_type='wo')
            self.structure.set_density('g/cm3',6.05) 
            # This density value is sus and needs a good source --jlball 
            # This value is from Metals Handbook, 9th ed, vol 2: "Properties and Selection: Nonferrous Alloys and Pure Metals" (1979) --ppark 2025-07-22
        
        elif self.breeder_name in ['LL']:
            # First wall front (F82H steel)
            self.firstwall_front = openmc.Material(name='firstwall_front', temperature=self.temp_k)
            self.firstwall_front.depletable = False
            self.firstwall_front.add_nuclide("C0",    0.00040000000, "ao")
            self.firstwall_front.add_nuclide("Si28",  0.00015679100, "ao")
            self.firstwall_front.add_nuclide("Si29",  0.00000795600, "ao")
            self.firstwall_front.add_nuclide("Si30",  0.00000525300, "ao")
            self.firstwall_front.add_nuclide("V51",   0.00018400000, "ao")
            self.firstwall_front.add_nuclide("Cr50",  0.00031327500, "ao")
            self.firstwall_front.add_nuclide("Cr52",  0.00604119000, "ao")
            self.firstwall_front.add_nuclide("Cr53",  0.00068502200, "ao")
            self.firstwall_front.add_nuclide("Cr54",  0.00017051700, "ao")
            self.firstwall_front.add_nuclide("Mn55",  0.00008600000, "ao")
            self.firstwall_front.add_nuclide("Fe54",  0.00438860000, "ao")
            self.firstwall_front.add_nuclide("Fe56",  0.06889170000, "ao")
            self.firstwall_front.add_nuclide("Fe57",  0.00159101000, "ao")
            self.firstwall_front.add_nuclide("Fe58",  0.00021173400, "ao")
            self.firstwall_front.add_nuclide("Ta181",  0.00001000000, "ao") 
            self.firstwall_front.add_nuclide("W182",  0.00013515000, "ao")
            self.firstwall_front.add_nuclide("W183",  0.00007298100, "ao")
            self.firstwall_front.add_nuclide("W184",  0.00015626400, "ao")
            self.firstwall_front.add_nuclide("W186",  0.00014499300, "ao")
            self.firstwall_front.set_density("atom/b-cm", 0.08365243600)  # from MCNP cell 112

            self.firstwall_cooling = openmc.Material(name="firstwall_cooling", temperature=self.temp_k)
            self.firstwall_cooling.depletable = False
            self.firstwall_cooling.add_nuclide("He4", 0.00049800000, "ao")
            self.firstwall_cooling.add_nuclide("C0",  0.00006800000, "ao")
            self.firstwall_cooling.add_nuclide("Si28",0.00002665450, "ao")
            self.firstwall_cooling.add_nuclide("Si29",0.00000135250, "ao")
            self.firstwall_cooling.add_nuclide("Si30",0.00000089300, "ao")
            self.firstwall_cooling.add_nuclide("V51", 0.00003128000, "ao")
            self.firstwall_cooling.add_nuclide("Cr50",0.00005325680, "ao")
            self.firstwall_cooling.add_nuclide("Cr52",0.00102700000, "ao")
            self.firstwall_cooling.add_nuclide("Cr53",0.00011645400, "ao")
            self.firstwall_cooling.add_nuclide("Cr54",0.00002898790, "ao")
            self.firstwall_cooling.add_nuclide("Mn55",0.00001462000, "ao")
            self.firstwall_cooling.add_nuclide("Fe54",0.00074606200, "ao")
            self.firstwall_cooling.add_nuclide("Fe56",0.01171160000, "ao")
            self.firstwall_cooling.add_nuclide("Fe57",0.00027047200, "ao")
            self.firstwall_cooling.add_nuclide("Fe58",0.00003599480, "ao")
            self.firstwall_cooling.add_nuclide("Ta181",0.00000170000, "ao")
            self.firstwall_cooling.add_nuclide("W182",0.00002297550, "ao")
            self.firstwall_cooling.add_nuclide("W183",0.00001240680, "ao")
            self.firstwall_cooling.add_nuclide("W184",0.00002656490, "ao")
            self.firstwall_cooling.add_nuclide("W186",0.00002464880, "ao")
            self.firstwall_cooling.set_density("atom/b-cm", 0.01471892350)

            self.back_plate = openmc.Material(name="back_plate", temperature=self.temp_k)
            self.back_plate.depletable = False
            self.back_plate.add_nuclide("C0",    0.00040000000, "ao")
            self.back_plate.add_nuclide("Si28",  0.00015679100, "ao")
            self.back_plate.add_nuclide("Si29",  0.00000795600, "ao")
            self.back_plate.add_nuclide("Si30",  0.00000525300, "ao")
            self.back_plate.add_nuclide("V51",   0.00018400000, "ao")
            self.back_plate.add_nuclide("Cr50",  0.00031327500, "ao")
            self.back_plate.add_nuclide("Cr52",  0.00604119000, "ao")
            self.back_plate.add_nuclide("Cr53",  0.00068502200, "ao")
            self.back_plate.add_nuclide("Cr54",  0.00017051700, "ao")
            self.back_plate.add_nuclide("Mn55",  0.00008600000, "ao")
            self.back_plate.add_nuclide("Fe54",  0.00438860000, "ao")
            self.back_plate.add_nuclide("Fe56",  0.06889170000, "ao")
            self.back_plate.add_nuclide("Fe57",  0.00159101000, "ao")
            self.back_plate.add_nuclide("Fe58",  0.00021173400, "ao")
            self.back_plate.add_nuclide("Ta181", 0.00001000000, "ao")
            self.back_plate.add_nuclide("W182",  0.00013515000, "ao")
            self.back_plate.add_nuclide("W183",  0.00007298100, "ao")
            self.back_plate.add_nuclide("W184",  0.00015626400, "ao")
            self.back_plate.add_nuclide("W186",  0.00014499300, "ao")
            self.back_plate.set_density("atom/b-cm", 0.08365243600)
        # ------------------------------------------------------------------
        # Breeder material
        # ------------------------------------------------------------------

        self.breeder = openmc.Material(name='breeder', temperature=self.temp_k)
        self.breeder.set_density('g/cm3', self.breeder_density)

        if self.breeder_name in ['FLiBe','ARC']:
            self.breeder.add_elements_from_formula('F4Li2Be', 'ao', enrichment_target='Li6', enrichment_type='wo', enrichment=self.breeder_enrich)

        elif self.breeder_name in ['LL']:
            # --emma: might need He4...
            self.breeder.add_element('Pb', 0.83, percent_type='wo')  # Pb isotopes expanded below
            self.breeder.add_element('Li', 0.17, percent_type='wo', enrichment=self.breeder_enrich, enrichment_target='Li6', enrichment_type='wo') #Li-6 enrichment to 90%
        
        # elif self.blanket.lower() == 'pb':
        #     self.blanket_density = DENSITY_PB
        #     self.blanket_enrich  = ENRICH_PB

        # ------------------------------------------------------------------
        # Divider material
        # ------------------------------------------------------------------
        
        if self.breeder_name in ['LL']:
            self.divider1 = openmc.Material(name="divider1", temperature=self.temp_k)
            self.divider1.depletable = False
            self.divider1.add_nuclide("He4", 0.00029280000, "ao")
            self.divider1.add_nuclide("C0", 0.00020480000, "ao")
            self.divider1.add_nuclide("Si28", 0.00008027700, "ao")
            self.divider1.add_nuclide("Si29", 0.00000407347, "ao")
            self.divider1.add_nuclide("Si30", 0.00000268954, "ao")
            self.divider1.add_nuclide("V51", 0.00009420800, "ao")
            self.divider1.add_nuclide("Cr50", 0.00016039700, "ao")
            self.divider1.add_nuclide("Cr52", 0.00309309000, "ao")
            self.divider1.add_nuclide("Cr53", 0.00035073100, "ao")
            self.divider1.add_nuclide("Cr54", 0.00008730470, "ao")
            self.divider1.add_nuclide("Mn55", 0.00004403200, "ao")
            self.divider1.add_nuclide("Fe54", 0.00224696000, "ao")
            self.divider1.add_nuclide("Fe56", 0.03527260000, "ao")
            self.divider1.add_nuclide("Fe57", 0.00081459700, "ao")
            self.divider1.add_nuclide("Fe58", 0.00010840800, "ao")
            self.divider1.add_nuclide("Ta181", 0.00000512000, "ao")
            self.divider1.add_nuclide("W182", 0.00006919680, "ao")
            self.divider1.add_nuclide("W183", 0.00003736630, "ao")
            self.divider1.add_nuclide("W184", 0.00008000720, "ao")
            self.divider1.add_nuclide("W186", 0.00007423640, "ao")
            self.divider1.set_density("atom/b-cm", 0.04312289441)

            self.divider2 = openmc.Material(name="divider2", temperature=self.temp_k)
            self.divider2.depletable = False
            self.divider2.add_nuclide("He4", 0.00029280000, "ao")
            self.divider2.add_nuclide("C0", 0.00020480000, "ao")
            self.divider2.add_nuclide("Si28", 0.00008027700, "ao")
            self.divider2.add_nuclide("Si29", 0.00000407347, "ao")
            self.divider2.add_nuclide("Si30", 0.00000268954, "ao")
            self.divider2.add_nuclide("V51", 0.00009420800, "ao")
            self.divider2.add_nuclide("Cr50", 0.00016039700, "ao")
            self.divider2.add_nuclide("Cr52", 0.00309309000, "ao")
            self.divider2.add_nuclide("Cr53", 0.00035073100, "ao")
            self.divider2.add_nuclide("Cr54", 0.00008730470, "ao")
            self.divider2.add_nuclide("Mn55", 0.00004403200, "ao")
            self.divider2.add_nuclide("Fe54", 0.00224696000, "ao")
            self.divider2.add_nuclide("Fe56", 0.03527260000, "ao")
            self.divider2.add_nuclide("Fe57", 0.00081459700, "ao")
            self.divider2.add_nuclide("Fe58", 0.00010840800, "ao")
            self.divider2.add_nuclide("Ta181", 0.00000512000, "ao")
            self.divider2.add_nuclide("W182", 0.00006919680, "ao")
            self.divider2.add_nuclide("W183", 0.00003736630, "ao")
            self.divider2.add_nuclide("W184", 0.00008000720, "ao")
            self.divider2.add_nuclide("W186", 0.00007423640, "ao")
            self.divider2.set_density("atom/b-cm", 0.04312289441)

            #--INNER MANIFOLD-- 
            self.inner_manifold = openmc.Material(name="inner_manifold", temperature=self.temp_k)
            self.inner_manifold.depletable = False
            self.inner_manifold.add_nuclide("He4", 0.00032820000, "ao")
            self.inner_manifold.add_nuclide("C0", 0.00018120000, "ao")
            self.inner_manifold.add_nuclide("Si28", 0.00007102630, "ao")
            self.inner_manifold.add_nuclide("Si29", 0.00000360407, "ao")
            self.inner_manifold.add_nuclide("Si30", 0.00000237961, "ao")
            self.inner_manifold.add_nuclide("V51", 0.00008335200, "ao")
            self.inner_manifold.add_nuclide("Cr50", 0.00014191400, "ao")
            self.inner_manifold.add_nuclide("Cr52", 0.00273666000, "ao")
            self.inner_manifold.add_nuclide("Cr53", 0.00031031500, "ao")
            self.inner_manifold.add_nuclide("Cr54", 0.00007724420, "ao")
            self.inner_manifold.add_nuclide("Mn55", 0.00003895800, "ao")
            self.inner_manifold.add_nuclide("Fe54", 0.00198804000, "ao")
            self.inner_manifold.add_nuclide("Fe56", 0.03120790000, "ao")
            self.inner_manifold.add_nuclide("Fe57", 0.00072072800, "ao")
            self.inner_manifold.add_nuclide("Fe58", 0.00009591550, "ao")
            self.inner_manifold.add_nuclide("Ta181", 0.00000453000, "ao")
            self.inner_manifold.add_nuclide("W182", 0.00006122300, "ao")
            self.inner_manifold.add_nuclide("W183", 0.00003306040, "ao")
            self.inner_manifold.add_nuclide("W184", 0.00007078760, "ao")
            self.inner_manifold.add_nuclide("W186", 0.00006568180, "ao")
            self.inner_manifold.set_density("atom/b-cm", 0.03822271948)
        # ------------------------------------------------------------------
        # Fertile material
        # ------------------------------------------------------------------

        self.fertile = openmc.Material(name='fertile', temperature=self.temp_k)

        if self.breeder_name in ['FLiBe','ARC']:
            if self.fertile_element == 'U':
                self.fertile.add_elements_from_formula('UF4','ao') # 'ao' here refers to 1:4 atomic ratio of U:F in UF4
                self.fertile.set_density('g/cm3', DENSITY_UF4) 
            elif self.fertile_element == 'Th':
                self.fertile.add_elements_from_formula('ThF4','ao') # 'ao' here refers to 1:4 atomic ratio of U:F in UF4
                self.fertile.set_density('g/cm3', DENSITY_ThF4) 

        elif self.breeder_name in ['LL', 'PB']:
           if self.fertile_element == 'U':
                # UO2 kernel (enriched uranium)
                uo2 = openmc.Material(name='UO2')
                uo2.add_elements_from_formula('UO2', enrichment=ENRICH_U,
                                            enrichment_target='U235',
                                            enrichment_type='wo')  
                uo2.set_density('g/cm3', 10.5)  # nominal UO2 density

                # SiC coating
                sic = openmc.Material(name='SiC')
                sic.add_elements_from_formula('SiC')
                sic.set_density('g/cm3', 3.2)

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
                biso.set_density('g/cm3', DENSITY_BISO)  

                self.fertile = biso

        #     elif self.fertile_element == 'Th':
        #         pass  # ThO2 pebbles


        # ------------------------------------------------------------------
        # Breeder and fertile material mixed in the blanket breeding regions 
        # ------------------------------------------------------------------

        if self.breeder_name  in ['FLiBe','ARC']:
            breeder_mass_frac, fertile_compound_mass_frac = calc_blanket_mass_fracs(self.fertile_bulk_density_kgm3, 
                                                                                    self.breeder_volume,
                                                                                    fertile_element=self.fertile_element, 
                                                                                    fertile_enrich=ENRICH_U, 
                                                                                    breeder_density_kgm3=DENSITY_FLIBE*1e3)
            self.blanket = openmc.Material.mix_materials([self.breeder, self.fertile], [breeder_mass_frac, fertile_compound_mass_frac], 'wo') # fractions in 'mix_materials' MUST add up to 1
            self.blanket.name, self.blanket.temperature = self.name, self.temp_k
            

        elif self.breeder_name in ['LL']:
            breeder_mass_frac, fertile_compound_mass_frac = calc_biso_blanket_mass_fracs(self.fertile_bulk_density_kgm3,
                                                                        self.breeder_volume,
                                                                        fertile_element=self.fertile_element,
                                                                        fertile_enrich=ENRICH_U,
                                                                        breeder_density_kgm3=DENSITY_LL*1e3)
            self.blanket = openmc.Material.mix_materials([self.breeder, self.fertile], [breeder_mass_frac, fertile_compound_mass_frac], 'wo')
            self.blanket.name, self.blanket.temperature = self.name, self.temp_k


        # ------------------------------------------------------------------
        # Add materials 
        # ------------------------------------------------------------------
        if self.breeder_name  in ['FLiBe','ARC']:
            self.materials = openmc.Materials([self.air, self.firstwall, self.firstwall_front, self.firstwall_cooling, self.backplate, self.blanket])
            # self.materials.export_to_xml(self.path)
        elif self.breeder_name == 'LL':
            self.materials = openmc.Materials([self.air, self.firstwall, self.firstwall_cooling, self.divider1, self.divider2, self.inner_manifold, self.blanket]) #check if need to add self.backplate

    @timer
    def geometry(self):

        # ------------------------------------------------------------------
        # Tokamak geometry parameters 
        # ------------------------------------------------------------------

        if self.breeder_name == 'FLiBe':
            self.R0, self.a, self.kappa, self.delta = FLIBE_R0, FLIBE_A, FLIBE_KAPPA, FLIBE_DELTA 
            d_fw  = FLIBE_FW_CM 
            d_st1 = d_fw  + FLIBE_ST1_CM
            d_br1 = d_st1 + FLIBE_BR1_CM
            d_st2 = d_br1 + FLIBE_ST2_CM
            d_br2 = d_st2 + FLIBE_BR2_CM
            d_st3 = d_br2 + FLIBE_ST3_CM

            self.extent_r = (self.R0 + self.a + d_st3)*1.1 # 110%
            self.extent_z = (self.kappa*self.a + d_st3)*1.1

        elif self.breeder_name == 'ARC':
            self.R0, self.a, self.kappa, self.delta = ARC_R0, ARC_A, ARC_KAPPA, ARC_DELTA 
            d_fw  = ARC_FW_CM 
            d_st1 = d_fw  + ARC_ST1_CM
            d_br1 = d_st1 + ARC_BR1_CM
            d_st2 = d_br1 + ARC_ST2_CM
            d_br2 = d_st2 + ARC_BR2_CM
            d_st3 = d_br2 + ARC_ST3_CM

            self.extent_r = (self.R0 + self.a + d_st3)*1.2 # 110%
            self.extent_z = (self.kappa*self.a + d_st3)*1.2
        
        elif self.breeder_name == 'LL':
            self.R0, self.a, self.kappa, self.delta = LL_R0, LL_A, LL_KAPPA, LL_DELTA

            #----Outboard thickness stack----
            d_fw_o   = LL_FW_O_CM
            d_fwf_o  = d_fw_o   + LL_FWF_O_CM
            d_fwc_o  = d_fwf_o  + LL_FWC_O_CM
            d_fwb_o  = d_fwc_o  + LL_FWB_O_CM
            d_br1_o  = d_fwb_o  + LL_BR1_O_CM
            d_d1_o   = d_br1_o  + LL_D1_O_CM
            d_br2_o  = d_d1_o   + LL_BR2_O_CM
            d_d2_o   = d_br2_o  + LL_D2_O_CM
            d_br3_o  = d_d2_o   + LL_BR3_O_CM
            d_im_o  = d_br3_o  + LL_IM_O_CM

            #----Inboard thickness stack----
            d_fw_i   = LL_FW_I_CM
            d_fwf_i  = d_fw_i   + LL_FWF_I_CM
            d_fwc_i  = d_fwf_i  + LL_FWC_I_CM
            d_fwb_i  = d_fwc_i  + LL_FWB_I_CM
            d_br1_i  = d_fwb_i  + LL_BR1_I_CM
            d_d1_i   = d_br1_i  + LL_D1_I_CM
            d_br2_i  = d_d1_i   + LL_BR2_I_CM
            d_im_i  = d_br2_i  + LL_IM_I_CM

            # Take the larger of outboard vs inboard stack for safety margin (check, i think it could be reversed -ezoccoli)
            d_outboard = d_im_o
            d_inboard  = d_im_i
            d_total    = max(d_outboard, d_inboard)

            self.extent_r = (self.R0 + self.a + d_total) * 1.2   # 110% radial extent
            self.extent_z = (self.kappa*self.a + d_total) * 1.2  # 110% vertical extent

        # ------------------------------------------------------------------
        # Surfaces 
        # ------------------------------------------------------------------
        if self.breeder_name  in ['FLiBe','ARC']:
            # Arrays of a (minor r) and z points
            points_vc  =  miller_model(self.R0, self.a, self.kappa, self.delta)         # coords around vacuum chamber
            points_fw  = miller_offset(self.R0, self.a, self.kappa, self.delta, d_fw)   # outer coords around first wall
            points_st1 = miller_offset(self.R0, self.a, self.kappa, self.delta, d_st1)  # outer coords around structural region 1
            points_br1 = miller_offset(self.R0, self.a, self.kappa, self.delta, d_br1)  # outer coords around breeding region 1
            points_st2 = miller_offset(self.R0, self.a, self.kappa, self.delta, d_st2)  # outer coords around structural region 2
            points_br2 = miller_offset(self.R0, self.a, self.kappa, self.delta, d_br2)  # outer coords around breeding region 2
            points_st3 = miller_offset(self.R0, self.a, self.kappa, self.delta, d_st3)  # outer coords around structural region 3

            # Create OpenMC surfaces and STORE THEM AS CLASS ATTRIBUTES for tally use
            self.surface_vc  = openmc.model.Polygon(points_vc , basis='rz')  # Plasma-facing surface (inner FW)
            self.surface_fw  = openmc.model.Polygon(points_fw , basis='rz')  # Outer surface of first wall
            self.surface_st1 = openmc.model.Polygon(points_st1, basis='rz')
            self.surface_br1 = openmc.model.Polygon(points_br1, basis='rz')
            self.surface_st2 = openmc.model.Polygon(points_st2, basis='rz')
            self.surface_br2 = openmc.model.Polygon(points_br2, basis='rz')
            self.surface_st3 = openmc.model.Polygon(points_st3, basis='rz')

            # Create OpenMC surfaces
            surfaces = [self.surface_vc, self.surface_fw, self.surface_st1, 
                        self.surface_br1, self.surface_st2, self.surface_br2, self.surface_st3]

            # Add boundary surfaces
            # r_max = max([point[0] for point in points_st3]) + 1.0  # 1 cm beyond outermost surface
            # z_max = max([point[1] for point in points_st3]) + 1.0
            # z_min = min([point[1] for point in points_st3]) - 1.0

            outer_cylinder = openmc.ZCylinder(r=self.extent_r, boundary_type='vacuum')
            top_plane      = openmc.ZPlane(z0=self.extent_z, boundary_type='vacuum')
            bottom_plane   = openmc.ZPlane(z0=-self.extent_z, boundary_type='vacuum')
        
        elif self.breeder_name == 'LL':
             # ---- Outboard offsets ----
            points_fw_o   = miller_offset_split(self.R0, self.a, self.kappa, self.delta, 0, d_fw_o)
            points_fwf_o  = miller_offset_split(self.R0, self.a, self.kappa, self.delta, 0, d_fwf_o)
            points_fwc_o  = miller_offset_split(self.R0, self.a, self.kappa, self.delta, 0, d_fwc_o)
            points_fwb_o  = miller_offset_split(self.R0, self.a, self.kappa, self.delta, 0, d_fwb_o)
            points_br1_o  = miller_offset_split(self.R0, self.a, self.kappa, self.delta, 0, d_br1_o)
            points_d1_o   = miller_offset_split(self.R0, self.a, self.kappa, self.delta, 0, d_d1_o)
            points_br2_o  = miller_offset_split(self.R0, self.a, self.kappa, self.delta, 0, d_br2_o)
            points_d2_o   = miller_offset_split(self.R0, self.a, self.kappa, self.delta, 0, d_d2_o)
            points_br3_o  = miller_offset_split(self.R0, self.a, self.kappa, self.delta, 0, d_br3_o)
            points_im_o   = miller_offset_split(self.R0, self.a, self.kappa, self.delta, 0, d_im_o)

            # ---- Inboard offsets ----
            points_fw_i   = miller_offset_split(self.R0, self.a, self.kappa, self.delta, d_fw_i, 0)
            points_fwf_i  = miller_offset_split(self.R0, self.a, self.kappa, self.delta, d_fwf_i, 0)
            points_fwc_i  = miller_offset_split(self.R0, self.a, self.kappa, self.delta, d_fwc_i, 0)
            points_fwb_i  = miller_offset_split(self.R0, self.a, self.kappa, self.delta, d_fwb_i, 0)
            points_br1_i  = miller_offset_split(self.R0, self.a, self.kappa, self.delta, d_br1_i, 0)
            points_d1_i   = miller_offset_split(self.R0, self.a, self.kappa, self.delta, d_d1_i, 0)
            points_br2_i  = miller_offset_split(self.R0, self.a, self.kappa, self.delta, d_br2_i, 0)
            points_im_i   = miller_offset_split(self.R0, self.a, self.kappa, self.delta, d_im_i, 0)

             # ---- Store as surfaces ----
            self.surface_vc = openmc.model.Polygon(points_vc , basis='rz')

            # Outboard surfaces
            self.surface_fw_o   = openmc.model.Polygon(points_fw_o, basis='rz')
            self.surface_fwf_o  = openmc.model.Polygon(points_fwf_o, basis='rz')
            self.surface_fwc_o  = openmc.model.Polygon(points_fwc_o, basis='rz')
            self.surface_fwb_o  = openmc.model.Polygon(points_fwb_o, basis='rz')
            self.surface_br1_o  = openmc.model.Polygon(points_br1_o, basis='rz')
            self.surface_d1_o   = openmc.model.Polygon(points_d1_o,  basis='rz')
            self.surface_br2_o  = openmc.model.Polygon(points_br2_o, basis='rz')
            self.surface_d2_o   = openmc.model.Polygon(points_d2_o,  basis='rz')
            self.surface_br3_o  = openmc.model.Polygon(points_br3_o, basis='rz')
            self.surface_im_o   = openmc.model.Polygon(points_im_o,  basis='rz')

            # Inboard surfaces
            self.surface_fw_i   = openmc.model.Polygon(points_fw_i, basis='rz')
            self.surface_fwf_i  = openmc.model.Polygon(points_fwf_i, basis='rz')
            self.surface_fwc_i  = openmc.model.Polygon(points_fwc_i, basis='rz')
            self.surface_fwb_i  = openmc.model.Polygon(points_fwb_i, basis='rz')
            self.surface_br1_i  = openmc.model.Polygon(points_br1_i, basis='rz')
            self.surface_d1_i   = openmc.model.Polygon(points_d1_i,  basis='rz')
            self.surface_br2_i  = openmc.model.Polygon(points_br2_i, basis='rz')
            self.surface_im_i   = openmc.model.Polygon(points_im_i,  basis='rz')

            
        # ------------------------------------------------------------------
        # Cells 
        # ------------------------------------------------------------------

        cell_vc   = openmc.Cell(cell_id=10, region= -surfaces[0])
        cell_vc.importance = {'neutron':1}
        cell_fw   = openmc.Cell(cell_id=11, region= +surfaces[0] & -surfaces[1] , fill=self.firstwall) 
        cell_st1  = openmc.Cell(cell_id=21, region= +surfaces[1] & -surfaces[2] , fill=self.structure)
        cell_br1  = openmc.Cell(cell_id=31, region= +surfaces[2] & -surfaces[3] , fill=self.blanket)
        cell_st2  = openmc.Cell(cell_id=22, region= +surfaces[3] & -surfaces[4] , fill=self.structure)
        cell_br2  = openmc.Cell(cell_id=32, region= +surfaces[4] & -surfaces[5] , fill=self.blanket)
        cell_st3  = openmc.Cell(cell_id=23, region= +surfaces[5] & -surfaces[6] , fill=self.structure) 

        # Surrounding air cell with proper boundaries (otherwise causes error with just Polygons)
        # cell_air = openmc.Cell(cell_id=99, region= +surfaces[6] & -outer_cylinder & +bottom_plane & -top_plane, fill=self.air)
        cell_void = openmc.Cell(cell_id=99, region= +surfaces[6] & -outer_cylinder & +bottom_plane & -top_plane)
        cell_void.importance = {'neutron':0}

        self.cells = [cell_vc, cell_fw, cell_st1, cell_br1, cell_st2, cell_br2, cell_st3, cell_void]
        self.geometry = openmc.Geometry(openmc.Universe(cells=self.cells))
        # self.geometry.export_to_xml(self.path)


    @timer
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
        flux_energy_tally = self.make_tally('flux', ['flux'], filters=[cell_filter, energy_filter])

        # Uranium reaction rates
        U_tally        = self.make_tally('U rxn rates',          ['(n,gamma)', 'fission', 'elastic'], filters=[cell_filter], nuclides=['U238', 'U235'])
        U_energy_tally = self.make_tally('U rxn rates spectrum', ['(n,gamma)', 'fission', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['U238', 'U235'])

        # Lithium reaction rates
        Li_tally        = self.make_tally('Li rxn rates',          ['(n,gamma)', '(n,Xt)', 'elastic'], filters=[cell_filter], nuclides=['Li6', 'Li7'])
        Li_energy_tally = self.make_tally('Li rxn rates spectrum', ['(n,gamma)', '(n,Xt)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['Li6', 'Li7'])

        # Fluorine reaction rates
        F_tally        = self.make_tally('F rxn rates', ['(n,gamma)', 'elastic'], filters=[cell_filter], nuclides=['F19'])
        F_energy_tally = self.make_tally('F rxn rates spectrum', ['(n,gamma)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['F19'])

        # Beryllium reaction rates
        Be_tally        = self.make_tally('Be rxn rates', ['(n,gamma)', '(n,2n)', 'elastic'], filters=[cell_filter], nuclides=['Be9'])
        Be_energy_tally = self.make_tally('Be rxn rates spectrum', ['(n,gamma)', '(n,2n)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['Be9'])

        self.tallies.extend([flux_tally, U_tally, Li_tally, F_tally, Be_tally])
        self.tallies.extend([flux_energy_tally, U_energy_tally, Li_energy_tally, F_energy_tally, Be_energy_tally])


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
        self.settings.run_mode = 'fixed source'
        self.settings.particles = int(1e5)
        self.settings.batches = 10


    # @timer
    def openmc(self):
        self.materials.cross_sections = set_xs_path()
        self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
        self.model.export_to_model_xml(self.path)
        self.model.run(cwd=self.path)


    # @timer
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
        

        # ---------------------------
        # Process volume calc results
        # ---------------------------
        vol_results = openmc.VolumeCalculation.from_hdf5(f"{self.path}/volume_1.h5")
        
        print(f"{Colors.GREEN}Stochastic Volume Calculation Results:{Colors.END}")

        vol_dict = {}
        for k, v in vol_results.volumes.items():
            vol_dict[k] = (v.nominal_value/1e6, v.std_dev/1e6)
        # vol_dict['sum'] = ( sum(v.nominal_value/1e6 for v in vol_results.volumes.values()),
        #                     sum(v.std_dev/1e6 for v in vol_results.volumes.values())       )

        df = pd.DataFrame.from_dict(vol_dict, orient='index', columns=['volume_m3', 'stdev_m3'])
        df.reset_index(inplace=True)
        df.rename(columns={'index': 'cells'}, inplace=True)
        df.to_csv(f'{self.path}/volume_1.csv',index=False)

        print("CSV file 'volume_1.csv' created successfully")
        print(df)

        """
        vol_results.volumes is a dictionary that looks like:
          {2: "909560.9858392784+/-140318.8308700002",
           3: "5543990.770829888+/-346055.6196989608"}
        so we will rewrite it to separate the value from the standev, like this:
          {2: (909560.9858392784, 140318.8308700002),
           3: (5543990.770829888, 346055.6196989608)}

        But also whoever designed the output of VolueCalculation.volumes() 
        honestly needs to WRITE in the documentation that the data is stored/written 
        using the uncertainties package. I was going crazy using all sorts of 
        regex match patterns to try to understand why splitting the vol_results.volumes.items()
        into k, v was causing issues vs. what v looked like when I was printing it. 
        It was the uncertainties package changing the formatting of v when you go print it. 
          --ppark  2025-09-20
        """


    # @timer
    def plot(self):
        # ---------- Plots.xml + generate images ----------
        # XY toroidal slice
        xy = openmc.Plot()
        xy.filename = "tokamak_xy"
        xy.basis = "xy"
        xy.width  = (1800, 1800)
        xy.pixels = (18000, 18000)
        xy.color_by = "material"
        xy.colors = {
            self.firstwall: (255, 0, 0),    # Red for first wall
            self.structure: (0, 255, 0),    # Green for structure
            self.blanket: (0, 0, 255),    # Blue for breeder
            # Void regions will be white by default
        }

        # XZ poloidal slice
        xz = openmc.Plot()
        xz.filename = "tokamak_xz"
        xz.basis = "xz"
        xz.width  = (1800, 1000)
        xz.pixels = (18000, 10000)
        xz.color_by = "material"
        xz.colors = {
            self.firstwall: (255, 0, 0),    # Red for first wall
            self.structure: (0, 255, 0),    # Green for structure
            self.blanket: (0, 0, 255),    # Blue for breeder
            # Void regions will be white by default
        }

        openmc.Plots([xy, xz]).export_to_xml(self.path)
        openmc.plot_geometry(path_input=self.path)  # writes tokamak_rz.ppm and tokamak_xz.ppm
        print("Done. Files: tokamak_rz.ppm, tokamak_xz.ppm")

                
