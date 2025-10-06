import numpy as np
import openmc

from Python.reactor    import *
from Python.parameters import *
from Python.utilities  import *


class LL(Reactor):

    def initialize(self):

        self.temp_k          = TEMP_K
        self.breeder_name    = 'PB'
        self.breeder_density = DENSITY_PB # g/cm^3
        self.breeder_enrich  = ENRICH_PB  # wt% 
        self.breeder_volume  = PB_BR_VOL  

        # Name file based on reactor config - should come out to smth like: tallies_FLiBe_U010kgm3_Li7.5_900K
        self.name = f"{self.run_type}_{self.breeder_name}_{self.fertile_element}{self.fertile_bulk_density_kgm3:06.2f}kgm3_Li{self.breeder_enrich:04.1f}_{self.temp_k}K"         
        self.path = f"./OpenMC/{self.name}"
        
        os.makedirs(self.path, exist_ok=True)

        start_msg = f"\n======== {self.breeder_name} reactor - {self.fertile_element} {self.fertile_bulk_density_kgm3:6.2f} kg/m3 - {self.breeder_enrich:4.1f}%-enriched Li - {self.temp_k} K ========"
        print(f"{Colors.MAGENTA}{start_msg}{Colors.END}")


    def materials(self):

        # ------------------------------------------------------------------
        # Air
        # ------------------------------------------------------------------
        # self.air = openmc.Material(name='air')
        # self.air.set_density('g/cm3', 0.001225)

        # # Atom fractions for dry air
        # self.air.add_element('N', 0.78084, percent_type='ao')   # Nitrogen
        # self.air.add_element('O', 0.20946, percent_type='ao')   # Oxygen
        # self.air.add_element('Ar', 0.00934, percent_type='ao')  # Argon
        # self.air.add_element('C', 0.00036, percent_type='ao')   # Carbon from CO2

        # ------------------------------------------------------------------
        # First wall 
        # ------------------------------------------------------------------

        self.firstwall = openmc.Material(name='firstwall', temperature=self.temp_k)
        self.firstwall.depletable = False
        self.firstwall.add_element('W',1)
        self.firstwall.set_density('g/cm3',19.3)
        # The original first wall specs we were using from Ball 25 is 99.9969 wt% W 
        # and the rest O, N ,C, Na, K, Al, Ca, Cr, Cu, Fe impurities...
        # I think we can just set it to be W to save compute time loading nuclide data
        

        # ------------------------------------------------------------------
        # Structure 
        # ------------------------------------------------------------------

        self.eurofer = openmc.Material(name='Eurofer', temperature=self.temp_k)
        self.eurofer.depletable = False
        self.eurofer.set_density('g/cm3', 7.8)
        self.eurofer.add_element('Fe', 89.0026, percent_type='wo')
        self.eurofer.add_element('B', 0.001, percent_type='wo')
        self.eurofer.add_element('C', 0.1049, percent_type='wo')
        self.eurofer.add_element('N', 0.04, percent_type='wo')
        self.eurofer.add_element('O', 0.001, percent_type='wo')
        self.eurofer.add_element('Al', 0.004, percent_type='wo')
        self.eurofer.add_element('Si', 0.026, percent_type='wo')
        self.eurofer.add_element('P', 0.002, percent_type='wo')
        self.eurofer.add_element('S', 0.003, percent_type='wo')
        self.eurofer.add_element('Ti', 0.001, percent_type='wo')
        self.eurofer.add_element('V', 0.01963, percent_type='wo')
        self.eurofer.add_element('Cr', 9.00, percent_type='wo')
        self.eurofer.add_element('Mn', 0.55, percent_type='wo')
        self.eurofer.add_element('Co', 0.005, percent_type='wo')
        self.eurofer.add_element('Ni', 0.01, percent_type='wo')
        self.eurofer.add_element('Cu', 0.003, percent_type='wo')
        self.eurofer.add_element('Nb', 0.005, percent_type='wo')
        self.eurofer.add_element('Mo', 0.003, percent_type='wo')
        self.eurofer.add_element('Ta', 0.12, percent_type='wo')
        self.eurofer.add_element('W', 1.0987, percent_type='wo')


        # ------------------------------------------------------------------
        # Breeder material
        # ------------------------------------------------------------------

        self.breeder = openmc.Material(name='breeder', temperature=self.temp_k)
        self.lithium_ceramic = openmc.Material(name='LithiumCeramic') 
        self.lithium_ceramic.set_density('g/cm3', 2.42) 
        self.lithium_ceramic.add_element('Li', 22.415, percent_type='wo', enrichment_target='Li6', enrichment_type='ao', enrichment=ENRICH_PB) 
        self.lithium_ceramic.add_element('Si', 24.077, percent_type='wo') 
        self.lithium_ceramic.add_element('O', 53.39, percent_type='wo') 
        self.lithium_ceramic.add_element('Al', 0.003, percent_type='wo') 
        self.lithium_ceramic.add_element('C', 0.1, percent_type='wo') 
        self.lithium_ceramic.add_element('Ca', 0.003, percent_type='wo') 
        self.lithium_ceramic.add_element('Co', 0.0002, percent_type='wo') 
        self.lithium_ceramic.add_element('Cr', 0.0001, percent_type='wo') 
        self.lithium_ceramic.add_element('Cu', 0.0001, percent_type='wo') 
        self.lithium_ceramic.add_element('Fe', 0.0005, percent_type='wo') 
        self.lithium_ceramic.add_element('K', 0.0001, percent_type='wo') 
        self.lithium_ceramic.add_element('Mg', 0.0005, percent_type='wo') 
        self.lithium_ceramic.add_element('Mn', 0.0001, percent_type='wo') 
        self.lithium_ceramic.add_element('Pt', 0.009, percent_type='wo') 
        self.lithium_ceramic.add_element('Na', 0.002, percent_type='wo') 
        self.lithium_ceramic.add_element('Ni', 0.0002, percent_type='wo') 
        self.lithium_ceramic.add_element('Ti', 0.0005, percent_type='wo') 
        self.lithium_ceramic.add_element('Zn', 0.0002, percent_type='wo') 
        self.lithium_ceramic.add_element('Zr', 0.0001, percent_type='wo') 
        # Beryllium 
        self.beryllium = openmc.Material(name='Beryllium') 
        self.beryllium.set_density('g/cm3', 1.85) 
        self.beryllium.add_element('Be', 98.749, percent_type='wo') 
        self.beryllium.add_element('O', 0.9, percent_type='wo') 
        self.beryllium.add_element('Al', 0.09, percent_type='wo') 
        self.beryllium.add_element('Fe', 0.1, percent_type='wo') 
        self.beryllium.add_element('Mg', 0.08, percent_type='wo') 
        self.beryllium.add_element('Si', 0.06, percent_type='wo') 
        self.beryllium.add_element('Mn', 0.01, percent_type='wo') 
        self.beryllium.add_element('U', 0.01, percent_type='wo') 
        self.beryllium.add_element('Co', 0.001, percent_type='wo') 
        self.beryllium.add_element('Cu', 0.001, percent_type='wo') 
        self.beryllium.add_element('K', 0.001, percent_type='wo') 
        self.beryllium.add_element('Mg', 0.0005, percent_type='wo') 
        self.beryllium.add_element('Mn', 0.0005, percent_type='wo') 
        self.beryllium.add_element('Na', 0.001, percent_type='wo') 
        self.beryllium.add_element('Nb', 0.001, percent_type='wo') 
        self.beryllium.add_element('Ni', 0.0005, percent_type='wo') 
        self.beryllium.add_element('Pb', 0.0005, percent_type='wo') 
        self.beryllium.add_element('Ta', 0.002, percent_type='wo') 

        self.he = openmc.Material(name='Helium') 
        self.he.set_density('g/cm3', 0.004279) #Helium density at 900k ~80bar 
        self.he.add_element('He', 1) 
        
        vf_lic = 0.1304 
        vf_be = 0.379 
        vf_euro = 0.1176 #volume fractions from EU Activation Analysis Table 2 
        vf_he = 1 - (vf_lic + vf_be + vf_euro) 
        self.breeder = openmc.Material.mix_materials([self.lithium_ceramic, self.he, self.eurofer, self.beryllium],[vf_lic, vf_he, vf_euro, vf_be],'vo') 
        self.breeder_density = self.breeder.get_density()

        # ------------------------------------------------------------------
        # Fertile material
        # ------------------------------------------------------------------

        self.fertile = openmc.Material(name='fertile', temperature=self.temp_k)

        if self.fertile_element == 'U':
            # UO2 kernel (enriched uranium)
            uo2 = openmc.Material(name='UO2')
            uo2.add_elements_from_formula('UO2', enrichment=ENRICH_U)
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

        elif self.fertile_element == 'Th': 
            #ThO2 kernel 
            tho2 = openmc.Material(name='ThO2') 
            tho2.add_elements_from_formula('ThO2') 
            tho2.set_density('g/cm3', 10.0) # nominal ThO2 density 

            # SiC coating 
            sic = openmc.Material(name='SiC') 
            sic.add_elements_from_formula('SiC') 
            sic.set_density('g/cm3', 3.2) 
            
            # Geometry for BISO 
            r_tho2 = 400e-4 # 400 µm kernel radius 
            r_sic = 500e-4 # 500 µm outer radius (100 µm coating) 
            V_biso_particle = (4 / 3) * np.pi * (r_sic)**3 
            V_tho2_in_biso = (4 / 3) * np.pi * (r_tho2)**3 
            Vf_tho2_in_biso = V_tho2_in_biso / V_biso_particle 
            Vf_sic_in_biso = 1.0 - Vf_tho2_in_biso 
            biso = openmc.Material.mix_materials([tho2, sic], [Vf_tho2_in_biso, Vf_sic_in_biso], 'vo') 
            biso.set_density('g/cm3', 10.0*Vf_tho2_in_biso + 3.2*Vf_sic_in_biso) 
            self.fertile = biso


        # ------------------------------------------------------------------
        # Breeder and fertile material mixed in the blanket breeding regions 
        # ------------------------------------------------------------------
        '''Emma changed fertile_enrich=ENRICH_U to 
        ENRICH_U if self.fertile_element == 'U' else 0.0 is this right?? 
        Did this even need to be fixed''' 
        breeder_mass_frac, fertile_compound_mass_frac = calc_biso_blanket_mass_fracs(self.fertile_bulk_density_kgm3,
                                                                    self.breeder_volume,
                                                                    fertile_element=self.fertile_element,
                                                                    fertile_enrich=ENRICH_U if self.fertile_element == 'U' else 0.0,
                                                                    breeder_density_kgm3=self.breeder_density*1e3)
        self.blanket = openmc.Material.mix_materials([self.breeder, self.fertile], [breeder_mass_frac, fertile_compound_mass_frac], 'wo')
        self.blanket.name, self.blanket.temperature = self.name, self.temp_k


        # ------------------------------------------------------------------
        # Add materials 
        # ------------------------------------------------------------------

        self.materials = openmc.Materials([self.firstwall, self.eurofer, self.blanket]) 


    def geometry(self):

        # ------------------------------------------------------------------
        # Tokamak geometry parameters 
        # ------------------------------------------------------------------

        self.R0, self.a, self.kappa, self.delta = PB_R0, PB_A, PB_KAPPA, PB_DELTA

        d_fw  = PB_FW_CM 
        d_st1 = d_fw  + PB_ST1_CM
        d_br1 = d_st1 + PB_BR1_CM
        d_st2 = d_br1 + PB_ST2_CM
        
        d_br1_o = d_st1    + PB_BR1_O_CM   # only on outboard blanket
        d_st2_o = d_br1_o  + PB_ST2_CM     # only on outboard blanket
        
        
        self.extent_r = (self.R0 + self.a + d_st2_o)*1.2 # 110%
        self.extent_z = (self.kappa*self.a + d_st2_o)*1.2

        # ------------------------------------------------------------------
        # Surfaces 
        # ------------------------------------------------------------------

        points_vc    =  miller_model(self.R0, self.a, self.kappa, self.delta)
        points_fw    = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_fw)
        points_st1   = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st1)
        points_br1   = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br1)
        points_st2   = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st2)
        points_br1_o = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br1_o)  # outboard blanket
        points_st2_o = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st2_o)  # outboard structure
        
        # Create OpenMC surfaces
        self.surface_vc    = openmc.model.Polygon(points_vc , basis='rz')
        self.surface_fw    = openmc.model.Polygon(points_fw, basis='rz')
        self.surface_st1   = openmc.model.Polygon(points_st1, basis='rz')
        self.surface_br1   = openmc.model.Polygon(points_br1, basis='rz')
        self.surface_st2   = openmc.model.Polygon(points_st2, basis='rz')
        self.surface_br1_o = openmc.model.Polygon(points_br1_o,  basis='rz') # outboard blanket
        self.surface_st2_o = openmc.model.Polygon(points_st2_o, basis='rz')  # outboard structure

        dividing_cylinder = openmc.ZCylinder(r=self.R0)
        outer_cylinder    = openmc.ZCylinder(r=self.extent_r, boundary_type='vacuum')
        top_plane         = openmc.ZPlane(z0=self.extent_z, boundary_type='vacuum')
        bottom_plane      = openmc.ZPlane(z0=-self.extent_z, boundary_type='vacuum')
        
        # ------------------------------------------------------------------
        # Cells | 10: vc, fw | 2X: structure | 3X: breeding
        # ------------------------------------------------------------------

        cell_vc   = openmc.Cell(cell_id=10, region= -self.surface_vc)
        cell_vc.importance = {'neutron':1}

        cell_fw    = openmc.Cell(cell_id=11, region= +self.surface_vc  & -self.surface_fw,  fill=self.firstwall)
        cell_st1   = openmc.Cell(cell_id=21, region= +self.surface_fw  & -self.surface_st1, fill=self.eurofer)
        cell_br1   = openmc.Cell(cell_id=31, region= -dividing_cylinder & +self.surface_st1 & -self.surface_br1, fill=self.blanket)
        cell_st2   = openmc.Cell(cell_id=22, region= -dividing_cylinder & +self.surface_br1 & -self.surface_st2, fill=self.eurofer)
        cell_br1_o = openmc.Cell(cell_id=32, region= +dividing_cylinder & +self.surface_st1 & -self.surface_br1_o, fill=self.blanket)
        cell_st2_o = openmc.Cell(cell_id=23, region= +dividing_cylinder & +self.surface_br1_o & -self.surface_st2_o,  fill=self.eurofer)
        
        # Void cells
        cell_void_i = openmc.Cell(cell_id=99, region= +self.surface_st2 & -dividing_cylinder & -outer_cylinder & +bottom_plane & -top_plane)
        cell_void_o = openmc.Cell(cell_id=98, region= +self.surface_st2_o & +dividing_cylinder & -outer_cylinder & +bottom_plane & -top_plane)
        cell_void_i.importance = {'neutron': 0}
        cell_void_o.importance = {'neutron': 0}

        self.cells = [
            cell_vc,
            cell_fw, cell_st1, cell_br1, cell_st2, cell_br1_o, cell_st2_o,
            cell_void_i, cell_void_o
        ]
    
        self.geometry = openmc.Geometry(openmc.Universe(cells=self.cells))