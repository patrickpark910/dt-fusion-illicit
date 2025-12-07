import numpy as np
import openmc

from Python.reactor    import *
from Python.parameters import *
from Python.utilities  import *


class HCPB(Reactor):

    def initialize(self):

        self.temp_k          = TEMP_K
        self.breeder_name    = 'HCPB'
        self.breeder_enrich  = ENRICH_HCPB  # at% 
        self.breeder_volume  = HCPB_BR_VOL  

        # Name file based on reactor config - should come out to smth like: tallies_FLiBe_U010kgm3_Li7.5_900K
        self.name = f"{self.run_type}_{self.breeder_name}_{self.temp_k}K_Li{self.breeder_enrich:04.1f}_{self.fertile_element}{self.fertile_bulk_density_kgm3:06.2f}kgm3"         
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
        self.firstwall.set_density('g/cm3',19.0)
        self.firstwall.depletable = False
        # self.firstwall.add_element('O',5/1e6,percent_type='wo')
        # self.firstwall.add_element('N',5/1e6,percent_type='wo')
        # self.firstwall.add_element('C',5/1e6,percent_type='wo')
        # self.firstwall.add_element('Na',4/1e6,percent_type='wo')
        # self.firstwall.add_element('K',2.5/1e6,percent_type='wo')
        # self.firstwall.add_element('Al',3/1e6,percent_type='wo')
        # self.firstwall.add_element('Ca',0.5/1e6,percent_type='wo')
        # self.firstwall.add_element('Cr',0.5/1e6,percent_type='wo')
        # self.firstwall.add_element('Cu',0.5/1e6,percent_type='wo')
        # self.firstwall.add_element('Fe',5/1e6,percent_type='wo')
        # self.firstwall.add_element('W',1-(5+5+5+4+2.5+3+0.5+0.5+0.5+5)/1e6,percent_type='wo')
        self.firstwall.add_element('W',1)


        # ------------------------------------------------------------------
        # Structure 
        # - Removed impurities to nominal specs per Gaganidze (2018) 
        #   "EUROFER 97 Database and Mat Prop Handbook"
        # ------------------------------------------------------------------

        self.structure = openmc.Material(name='Eurofer', temperature=self.temp_k)
        self.structure.depletable = False
        self.structure.set_density('g/cm3', 7.8)
        self.structure.add_element('Fe', 89.36, percent_type='wo')
        self.structure.add_element('C' ,  0.11, percent_type='wo')
        self.structure.add_element('Cr',  9.00, percent_type='wo')
        self.structure.add_element('W' ,  1.10, percent_type='wo')
        self.structure.add_element('Mn',  0.40, percent_type='wo')
        self.structure.add_element('N' ,  0.03, percent_type='wo')

        # Original Eurofer specs from Lu (2017) "HCPB Analysis"
        # self.structure.add_element('Fe', 89.0026, percent_type='wo')
        # self.structure.add_element('B', 0.001, percent_type='wo')
        # self.structure.add_element('C', 0.1049, percent_type='wo')
        # self.structure.add_element('N', 0.04, percent_type='wo')
        # self.structure.add_element('O', 0.001, percent_type='wo')
        # self.structure.add_element('Al', 0.004, percent_type='wo')
        # self.structure.add_element('Si', 0.026, percent_type='wo')
        # self.structure.add_element('P', 0.002, percent_type='wo')
        # self.structure.add_element('S', 0.003, percent_type='wo')
        # self.structure.add_element('Ti', 0.001, percent_type='wo')
        # self.structure.add_element('V', 0.01963, percent_type='wo')
        # self.structure.add_element('Cr', 9.00, percent_type='wo')
        # self.structure.add_element('Mn', 0.55, percent_type='wo')
        # self.structure.add_element('Co', 0.005, percent_type='wo')
        # self.structure.add_element('Ni', 0.01, percent_type='wo')
        # self.structure.add_element('Cu', 0.003, percent_type='wo')
        # self.structure.add_element('Nb', 0.005, percent_type='wo')
        # self.structure.add_element('Mo', 0.003, percent_type='wo')
        # # self.structure.add_element('Ta', 0.12, percent_type='wo') # no cross sections
        # self.structure.add_element('W', 1.0987, percent_type='wo')


        # ------------------------------------------------------------------
        # Breeder material
        # ------------------------------------------------------------------

        li4sio4 = openmc.Material(name='Li4SiO4', temperature=self.temp_k) 
        li4sio4.set_density('g/cm3', 2.17)  # normalized for ceramic porosity and 900 K (pure, room temp g/cm3 = 2.42)
        li4sio4.add_elements_from_formula('Li4SiO4', enrichment_target='Li6', enrichment_type='ao', enrichment=ENRICH_HCPB) 
        # self.lithium_ceramic.add_elements('Li', 22.415, percent_type='wo', enrichment_target='Li6', enrichment_type='ao', enrichment=ENRICH_HCPB) 
        # self.lithium_ceramic.add_element('Si', 24.077, percent_type='wo') 
        # self.lithium_ceramic.add_element('O', 53.39, percent_type='wo') 
        # self.lithium_ceramic.add_element('Al', 0.003, percent_type='wo') 
        # self.lithium_ceramic.add_element('C', 0.1, percent_type='wo') 
        # self.lithium_ceramic.add_element('Ca', 0.003, percent_type='wo') 
        # self.lithium_ceramic.add_element('Co', 0.0002, percent_type='wo') 
        # self.lithium_ceramic.add_element('Cr', 0.0001, percent_type='wo') 
        # self.lithium_ceramic.add_element('Cu', 0.0001, percent_type='wo') 
        # self.lithium_ceramic.add_element('Fe', 0.0005, percent_type='wo') 
        # self.lithium_ceramic.add_element('K', 0.0001, percent_type='wo') 
        # self.lithium_ceramic.add_element('Mg', 0.0005, percent_type='wo') 
        # self.lithium_ceramic.add_element('Mn', 0.0001, percent_type='wo') 
        # self.lithium_ceramic.add_element('Pt', 0.009, percent_type='wo') 
        # self.lithium_ceramic.add_element('Na', 0.002, percent_type='wo') 
        # self.lithium_ceramic.add_element('Ni', 0.0002, percent_type='wo') 
        # self.lithium_ceramic.add_element('Ti', 0.0005, percent_type='wo') 
        # self.lithium_ceramic.add_element('Zn', 0.0002, percent_type='wo') 
        # self.lithium_ceramic.add_element('Zr', 0.0001, percent_type='wo') 
        
        # Beryllium 
        be = openmc.Material(name='Beryllium') 
        be.set_density('g/cm3', 1.80)  # normalized from 1.85 for 900 K
        # self.beryllium.add_element('Be', 98.749, percent_type='wo') 
        # self.beryllium.add_element('O', 0.9, percent_type='wo') 
        # self.beryllium.add_element('Al', 0.09, percent_type='wo') 
        # self.beryllium.add_element('Fe', 0.1, percent_type='wo') 
        # self.beryllium.add_element('Mg', 0.08, percent_type='wo') 
        # self.beryllium.add_element('Si', 0.06, percent_type='wo') 
        # self.beryllium.add_element('Mn', 0.01, percent_type='wo') 
        # self.beryllium.add_element('U', 0.01, percent_type='wo') 
        # self.beryllium.add_element('Co', 0.001, percent_type='wo') 
        # self.beryllium.add_element('Cu', 0.001, percent_type='wo') 
        # self.beryllium.add_element('K', 0.001, percent_type='wo') 
        # self.beryllium.add_element('Mg', 0.0005, percent_type='wo') 
        # self.beryllium.add_element('Mn', 0.0005, percent_type='wo') 
        # self.beryllium.add_element('Na', 0.001, percent_type='wo') 
        # self.beryllium.add_element('Nb', 0.001, percent_type='wo') 
        # self.beryllium.add_element('Ni', 0.0005, percent_type='wo') 
        # self.beryllium.add_element('Pb', 0.0005, percent_type='wo') 
        # self.beryllium.add_element('Ta', 0.002, percent_type='wo') 

        # Helium-4
        he = openmc.Material(name='Helium') 
        he.set_density('g/cm3', 0.004279) # Helium density at 900 K ~80 bar 
        he.add_element('He', 1) 
        

        # ------------------------------------------------------------------
        # Fertile material - BISO Particle
        # ------------------------------------------------------------------

        if self.fertile_element == 'U':
            kernel = openmc.Material(name='UO2')
            kernel.add_elements_from_formula('UO2', enrichment=ENRICH_U)
            kernel.set_density('g/cm3', DENSITY_UO2)  

        elif self.fertile_element == 'Th': 
            kernel = openmc.Material(name='ThO2') 
            kernel.add_elements_from_formula('ThO2') 
            kernel.set_density('g/cm3', DENSITY_ThO2) 

        # SiC coating
        sic = openmc.Material(name='SiC')
        sic.add_elements_from_formula('SiC')
        sic.set_density('g/cm3', 3.2)

        # BISO particle
        biso = openmc.Material.mix_materials([kernel, sic], [BISO_KERNEL_VOL_FRAC, BISO_COAT_VOL_FRAC], 'vo') 
        # biso.set_density( )  # get BISO density from mix_materials


        # ------------------------------------------------------------------
        # Breeder and fertile material mixed in the blanket breeding regions 
        # - changed old mass frac function with new/simpler vol frac function --ppark 2025-11-07
        # - changed fertile kg/m³ to be per Li+Be vol, not whole breeder vol
        # 
        # We want "fertile bulk density" to be kg of U-238 per m³ of *breeding* material
        # so we mix Li4SiO4 + Be first, and then mix Li4SiO4-Be + BISO,
        # and then mix Li4SiO4-Be-BISO with the Eurofer and He coolant
        # ------------------------------------------------------------------

        # Volume fractions from Lu (2017) Table 2 
        vf_li4sio4 = 0.1304 ; vf_be = 0.3790 ; vf_eurofer = 0.1176 ; vf_he = 1 - (vf_li4sio4 + vf_be + vf_eurofer) 
        

        # Mix Li4SiO4 and Be (should be 25.6, 74.4 vol% respectively)
        li4sio4_be = openmc.Material.mix_materials([li4sio4, be], [vf_li4sio4/(vf_li4sio4+vf_be), vf_be/(vf_li4sio4+vf_be)], 'vo') 

        # Mix Li4SiO4-Be with BISO
        vol_li4sio4_be = self.breeder_volume*(vf_li4sio4+vf_be)
        
        vf_li4sio4_be, vf_biso = calc_biso_blanket_vol_fracs(self.fertile_bulk_density_kgm3, vol_li4sio4_be, fertile_element=self.fertile_element)
        
        breeder = openmc.Material.mix_materials([li4sio4_be, biso], [vf_li4sio4_be, vf_biso], 'vo')
        
        # So now we mix in Li4SiO4-Be-BISO with the Eurofer structure and He coolant materials
        # Li4SiO4-Be-BISO should be the volume fractions of Li4SiO4 (0.1304) + Be (0.3790) = 0.5094
        self.blanket = openmc.Material.mix_materials([breeder, self.structure, he], [(vf_li4sio4+vf_be), vf_eurofer, vf_he], 'vo') 
        self.blanket.name, self.blanket.temperature = self.name, self.temp_k    


        # ------------------------------------------------------------------
        # Add materials 
        # ------------------------------------------------------------------

        self.materials = openmc.Materials([self.firstwall, self.structure, self.blanket]) 


    def geometry(self):

        # ------------------------------------------------------------------
        # Tokamak geometry parameters 
        # ------------------------------------------------------------------

        self.R0, self.a, self.kappa, self.delta = HCPB_R0, HCPB_A, HCPB_KAPPA, HCPB_DELTA

        d_fw  = HCPB_FW_CM 
        d_st1 = d_fw  + HCPB_ST1_CM
        d_br1 = d_st1 + HCPB_BR1_I_CM
        d_st2 = d_br1 + HCPB_ST2_CM
        
        d_br1_o = d_st1    + HCPB_BR1_O_CM   # only on outboard blanket
        d_st2_o = d_br1_o  + HCPB_ST2_CM     # only on outboard blanket
        
        self.extent_r = (self.R0 + self.a + d_st2_o)*1.2 # 110%
        self.extent_z = (self.kappa*self.a + d_st2_o)*1.2


        # ------------------------------------------------------------------
        # Surfaces 
        # ------------------------------------------------------------------

        points_vc    = miller_model(self.R0, self.a, self.kappa, self.delta)
        points_fw    = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_fw)
        points_st1   = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st1)
        points_br1   = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br1)
        points_st2   = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st2)
        points_br1_o = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br1_o)  # outboard blanket
        points_st2_o = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st2_o)  # outboard structure
        
        # Create OpenMC surfaces
        self.surface_vc    = openmc.model.Polygon(points_vc,    basis='rz')
        self.surface_fw    = openmc.model.Polygon(points_fw,    basis='rz')
        self.surface_st1   = openmc.model.Polygon(points_st1,   basis='rz')
        self.surface_br1   = openmc.model.Polygon(points_br1,   basis='rz')
        self.surface_st2   = openmc.model.Polygon(points_st2,   basis='rz')
        self.surface_br1_o = openmc.model.Polygon(points_br1_o, basis='rz') # outboard blanket
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
        cell_st1   = openmc.Cell(cell_id=21, region= +self.surface_fw  & -self.surface_st1, fill=self.structure)
        cell_br1   = openmc.Cell(cell_id=31, region= -dividing_cylinder & +self.surface_st1 & -self.surface_br1, fill=self.blanket)
        cell_st2   = openmc.Cell(cell_id=22, region= -dividing_cylinder & +self.surface_br1 & -self.surface_st2, fill=self.structure)
        cell_br1_o = openmc.Cell(cell_id=32, region= +dividing_cylinder & +self.surface_st1 & -self.surface_br1_o, fill=self.blanket)
        cell_st2_o = openmc.Cell(cell_id=23, region= +dividing_cylinder & +self.surface_br1_o & -self.surface_st2_o,  fill=self.structure)
        
        # Void cells
        cell_void_i = openmc.Cell(cell_id=98, region= +self.surface_st2 & -dividing_cylinder & -outer_cylinder & +bottom_plane & -top_plane)
        cell_void_o = openmc.Cell(cell_id=97, region= +self.surface_st2_o & +dividing_cylinder & -outer_cylinder & +bottom_plane & -top_plane)
        cell_void_i.importance = {'neutron': 0}
        cell_void_o.importance = {'neutron': 0}

        self.cells = [cell_vc,
                      cell_fw, cell_st1, cell_br1, cell_st2, 
                      cell_br1_o, cell_st2_o,
                      cell_void_i, cell_void_o,]
    
        self.geometry = openmc.Geometry(openmc.Universe(cells=self.cells))
