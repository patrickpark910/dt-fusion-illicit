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
        self.name = f"{self.run_type}_{self.breeder_name}_{self.temp_k}K_Li{self.breeder_enrich:04.1f}_{self.fertile_isotope}_{self.fertile_str}kgm3"         
        self.path = f"./OpenMC/{self.name}"
        
        os.makedirs(self.path, exist_ok=True)

        start_msg = f"\n======== {self.breeder_name} reactor - {self.fertile_isotope} {self.fertile_str} kg/m3 - {self.breeder_enrich:4.1f}%-enriched Li - {self.temp_k} K ========"
        print(f"{Colors.MAGENTA}{start_msg}{Colors.END}")

        self.summary = ""

        self.helpers()


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
        li4sio4.set_density('g/cm3', DENSITY_LI4SIO4)  
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
        be.set_density('g/cm3', DENSITY_BE)  
        be.add_element('Be', 1, percent_type='wo') 
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

        if self.fertile_isotope == 'U238':
            kernel = openmc.Material(name='UO2')
            kernel.add_elements_from_formula('UO2', enrichment=ENRICH_U)
            kernel.set_density('g/cm3', DENSITY_UO2)  

        elif self.fertile_isotope == 'Th232': 
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

        vf_li_bv = self.vf_breeder_bv * VF_LI_NOM / (VF_LI_NOM + VF_BE_NOM)
        vf_be_bv = self.vf_breeder_bv * VF_BE_NOM / (VF_LI_NOM + VF_BE_NOM)

        breeder = openmc.Material.mix_materials([li4sio4, be, biso], [vf_li_bv, vf_be_bv, self.vf_biso_bv], 'vo')
        
        # So now we mix in Li4SiO4-Be-BISO with the Eurofer structure and He coolant materials
        self.blanket = openmc.Material.mix_materials([breeder, self.structure, he], [(VF_LI_NOM+VF_BE_NOM), VF_EU_NOM, VF_HE_NOM], 'vo') 
        self.blanket.name = (f"{self.fertile_str} kg/m3"
                             f" | {self.biso_per_cc_bv:.4f} spheres/cc = {(self.vf_biso_bv*100):.4f} vol% in breeder volume"
                             f" | {self.biso_per_cc_br:.4f} spheres/cc = {(self.vf_biso_br*100):.4f} vol% in breeding region")
        self.blanket.temperature = self.temp_k    


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
        self.surface_br1_o = openmc.model.Polygon(points_br1_o, basis='rz')  # outboard blanket
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


    def helpers(self):

        ''' 
        Nominally, the volume fracs are: Li4SiO4 = 13.04%, Be = 37.90%, Eurofer = 11.76%, He = 37.30%
        We want to deduct the BISO volume we add from the Li4SiO4.

        Here I use: "breeder" = Li4SiO4, Be
                    "breeding region" (br) = the physical blanket layer that contains the breeder
                                             this layer usually contains breeder + structure + coolant
                    "breeding volume" (bv) = amount of volume nominally occupied in blanket by breeder

        That is, if the breeding region has volume 100 m³ of which 60% is (Li4SiO4 + Be) and the other 
        40% is (Eurofer + Be), then 60 m³ is the breeding volume.
            
        Given fertile_kgm3 [kg U-238 / m³ breeder] we convert this to m³ of BISO.

        **Be careful NOT to renormalize everything. DEDUCT the volume fraction of BISO from 1, and then
        set THAT as the volume fraction of breeder, which then subdivides into Li4SiO4 and Be in a 
        13.04:37.90 ratio.

        We use the fact that 13.04 vol% Li4SiO4 = 0.1304 m³ of Li4SiO4 per m³ of blanket to do our work 
        below in terms of volume fractions.
        '''

        self.vf_biso_bv, self.vf_breeder_bv, self.biso_per_cc_bv = calc_biso_breeder_vol_fracs(self.fertile_kgm3, fertile_isotope=self.fertile_isotope)

        # Multiply nominal volume fractions of Li4SiO4, Be in bv by the new volume fraction of breeders in bv
        vf_li_bv     = VF_LI_NOM / (VF_LI_NOM + VF_BE_NOM) * self.vf_breeder_bv
        vf_be_bv     = VF_BE_NOM / (VF_LI_NOM + VF_BE_NOM) * self.vf_breeder_bv
        # vf_checksum1 = self.vf_biso_bv + vf_li_bv + vf_be_bv  # should equal 1

        # New volume ratios of everyting w.r.t. everything else in the breeding region
        self.vf_biso_br = self.vf_biso_bv * (VF_LI_NOM + VF_BE_NOM)
        vf_li_br   = vf_li_bv * (VF_LI_NOM + VF_BE_NOM)
        vf_be_br   = vf_be_bv * (VF_LI_NOM + VF_BE_NOM)
        # vol frac of Eurofer and He-4 doesn't change
        vf_checksum2 = self.vf_biso_br + vf_li_br + vf_be_br + VF_EU_NOM + VF_HE_NOM  # should equal 1

        # Number of BISO spheres per cc of breeding region
        self.biso_per_cc_br = self.vf_biso_br / BISO_VOLUME
