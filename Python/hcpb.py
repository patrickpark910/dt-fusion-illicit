import numpy as np
import openmc

from Python.reactor    import *
from Python.parameters import *
from Python.utilities  import *


class HCPB(Reactor):

    def initialize(self):

        self.temp_k          = TEMP_K
        self.blanket_name    = 'HCPB'
        self.blanket_volume  = HCPB_BL_VOL   # [m³]
        self.breeder_volume  = HCPB_BR_VOL   # [m³]
        self.breeder_enrich  = ENRICH_HCPB   # at% 

        # Name file based on reactor config - should come out to smth like: tallies_FLiBe_U010kgm3_Li7.5_900K
        self.name = f"{self.run_type}_{self.blanket_name}_{self.temp_k}K_Li{self.breeder_enrich:04.1f}_{self.fertile_isotope}_{self.fertile_kgm3:06.2f}kgm3"         
        self.path = f"./OpenMC/{self.name}"
        
        os.makedirs(self.path, exist_ok=True)

        start_msg = f"\n======== {self.blanket_name} blanket - {self.fertile_kgm3:06.2f} kg({self.fertile_isotope})/m³ - {self.breeder_enrich:4.1f}at%-enriched Li - {self.temp_k} K ========"
        print(f"{C.MAGENTA}{start_msg}{C.END}")


    def materials(self):

        # ------------------------------------------------------------------
        # Helium-4 gas
        # ------------------------------------------------------------------

        he = openmc.Material(name='helium') 
        he.set_density('atom/b-cm', 0.00049800000)  # from Glaser et al. (2025) MCNP -- 0.0033 g/cm3 @ 293 K
        he.add_element('He', 1) 


        # ------------------------------------------------------------------
        # First wall 
        # ------------------------------------------------------------------

        self.firstwall = openmc.Material(name='firstwall', temperature=self.temp_k)
        self.firstwall.set_density('atom/b-cm', 0.06322200000)  # from Glaser et al. (2025) MCNP = 19.3 g/cm3
        self.firstwall.add_element('W',1)


        # ------------------------------------------------------------------
        # Structure 
        # - Removed impurities to nominal specs per Gaganidze (2018) 
        #   "EUROFER 97 Database and Mat Prop Handbook"
        # ------------------------------------------------------------------

        self.eurofer = openmc.Material(name='Eurofer', temperature=self.temp_k)
        self.eurofer.set_density('g/cm3', 7.8)
        self.eurofer.add_element('Fe', 89.36, percent_type='wo')
        self.eurofer.add_element('C' ,  0.11, percent_type='wo')
        self.eurofer.add_element('Cr',  9.00, percent_type='wo')
        self.eurofer.add_element('W' ,  1.10, percent_type='wo')
        self.eurofer.add_element('Mn',  0.40, percent_type='wo')
        self.eurofer.add_element('N' ,  0.03, percent_type='wo')

        # Original Eurofer specs from Lu (2017) "HCPB Analysis"
        # self.eurofer.add_element('Fe', 89.0026, percent_type='wo')
        # self.eurofer.add_element('B', 0.001, percent_type='wo')
        # self.eurofer.add_element('C', 0.1049, percent_type='wo')
        # self.eurofer.add_element('N', 0.04, percent_type='wo')
        # self.eurofer.add_element('O', 0.001, percent_type='wo')
        # self.eurofer.add_element('Al', 0.004, percent_type='wo')
        # self.eurofer.add_element('Si', 0.026, percent_type='wo')
        # self.eurofer.add_element('P', 0.002, percent_type='wo')
        # self.eurofer.add_element('S', 0.003, percent_type='wo')
        # self.eurofer.add_element('Ti', 0.001, percent_type='wo')
        # self.eurofer.add_element('V', 0.01963, percent_type='wo')
        # self.eurofer.add_element('Cr', 9.00, percent_type='wo')
        # self.eurofer.add_element('Mn', 0.55, percent_type='wo')
        # self.eurofer.add_element('Co', 0.005, percent_type='wo')
        # self.eurofer.add_element('Ni', 0.01, percent_type='wo')
        # self.eurofer.add_element('Cu', 0.003, percent_type='wo')
        # self.eurofer.add_element('Nb', 0.005, percent_type='wo')
        # self.eurofer.add_element('Mo', 0.003, percent_type='wo')
        # # self.eurofer.add_element('Ta', 0.12, percent_type='wo') # no cross sections
        # self.eurofer.add_element('W', 1.0987, percent_type='wo')


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


        # ------------------------------------------------------------------
        # Mix blanket materials
        # ------------------------------------------------------------------

        # BISO and (Li4SiO4 + Be) volume fractions relative to BREEDER (Li4SiO4 + Be)
        vf_biso_br, vf_libe_br, biso_per_cc_br = calc_biso_breeder_vol_fracs(self.fertile_kgm3, fertile_isotope=self.fertile_isotope)

        # New volume ratios of everyting relative to BLANKET
        vf_biso_bl = vf_biso_br * (HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)
        vf_libe_bl = vf_libe_br * (HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)
        vf_li_bl   = vf_libe_bl * HCPB_VF_LI_NOM / (HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)
        vf_be_bl   = vf_libe_bl * HCPB_VF_BE_NOM / (HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)

        # Number of BISO spheres per cm³ of BLANKET
        biso_per_cc_bl = vf_biso_bl / BISO_VOLUME

        # Checksums
        checksum1  = vf_biso_br + vf_libe_br  # should equal 1
        checksum2  = vf_biso_bl + vf_libe_bl + HCPB_VF_EU_NOM + HCPB_VF_HE_NOM  # should equal 1
        checksum3  = vf_biso_bl + vf_li_bl + vf_be_bl  # should equal HCPB_VF_LI_NOM + HCPB_VF_BE_NOM = 0.5094
        
        # Blanket material
        self.blanket = openmc.Material.mix_materials([biso,       li4sio4,  be,       self.eurofer, he], 
                                                     [vf_biso_bl, vf_li_bl, vf_be_bl, HCPB_VF_EU_NOM, HCPB_VF_HE_NOM], 'vo') 
        # self.blanket.set_density('atom/b-cm', _)  # Compute from OpenMC
        self.blanket.temperature = self.temp_k 
        self.blanket.name = (f"{self.fertile_kgm3:06.2f} kg/m3"
                             f" | {biso_per_cc_br:.4f} spheres/cc = {(vf_biso_br*100):.4f} vol% in breeder"
                             f" | {biso_per_cc_bl:.4f} spheres/cc = {(vf_biso_bl*100):.4f} vol% in blanket")
        

        # ------------------------------------------------------------------
        # Add materials 
        # ------------------------------------------------------------------

        self.materials = openmc.Materials([self.firstwall, self.eurofer, self.blanket]) 


        # ------------------------------------------------------------------
        # Debugging printouts
        # ------------------------------------------------------------------

        uh_oh_did_i_make_a_fucky_wucky = True

        if self.run_debug and uh_oh_did_i_make_a_fucky_wucky:
            print(f"")
            print(f"| DEBUG PRINTOUT - MATERIALS")
            print(f"| ")
            print(f"| Breeder = Li4SiO4 (ceramic) + Be (metal)")
            print(f"| Breeder vol: {self.breeder_volume:.6f} [cm³] (breeder in blanket)")
            print(f"| Blanket vol: {self.blanket_volume:.6f} [cm³] (breeder + structure + coolant)")
            print(f"| ")
            print(f"| BISO/cm³ of breeder: {biso_per_cc_br:.6f} spheres/cm³")
            print(f"|     /cm³ of blanket: {biso_per_cc_bl:.6f} spheres/cm³")
            print(f"| ")
            print(f"| With respect to the BREEDER, i.e., per 1 m³ of breeder volume, we have these new volume fractions:")
            print(f"|   vf_biso_br =  {(vf_biso_br*100):.6f} vol%")
            print(f"|   vf_libe_br =  {(vf_libe_br*100):.6f} vol%")
            print(f"|   check they add up = {(checksum1*100):.6f} vol%")
            print(f"| ")
            print(f"| With respect to the BLANKET, we have these new volume fractions:")
            print(f"|   vf_biso_bl        =  {(vf_biso_bl*100):.6f} vol%")
            print(f"|   vf_li_bl          =  {(vf_li_bl*100):.6f} vol%")
            print(f"|   vf_be_bl          =  {(vf_be_bl*100):.6f} vol%")
            print(f"|   eurofer           =  {(HCPB_VF_EU_NOM*100):.6f} vol%")
            print(f"|   helium gas        =   {(DCLL_VF_HE_NOM*100):.6f} vol%")
            print(f"|   check they add up = {(checksum2*100):.6f} vol%")
            print(f"|   check that BISO + Li4SiO4 + Be adds up to the nominal Li4SiO4 + Be fraction")
            print(f"|     vf_biso_bl + vf_li_bl + vf_be_bl = {(checksum3*100):.6f} ")
            print(f"|     HCPB_VF_LI_NOM + HCPB_VF_BE_NOM  = {((HCPB_VF_LI_NOM + HCPB_VF_BE_NOM)*100):.6f}")
            print(f"")


    def geometry(self):

        # ------------------------------------------------------------------
        # Tokamak geometry parameters 
        # ------------------------------------------------------------------

        self.R0, self.a, self.kappa, self.delta = HCPB_R0, HCPB_A, HCPB_KAPPA, HCPB_DELTA

        d_fw    = HCPB_FW_CM 
        d_st1   = d_fw     + HCPB_ST1_CM
        d_br1_i = d_st1    + HCPB_BR1_I_CM   # inboard blanket
        d_st2_i = d_br1_i  + HCPB_ST2_CM     # inboard structure
        d_br1_o = d_st1    + HCPB_BR1_O_CM   # outboard blanket
        d_st2_o = d_br1_o  + HCPB_ST2_CM     # outboard blanket
        
        self.extent_r = (self.R0 + self.a + d_st2_o)*1.2 # 110%
        self.extent_z = (self.kappa*self.a + d_st2_o)*1.2


        # ------------------------------------------------------------------
        # Surfaces 
        # ------------------------------------------------------------------

        points_vc    = miller_model(self.R0, self.a, self.kappa, self.delta)
        points_fw    = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_fw)
        points_st1   = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st1)
        points_br1_i = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br1_i)
        points_st2_i = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st2_i)
        points_br1_o = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br1_o)  # outboard blanket
        points_st2_o = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st2_o)  # outboard structure
        
        # Create OpenMC surfaces
        surface_vc    = openmc.model.Polygon(points_vc,    basis='rz')
        surface_fw    = openmc.model.Polygon(points_fw,    basis='rz')
        surface_st1   = openmc.model.Polygon(points_st1,   basis='rz')
        surface_br1_i = openmc.model.Polygon(points_br1_i, basis='rz')  # inboard blanket
        surface_st2_i = openmc.model.Polygon(points_st2_i, basis='rz')  # inboard structure
        surface_br1_o = openmc.model.Polygon(points_br1_o, basis='rz')  # outboard blanket
        surface_st2_o = openmc.model.Polygon(points_st2_o, basis='rz')  # outboard structure

        dividing_cylinder = openmc.ZCylinder(r=self.R0)
        outer_cylinder    = openmc.ZCylinder(r=self.extent_r, boundary_type='vacuum')
        top_plane         = openmc.ZPlane(z0=self.extent_z, boundary_type='vacuum')
        bottom_plane      = openmc.ZPlane(z0=-self.extent_z, boundary_type='vacuum')
        

        # ------------------------------------------------------------------
        # Cells | 10: vc, fw | 2X: structure | 3X: breeding
        # ------------------------------------------------------------------

        cell_vc = openmc.Cell(cell_id=10, region= -surface_vc)
        cell_vc.importance = {'neutron':1}

        cell_fw    = openmc.Cell(cell_id=11, region= +surface_vc  & -surface_fw,  fill=self.firstwall)
        cell_st1   = openmc.Cell(cell_id=21, region= +surface_fw  & -surface_st1, fill=self.eurofer)
        cell_br1_i = openmc.Cell(cell_id=31, region= -dividing_cylinder & +surface_st1   & -surface_br1_i, fill=self.blanket)  # inboard blanket
        cell_st2_i = openmc.Cell(cell_id=22, region= -dividing_cylinder & +surface_br1_i & -surface_st2_i, fill=self.eurofer)  # inboard structure
        cell_br1_o = openmc.Cell(cell_id=32, region= +dividing_cylinder & +surface_st1   & -surface_br1_o, fill=self.blanket)  # outboard blanket
        cell_st2_o = openmc.Cell(cell_id=23, region= +dividing_cylinder & +surface_br1_o & -surface_st2_o, fill=self.eurofer)  # outboard structure
        
        # Void cells
        cell_void_i = openmc.Cell(cell_id=98, region= +surface_st2_i & -dividing_cylinder & -outer_cylinder & +bottom_plane & -top_plane)
        cell_void_o = openmc.Cell(cell_id=97, region= +surface_st2_o & +dividing_cylinder & -outer_cylinder & +bottom_plane & -top_plane)
        cell_void_i.importance = {'neutron': 0}
        cell_void_o.importance = {'neutron': 0}

        self.cells = [cell_vc,
                      cell_fw, cell_st1, 
                      cell_br1_i, cell_st2_i, 
                      cell_br1_o, cell_st2_o,
                      cell_void_i, cell_void_o,]
    
        self.geometry = openmc.Geometry(openmc.Universe(cells=self.cells))
