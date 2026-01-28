import numpy as np
import openmc

from Python.reactor    import *
from Python.parameters import *
from Python.utilities  import *


class DCLL(Reactor):

    def initialize(self):

        self.temp_k          = TEMP_K
        self.blanket_name    = 'DCLL'
        self.blanket_volume  = DCLL_BL_VOL   # [m³]
        self.breeder_volume  = DCLL_BR_VOL   # [m³]
        self.breeder_enrich  = ENRICH_DCLL   # [at%]

        # Name file based on reactor config - should come out to smth like: tallies_FLiBe_U010kgm3_Li7.5_900K
        self.name = f"{self.run_type}_{self.blanket_name}_{self.temp_k}K_Li{self.breeder_enrich:04.1f}_{self.fertile_isotope}_{self.fertile_kgm3:06.2f}kgm3"         
        self.path = f"./OpenMC/{self.name}"
        
        os.makedirs(self.path, exist_ok=True)

        start_msg = f"\n======== {self.blanket_name} blanket - {self.fertile_kgm3:6.2f} kg({self.fertile_isotope})/m³ - {self.breeder_enrich:4.1f} at%-enriched Li - {self.temp_k} K ========"
        print(f"{C.MAGENTA}{start_msg}{C.END}")


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
        # - The original first wall specs we were using from JL Ball (2025) 
        #   is 99.9969 wt% W and the rest impurities. Changing to nominal W 
        #   to isolate fertile reaction effects and save compute time   
        # ------------------------------------------------------------------

        self.firstwall = openmc.Material(name='firstwall', temperature=self.temp_k)
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
        self.firstwall.set_density('g/cm3',19.3)
        # self.firstwall.set_density('atom/b-cm', 0.06322200000)  # from Glaser et al. (2025) MCNP = 19.3 g/cm3
        self.firstwall.add_element('W',1)
        

        # ------------------------------------------------------------------
        # Structure 
        # - F82H steel from T. Hirose (2014) Fus Engr Des
        # ------------------------------------------------------------------

        # F82H steel (first wall front, back + part of other structural components)
        self.f82h = openmc.Material(name='f82h', temperature=self.temp_k)
        self.f82h.add_element('Fe', 89.3686, percent_type='wo')
        self.f82h.add_element('C' ,  0.1000, percent_type='wo')
        self.f82h.add_element('Si',  0.1000, percent_type='wo')
        self.f82h.add_element('Mn',  0.1300, percent_type='wo')
        self.f82h.add_element('Cr',  8.1600, percent_type='wo')
        self.f82h.add_element('W' ,  1.9400, percent_type='wo')
        self.f82h.add_element('V' ,  0.2000, percent_type='wo')
        self.f82h.add_element('N' ,  0.0014, percent_type='wo')
        self.f82h.set_density('g/cm3', 7.78)
        # self.f82h.set_density('atom/b-cm', 0.08365243600)  # from Glaser et al. (2025) MCNP -- 7.887 g/cm3 @ 293 K


        # ------------------------------------------------------------------
        # Coolant (17 vol% F82H + 83 vol% He-4)
        # ------------------------------------------------------------------

        # Helium-4 gas 
        he = openmc.Material(name='helium') 
        he.set_density('g/cm3', 0.004)
        # he.set_density('atom/b-cm', 0.00049800000)  # from Glaser et al. (2025) MCNP -- 0.0033 g/cm3 @ 293 K
        he.add_element('He', 1) 

        self.coolant = openmc.Material.mix_materials([self.f82h, he], [0.170, 0.830], 'vo')
        self.coolant.set_density('atom/b-cm', 0.01471892350) # from Glaser et al. (2025) MCNP
        self.coolant.temperature = self.temp_k
        self.coolant.name = "coolant (17.0 vol% F82H, 83.0 vol% He-4)" 


        # ------------------------------------------------------------------
        # Divider material
        # ------------------------------------------------------------------
        
        self.divider = openmc.Material.mix_materials([self.f82h, he], [0.512, 0.488], 'vo')
        # self.divider.set_density('atom/b-cm', 0.04312289441)
        self.divider.temperature = self.temp_k
        self.divider.name        = "divider (51.2 vol% F82H, 48.8 vol% He-4)"


        # ------------------------------------------------------------------
        # Manifold material
        # ------------------------------------------------------------------

        self.manifold = openmc.Material.mix_materials([self.f82h, he], [0.453, 0.547], 'vo')
        # self.manifold.set_density('atom/b-cm', 0.03822271948)
        self.manifold.temperature = self.temp_k
        self.manifold.name        = "inner manifold (45.3 vol% F82H, 54.7 vol% He-4)"


        # ------------------------------------------------------------------
        # Steel shield material
        # ------------------------------------------------------------------

        self.shield = openmc.Material.mix_materials([self.f82h, he], [0.80, 0.20], 'vo')
        # self.shield.set_density('atom/b-cm', 0.06704198900)  # from Glaser et al. (2025) MCNP
        self.shield.temperature = self.temp_k
        self.shield.name = "steel shield (80 vol% F82H, 20 vol% He-4)"


        # ------------------------------------------------------------------
        # Breeder material
        # ------------------------------------------------------------------

        pbli = openmc.Material(name='breeder', temperature=self.temp_k)
        pbli.set_density('g/cm3', DENSITY_DCLL)
        pbli.add_element('Pb', 0.83, percent_type='ao') 
        pbli.add_element('Li', 0.17, percent_type='ao', enrichment_target='Li6', enrichment_type='ao', enrichment=self.breeder_enrich) # Li-6 enrichment to 90 at% 


        # ------------------------------------------------------------------
        # Fertile material - BISO Particle
        # ------------------------------------------------------------------

        if self.fertile_isotope == 'U238':
            kernel = openmc.Material(name='UO2', temperature=self.temp_k)
            kernel.add_elements_from_formula('UO2', enrichment=ENRICH_U)
            kernel.set_density('g/cm3', DENSITY_UO2)  

        elif self.fertile_isotope == 'Th232': 
            kernel = openmc.Material(name='ThO2', temperature=self.temp_k) 
            kernel.add_elements_from_formula('ThO2') 
            kernel.set_density('g/cm3', DENSITY_ThO2) 

        # SiC coating
        sic = openmc.Material(name='SiC', temperature=self.temp_k)
        sic.add_elements_from_formula('SiC')
        sic.set_density('g/cm3', 3.2)

        # BISO particle
        biso = openmc.Material.mix_materials([kernel, sic], [BISO_KERNEL_VOL_FRAC, BISO_COAT_VOL_FRAC], 'vo') 


        # ------------------------------------------------------------------
        # Mix blanket materials
        # ------------------------------------------------------------------

        # BISO and Pb-Li volume fractions relative to BREEDER (Pb-Li)
        vf_biso_br, vf_pbli_br, biso_per_cc_br = calc_biso_breeder_vol_fracs(self.fertile_kgm3, fertile_isotope=self.fertile_isotope)

        # New volume ratios of everyting relative to BLANKET
        vf_biso_bl = vf_biso_br * DCLL_VF_LL_NOM
        vf_pbli_bl = vf_pbli_br * DCLL_VF_LL_NOM

        # Number of BISO spheres per cm³ of BLANKET
        biso_per_cc_bl = vf_biso_bl / BISO_VOLUME

        # Checksums
        checksum1  = vf_biso_br + vf_pbli_br  # should equal 1
        checksum2  = vf_biso_bl + vf_pbli_bl + DCLL_VF_FS_NOM + DCLL_VF_SI_NOM + DCLL_VF_HE_NOM  # should equal 1
        checksum3  = vf_biso_bl + vf_pbli_bl  # should equal DCLL_VF_LL_NOM = 0.808 
        
        # Blanket material
        self.blanket = openmc.Material.mix_materials([biso, pbli, self.f82h, sic, he], [vf_biso_bl, vf_pbli_bl, DCLL_VF_FS_NOM, DCLL_VF_SI_NOM, DCLL_VF_HE_NOM], 'vo') 
        self.blanket.set_density('atom/b-cm', 0.03541604638)  # from Glaser et al. (2025) MCNP
        self.blanket.temperature = self.temp_k
        self.blanket.name = (f"{self.fertile_kgm3:06.2f} kg/m3"
                             f" | {biso_per_cc_br:.4f} spheres/cc = {(vf_biso_br*100):.4f} vol% in breeder"
                             f" | {biso_per_cc_bl:.4f} spheres/cc = {(vf_biso_bl*100):.4f} vol% in blanket")
        

        # ------------------------------------------------------------------
        # Add materials 
        # ------------------------------------------------------------------

        self.materials = openmc.Materials([self.firstwall, self.f82h, self.coolant, self.blanket, self.divider, self.manifold, self.shield])


        # ------------------------------------------------------------------
        # Debugging printouts
        # ------------------------------------------------------------------

        uh_oh_did_i_make_a_fucky_wucky = True

        if self.run_debug and uh_oh_did_i_make_a_fucky_wucky:
            print(f"")
            print(f"| DEBUG PRINTOUT - MATERIALS")
            print(f"| ")
            print(f"| Breeder = Pb-17Li")
            print(f"| Breeder vol: {self.breeder_volume:.6f} [cm³] (breeder in blanket)")
            print(f"| Blanket vol: {self.blanket_volume:.6f} [cm³] (breeder + structure + coolant)")
            print(f"| ")
            print(f"| BISO/cm³ of breeder: {biso_per_cc_br:.6f} spheres/cm³")
            print(f"|     /cm³ of blanket: {biso_per_cc_bl:.6f} spheres/cm³")
            print(f"| ")
            print(f"| With respect to the BREEDER, i.e., per 1 m³ of breeder volume, we have these new volume fractions:")
            print(f"|   vf_biso_br =  {(vf_biso_br*100):.6f} vol%")
            print(f"|   vf_pbli_br =  {(vf_pbli_br*100):.6f} vol%")
            print(f"|   check they add up = {(checksum1*100):.6f} vol%")
            print(f"| ")
            print(f"| With respect to the BLANKET, we have these new volume fractions:")
            print(f"|   vf_biso_bl        =  {(vf_biso_bl*100):.6f} vol%")
            print(f"|   vf_pbli_bl        =  {(vf_pbli_bl*100):.6f} vol%")
            print(f"|   ferritic steel    =   {(DCLL_VF_FS_NOM*100):.6f} vol%")
            print(f"|   silicon carbide   =   {(DCLL_VF_SI_NOM*100):.6f} vol%")
            print(f"|   helium coolant    =   {(DCLL_VF_HE_NOM*100):.6f} vol%")
            print(f"|   check they add up = {(checksum2*100):.6f} vol%")
            print(f"|   check that BISO + PbLi adds up to the nominal PbLi fraction")
            print(f"|   vf_biso_bl + vf_pbli_bl = {(checksum3*100):.6f} : DCLL_VF_LL_NOM = {(DCLL_VF_LL_NOM*100):.6f}")
            print(f"")


    def geometry(self):

        # ------------------------------------------------------------------
        # Tokamak geometry parameters 
        # ------------------------------------------------------------------

        self.R0, self.a, self.kappa, self.delta = DCLL_R0, DCLL_A, DCLL_KAPPA, DCLL_DELTA

        d_fw   = DCLL_FW_CM
        d_fwf  = d_fw   + DCLL_FWF_CM  # front wall front
        d_fwc  = d_fwf  + DCLL_FWC_CM  # front wall cooling channel
        d_fwb  = d_fwc  + DCLL_FWB_CM  # front wall back
        d_br1  = d_fwb  + DCLL_BR1_CM  # breeding region 1
        d_d1   = d_br1  + DCLL_D1_CM   # divider 1 
        d_br2  = d_d1   + DCLL_BR2_CM  # breeding region 2

        d_im_i = d_br2  + DCLL_IM_CM   # inner manifold - inboard
        d_bp_i = d_im_i + DCLL_BP_CM   # back plate     - inboard
        d_ss_i = d_bp_i + DCLL_SS_CM   # steel shield   - inboard
        
        d_d2_o  = d_br2   + DCLL_D2_CM   # divider 2         - only on outboard 
        d_br3_o = d_d2_o  + DCLL_BR3_CM  # breeding region 3 - only on outboard 
        d_im_o  = d_br3_o + DCLL_IM_CM   # inner manifold    - outboard  
        d_bp_o = d_im_o   + DCLL_BP_CM   # back plate        - outboard
        d_ss_o = d_bp_o   + DCLL_SS_CM   # steel shield      - outboard
        

        self.extent_r = (self.R0 + self.a + d_im_o)  * 1.05  # 105% radial extent
        self.extent_z = (self.kappa*self.a + d_im_o) * 1.05  # 105% vertical extent


        # ------------------------------------------------------------------
        # Surfaces 
        # ------------------------------------------------------------------

        # Draw toroidal D-shapes (points_) -- Returns list of points 
        points_vc   = miller_model(self.R0, self.a, self.kappa, self.delta)
        points_fw   = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_fw)
        points_fwf  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_fwf)
        points_fwc  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_fwc)
        points_fwb  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_fwb)
        points_br1  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br1)
        points_d1   = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_d1)
        points_br2  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br2)

        points_im_i = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_im_i)   # inner manifold - inboard
        points_bp_i = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_bp_i)   # back plate     - inboard
        points_ss_i = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_ss_i)   # steel shield   - inboard
        
        points_d2_o  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_d2_o)  # divider 2         - only on outboard 
        points_br3_o = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br3_o) # breeding region 3 - only on outboard 
        points_im_o  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_im_o)  # inner manifold    - outboard  
        points_bp_o  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_bp_o)  # back plate        - outboard
        points_ss_o  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_ss_o)  # steel shield      - outboard
        

        # Create OpenMC surfaces (s_)
        surface_vc    = openmc.model.Polygon(points_vc,    basis='rz')
        surface_fw    = openmc.model.Polygon(points_fw,    basis='rz')
        surface_fwf   = openmc.model.Polygon(points_fwf,   basis='rz')
        surface_fwc   = openmc.model.Polygon(points_fwc,   basis='rz')
        surface_fwb   = openmc.model.Polygon(points_fwb,   basis='rz')
        surface_br1   = openmc.model.Polygon(points_br1,   basis='rz')
        surface_d1    = openmc.model.Polygon(points_d1,    basis='rz')
        surface_br2   = openmc.model.Polygon(points_br2,   basis='rz')

        surface_im_i  = openmc.model.Polygon(points_im_i,  basis='rz')  # inner manifold - inboard
        surface_bp_i  = openmc.model.Polygon(points_bp_i,  basis='rz')  # back plate     - inboard
        surface_ss_i  = openmc.model.Polygon(points_ss_i,  basis='rz')  # steel shield   - inboard

        surface_d2_o  = openmc.model.Polygon(points_d2_o,  basis='rz')  # divider 2         - only on outboard 
        surface_br3_o = openmc.model.Polygon(points_br3_o, basis='rz')  # breeding region 3 - only on outboard 
        surface_im_o  = openmc.model.Polygon(points_im_o,  basis='rz')  # inner manifold    - outboard  
        surface_bp_o  = openmc.model.Polygon(points_bp_o,  basis='rz')  # back plate        - outboard
        surface_ss_o  = openmc.model.Polygon(points_ss_o,  basis='rz')  # steel shield      - outboard

        dividing_cylinder = openmc.ZCylinder(r=self.R0)
        outer_cylinder    = openmc.ZCylinder(r=self.extent_r, boundary_type='vacuum')
        top_plane         = openmc.ZPlane(z0=self.extent_z, boundary_type='vacuum')
        bottom_plane      = openmc.ZPlane(z0=-self.extent_z, boundary_type='vacuum')
        

        # ------------------------------------------------------------------
        # Cells | 10: vc, fw | 2X: structure | 3X: breeding
        # ------------------------------------------------------------------

        cell_vc   = openmc.Cell(cell_id=10, region= -surface_vc)
        cell_vc.importance = {'neutron':1}

        cell_fw   = openmc.Cell(cell_id=11, region= +surface_vc  & -surface_fw,  fill=self.firstwall)
        cell_fwf  = openmc.Cell(cell_id=21, region= +surface_fw  & -surface_fwf, fill=self.f82h)
        cell_fwc  = openmc.Cell(cell_id=12, region= +surface_fwf & -surface_fwc, fill=self.coolant)
        cell_fwb  = openmc.Cell(cell_id=22, region= +surface_fwc & -surface_fwb, fill=self.f82h)
        cell_br1  = openmc.Cell(cell_id=31, region= +surface_fwb & -surface_br1, fill=self.blanket)
        cell_d1   = openmc.Cell(cell_id=23, region= +surface_br1 & -surface_d1,  fill=self.divider)
        cell_br2  = openmc.Cell(cell_id=32, region= +surface_d1  & -surface_br2, fill=self.blanket)

        cell_im_i = openmc.Cell(cell_id=24, region= -dividing_cylinder & +surface_br2  & -surface_im_i, fill=self.manifold)    # inner manifold - inboard
        cell_bp_i = openmc.Cell(cell_id=41, region= -dividing_cylinder & +surface_im_i & -surface_bp_i, fill=self.f82h)        # back plate     - inboard
        cell_ss_i = openmc.Cell(cell_id=42, region= -dividing_cylinder & +surface_bp_i & -surface_ss_i, fill=self.shield)      # steel shield   - inboard

        cell_d2_o  = openmc.Cell(cell_id=25, region= +dividing_cylinder & +surface_br2   & -surface_d2_o,  fill=self.divider)  # divider 2         - only on outboard 
        cell_br3_o = openmc.Cell(cell_id=33, region= +dividing_cylinder & +surface_d2_o  & -surface_br3_o, fill=self.blanket)  # breeding region 3 - only on outboard  
        cell_im_o  = openmc.Cell(cell_id=26, region= +dividing_cylinder & +surface_br3_o & -surface_im_o,  fill=self.manifold) # inner manifold    - outboard  
        cell_bp_o  = openmc.Cell(cell_id=43, region= +dividing_cylinder & +surface_im_o  & -surface_bp_o,  fill=self.f82h)     # back plate        - outboard
        cell_ss_o  = openmc.Cell(cell_id=44, region= +dividing_cylinder & +surface_bp_o  & -surface_ss_o,  fill=self.shield)   # steel shield      - outboard


        # Void cells
        cell_void_i = openmc.Cell(cell_id=99, region= +surface_im_i & -dividing_cylinder & -outer_cylinder & +bottom_plane & -top_plane)
        cell_void_o = openmc.Cell(cell_id=98, region= +surface_im_o & +dividing_cylinder & -outer_cylinder & +bottom_plane & -top_plane)
        cell_void_i.importance = {'neutron': 0}
        cell_void_o.importance = {'neutron': 0}


        # -------------------------------------------------------
        # Collect all cells
        # -------------------------------------------------------
        self.cells = [
            cell_vc,
            cell_fw, cell_fwf, cell_fwc, cell_fwb, cell_br1, cell_d1, cell_br2, 
            cell_im_i, cell_bp_i, cell_ss_i,
            cell_d2_o, cell_br3_o, cell_im_o, cell_bp_o, cell_ss_o,
            cell_void_i, cell_void_o,
        ]
        self.geometry = openmc.Geometry(openmc.Universe(cells=self.cells))
