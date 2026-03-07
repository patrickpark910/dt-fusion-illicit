import numpy as np
import openmc

from Python.reactor    import *
from Python.parameters import *
from Python.utilities  import *


class FLiBe(Reactor):

    def initialize(self):

        self.temp_k          = TEMP_K
        self.blanket_name    = 'FLiBe'
        self.blanket_volume  = FLIBE_BL_VOL         # m³ total breeding region (FLiBe + XF4)
        self.breeder_density = DENSITY_FLIBE        # g/cm³
        self.breeder_enrich  = 0 # ENRICH_FLIBE     # wt% 
        self.breeder_volume  = self.blanket_volume  # set in materials() after XF4 deduction


        # Name file based on reactor config - should come out to smth like: tallies_FLiBe_U010kgm3_Li7.5_900K
        self.name = f"{self.run_type}_{self.blanket_name}_{self.temp_k}K_Li{self.breeder_enrich:04.1f}_{self.fertile_isotope}_{self.fertile_kgm3:06.2f}kgm3"         
        self.path = f"./OpenMC/{self.name}"
        
        os.makedirs(self.path, exist_ok=True)

        start_msg = f"\n======== {self.blanket_name} blanket - {self.fertile_isotope} {self.fertile_kgm3:6.2f} kg/m3 - {self.breeder_enrich:4.1f}%-enriched Li - {self.temp_k} K ========"
        print(f"{C.CYAN}{start_msg}{C.END}")


    def materials(self):

        # ------------------------------------------------------------------
        # First wall 
        # ------------------------------------------------------------------

        self.firstwall = openmc.Material(name='firstwall', temperature=self.temp_k)
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
        self.firstwall.set_density('g/cm3',19.0)
        # The original first wall specs we were using from Ball 25 is 99.9969 wt% W 
        # and the rest O, N ,C, Na, K, Al, Ca, Cr, Cu, Fe impurities...
        # I think we can just set it to be W to save compute time loading nuclide data
        
        
        # ------------------------------------------------------------------
        # Beryllium    
        # ------------------------------------------------------------------

        self.beryllium = openmc.Material(name='beryllium', temperature=self.temp_k)
        self.beryllium.depletable = False
        self.beryllium.add_element('Be', 1) 
        self.beryllium.set_density('g/cm3', DENSITY_BE) 


        # ------------------------------------------------------------------
        # Structure (V-4Cr-4Ti)
        # ------------------------------------------------------------------

        self.structure = openmc.Material(name='structure', temperature=self.temp_k)
        self.structure.depletable = False
        # This code is from jlball but V-4Cr-4Ti specs also can be found
        # from ANL material id BL-47, p.3 (1994) altho these are slightly different
        # cf. osti.gov/servlets/purl/10194461 --ppark 2025-07-22
        self.structure.add_element('Cr',0.04,percent_type='wo')
        self.structure.add_element('Ti',0.04,percent_type='wo')
        # self.structure.add_element('C',56/1e6,percent_type='wo')
        # self.structure.add_element('O',181/1e6,percent_type='wo')
        # self.structure.add_element('N',103/1e6,percent_type='wo')
        # self.structure.add_element('B',7/1e6,percent_type='wo')
        # self.structure.add_element('Na',17/1e6,percent_type='wo')
        # self.structure.add_element('Mg',0.5/1e6,percent_type='wo')
        # self.structure.add_element('Al',119/1e6,percent_type='wo')
        # self.structure.add_element('Si',280/1e6,percent_type='wo')
        # self.structure.add_element('Mn',0.5/1e6,percent_type='wo')
        # self.structure.add_element('Fe',80/1e6,percent_type='wo')
        # self.structure.add_element('Ni',13/1e6,percent_type='wo')
        # self.structure.add_element('Cu',4/1e6,percent_type='wo')
        self.structure.add_element('V',0.92,percent_type='wo') # 0.919139
        # 1-0.04-0.04-(56+181+103+7+17+0.5+119+280+0.5+80+13+4)/1e6 = 0.919139
        self.structure.set_density('g/cm3',6.05) 
        # This density value is sus and needs a good source --jlball 
        # This value is from Metals Handbook, 9th ed, vol 2: "Properties and Selection: Nonferrous Alloys and Pure Metals" (1979) --ppark 2025-07-22
        

        # ------------------------------------------------------------------
        # Breeder material
        # ------------------------------------------------------------------

        self.breeder = openmc.Material(name='breeder', temperature=self.temp_k)
        self.breeder.set_density('g/cm3', self.breeder_density)
        self.breeder.add_elements_from_formula('F4Li2Be', 'ao', 
                                               enrichment_target='Li6', 
                                               enrichment_type='ao',  # I don't know why Alex says this should be wt%. I'm deferring to the MIT papers that say 7.5 at% here. --ppark
                                               enrichment=self.breeder_enrich)

        
        # ------------------------------------------------------------------
        # Breeder and fertile material mixed in the blanket breeding regions 
        # ------------------------------------------------------------------        

        if self.fertile_isotope == 'U238':

            self.fertile = openmc.Material(name='fertile', temperature=self.temp_k)
            self.fertile.add_elements_from_formula('UF4','ao') # 'ao' here refers to 1:4 atomic ratio of U:F in UF4
            self.fertile.set_density('g/cm3', DENSITY_UF4) 
            xf4_x_ratio = AMU_UF4 / AMU_U238

        elif self.fertile_isotope == 'Th232':

            self.fertile = openmc.Material(name='fertile', temperature=self.temp_k)
            self.fertile.add_elements_from_formula('ThF4','ao') # 'ao' here refers to 1:4 atomic ratio of Th:F in ThF4
            self.fertile.set_density('g/cm3', DENSITY_ThF4) 
            xf4_x_ratio = AMU_ThF4 / AMU_Th232


        # Per 1 m³ blanket, we deduct the volume of XF4 based on its mass density and specified kg/m³, then we compute the new volume fractions of FLiBe and XF4.
        # ORNL data show mixing UF4/ThF4 in 2(LiF)-BeF2 is additive in molar volume, i.e., XF4 takes up its own volume in the mixture, up to c.1000 C.
        # cf. Cantor 73 (ORNL-TM-4308) p.19, Cantor 68 (ORNL-TM-2316) p.34

        vf_xf4   = self.fertile_kgm3 * xf4_x_ratio / (self.fertile.density*1000)  # 1 g/cm³ = 1000 kg/m³
        vf_flibe = 1.0 - vf_xf4

        # FLiBe volume = blanket volume minus XF4 volume (deduct as U/Th loading increases)
        self.breeder_volume = self.blanket_volume * vf_flibe

        self.blanket = openmc.Material.mix_materials([self.breeder, self.fertile], [vf_flibe, vf_xf4], 'vo') # fractions in 'mix_materials' MUST add up to 1
        self.blanket.name, self.blanket.temperature = self.name, self.temp_k


        # ------------------------------------------------------------------
        # Add materials 
        # ------------------------------------------------------------------

        self.materials = openmc.Materials([self.firstwall, self.beryllium, self.structure, self.blanket]) # self.air, 
        # self.materials.export_to_xml(self.path)


        # ------------------------------------------------------------------
        # Debugging printouts
        # ------------------------------------------------------------------

        uh_oh_did_i_make_a_fucky_wucky = True

        if self.run_debug and uh_oh_did_i_make_a_fucky_wucky:
            print(f"")
            print(f"| DEBUG PRINTOUT - MATERIALS")
            print(f"| ")
            print(f"| Breeder = 2(LiF)-BeF2")
            print(f"| Blanket vol: {self.blanket_volume:.6f} [m³] (total blanket volume)")
            print(f"| Breeder vol: {self.breeder_volume:.6f} [m³] (breeder in blanket)")
            print(f"| XF4     vol: {(self.blanket_volume-self.breeder_volume):.6f} [m³] (XF4 in blanket)")
            print(f"| ")
            print(f"| With respect to the BREEDER, i.e., per 1 m³ of breeder volume, we have these new volume fractions:")
            print(f"|   vf_flibe =  {(vf_flibe*100):.6f} vol%")
            print(f"|   vf_xf4   =  {(vf_xf4*100):.6f} vol%")
            print(f"|   check they add up = {(vf_xf4+vf_flibe)*100:.6f} vol%")
            print(f"")


    def geometry(self):

        # ------------------------------------------------------------------
        # Tokamak geometry parameters 
        # ------------------------------------------------------------------

        self.R0, self.a, self.kappa, self.delta = FLIBE_R0, FLIBE_A, FLIBE_KAPPA, FLIBE_DELTA 

        d_fw  = FLIBE_FW_CM 
        d_be  = d_fw  + FLIBE_BE_CM
        d_st0 = d_be  + FLIBE_ST1_CM
        d_br0 = d_st0 + FLIBE_BR1_CM
        d_st1 = d_br0 + FLIBE_ST2_CM
        d_br1 = d_st1 + FLIBE_BR2_CM
        d_st2 = d_br1 + FLIBE_ST3_CM

        self.extent_r = (self.R0 + self.a + d_st2)*1.05 # 105%
        self.extent_z = (self.kappa*self.a + d_st2)*1.05


        # ------------------------------------------------------------------
        # Surfaces 
        # ------------------------------------------------------------------

        # Arrays of a (minor r) and z points
        points_vc  = miller_model(self.R0, self.a, self.kappa, self.delta)                 # coords around vacuum chamber
        points_fw  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_fw)   # outer coords around first wall
        points_be  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_be)   # outer coords around beryllium layer
        points_st0 = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st0)  # outer coords around structural region 0
        points_br0 = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br0)  # outer coords around breeding channel
        points_st1 = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st1)  # outer coords around structural region 1
        points_br1 = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br1)  # outer coords around breeding region 1
        points_st2 = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st2)  # outer coords around structural region 2

        # Create OpenMC surfaces and STORE THEM AS CLASS ATTRIBUTES
        self.surface_vc  = openmc.model.Polygon(points_vc , basis='rz')  # Plasma-facing surface (inner FW)
        self.surface_fw  = openmc.model.Polygon(points_fw , basis='rz')  # Outer surface of first wall
        self.surface_be  = openmc.model.Polygon(points_be , basis='rz')  # Outer surface of beryllium layer
        self.surface_st0 = openmc.model.Polygon(points_st0, basis='rz')
        self.surface_br0 = openmc.model.Polygon(points_br0, basis='rz')
        self.surface_st1 = openmc.model.Polygon(points_st1, basis='rz')
        self.surface_br1 = openmc.model.Polygon(points_br1, basis='rz')
        self.surface_st2 = openmc.model.Polygon(points_st2, basis='rz')

        # Add boundary surfaces
        outer_cylinder = openmc.ZCylinder(r=self.extent_r, boundary_type='vacuum')
        top_plane      = openmc.ZPlane(z0=self.extent_z, boundary_type='vacuum')
        bottom_plane   = openmc.ZPlane(z0=-self.extent_z, boundary_type='vacuum')
        

        # ------------------------------------------------------------------
        # Cells | 10: vc | 11: fw | 41: be | 2X: structure | 3X: breeding
        # ------------------------------------------------------------------
        cell_vc   = openmc.Cell(cell_id=10, region= -self.surface_vc)
        cell_vc.importance = {'neutron':1}
        cell_fw   = openmc.Cell(cell_id=11, region= +self.surface_vc  & -self.surface_fw  , fill=self.firstwall) 
        cell_be   = openmc.Cell(cell_id=41, region= +self.surface_fw  & -self.surface_be  , fill=self.beryllium)
        cell_st0  = openmc.Cell(cell_id=21, region= +self.surface_be  & -self.surface_st0 , fill=self.structure)
        cell_br0  = openmc.Cell(cell_id=31, region= +self.surface_st0 & -self.surface_br0 , fill=self.blanket)
        cell_st1  = openmc.Cell(cell_id=22, region= +self.surface_br0 & -self.surface_st1 , fill=self.structure)
        cell_br1  = openmc.Cell(cell_id=32, region= +self.surface_st1 & -self.surface_br1 , fill=self.blanket)
        cell_st2  = openmc.Cell(cell_id=23, region= +self.surface_br1 & -self.surface_st2 , fill=self.structure) 

        # Surrounding air cell with proper boundaries (otherwise causes error with just Polygons)
        cell_void = openmc.Cell(cell_id=99, region= +self.surface_st2 & -outer_cylinder & +bottom_plane & -top_plane) # fill=self.air
        cell_void.importance = {'neutron':0}

        self.cells = [cell_vc, cell_fw, cell_be, cell_st0, cell_br0, cell_st1, cell_br1, cell_st2, cell_void]
        self.geometry = openmc.Geometry(openmc.Universe(cells=self.cells))
        # self.geometry.export_to_xml(self.path)
