import numpy as np
import openmc

from Python.reactor import * #     import 
from Python.parameters import *
from Python.utilities  import *


class ARC(Reactor):

    def initialize(self):

        self.temp_k          = TEMP_K
        self.breeder_name    = 'ARC'
        self.breeder_density = DENSITY_FLIBE # g/cm³
        self.breeder_enrich  = ENRICH_FLIBE  # wt%
        self.breeder_volume  = ARC_BR_VOL    # m³
        self.R0, self.a, self.kappa, self.delta = ARC_R0, ARC_A, ARC_KAPPA, ARC_DELTA 

        # Name file based on reactor config - should come out to smth like: tallies_FLiBe_U010kgm3_Li7.5_900K
        self.name = f"{self.run_type}_{self.breeder_name}_{self.temp_k}K_Li{self.breeder_enrich:04.1f}_{self.fertile_isotope}_{self.fertile_kgm3:06.2f}kgm3"         
        self.path = f"./OpenMC/{self.name}"
        
        os.makedirs(self.path, exist_ok=True)

        start_msg = f"\n======== {self.breeder_name} reactor - {self.fertile_isotope} {self.fertile_kgm3:6.2f} kg/m3 - {self.breeder_enrich:4.1f}%-enriched Li - {self.temp_k} K ========"
        print(f"{Colors.BLUE}{start_msg}{Colors.END}")


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
                                               enrichment_type='wo', 
                                               enrichment=self.breeder_enrich)


        # ------------------------------------------------------------------
        # Breeder and fertile material mixed in the blanket breeding regions 
        # ------------------------------------------------------------------        

        if self.fertile_isotope == 'U238':

            self.fertile = openmc.Material(name='fertile', temperature=self.temp_k)
            self.fertile.add_elements_from_formula('UF4','ao') # 'ao' here refers to 1:4 atomic ratio of U:F in UF4
            self.fertile.set_density('g/cm3', DENSITY_UF4) 

            breeder_mass_frac, fertile_compound_mass_frac = calc_blanket_mass_fracs(self.fertile_kgm3, 
                                                                                    self.breeder_volume,
                                                                                    fertile_isotope=self.fertile_isotope, 
                                                                                    fertile_enrich=ENRICH_U, 
                                                                                    breeder_density_kgm3=DENSITY_FLIBE*1e3)
            self.blanket = openmc.Material.mix_materials([self.breeder, self.fertile], [breeder_mass_frac, fertile_compound_mass_frac], 'wo') # fractions in 'mix_materials' MUST add up to 1
            self.blanket.name, self.blanket.temperature = self.name, self.temp_k

        elif self.fertile_isotope == 'Th232':

            self.fertile = openmc.Material(name='fertile', temperature=self.temp_k)
            self.fertile.add_elements_from_formula('ThF4','ao') # 'ao' here refers to 1:4 atomic ratio of U:F in UF4
            self.fertile.set_density('g/cm3', DENSITY_ThF4) 

            breeder_mass_frac, fertile_compound_mass_frac = calc_blanket_mass_fracs(self.fertile_kgm3, 
                                                                                    self.breeder_volume,
                                                                                    fertile_isotope=self.fertile_isotope, 
                                                                                    fertile_enrich=100, 
                                                                                    breeder_density_kgm3=DENSITY_FLIBE*1e3)
            print(breeder_mass_frac, fertile_compound_mass_frac)
            self.blanket = openmc.Material.mix_materials([self.breeder, self.fertile], [breeder_mass_frac, fertile_compound_mass_frac], 'wo') # fractions in 'mix_materials' MUST add up to 1
            self.blanket.name, self.blanket.temperature = self.name, self.temp_k



        # ------------------------------------------------------------------
        # Add materials 
        # ------------------------------------------------------------------
        
        self.materials = openmc.Materials([self.firstwall, self.structure, self.blanket]) # self.air, 
        # self.materials.export_to_xml(self.path)


    def geometry(self):

        # ------------------------------------------------------------------
        # Tokamak geometry parameters 
        # ------------------------------------------------------------------

        d_fw  = ARC_FW_CM 
        d_st0 = d_fw  + ARC_ST1_CM
        d_br0 = d_st0 + ARC_BR1_CM
        d_st1 = d_br0 + ARC_ST2_CM
        d_br1 = d_st1 + ARC_BR2_CM
        d_st2 = d_br1 + ARC_ST3_CM

        self.extent_r = (self.R0 + self.a + d_st2)*1.1  # 110%
        self.extent_z = (self.kappa*self.a + d_st2)*1.2 # 120%


        # ------------------------------------------------------------------
        # Surfaces 
        # ------------------------------------------------------------------

        # Arrays of a (minor r) and z points
        points_vc  = miller_model(self.R0, self.a, self.kappa, self.delta)                 # coords around vacuum chamber
        points_fw  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_fw)   # outer coords around first wall
        points_st0 = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st0)  # outer coords around structural region 0
        points_br0 = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br0)  # outer coords around breeding channel
        points_st1 = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st1)  # outer coords around structural region 1
        points_br1 = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br1)  # outer coords around breeding region 1
        points_st2 = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_st2)  # outer coords around structural region 2

        # Create OpenMC surfaces and STORE THEM AS CLASS ATTRIBUTES
        self.surface_vc  = openmc.model.Polygon(points_vc , basis='rz')  # Plasma-facing surface (inner FW)
        self.surface_fw  = openmc.model.Polygon(points_fw , basis='rz')  # Outer surface of first wall
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
        # Cells | 10: vc, fw | 2X: structure | 3X: breeding
        # ------------------------------------------------------------------
        cell_vc   = openmc.Cell(cell_id=10, region= -self.surface_vc)
        cell_vc.importance = {'neutron':1}
        cell_fw   = openmc.Cell(cell_id=11, region= +self.surface_vc  & -self.surface_fw  , fill=self.firstwall) 
        cell_st0  = openmc.Cell(cell_id=21, region= +self.surface_fw  & -self.surface_st0 , fill=self.structure)
        cell_br0  = openmc.Cell(cell_id=31, region= +self.surface_st0 & -self.surface_br0 , fill=self.blanket)
        cell_st1  = openmc.Cell(cell_id=22, region= +self.surface_br0 & -self.surface_st1 , fill=self.structure)
        cell_br1  = openmc.Cell(cell_id=32, region= +self.surface_st1 & -self.surface_br1 , fill=self.blanket)
        cell_st2  = openmc.Cell(cell_id=23, region= +self.surface_br1 & -self.surface_st2 , fill=self.structure) 

        # Surrounding air cell with proper boundaries (otherwise causes error with just Polygons)
        cell_void = openmc.Cell(cell_id=99, region= +self.surface_st2 & -outer_cylinder & +bottom_plane & -top_plane) # fill=self.air
        cell_void.importance = {'neutron':0}

        self.cells = [cell_vc, cell_fw, cell_st0, cell_br0, cell_st1, cell_br1, cell_st2, cell_void]
        self.geometry = openmc.Geometry(openmc.Universe(cells=self.cells))
        # self.geometry.export_to_xml(self.path)
