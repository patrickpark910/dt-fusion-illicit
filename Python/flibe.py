import numpy as np
import openmc

from Python.reactor    import *
from Python.parameters import *
from Python.utilities  import *


class FLiBe(Reactor):

    def initialize(self):

        self.temp_k          = TEMP_K
        self.breeder_name    = 'FLiBe'
        self.breeder_density = DENSITY_FLIBE # g/cm^3
        self.breeder_enrich  = ENRICH_FLIBE  # wt% 
        self.breeder_volume  = FLIBE_BR_VOL  # m^3
        self.R0, self.a, self.kappa, self.delta = FLIBE_R0, FLIBE_A, FLIBE_KAPPA, FLIBE_DELTA 

        # Name file based on reactor config - should come out to smth like: tallies_FLiBe_U010kgm3_Li7.5_900K
        self.name = f"{self.run_type}_{self.breeder_name}_{self.temp_k}K_Li{self.breeder_enrich:04.1f}_{self.fertile_element}{self.fertile_bulk_density_kgm3:06.2f}kgm3"         
        self.path = f"./OpenMC/{self.name}"
        
        os.makedirs(self.path, exist_ok=True)

        start_msg = f"\n======== {self.breeder_name} reactor - {self.fertile_element} {self.fertile_bulk_density_kgm3:6.2f} kg/m3 - {self.breeder_enrich:4.1f}%-enriched Li - {self.temp_k} K ========"
        print(f"{Colors.CYAN}{start_msg}{Colors.END}")


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
        # Structure (V-4Cr-4Ti)
        # ------------------------------------------------------------------

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
        self.structure.add_element('V',0.919139,percent_type='wo')
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
                                               enrichment_type='ao', 
                                               enrichment=self.breeder_enrich)

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

        self.materials = openmc.Materials([self.firstwall, self.structure, self.blanket]) # self.air, 
        # self.materials.export_to_xml(self.path)


    def geometry(self):

        # ------------------------------------------------------------------
        # Tokamak geometry parameters 
        # ------------------------------------------------------------------

        d_fw  = FLIBE_FW_CM 
        d_st0 = d_fw  + FLIBE_ST1_CM
        d_br0 = d_st0 + FLIBE_BR1_CM
        d_st1 = d_br0 + FLIBE_ST2_CM
        d_br1 = d_st1 + FLIBE_BR2_CM
        d_st2 = d_br1 + FLIBE_ST3_CM

        self.extent_r = (self.R0 + self.a + d_st2)*1.2 # 110%
        self.extent_z = (self.kappa*self.a + d_st2)*1.2


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