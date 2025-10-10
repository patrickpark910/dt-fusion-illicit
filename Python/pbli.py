import numpy as np
import openmc

from Python.reactor    import *
from Python.parameters import *
from Python.utilities  import *


class LL(Reactor):

    def initialize(self):

        self.temp_k          = TEMP_K
        self.breeder_name    = 'LL'
        self.breeder_density = DENSITY_LL # g/cm^3
        self.breeder_enrich  = ENRICH_LL  # wt% 
        self.breeder_volume  = LL_BR_VOL  # m^3 must check could be fishy -ezoccoli
        self.R0, self.a, self.kappa, self.delta = LL_R0, LL_A, LL_KAPPA, LL_DELTA

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
        self.firstwall.depletable = False
        self.firstwall.add_element('W',1)
        self.firstwall.set_density('g/cm3',19.3)
        # The original first wall specs we were using from Ball 25 is 99.9969 wt% W 
        # and the rest O, N ,C, Na, K, Al, Ca, Cr, Cu, Fe impurities...
        # I think we can just set it to be W to save compute time loading nuclide data
        

        # ------------------------------------------------------------------
        # Structure 
        # ------------------------------------------------------------------

        # First wall front (F82H steel)
        self.structure = openmc.Material(name='firstwall_front', temperature=self.temp_k)
        self.structure.depletable = False
        self.structure.add_nuclide("C12",    0.00040000000, "ao")
        self.structure.add_nuclide("Si28",  0.00015679100, "ao")
        self.structure.add_nuclide("Si29",  0.00000795600, "ao")
        self.structure.add_nuclide("Si30",  0.00000525300, "ao")
        self.structure.add_nuclide("V51",   0.00018400000, "ao")
        self.structure.add_nuclide("Cr50",  0.00031327500, "ao")
        self.structure.add_nuclide("Cr52",  0.00604119000, "ao")
        self.structure.add_nuclide("Cr53",  0.00068502200, "ao")
        self.structure.add_nuclide("Cr54",  0.00017051700, "ao")
        self.structure.add_nuclide("Mn55",  0.00008600000, "ao")
        self.structure.add_nuclide("Fe54",  0.00438860000, "ao")
        self.structure.add_nuclide("Fe56",  0.06889170000, "ao")
        self.structure.add_nuclide("Fe57",  0.00159101000, "ao")
        self.structure.add_nuclide("Fe58",  0.00021173400, "ao")
        self.structure.add_nuclide("Ta181",  0.00001000000, "ao") 
        self.structure.add_nuclide("W182",  0.00013515000, "ao")
        self.structure.add_nuclide("W183",  0.00007298100, "ao")
        self.structure.add_nuclide("W184",  0.00015626400, "ao")
        self.structure.add_nuclide("W186",  0.00014499300, "ao")
        self.structure.set_density("atom/b-cm", 0.08365243600)  # from MCNP cell 112

        self.coolant = openmc.Material(name="coolant", temperature=self.temp_k)
        self.coolant.depletable = False
        self.coolant.add_nuclide("He4", 0.00049800000, "ao")
        self.coolant.add_nuclide("C12",  0.00006800000, "ao")
        self.coolant.add_nuclide("Si28",0.00002665450, "ao")
        self.coolant.add_nuclide("Si29",0.00000135250, "ao")
        self.coolant.add_nuclide("Si30",0.00000089300, "ao")
        self.coolant.add_nuclide("V51", 0.00003128000, "ao")
        self.coolant.add_nuclide("Cr50",0.00005325680, "ao")
        self.coolant.add_nuclide("Cr52",0.00102700000, "ao")
        self.coolant.add_nuclide("Cr53",0.00011645400, "ao")
        self.coolant.add_nuclide("Cr54",0.00002898790, "ao")
        self.coolant.add_nuclide("Mn55",0.00001462000, "ao")
        self.coolant.add_nuclide("Fe54",0.00074606200, "ao")
        self.coolant.add_nuclide("Fe56",0.01171160000, "ao")
        self.coolant.add_nuclide("Fe57",0.00027047200, "ao")
        self.coolant.add_nuclide("Fe58",0.00003599480, "ao")
        self.coolant.add_nuclide("Ta181",0.00000170000, "ao")
        self.coolant.add_nuclide("W182",0.00002297550, "ao")
        self.coolant.add_nuclide("W183",0.00001240680, "ao")
        self.coolant.add_nuclide("W184",0.00002656490, "ao")
        self.coolant.add_nuclide("W186",0.00002464880, "ao")
        self.coolant.set_density("atom/b-cm", 0.01471892350)


        # ------------------------------------------------------------------
        # Breeder material
        # ------------------------------------------------------------------

        self.breeder = openmc.Material(name='breeder', temperature=self.temp_k)
        self.breeder.set_density('g/cm3', self.breeder_density)
        self.breeder.add_element('Pb', 0.83, percent_type='ao') 
        self.breeder.add_element('Li', 0.17, percent_type='ao', enrichment_target='Li6', enrichment_type='ao', enrichment=self.breeder_enrich) # Li-6 enrichment to 90% 

        # ------------------------------------------------------------------
        # Divider material
        # ------------------------------------------------------------------
        
        self.divider = openmc.Material(name="divider1", temperature=self.temp_k)
        self.divider.depletable = False
        self.divider.add_nuclide("He4", 0.00029280000, "ao")
        self.divider.add_nuclide("C12", 0.00020480000, "ao")
        self.divider.add_nuclide("Si28", 0.00008027700, "ao")
        self.divider.add_nuclide("Si29", 0.00000407347, "ao")
        self.divider.add_nuclide("Si30", 0.00000268954, "ao")
        self.divider.add_nuclide("V51", 0.00009420800, "ao")
        self.divider.add_nuclide("Cr50", 0.00016039700, "ao")
        self.divider.add_nuclide("Cr52", 0.00309309000, "ao")
        self.divider.add_nuclide("Cr53", 0.00035073100, "ao")
        self.divider.add_nuclide("Cr54", 0.00008730470, "ao")
        self.divider.add_nuclide("Mn55", 0.00004403200, "ao")
        self.divider.add_nuclide("Fe54", 0.00224696000, "ao")
        self.divider.add_nuclide("Fe56", 0.03527260000, "ao")
        self.divider.add_nuclide("Fe57", 0.00081459700, "ao")
        self.divider.add_nuclide("Fe58", 0.00010840800, "ao")
        self.divider.add_nuclide("Ta181", 0.00000512000, "ao")
        self.divider.add_nuclide("W182", 0.00006919680, "ao")
        self.divider.add_nuclide("W183", 0.00003736630, "ao")
        self.divider.add_nuclide("W184", 0.00008000720, "ao")
        self.divider.add_nuclide("W186", 0.00007423640, "ao")
        self.divider.set_density("atom/b-cm", 0.04312289441)

        #--INNER MANIFOLD-- 
        self.inner_manifold = openmc.Material(name="inner_manifold", temperature=self.temp_k)
        self.inner_manifold.depletable = False
        self.inner_manifold.add_nuclide("He4", 0.00032820000, "ao")
        self.inner_manifold.add_nuclide("C12", 0.00018120000, "ao")
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

        breeder_mass_frac, fertile_compound_mass_frac = calc_biso_blanket_mass_fracs(self.fertile_bulk_density_kgm3,
                                                                    self.breeder_volume,
                                                                    fertile_element=self.fertile_element,
                                                                    fertile_enrich=ENRICH_U if self.fertile_element == 'U' else 0.0,
                                                                    breeder_density_kgm3=DENSITY_LL*1e3)
        self.blanket = openmc.Material.mix_materials([self.breeder, self.fertile], [breeder_mass_frac, fertile_compound_mass_frac], 'wo')
        self.blanket.name, self.blanket.temperature = self.name, self.temp_k


        # ------------------------------------------------------------------
        # Add materials 
        # ------------------------------------------------------------------

        self.materials = openmc.Materials([self.firstwall, self.structure, self.coolant, self.divider, self.inner_manifold, self.blanket]) # check if need to add self.backplate


    def geometry(self):

        # ------------------------------------------------------------------
        # Tokamak geometry parameters 
        # ------------------------------------------------------------------

        d_fw   = LL_FW_O_CM
        d_fwf  = d_fw   + LL_FWF_O_CM  # front wall front
        d_fwc  = d_fwf  + LL_FWC_O_CM  # front wall channel
        d_fwb  = d_fwc  + LL_FWB_O_CM  # front wall back
        d_br1  = d_fwb  + LL_BR1_O_CM  # breeding region 1
        d_d1   = d_br1  + LL_D1_O_CM   # divider 1 
        d_br2  = d_d1   + LL_BR2_O_CM  # breeding region 2
        
        d_d2_o  = d_br2    + LL_D2_O_CM    # only on outboard blanket
        d_br3_o = d_d2_o   + LL_BR3_O_CM   # only on outboard blanket
        d_im_o  = d_br3_o  + LL_IM_O_CM    # only on outboard blanket
        
        d_im_i  = d_br2  + LL_IM_I_CM   # inboard blanket

        self.extent_r = (self.R0 + self.a + d_im_o)  * 1.2  # 120% radial extent
        self.extent_z = (self.kappa*self.a + d_im_o) * 1.2  # 120% vertical extent


        # ------------------------------------------------------------------
        # Surfaces 
        # ------------------------------------------------------------------

        points_vc   =  miller_model(self.R0, self.a, self.kappa, self.delta)
        points_fw   = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_fw)
        points_fwf  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_fwf)
        points_fwc  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_fwc)
        points_fwb  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_fwb)
        points_br1  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br1)
        points_d1   = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_d1)
        points_br2  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br2)
        points_im_i = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_im_i)   # inboard manifold
        points_d2_o  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_d2_o)  # outboard divider
        points_br3_o = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_br3_o) # outboard blanket module
        points_im_o  = miller_model(self.R0, self.a, self.kappa, self.delta, extrude=d_im_o)  # outboard manifold
        
        # Create OpenMC surfaces
        self.surface_vc    = openmc.model.Polygon(points_vc , basis='rz')
        self.surface_fw    = openmc.model.Polygon(points_fw, basis='rz')
        self.surface_fwf   = openmc.model.Polygon(points_fwf, basis='rz')
        self.surface_fwc   = openmc.model.Polygon(points_fwc, basis='rz')
        self.surface_fwb   = openmc.model.Polygon(points_fwb, basis='rz')
        self.surface_br1   = openmc.model.Polygon(points_br1, basis='rz')
        self.surface_d1    = openmc.model.Polygon(points_d1,  basis='rz')
        self.surface_br2   = openmc.model.Polygon(points_br2, basis='rz')
        self.surface_im_i  = openmc.model.Polygon(points_im_i,  basis='rz') # inboard manifold
        self.surface_d2_o  = openmc.model.Polygon(points_d2_o,  basis='rz') # outboard divider
        self.surface_br3_o = openmc.model.Polygon(points_br3_o, basis='rz') # outboard blanket module
        self.surface_im_o  = openmc.model.Polygon(points_im_o,  basis='rz') # outboard manifold

        dividing_cylinder = openmc.ZCylinder(r=self.R0)
        outer_cylinder    = openmc.ZCylinder(r=self.extent_r, boundary_type='vacuum')
        top_plane         = openmc.ZPlane(z0=self.extent_z, boundary_type='vacuum')
        bottom_plane      = openmc.ZPlane(z0=-self.extent_z, boundary_type='vacuum')
        

        # ------------------------------------------------------------------
        # Cells | 10: vc, fw | 2X: structure | 3X: breeding
        # ------------------------------------------------------------------

        cell_vc   = openmc.Cell(cell_id=10, region= -self.surface_vc)
        cell_vc.importance = {'neutron':1}

        cell_fw   = openmc.Cell(cell_id=11, region= +self.surface_vc  & -self.surface_fw,  fill=self.firstwall)
        cell_fwf  = openmc.Cell(cell_id=21, region= +self.surface_fw  & -self.surface_fwf, fill=self.structure)
        cell_fwc  = openmc.Cell(cell_id=12, region= +self.surface_fwf & -self.surface_fwc, fill=self.coolant)
        cell_fwb  = openmc.Cell(cell_id=22, region= +self.surface_fwc & -self.surface_fwb, fill=self.structure)
        cell_br1  = openmc.Cell(cell_id=31, region= +self.surface_fwb & -self.surface_br1, fill=self.blanket)
        cell_d1   = openmc.Cell(cell_id=23, region= +self.surface_br1 & -self.surface_d1,  fill=self.divider)
        cell_br2  = openmc.Cell(cell_id=32, region= +self.surface_d1  & -self.surface_br2, fill=self.blanket)
        cell_im_i   = openmc.Cell(cell_id=24, region= -dividing_cylinder & +self.surface_br2 & -self.surface_im_i, fill=self.inner_manifold)    # inboard manifold
        cell_d2_o   = openmc.Cell(cell_id=25, region= +dividing_cylinder & +self.surface_br2 & -self.surface_d2_o,   fill=self.divider)         # outboard divider
        cell_br3_o  = openmc.Cell(cell_id=33, region= +dividing_cylinder & +self.surface_d2_o & -self.surface_br3_o, fill=self.blanket)         # outboard blanket module
        cell_im_o   = openmc.Cell(cell_id=26, region= +dividing_cylinder & +self.surface_br3_o & -self.surface_im_o, fill=self.inner_manifold)  # outboard manifold

        # Void cells
        cell_void_i = openmc.Cell(cell_id=99, region= +self.surface_im_i & -dividing_cylinder & -outer_cylinder & +bottom_plane & -top_plane)
        cell_void_o = openmc.Cell(cell_id=98, region= +self.surface_im_o & +dividing_cylinder & -outer_cylinder & +bottom_plane & -top_plane)
        cell_void_i.importance = {'neutron': 0}
        cell_void_o.importance = {'neutron': 0}


        # -------------------------------------------------------
        # Collect all cells
        # -------------------------------------------------------
        self.cells = [
            cell_vc,
            cell_fw, cell_fwf, cell_fwc, cell_fwb, cell_br1, cell_d1, cell_br2, 
            cell_im_i,
            cell_d2_o, cell_br3_o, cell_im_o,
            cell_void_i, cell_void_o,
        ]
        self.geometry = openmc.Geometry(openmc.Universe(cells=self.cells))