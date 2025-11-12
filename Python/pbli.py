import numpy as np
import openmc

from Python.reactor    import *
from Python.parameters import *
from Python.utilities  import *


class LL(Reactor):

    def initialize(self):

        self.temp_k          = TEMP_K
        self.breeder_name    = 'LL'
        self.breeder_density = DENSITY_LL # [g/cm³]
        self.breeder_enrich  = ENRICH_LL  # [at%]
        self.breeder_volume  = LL_BR_VOL  # [m³]
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
        # - The original first wall specs we were using from JL Ball (2025) 
        #   is 99.9969 wt% W and the rest impurities. Changing to nominal W 
        #   to isolate fertile reaction effects and save compute time   
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
        # - F82H steel from T. Hirose (2014) Fus Engr Des
        # ------------------------------------------------------------------

        # F82H steel (first wall front, back + part of other structural components)
        self.f82h = openmc.Material(name='f82h', temperature=self.temp_k)
        self.f82h.depletable = False
        self.f82h.add_element('Fe', 89.3686, percent_type='wo')
        self.f82h.add_element('C' , 0.10, percent_type='wo')
        self.f82h.add_element('Si', 0.10, percent_type='wo')
        self.f82h.add_element('Mn', 0.13, percent_type='wo')
        self.f82h.add_element('Cr', 8.16, percent_type='wo')
        self.f82h.add_element('W' , 1.94, percent_type='wo')
        self.f82h.add_element('V' , 0.20, percent_type='wo')
        self.f82h.add_element('N' , 0.0014, percent_type='wo')
        self.f82h.set_density("g/cm3", 7.6)  # from MCNP cell 112


        # ------------------------------------------------------------------
        # Coolant (17 vol% F82H + 83 vol% He-4)
        # ------------------------------------------------------------------

        # Helium-4 gas 
        he = openmc.Material(name='helium') 
        he.set_density('g/cm3', 0.004279) # Helium density at 900 K ~80 bar 
        he.add_element('He', 1) 

        self.coolant = openmc.Material.mix_materials([self.f82h, he], [0.170, 0.830], 'vo')
        # self.coolant.set_density("atom/b-cm", 0.01471892350) # from MCNP
        self.coolant.name, self.coolant.temperature = "coolant (17.0 vol% F82H, 83.0 vol% He-4)", self.temp_k


        # ------------------------------------------------------------------
        # Divider material
        # ------------------------------------------------------------------
        
        self.divider = openmc.Material.mix_materials([self.f82h, he], [0.512, 0.488], 'vo')
        self.divider.depletable = False
        self.divider.name, self.divider.temperature = "divider (51.2 vol% F82H, 48.8 vol% He-4)", self.temp_k
        # self.divider.set_density("atom/b-cm", 0.04312289441)


        # ------------------------------------------------------------------
        # Manifold material
        # ------------------------------------------------------------------

        self.manifold = openmc.Material.mix_materials([self.f82h, he], [0.453, 0.547], 'vo')
        self.manifold.depletable = False
        self.manifold.name, self.manifold.temperature = "inner manifold (45.3 vol% F82H, 54.7 vol% He-4)", self.temp_k
        # self.manifold.set_density("atom/b-cm", 0.03822271948)


        # ------------------------------------------------------------------
        # Breeder material
        # ------------------------------------------------------------------

        pbli = openmc.Material(name='breeder', temperature=self.temp_k)
        pbli.set_density('g/cm3', self.breeder_density)
        pbli.add_element('Pb', 0.83, percent_type='ao') 
        pbli.add_element('Li', 0.17, percent_type='ao', enrichment_target='Li6', enrichment_type='ao', enrichment=self.breeder_enrich) # Li-6 enrichment to 90 at% 


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
        # - changed fertile kg/m³ to be per Pb-Li vol, not whole breeder region vol
        # 
        # We want "fertile bulk density" to be kg of U-238 per m³ of *breeding* material
        # so we mix Pb-Li + BISO, and then mix Pb-Li-BISO with F82H, SiC, and He-4
        # ------------------------------------------------------------------

        # Volume fractions from Glaser & Goldston (2012) - Table 1, Breeding Ch.2
        vf_fs = 0.019; vf_ll = 0.805; vf_sic = 0.079; vf_he = 0.097 # til ; can separate statements in python --ppark
        
        # Mix Pb-Li with BISO
        vol_ll = self.breeder_volume*vf_ll
        
        vf_pbli, vf_biso = calc_biso_blanket_vol_fracs(self.fertile_bulk_density_kgm3, vol_ll, fertile_element=self.fertile_element)
        
        breeder = openmc.Material.mix_materials([pbli, biso], [vf_pbli, vf_biso], 'vo')
        # DO NOT CONFUSE vf_pbli WITH vf_ll !!! -- ppark 2025-11-11
        
        # So now we mix in Li4SiO4-Be-BISO with the Eurofer structure and He coolant materials
        # Li4SiO4-Be-BISO should be the volume fractions of Li4SiO4 (0.1304) + Be (0.3790) = 0.5094
        self.blanket = openmc.Material.mix_materials([self.f82h, breeder, sic, he], [vf_fs, vf_ll, vf_sic, vf_he], 'vo') 
        self.blanket.name, self.blanket.temperature = self.name, self.temp_k    


        # ------------------------------------------------------------------
        # Add materials 
        # ------------------------------------------------------------------

        self.materials = openmc.Materials([self.firstwall, self.f82h, self.coolant, self.blanket, self.divider, self.manifold])


    def geometry(self):

        # ------------------------------------------------------------------
        # Tokamak geometry parameters 
        # ------------------------------------------------------------------

        d_fw   = LL_FW_CM
        d_fwf  = d_fw   + LL_FWF_CM  # front wall front
        d_fwc  = d_fwf  + LL_FWC_CM  # front wall cooling channel
        d_fwb  = d_fwc  + LL_FWB_CM  # front wall back
        d_br1  = d_fwb  + LL_BR1_CM  # breeding region 1
        d_d1   = d_br1  + LL_D1_CM   # divider 1 
        d_br2  = d_d1   + LL_BR2_CM  # breeding region 2
        
        d_d2_o  = d_br2    + LL_D2_O_CM    # only on outboard blanket
        d_br3_o = d_d2_o   + LL_BR3_O_CM   # only on outboard blanket
        d_im_o  = d_br3_o  + LL_IM_CM    # only on outboard blanket
        
        d_im_i  = d_br2  + LL_IM_CM   # inboard blanket

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
        cell_fwf  = openmc.Cell(cell_id=21, region= +self.surface_fw  & -self.surface_fwf, fill=self.f82h)
        cell_fwc  = openmc.Cell(cell_id=12, region= +self.surface_fwf & -self.surface_fwc, fill=self.coolant)
        cell_fwb  = openmc.Cell(cell_id=22, region= +self.surface_fwc & -self.surface_fwb, fill=self.f82h)
        cell_br1  = openmc.Cell(cell_id=31, region= +self.surface_fwb & -self.surface_br1, fill=self.blanket)
        cell_d1   = openmc.Cell(cell_id=23, region= +self.surface_br1 & -self.surface_d1,  fill=self.divider)
        cell_br2  = openmc.Cell(cell_id=32, region= +self.surface_d1  & -self.surface_br2, fill=self.blanket)
        cell_im_i   = openmc.Cell(cell_id=24, region= -dividing_cylinder & +self.surface_br2 & -self.surface_im_i, fill=self.manifold)    # inboard manifold
        cell_d2_o   = openmc.Cell(cell_id=25, region= +dividing_cylinder & +self.surface_br2 & -self.surface_d2_o,   fill=self.divider)   # outboard divider
        cell_br3_o  = openmc.Cell(cell_id=33, region= +dividing_cylinder & +self.surface_d2_o & -self.surface_br3_o, fill=self.blanket)   # outboard blanket module
        cell_im_o   = openmc.Cell(cell_id=26, region= +dividing_cylinder & +self.surface_br3_o & -self.surface_im_o, fill=self.manifold)  # outboard manifold

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