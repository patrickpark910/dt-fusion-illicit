import numpy as np
import openmc

from Python.reactor    import *
from Python.parameters import *
from Python.utilities  import *


class PBHet(Reactor):

    def initialize(self):

        self.temp_k          = TEMP_K
        self.breeder_name    = 'PBHet'  # Changed from 'PB' to distinguish from PBHomo
        self.breeder_enrich  = ENRICH_PB  # at% 
        self.breeder_volume  = PBCoupon_BR_VOL  

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
        li4sio4.add_elements_from_formula('Li4SiO4', enrichment_target='Li6', enrichment_type='ao', enrichment=self.breeder_enrich) 
        # self.lithium_ceramic.add_elements('Li', 22.415, percent_type='wo', enrichment_target='Li6', enrichment_type='ao', enrichment=ENRICH_PB) 
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
        be.set_density('g/cm3', 1.80) 
        be.add_element('Be', 1.0, percent_type='ao')
        # normalized from 1.85 for 900 K
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
        # --- Kernel volume [cm^3]
        V_k_cm3 = (4.0/3.0) * np.pi * (BISO_KERNEL_RADIUS**3)

        # --- Outer particle volume [cm^3] (for packing fraction)
        V_p_cm3 = (4.0/3.0) * np.pi * (BISO_RADIUS**3)
       
        # --- Kernel mass [g]
        if self.fertile_element == 'U':
            rho_k = DENSITY_UO2  # g/cm3
            m_k_g = rho_k * V_k_cm3

            # mass fraction of U in UO2 (enrichment only changes U molar mass slightly)
            fert_enrich = ENRICH_U * 0.01
            M_U = (1-fert_enrich)*AMU_U238 + fert_enrich*AMU_U235
            M_UO2 = M_U + 2.0*AMU_O
            w_U_in_UO2 = M_U / M_UO2

            m_fertile_per_particle_g = m_k_g * w_U_in_UO2

        elif self.fertile_element == 'Th':
            rho_k = DENSITY_ThO2
            m_k_g = rho_k * V_k_cm3
            M_ThO2 = AMU_Th + 2.0*AMU_O
            w_Th_in_ThO2 = AMU_Th / M_ThO2
            m_fertile_per_particle_g = m_k_g * w_Th_in_ThO2
        
        self.V_in  = 151.678656e-6    # m3 volume of inner blanket frustum
        self.V_out = 368.481344e-6    # m3 vol of outer blanket frustum (new fave vocab word)
    
        rho_fert_kg_m3 = self.fertile_bulk_density_kgm3
        m_in_g  = rho_fert_kg_m3 * (self.V_in * 1e3)
        m_out_g = rho_fert_kg_m3 * (self.V_out * 1e3)
        
        # mass of fertile needed into particle count
        if m_in_g <= 0.0 or m_fertile_per_particle_g <= 0.0:
            self.N_in = 0
        else:
            self.N_in = int(np.rint(m_in_g / m_fertile_per_particle_g))

        if m_out_g <= 0.0 or m_fertile_per_particle_g <= 0.0:
            self.N_out = 0
        else:
            self.N_out = int(np.rint(m_out_g / m_fertile_per_particle_g))
        
        # packing fractions
        pf_in  = (self.N_in  * V_p_cm3) / self.V_in 
        pf_out = (self.N_out * V_p_cm3) / self.V_out
        print("N_in =", self.N_in, " N_out =", self.N_out)
        print("pf_in =", pf_in, " pf_out =", pf_out)

        if self.fertile_element == 'U':
            self.kernel = openmc.Material(name='UO2', temperature=self.temp_k)
            self.kernel.add_elements_from_formula('UO2', enrichment=ENRICH_U)
            self.kernel.set_density('g/cm3', DENSITY_UO2) 
            self.kernel.depletable = False 

        elif self.fertile_element == 'Th': 
            self.kernel = openmc.Material(name='ThO2', temperature=self.temp_k)
            self.kernel.add_elements_from_formula('ThO2', enrichment=ENRICH_U)
            self.kernel.set_density('g/cm3', DENSITY_ThO2) 
            self.kernel.depletable = False 

        # SiC coating
        self.sic = openmc.Material(name='SiC', temperature=self.temp_k)
        self.sic.add_elements_from_formula('SiC')
        self.sic.set_density('g/cm3', 3.2)
        self.sic.depletable = False 
        
        # Volume fractions from Lu (2017) Table 2 
        vf_li4sio4 = 0.1304 ; vf_be = 0.3790 ; vf_eurofer = 0.1176 ; vf_he = 1 - (vf_li4sio4 + vf_be + vf_eurofer) 
        print("Li4SiO4 nuclides:", len(li4sio4.nuclides), li4sio4.nuclides[:3])
        print("Be nuclides:", len(be.nuclides), be.nuclides[:3])
        print("vf sum:", vf_li4sio4 + vf_be)

        # Mix Li4SiO4 and Be (should be 25.6, 74.4 vol% respectively)
        breeder = openmc.Material.mix_materials([li4sio4, be], [vf_li4sio4/(vf_li4sio4+vf_be), vf_be/(vf_li4sio4+vf_be)], 'vo') 

        # So now we mix in Li4SiO4-Be-BISO with the Eurofer structure and He coolant materials
        # Li4SiO4-Be-BISO should be the volume fractions of Li4SiO4 (0.1304) + Be (0.3790) = 0.5094
        self.blanket = openmc.Material.mix_materials([breeder, self.structure, he], [(vf_li4sio4+vf_be), vf_eurofer, vf_he], 'vo') 
        self.blanket.name, self.blanket.temperature = self.name, self.temp_k


        # ------------------------------------------------------------------
        # Add materials 
        # ------------------------------------------------------------------

        self.materials = openmc.Materials([self.firstwall, self.structure, self.blanket, self.kernel, self.sic]) 


    def geometry(self):

        #  geometry parameters -----------------------
        d1 = 365.0   # cm
        h1 = 45.0    # cm

        #total surface areas of inner and outer MILLER blanket volumes
        # inner blanket (both sides) = 6.850857e6 cm2
        # outer blanket = 9.135564e6 cm2
        # V(cuboid) = h * a^2 
        # based on volumes a1 = 54.86cm; a2 = 68.54cm
        # based on surface area ratios a1 = 58.05cm; 67.04cm
        #    based on SA Vinner=29.16%; Vouter=70.80% of 520.16cm2
        a1 = 58.05 #cm inboard blanket cuboid side length
        a2 = 67.04 #cm outboard blanket cuboid side length

        d_out = d1 + 465.0 
        h2 = 82.0


        PB_st1_ex            = d1 - PB_ST2_CM
        PB_BR1_exterior      = d1
        PB_BR1_plasmafacing  = d1 + h1
        PB_st1_pf            = PB_BR1_plasmafacing + PB_ST1_CM
        PB_fw1               = PB_st1_pf + PB_FW_CM

        PB_fw2               = PB_st1_pf - PB_FW_CM
        PB_st2_pf            = PB_BR1_plasmafacing - PB_ST1_CM
        PB_BR2_plasmafacing  = d_out
        PB_BR2_exterior      = d_out + h2
        PB_st2_ex            = PB_BR2_exterior + PB_ST2_CM
    
        # ------------------------------------------------------------------
        # Surfaces 
        # ------------------------------------------------------------------
        # X-direction: transmissive (neutrons can enter/exit front and back)
        x0i = openmc.XPlane(x0=PB_st1_ex, boundary_type='transmission')
        x1i = openmc.XPlane(x0=PB_fw1,   boundary_type='transmission')
        # Y and Z directions: REFLECTIVE (simulates infinite extent, like coupon test)
        y0i = openmc.YPlane(y0=-a1/2,    boundary_type='reflective')
        y1i = openmc.YPlane(y0=+a1/2,    boundary_type='reflective')
        z0i = openmc.ZPlane(z0=-a1/2,    boundary_type='reflective')
        z1i = openmc.ZPlane(z0=+a1/2,    boundary_type='reflective')

        self.inner_cuboid = (+x0i & -x1i & +y0i & -y1i & +z0i & -z1i)

        # X-direction: transmissive (neutrons can enter/exit front and back)
        x0o = openmc.XPlane(x0=PB_fw2,     boundary_type='transmission')
        x1o = openmc.XPlane(x0=PB_st2_ex,  boundary_type='transmission')
        # Y and Z directions: REFLECTIVE (simulates infinite extent, like coupon test)
        y0o = openmc.YPlane(y0=-a2/2,      boundary_type='reflective')
        y1o = openmc.YPlane(y0=+a2/2,      boundary_type='reflective')
        z0o = openmc.ZPlane(z0=-a2/2,      boundary_type='reflective')
        z1o = openmc.ZPlane(z0=+a2/2,      boundary_type='reflective')

        self.outer_cuboid = (+x0o & -x1o & +y0o & -y1o & +z0o & -z1o)

        p2_BL = (PB_fw2, -a2/2 , -a2/2)
        p2_BR = (PB_fw2, +a2/2 , -a2/2)
        p2_TL = (PB_fw2, -a2/2 , +a2/2)
        p2_TR = (PB_fw2, +a2/2 , +a2/2)
        p1_BL = (PB_fw1, -a2/2 , -a2/2)
        p1_BR = (PB_fw1, +a2/2 , -a2/2)
        p1_TL = (PB_fw1, -a2/2 , +a2/2)
        p1_TR = (PB_fw1, +a2/2 , +a2/2)
        y0p_B = openmc.Plane.from_points(p2_BL, p2_BR, p1_BL, boundary_type='reflective')
        y0p_L = openmc.Plane.from_points(p2_BL, p2_TL, p1_BL, boundary_type='reflective')
        y0p_T = openmc.Plane.from_points(p2_TL, p2_TR, p1_TL, boundary_type='reflective')
        y0p_R = openmc.Plane.from_points(p2_BR, p2_TR, p1_BR, boundary_type='reflective')

        self.plasma_region = (+y0p_B & -y0p_R & +y0p_L & -y0p_T & +x1i & -x0o)


        # --- Inner breeder packing box planes
        x0_in = openmc.XPlane(x0=PB_BR1_exterior)       # inner breeder exterior
        x1_in = openmc.XPlane(x0=PB_BR1_plasmafacing)   # inner breeder plasma-facing
        y0_in = openmc.YPlane(y0=-a1/2)
        y1_in = openmc.YPlane(y0=+a1/2)
        z0_in = openmc.ZPlane(z0=-a1/2)
        z1_in = openmc.ZPlane(z0=+a1/2)

        self.inner_blanket = (+x0_in & -x1_in & +y0_in & -y1_in & +z0_in & -z1_in)
        
        x0_out = openmc.XPlane(x0=PB_BR2_plasmafacing)
        x1_out = openmc.XPlane(x0=PB_BR2_exterior)
        y0_out = openmc.YPlane(y0=-a2/2)
        y1_out = openmc.YPlane(y0=+a2/2)
        z0_out = openmc.ZPlane(z0=-a2/2)
        z1_out = openmc.ZPlane(z0=+a2/2)

        self.outer_blanket = (+x0_out & -x1_out & +y0_out & -y1_out & +z0_out & -z1_out)

        # --- Plasma cell
        cell_vc = openmc.Cell(cell_id=60, region=self.plasma_region, fill=None)
        cell_vc.importance = {'neutron': 1}
        
        self.surface_st1_ex         = openmc.XPlane(x0=PB_st1_ex)
        self.surface_br1_exterior   = openmc.XPlane(x0=PB_BR1_exterior)
        self.surface_br1_pf         = openmc.XPlane(x0=PB_BR1_plasmafacing)
        self.surface_st1_pf         = openmc.XPlane(x0=PB_st1_pf)
        self.surface_fw1            = openmc.XPlane(x0=PB_fw1)

        self.surface_fw2            = openmc.XPlane(x0=PB_fw2)
        self.surface_st2_pf         = openmc.XPlane(x0=PB_st2_pf)
        self.surface_br2_pf         = openmc.XPlane(x0=PB_BR2_plasmafacing)
        self.surface_br2_exterior   = openmc.XPlane(x0=PB_BR2_exterior)
        self.surface_st2_ex         = openmc.XPlane(x0=PB_st2_ex)
      
        # ----------Cells----------------

        cell_fw1   = openmc.Cell(cell_id=11, region=self.inner_cuboid & -self.surface_fw1 & +self.surface_st1_pf,  fill=self.firstwall)
        cell_st1_pf   = openmc.Cell(cell_id=21, region=self.inner_cuboid & -self.surface_st1_pf & +self.surface_br1_pf, fill=self.structure)
        cell_st1_ex   = openmc.Cell(cell_id=41, region=self.inner_cuboid & -self.surface_br1_exterior & +self.surface_st1_ex, fill=self.structure)

        cell_fw2    = openmc.Cell(cell_id=12, region=self.outer_cuboid & +self.surface_fw2 & -self.surface_st2_pf,  fill=self.firstwall)
        cell_st2_pf   = openmc.Cell(cell_id=22, region=self.outer_cuboid & +self.surface_st2_pf  & -self.surface_br2_pf, fill=self.structure)
        cell_st2_ex   = openmc.Cell(cell_id=42, region=self.outer_cuboid & +self.surface_br2_exterior & -self.surface_st2_ex, fill=self.structure)
     

        # ---------BISO Packing--------
        if self.N_in <=0:
            cell_insert_in  = openmc.Cell(cell_id=31,  fill=self.blanket, region=self.inner_blanket)
            cell_insert_out = openmc.Cell(cell_id=32, fill=self.blanket, region=self.outer_blanket)

        else:
            s_k = openmc.Sphere(r=BISO_KERNEL_RADIUS)
            s_o = openmc.Sphere(r=BISO_RADIUS)

            c_kernel  = openmc.Cell(fill=self.kernel, region=-s_k)
            c_coating = openmc.Cell(fill=self.sic,    region=+s_k & -s_o)

            self.u_biso = openmc.Universe(name='BISO_universe', cells=[c_kernel, c_coating])
      
            inner_centers = openmc.model.pack_spheres(radius=BISO_RADIUS, region=self.inner_blanket, num_spheres=self.N_in, seed=1)
            outer_centers = openmc.model.pack_spheres(radius=BISO_RADIUS, region=self.outer_blanket, num_spheres=self.N_out, seed=2)
            # Background cells in the insert regions (homogenized blanket, but with sphere voids carved out)
            cell_bg_in  = openmc.Cell(name="bg_in",  fill=self.blanket, region=self.inner_blanket)
            cell_bg_out = openmc.Cell(name="bg_out", fill=self.blanket, region=self.outer_blanket)
            biso_cells_in = []
            for i, (x, y, z) in enumerate(inner_centers):
                s = openmc.Sphere(x0=float(x), y0=float(y), z0=float(z), r=BISO_RADIUS)
                c = openmc.Cell(name=f"biso_in_{i}", fill=self.u_biso, region=-s)
                biso_cells_in.append(c)
                cell_bg_in.region &= +s   # keep only outside each sphere in the background

            biso_cells_out = []
            for i, (x, y, z) in enumerate(outer_centers):
                s = openmc.Sphere(x0=float(x), y0=float(y), z0=float(z), r=BISO_RADIUS)
                c = openmc.Cell(name=f"biso_out_{i}", fill=self.u_biso, region=-s)
                biso_cells_out.append(c)
                cell_bg_out.region &= +s

            # Universes that contain background + explicit BISO spheres
            u_insert_in  = openmc.Universe(name="insert_in",  cells=[cell_bg_in]  + biso_cells_in)
            u_insert_out = openmc.Universe(name="insert_out", cells=[cell_bg_out] + biso_cells_out)

            # Cells that place those universes into the model
            cell_insert_in  = openmc.Cell(cell_id=31,  fill=u_insert_in,  region=self.inner_blanket)
            cell_insert_out = openmc.Cell(cell_id=32, fill=u_insert_out, region=self.outer_blanket)


        self.cells = [cell_vc,
                      cell_fw1, cell_st1_pf, cell_insert_in, cell_st1_ex, 
                      cell_fw2, cell_st2_pf, cell_insert_out, cell_st2_ex]
    
        self.geometry = openmc.Geometry(openmc.Universe(cells=self.cells))
