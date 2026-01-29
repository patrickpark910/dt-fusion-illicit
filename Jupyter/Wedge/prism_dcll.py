import os
import numpy as np
import openmc


class Prism():

    def __init__(self, case, fertile_kgm3, isotope='U238', write_openmc=True, run_openmc=False,):

        self.case = case
        self.fertile_kgm3 = fertile_kgm3
        self.fertile_str = f"{fertile_kgm3:06.2f}"
        self.isotope = isotope 
        self.name = f"prism_{self.case}_{self.isotope}_{self.fertile_str}kgm3"         
        self.path = f"./OpenMC/{self.name}"

        os.makedirs(self.path, exist_ok=True)


    def materials(self, debug=False):

        # Tungsten | First wall
        self.firstwall = openmc.Material(name='firstwall', temperature=self.temp_k)
        self.firstwall.set_density('g/cm3',19.3)
        self.tungsten.add_element('W',1)

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

        # He-4 (gas) 
        he = openmc.Material(name='helium') 
        he.set_density('g/cm3', 0.004) # Helium density at 900 K ~80 bar 
        he.add_element('He', 1) 

        # Breeder material
        pbli = openmc.Material(name='breeder', temperature=self.temp_k)
        pbli.set_density('g/cm3', DENSITY_DCLL)
        pbli.add_element('Pb', 0.83, percent_type='ao') 
        pbli.add_element('Li', 0.17, percent_type='ao', enrichment_target='Li6', enrichment_type='ao', enrichment=self.breeder_enrich) # Li-6 enrichment to 90 at% 

        # Coolant 
        self.coolant = openmc.Material.mix_materials([self.f82h, he], [0.170, 0.830], 'vo')
        self.coolant.set_density('atom/b-cm', 0.01471892350) # from Glaser et al. (2025) MCNP
        self.coolant.temperature = self.temp_k
        self.coolant.name = "coolant (17.0 vol% F82H, 83.0 vol% He-4)" 

        self.divider = openmc.Material.mix_materials([self.f82h, he], [0.512, 0.488], 'vo')
        self.divider.temperature = self.temp_k
        self.divider.name        = "divider (51.2 vol% F82H, 48.8 vol% He-4)"

        self.manifold = openmc.Material.mix_materials([self.f82h, he], [0.453, 0.547], 'vo')
        self.manifold.temperature = self.temp_k
        self.manifold.name        = "inner manifold (45.3 vol% F82H, 54.7 vol% He-4)"

        self.shield = openmc.Material.mix_materials([self.f82h, he], [0.80, 0.20], 'vo')
        self.shield.temperature = self.temp_k
        self.shield.name = "steel shield (80 vol% F82H, 20 vol% He-4)"


        ''' BISO particle '''
        # Fertile material
        if self.isotope == 'U238':
            self.kernel = openmc.Material(name='UO2', temperature=TEMP_K)
            self.kernel.add_elements_from_formula('UO2')
            self.kernel.set_density('g/cm3', DENSITY_UO2)  
            # self.kernel.add_s_alpha_beta('c_U_in_UO2')  # NotImplementedError
            # self.kernel.add_s_alpha_beta('c_O_in_UO2')  # NotImplementedError

        elif self.isotope == 'Th232': 
            self.kernel = openmc.Material(name='ThO2', temperature=TEMP_K) 
            self.kernel.add_elements_from_formula('ThO2') 
            self.kernel.set_density('g/cm3', DENSITY_ThO2)   

        # SiC coating
        self.sic = openmc.Material(name='SiC', temperature=TEMP_K)
        self.sic.add_elements_from_formula('SiC')
        self.sic.set_density('g/cm3', 3.2)
        # self.sic.add_s_alpha_beta('c_Si_in_SiC')  # NotImplementedError


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


    def geometry(self, debug=False):

        d_0    = DCLL_A
        d_fw   = d_0    + DCLL_FW_CM   # first wall 
        d_fwf  = d_fw   + DCLL_FWF_CM  # first wall front
        d_fwc  = d_fwf  + DCLL_FWC_CM  # first wall cooling channel
        d_fwb  = d_fwc  + DCLL_FWB_CM  # first wall back
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

        # -----------------------------------
        # SURFACES
        # -----------------------------------

        # X-planes
        s111 = openmc.XPlane(surface_id=111, x0= -self.geom['c1'], boundary_type="reflective")
        s112 = openmc.XPlane(surface_id=112, x0= +self.geom['c1'], boundary_type="reflective")
        s113 = openmc.XPlane(surface_id=113, x0= -self.geom['c2'], boundary_type="reflective")
        s114 = openmc.XPlane(surface_id=114, x0= +self.geom['c2'], boundary_type="reflective")

        # Y-planes
        s121 = openmc.YPlane(surface_id=121, y0= -self.geom['c1'], boundary_type="reflective")
        s122 = openmc.YPlane(surface_id=122, y0= +self.geom['c1'], boundary_type="reflective")
        s123 = openmc.YPlane(surface_id=123, y0= -self.geom['c2'], boundary_type="reflective")
        s124 = openmc.YPlane(surface_id=124, y0= +self.geom['c2'], boundary_type="reflective")

        # Z-planes - inboard
        s230 = openmc.ZPlane(surface_id=230, z0=d_fw  )  # first wall 
        s231 = openmc.ZPlane(surface_id=231, z0=d_fw  )  # first wall 
        s232 = openmc.ZPlane(surface_id=232, z0=d_fwf )  # first wall front
        s233 = openmc.ZPlane(surface_id=233, z0=d_fwc )  # first wall cooling channel
        s234 = openmc.ZPlane(surface_id=234, z0=d_fwb )  # first wall back
        s235 = openmc.ZPlane(surface_id=235, z0=d_br1 )  # breeding region 1
        s236 = openmc.ZPlane(surface_id=236, z0=d_d1  )  # divider 1 
        s237 = openmc.ZPlane(surface_id=237, z0=d_br2 )  # breeding region 2
        s238 = openmc.ZPlane(surface_id=238, z0=d_im_i)  # inner manifold - inboard
        s239 = openmc.ZPlane(surface_id=239, z0=d_bp_i)  # back plate     - inboard
        s240 = openmc.ZPlane(surface_id=240, z0=d_ss_i)  # steel shield   - inboard

        # Z-planes - outboard
        s330 = openmc.ZPlane(surface_id=330, z0=d_fw  )  # first wall 
        s331 = openmc.ZPlane(surface_id=331, z0=d_fw  )  # first wall 
        s332 = openmc.ZPlane(surface_id=332, z0=d_fwf )  # first wall front
        s333 = openmc.ZPlane(surface_id=333, z0=d_fwc )  # first wall cooling channel
        s334 = openmc.ZPlane(surface_id=334, z0=d_fwb )  # first wall back
        s335 = openmc.ZPlane(surface_id=335, z0=d_br1 )  # breeding region 1
        s336 = openmc.ZPlane(surface_id=336, z0=d_d1  )  # divider 1 
        s337 = openmc.ZPlane(surface_id=337, z0=d_br2 )  # breeding region 2
        s338 = openmc.ZPlane(surface_id=338, z0=d_d2_o)  # divider 2         - only on outboard 
        s339 = openmc.ZPlane(surface_id=339, z0=d_br3_o) # breeding region 3 - only on outboard 
        s330 = openmc.ZPlane(surface_id=340, z0=d_im_o)  # inner manifold    - outboard  
        s341 = openmc.ZPlane(surface_id=341, z0=d_bp_o)  # back plate        - outboard
        s342 = openmc.ZPlane(surface_id=342, z0=d_ss_o)  # steel shield      - outboard

        # Plasma chamber boundaries
        # Plane normal card entries: Ax + By + Cz - D = 0
        s141 = openmc.Plane(surface_id=141, a= +1.0, b=  0.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s142 = openmc.Plane(surface_id=142, a= -1.0, b=  0.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s143 = openmc.Plane(surface_id=143, a=  0.0, b= +1.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")
        s144 = openmc.Plane(surface_id=144, a=  0.0, b= -1.0, c=self.geom['pC'], d=self.geom['pD'], boundary_type="white")

        # BISO universe
        s150 = openmc.Sphere(surface_id=150, r=0.040) # Kernel
        s151 = openmc.Sphere(surface_id=151, r=0.050) # Outer Coating


        # -----------------------------------
        # CELLS
        # -----------------------------------

        c10 = openmc.Cell(cell_id=10, region= -s142 & -s141 & -s144 & -s143 & +s134 & -s135)  # plasma chamber void
        c10.importance = {'neutron':1}

        c11 = openmc.Cell(cell_id=11, region= +s111 & -s112 & +s121 & -s122 & +s133 & -s134, fill=self.tungsten) # first wall 
        c12 = openmc.Cell(cell_id=12, region= +s111 & -s112 & +s121 & -s122 & +s132 & -s133, fill=self.eurofer)  # after first wall
        c14 = openmc.Cell(cell_id=14, region= +s111 & -s112 & +s121 & -s122 & +s130 & -s131, fill=self.eurofer)  # manifold

        c21 = openmc.Cell(cell_id=21, region= +s113 & -s114 & +s123 & -s124 & +s135 & -s136, fill=self.tungsten)
        c22 = openmc.Cell(cell_id=22, region= +s113 & -s114 & +s123 & -s124 & +s136 & -s137, fill=self.eurofer)
        c24 = openmc.Cell(cell_id=24, region= +s113 & -s114 & +s123 & -s124 & +s138 & -s139, fill=self.eurofer)

            c13 = openmc.Cell(cell_id=13, region= +s111 & -s112 & +s121 & -s122 & +s131 & -s132)
            c23 = openmc.Cell(cell_id=23, region= +s113 & -s114 & +s123 & -s124 & +s137 & -s138)

            c31 = openmc.Cell(cell_id=31, fill=self.kernel,  region= -s150)
            c32 = openmc.Cell(cell_id=32, fill=self.sic,     region= +s150 & -s151) # 
            # c33 = openmc.Cell(cell_id=33, fill=self.blanket, region= +s151)
            u50 = openmc.Universe(name='BISO Universe', cells=[c31, c32]) # c33

            # Lower left (ll), upper right (ur) of inboard, outboard bounding boxes
            # Example values:
            #   Lower left : [-12.01134016 -12.01134016   0.01    ]
            #   Upper right: [ 12.01134016  12.01134016  45.45    ]
            
            ll_in, ur_in   = c13.region.bounding_box
            ll_out, ur_out = c23.region.bounding_box

            # ll_in, ur_in   = ll_in*np.array([1.01, 1.01, 0.99]), ur_in*np.array([1.01, 1.01, 1.01])
            # ll_out, ur_out = ll_out*np.array([1.01, 1.01, 0.99]), ur_out*np.array([1.01, 1.01, 1.01])

            shape_in  = (2, 2, 128)  # (Nx, Ny, Nz)
            shape_out = (2, 2, 375)

            pitch_in  = (ur_in - ll_in) / shape_in
            pitch_out = (ur_out - ll_out) / shape_out
            
                # Rectangular lattice
                centers_in  = lattice_coords(ll_in, shape_in, pitch_in)
                centers_out = lattice_coords(ll_out, shape_out, pitch_out)
                
                # openmc.model.TRISO(outer_radius, fill, center=(0.0, 0.0, 0.0)) -- each model.TRISO instance IS a cell
                biso_in  = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_in]
                biso_out = [openmc.model.TRISO(0.05, u50, center=c) for c in centers_out]

                c13.fill = openmc.model.create_triso_lattice(biso_in, ll_in, pitch_in, shape_in, self.blanket)
                c23.fill = openmc.model.create_triso_lattice(biso_out, ll_out, pitch_out, shape_out, self.blanket)



                if debug:
                    print(f"DCLL: Fixed lattice BISO placement")
                    print(f"  Inboard lattice:")
                    print(f"    Shape: {shape_in}")
                    print(f"    Pitch: {pitch_in} cm")
                    print(f"    Lower left: {ll_in}")
                    print(f"    Upper right: {ur_in}")
                    print(f"    Subtotal BISO particles: {np.prod(shape_in)}")
                    print(f"  Outboard lattice:")
                    print(f"    Shape: {shape_out}")
                    print(f"    Pitch: {pitch_out} cm")
                    print(f"    Lower left: {ll_out}")
                    print(f"    Upper right: {ur_out}")
                    print(f"    Subtotal BISO particles: {np.prod(shape_out)}")
                    print(f"  Total BISO particles: {np.prod(shape_in) + np.prod(shape_out)}")


        self.cells = [c10, c11, c12, c13, c14, c21, c22, c23, c24]
        self.geometry = openmc.Geometry(openmc.Universe(cells=self.cells))


    def tallies(self, debug=False):

        self.tallies = openmc.Tallies() # initialize
        cell_filter = openmc.CellFilter([cell for cell in self.cells if cell.fill is not None])
        energy_filter = openmc.EnergyFilter(logspace_per_decade(1e-5, 20e6, 100)) # './helpers/utilities.py'

        # Flux tally 
        flux_tally        = self.make_tally('flux', ['flux'], filters=[cell_filter])
        flux_energy_tally = self.make_tally('flux spectrum', ['flux'], filters=[cell_filter, energy_filter])

        # Fertile element reaction rates
        fertile_tally_tot    = self.make_tally(f'Total fertile rxn rate',     ['(n,gamma)', 'fission', 'elastic'], nuclides=['U238', 'U235', 'Th232'])
        fertile_tally        = self.make_tally(f'Fertile rxn rates',          ['(n,gamma)', 'fission', 'elastic'], filters=[cell_filter], nuclides=['U238', 'U235', 'Th232'])
        fertile_energy_tally = self.make_tally(f'Fertile rxn rates spectrum', ['(n,gamma)', 'fission', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['U238', 'U235', 'Th232'])

        # Lithium reaction rates
        Li_tally_tot    = self.make_tally('Total Li rxn rate',     ['(n,gamma)', '(n,Xt)', 'elastic'], nuclides=['Li6', 'Li7'])
        Li_tally        = self.make_tally('Li rxn rates',          ['(n,gamma)', '(n,Xt)', 'elastic'], filters=[cell_filter], nuclides=['Li6', 'Li7'])
        Li_energy_tally = self.make_tally('Li rxn rates spectrum', ['(n,gamma)', '(n,Xt)', 'elastic'], filters=[cell_filter, energy_filter], nuclides=['Li6', 'Li7'])

        # Put in order you want it to print in... recommend fertile_tally's first
        self.tallies.extend([fertile_tally_tot, Li_tally_tot, fertile_tally, Li_tally, flux_tally, fertile_energy_tally, flux_energy_tally, Li_energy_tally]) 


    def make_tally(self, name, scores, filters:list=None, nuclides:list=None):
        tally = openmc.Tally(name=name)
        tally.scores = scores
        if filters is not None:
            tally.filters = filters
        if nuclides is not None:
            tally.nuclides = nuclides
        return tally


    def settings(self, debug=False):

        # 14-MeV point source at center of plasma chamber
        source = openmc.IndependentSource()
        # source.space = openmc.stats.Point((0,0,247.7))                                                           
        source.space = openmc.stats.CartesianIndependent(x=openmc.stats.Uniform(a= -self.geom['c1'], b= +self.geom['c1']),
                                                         y=openmc.stats.Uniform(a= -self.geom['c1'], b= +self.geom['c1']),
                                                         z=openmc.stats.Discrete([247.7], [1.0]) )
        source.particle = 'neutron'
        source.energy   = openmc.stats.Discrete([14.0e6], [1.0])  # 14 MeV
        
        # Define the angular distribution -- mu=1 is straight up (+z), mu=-1 is straight down (-z)
        # We give each a 50% (0.5) probability.
        source.angle = openmc.stats.PolarAzimuthal(
                           mu=openmc.stats.Discrete([1.0, -1.0], [0.5, 0.5]),
                           phi=openmc.stats.Uniform(0, 2*np.pi) # azimuthal angle doesn't matter for mu=+/-1
                           )
        # source.angle = openmc.stats.Isotropic()

        # Create settings and assign source
        self.settings = openmc.Settings()
        self.settings.source = source

        # Run type
        self.settings.run_mode  = 'fixed source'
        self.settings.particles = self.nps
        self.settings.batches   = self.batches

        # note to self 2026-01-22: change to 20 or 50 batches

        # self.settings.trace = (1,1,1)
        self.settings.max_tracks = 4


    def prism_helpers(self, debug=False):

        fertile_kgm3 = self.fertile_kgm3

        if fertile_kgm3 < 10.0:
            self.nps = int(1e7)
            self.batches = int(10)
        else:
            self.nps = int(1e5)
            self.batches = int(10)

        target_total_biso = 2012 
        N1, N2 = 512, 1500

        # HCPB tokamak dimensions (along z in model)
        h1, h2 = 45.0, 82.0     # inboard, outboard thicknesses
        d1 = 2.5 + 0.2          # Eurofer + first wall 
        d0 = 400.0              # plasma chamber
        oo = h1 + d0 + 2 * d1   # outboard offset (where outboard starts in z)
        f1, f2 = 0.39, 0.61     # fraction of fluence into inboard vs. outboard in full OpenMC HCPB model


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
        below in terms of volume fractions.  --ppark
        '''

        # Nominal fraction the breeding volume is in the breeding region
        VF_BV_BR_NOM = VF_LI_NOM + VF_BE_NOM

        # Nominal volume fractions of Li4SiO4 vs. Be in the breeding volume (bv)
        vf_li_bv_nom = VF_LI_NOM / (VF_LI_NOM + VF_BE_NOM)
        vf_be_bv_nom = VF_BE_NOM / (VF_LI_NOM + VF_BE_NOM)

        # Number of BISO spheres per cc of breeding volume
        biso_per_cc_bv = fertile_kgm3_to_biso_per_cc(fertile_kgm3)
        vf_biso_bv     = biso_per_cc_bv * BISO_VOLUME

        if vf_biso_bv > 1.0:
            print(f"Fatal. Your fertile kg/m³ exceeds what can physically fit in the breeding volume!")
            print(f"Fatal. That is, your volume of BISO per cm³ of breeding volume exceeds 1.")

        # New volume ratios of everything w.r.t. breeders in the breeder volume
        vf_libe_bv   = 1 - vf_biso_bv
        vf_li_bv     = vf_li_bv_nom * vf_libe_bv
        vf_be_bv     = vf_be_bv_nom * vf_libe_bv
        vf_checksum1 = vf_biso_bv + vf_li_bv + vf_be_bv  # should equal 1

        # New volume ratios of everyting w.r.t. everything else in the breeding region
        vf_biso_br = vf_biso_bv * (VF_LI_NOM + VF_BE_NOM)
        vf_li_br   = vf_li_bv * (VF_LI_NOM + VF_BE_NOM)
        vf_be_br   = vf_be_bv * (VF_LI_NOM + VF_BE_NOM)
        # vol frac of Eurofer and He-4 doesn't change
        vf_checksum2 = vf_biso_br + vf_li_br + vf_be_br + VF_EU_NOM + VF_HE_NOM  # should equal 1

        # Number of BISO spheres per cc of breeding region
        biso_per_cc_br = vf_biso_br / BISO_VOLUME

        '''
        NB. Maximum possible value of vf_biso_br is (VF_LI_NOM + VF_BE_NOM) = 50.94%. 
        This means that at most 50.94% of the breeding region can be occupied by BISO
        Per 1 cm³ of breeding region, at most 0.5094 cm³ can be BISO, or 973 spheres/cm³.
        Practically, 400 spheres/cm³ = 1000 kg/m³ so the max possible fertile kg/m³ is 2432.50
        '''

        A_ref = target_total_biso / biso_per_cc_br / (h1+h2) 
        A2 = A_ref * (h1+h2) / ( f1/f2 * h1 + h2)                    # outboard XY face area 
        A1 = A2 * f1/f2                                              # inboard XY face area
        b1, b2 = np.round(np.sqrt(A1), 6), np.round(np.sqrt(A2), 6)  # face lengths
        c1, c2 = b1/2, b2/2                                # face HALF-lengths

        # Volumes of the breeding regions (i.e. area*height)
        V1, V2 = A1*h1, A2*h2


        if debug:
            print(40*f"=")
            print(f"Case: {fertile_kgm3} kg/cm³")
            print(f"")
            print(f"Dimensions (inboard, outboard):")
            print(f" XY face half-lengths c1, c2: {c1:.6f}, {c2:.6f} [cm]")
            print(f" XY face lengths      b1, b2: {b1:.6f}, {b2:.6f} [cm]")
            print(f" XY face areas        A1, A2: {A1:.6f}, {A2:.6f} [cm²]")
            print(f" Breeding region vols V1, V2: {V1:.6f}, {V2:.6f} [cm³]")
            print(f" Fraction of fluence  f1, f2: {f1:.6f}, {f2:.6f} ")
            print(f"")
            print(f"'breeder' = Li4SiO4, Be ")
            print(f"'breeding region' (br) = the physical blanket layer that contains the breeder ")
            print(f"                         this layer usually contains breeder + structure + coolant ")
            print(f"'breeding volume' (bv) = volume nominally occupied in breeding region by breeder ")
            print(f"")
            print(f"With respect to the BREEDER material, i.e., per 1 m³ of breeder volume, we have these new volume fractions:")
            print(f"  vf_biso_breeder_new =  {(vf_biso_bv*100):.6f} vol%")
            print(f"  vf_li_breeder_new   =  {(vf_li_bv*100):.6f} vol%")
            print(f"  vf_be_breeder_new   =  {(vf_be_bv*100):.6f} vol%")
            print(f"  check they add up   = {(vf_checksum1*100):.6f} vol%")
            print(f"")
            print(f"With respect to the whole BREEDING REGION, we have these new volume fractions:")
            print(f"  vf_biso_br =  {(vf_biso_br*100):.6f} vol%")
            print(f"  vf_li_br   =  {(vf_li_br*100):.6f} vol%")
            print(f"  vf_be_br   =  {(vf_be_br*100):.6f} vol%")
            print(f"  VF_EU_NOM  =  {(VF_EU_NOM*100):.6f} vol%")
            print(f"  VF_HE_NOM  =  {(VF_HE_NOM*100):.6f} vol%")
            print(f"  check they add up   = {(vf_checksum2*100):.6f} vol%")
            print(f"  check that BISO + Li4SiO4 + Be adds up to the nominal Li4SiO4 + Be fraction")
            print(f"  vf_biso_br + vf_li_br + vf_be_br = {((vf_biso_br+vf_li_br+vf_be_br))*100:.6f} : VF_LI_NOM + VF_BE_NOM = {((VF_LI_NOM+VF_BE_NOM))*100:.6f}")
            print(f"")
            print(f"Check BISO volume fraction wrt breeding region is correct:")
            print(f"  Inboard : N1 * BISO_VOLUME / (b1**2 * h1) = {(N1 * BISO_VOLUME / (b1**2 * h1) * 100):.6f} vol%")
            print(f"  Outboard: N2 * BISO_VOLUME / (b2**2 * h2) = {(N2 * BISO_VOLUME / (b2**2 * h2) * 100):.6f} vol%")
            print(f"  Total   : (N1+N2) * BISO_VOLUME / (b1**2 * h1 + b2**2 * h2) = {((N1+N2) * BISO_VOLUME / (b1**2 * h1 + b2**2 * h2) * 100):.6f} vol%")
            print(f"  ...should match                                  vf_biso_br = {(vf_biso_br*100):.6f} vol%")
            print(f"")
            print(f"BISO spheres per cc of breeding volume: {biso_per_cc_bv:.6f} spheres/cm³")
            print(f"             per cc of breeding region: {biso_per_cc_br:.6f} spheres/cm³")
            print(f"")


        # Calculate surfaces for trapezoidal plasma chamber
        # NB. MCNP plane: 1.0*x + 0.0*y + plane_C*z - plane_D = 0
        z_in_end    = 47.7
        z_out_start = 447.7
        dz = z_out_start - z_in_end
        dx = c2 - c1
        slope = dx / dz
        
        plane_C = -slope
        plane_D = c1 - (slope * z_in_end)

        # I moved the needed variables to these class dictionaries bc the amount of 'self.' 
        # in the arithmetic above was really screwing readability --ppark 2026-01-02

        self.biso_per_cc_bv = biso_per_cc_bv
        self.biso_per_cc_br = biso_per_cc_br
        self.vf_biso_bv = vf_biso_bv
        self.vf_biso_br = vf_biso_br
        self.vf_li_bv = vf_li_bv
        self.vf_be_bv = vf_be_bv

        self.geom = {'N1':N1, 'N2':N2, 'b1':b1, 'b2':b2, 'c1':c1, 'c2':c2, 'h1':h1, 'h2':h2, 'pC':plane_C, 'pD':plane_D}


    def openmc(self, debug=False, write=True, run=False):

        if write:

            self.prism_helpers()
            self.materials()
            self.geometry(debug=debug)
            self.settings()
            self.tallies()

            self.model = openmc.model.Model(self.geometry, self.materials, self.settings, self.tallies)
            self.model.export_to_model_xml(self.path)

        if run:
            if has_statepoint(self.path):
                print(f"Warning. File {self.path}/statepoint.h5 already exists, so this OpenMC run will be aborted...")
            else:
                self.model.run(cwd=self.path, tracks=True)


if __name__ == '__main__':

    os.makedirs(f"./OpenMC/", exist_ok=True)

    for case in ['C',]:
        for fertile_kgm3 in [0.10, 0.50, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99] : #  [0.10, 0.50, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]   # [0.10,60,150]
            current_run = Prism(case, fertile_kgm3, isotope='Th232')
            current_run.openmc(debug=True, write=True, run=True)