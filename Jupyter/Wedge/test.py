import numpy as np

# GLOBAL PARAMETERS
TEMP_K = 900
BISO_VOLUME = 4/3*np.pi*(0.05)**3
KERNEL_VOLUME = 4/3*np.pi*(0.04)**3

DENSITY_UO2  = 10.50 # [g/cm³]
DENSITY_ThO2 = 10.00 # [g/cm³]

AMU_U235 = 235.0439299
AMU_U238 = 238.05078826
AMU_U = 238.02891 # for natural enrichment
AMU_O = 15.999
AMU_UO2 = AMU_U + 2 * AMU_O
AMU_Th232 = 232.0381

# Volume fractions in HCPB blanket
VF_LI_NOM = 0.1304
VF_BE_NOM = 0.3790
VF_EU_NOM = 0.1176
VF_HE_NOM = 0.3730


def fertile_kgm3_to_biso_per_cc(fertile_kgm3, isotope='U238'):
    if isotope == 'U238':
        biso_per_cc = fertile_kgm3 * AMU_UO2 / AMU_U238 / KERNEL_VOLUME / DENSITY_UO2 / 100**3 * 1000
    return biso_per_cc


for fertile_kgm3 in [0.10, 0.50, 1.5, 15, 30, 60, 90, 120, 150, 250, 500, 750, 999.99]:

    VF_LI_NOM = 0.1304
    VF_BE_NOM = 0.3790
    VF_EU_NOM = 0.1176
    VF_HE_NOM = 0.3730


    target_total_biso = 2012
    N1, N2 =  512, 1500

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
    below in terms of volume fractions.
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
    This means that at most 50.94% of the breeding region can be occupied by BISO.
    Per 1 cm³ of breeding region, at most 0.5094 cm³ can be BISO, or 973 spheres/cm³ of breeding region.

    Of course, per 1 cm³ of breeder volume, BISO can be up to 1 cm³. 
    This equals 1910 spheres/cm³ of breeding volume. (And 1910*0.5094=973)

    Practically, 403 spheres/cm³ = 1000 kg/m³ = 21 vol% of breeding volume = 11 vol% of breeding region
    so the max possible fertile density is 4739 kg/m³ of breeding volume.
    '''

    A_ref = target_total_biso / biso_per_cc_br / (h1+h2) 
    A2 = A_ref * (h1+h2) / ( f1/f2 * h1 + h2)                    # outboard XY face area 
    A1 = A2 * f1/f2                                              # inboard XY face area
    b1, b2 = np.round(np.sqrt(A1), 6), np.round(np.sqrt(A2), 6)  # face lengths
    c1, c2 = b1/2, b2/2                                # face HALF-lengths

    # Volumes of the breeding regions (i.e. area*height)
    V1, V2 = A1*h1, A2*h2


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
    print(f"BISO spheres per cc of breeding volume: {biso_per_cc_bv:.6f} spheres/cm³ : per cc of breeding region : {biso_per_cc_br:.6f} spheres/cm³")
    print(f"")

