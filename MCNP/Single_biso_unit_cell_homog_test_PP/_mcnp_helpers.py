import numpy as np

"""
Nuclear constants -- from atom.kaeri.re.kr/nuchart/
"""
AVO = 6.022141076e+23
AMU_LI6, AMU_LI7 = 6.0150, 7.0160 # 6.01512288742, 7.01600343426 # amu = g/mol
AMU_F19 = 18.9984 # 18.99840316207
AMU_Be9 =  9.0120 # 9.012183062
AMU_U = 238.02891 # for natural enrichment
AMU_U234 = 234.0409456
AMU_U235 = 235.0439299
AMU_U238 = 238.05078826
AMU_He4 = 4.0026 # [amu]
AMU_O   = 15.999
AMU_Th = 232.0381
AMU_UF4 = AMU_U + 4 * AMU_F19
AMU_ThF4 = AMU_Th + 4 * AMU_F19
AMU_UO2 = AMU_U + 2 * AMU_O
AMU_ThO2 = AMU_Th + 2 * AMU_O
AMU_FLIBE = 98.89 # g/mol
DENSITY_UF4 = 6.7 # g/cm³ - no this ³ isn't chatgpt, yes i was petty enough to add it
DENSITY_ThF4 = 6.3 # g/cm³
DENSITY_UO2 = 10.5 # g/cm³
DENSITY_ThO2 = 10.0 # g/cm³
DENSITY_SIC = 3.2 # g/cm³
DENSITY_LI4SIO4 = 2.17 # g/cm3 at 900 K
DENSITY_He  = 0.00428 # g/cm3 at 900 K, 8 MPa
DENSITY_Be  = 1.80 # g/cm3 at 900 K
AMU_PU239 = 239.0521634
AMU_U233 = 233.039635207
SEC_PER_YR = 3600 * 24 * 365


BISO_KERNEL_RADIUS   = 0.04 # cm  # r = 400 μm = 0.0400 cm // "800 μm kernel"
BISO_RADIUS          = 0.05 # cm  # r = 500 μm = 0.0500 cm // "100 μm thickness"
BISO_VOLUME          = (4 / 3) * np.pi * (BISO_RADIUS)**3            # volume of single BISO particle
KERNEL_VOLUME        = (4 / 3) * np.pi * (BISO_KERNEL_RADIUS)**3     # volume of UO2/ThO2 kernel in single BISO particle
COAT_VOLUME          = BISO_VOLUME - KERNEL_VOLUME
BISO_KERNEL_VOL_FRAC = KERNEL_VOLUME / BISO_VOLUME  # vol frac kernel in single BISO
BISO_COAT_VOL_FRAC   = 1.0 - BISO_KERNEL_VOL_FRAC

def cube_dims(biso_vol_frac):
    cube_length = ((4/3 * np.pi * BISO_RADIUS**3)/biso_vol_frac)**(1/3)
    cube_half   = cube_length / 2
    return cube_length, cube_half

# Carbon isotope abundances (natural)
C12_ABUNDANCE = 0.9893
C13_ABUNDANCE = 0.0107
AMU_C12 = 12.0000000
AMU_C13 = 13.0033548

# Oxygen isotope abundances (natural)
O16_ABUNDANCE = 0.99757
O17_ABUNDANCE = 0.00038
O18_ABUNDANCE = 0.00205
AMU_O16 = 15.9949146
AMU_O17 = 16.9991317
AMU_O18 = 17.9991610

U234_ABUNDANCE = 0.000054
U235_ABUNDANCE = 0.007204
U238_ABUNDANCE = 0.992742

# Silicon isotope abundances (natural)
SI28_ABUNDANCE = 0.92223
SI29_ABUNDANCE = 0.04685
SI30_ABUNDANCE = 0.03092
AMU_SI28 = 27.9769265
AMU_SI29 = 28.9764947
AMU_SI30 = 29.9737702

# Additional constants for SiC
AMU_C  = 12.0107  # Natural carbon (mostly C-12)
AMU_SI = 28.0855  # Natural silicon
AMU_SIC = AMU_SI + AMU_C12  # SiC molecular weight

# Lead isotope abundances
PB204_ABUNDANCE = 0.0140
PB206_ABUNDANCE = 0.2410
PB207_ABUNDANCE = 0.2210
PB208_ABUNDANCE = 0.5240
AMU_PB204 = 203.9730435
AMU_PB206 = 205.9744652
AMU_PB207 = 206.9758968
AMU_PB208 = 207.9766520
AMU_PB    = AMU_PB204*PB204_ABUNDANCE + AMU_PB206*PB206_ABUNDANCE + AMU_PB207*PB207_ABUNDANCE + AMU_PB208*PB208_ABUNDANCE


def calculate_atomic_fractions(breeder, biso_vol_frac, li6_enrichment_at_percent):
    """
    Calculate atomic fractions of all isotopes in BISO-FLiBe system
    
    Parameters:
    -----------
    biso_vol_frac : float
        Volume fraction of BISO particles in the system
    li6_enrichment_at_percent : float
        Li-6 enrichment in atomic percent (0-100)
    
    Returns:
    --------
    dict : Dictionary of ZAID keys and atomic fraction values
    """
    
    print("Breeder-BISO Atomic Fraction Calculator")
    print("=" * 40)
    print("BISO composition:")
    print(f"  - UO2 kernel: radius = {BISO_KERNEL_RADIUS} cm, density = {DENSITY_UO2} g/cm³")
    print(f"  - U enrichment: 0.71 wt% U-235 (natural)")
    print(f"  - SiC coating: thickness = {(BISO_RADIUS - BISO_KERNEL_RADIUS):.2f} cm, density = {DENSITY_SIC} g/cm³")

    # Convert Li-6 enrichment from at% to fraction
    li6_frac = li6_enrichment_at_percent / 100.0
    li7_frac = 1.0 - li6_frac
    
    # FLiBe properties
    breeder_densities = {'flibe':1.94, 'pbli':9.45, 'pebble':0} # g/cm3 # pb-li at 900 K = 9448.891 kg/m3
    density_breeder = breeder_densities[breeder]

    # Calculate masses
    mass_uo2     = KERNEL_VOLUME * DENSITY_UO2  # g
    mass_sic     = COAT_VOLUME * DENSITY_SIC    # g
    mass_breeder = BISO_VOLUME / biso_vol_frac * (1 - biso_vol_frac) * density_breeder  # g

    mass_total = mass_uo2 + mass_sic + mass_breeder
    vol_total  = BISO_VOLUME / biso_vol_frac # = KERNEL_VOLUME + COAT_VOLUME + BISO_VOLUME / biso_vol_frac * (1 - biso_vol_frac)
    rho_total  = mass_total / vol_total 

    # Natural uranium enrichment (0.71 wt% U-235)
    u234_wt_frac = 0.000053
    u235_wt_frac = 0.007114
    u238_wt_frac = 0.992833
    
    # Calculate mols of each component
    
    # UO2 kernel
    mol_uo2  = mass_uo2 / AMU_UO2
    mol_o    = 2 * mol_uo2
    
    # SiC coating
    mol_sic = mass_sic / AMU_SIC
    
    # FLiBe (F4Li2Be = 2LiF + BeF2)
    # Calculate average molecular weight considering Li isotopes
    amu_li_avg = li6_frac * AMU_LI6 + li7_frac * AMU_LI7

    # Break down into isotopes
    mols = {'14028':0, '14029':0, '14030':0, # need these isotopes to be initialized here!
             '8016':0, '8017':0, '8018':0,}
    
    if breeder == 'flibe':
        amu_flibe = 4 * AMU_F19 + 2 * amu_li_avg + AMU_BE9
        mol_flibe = mass_breeder / amu_flibe
        mol_f   = 4 * mol_flibe
        mol_li  = 2 * mol_flibe
        mol_be  = mol_flibe
    
        # Fluorine (F-19, 100% natural abundance)
        mols['9019'] = mol_f 
        
        # Lithium isotopes
        mols['3006'] = mol_li * li6_frac
        mols['3007'] = mol_li * li7_frac
        
        # Beryllium (Be-9, 100% natural abundance)
        mols['4009'] = mol_be 

    elif breeder == 'pbli':
        amu_pbli  = 0.83*AMU_PB + 0.17*amu_li_avg
        mol_pbli = mass_breeder / amu_pbli
        mol_83pb = 0.83 * mol_pbli
        mol_17li = 0.17 * mol_pbli
    
        # Lead (F-19, 100% natural abundance)
        mols['82204'] = mol_83pb * PB204_ABUNDANCE
        mols['82206'] = mol_83pb * PB206_ABUNDANCE
        mols['82207'] = mol_83pb * PB207_ABUNDANCE
        mols['82208'] = mol_83pb * PB208_ABUNDANCE
        
        # Lithium isotopes
        mols['3006']  = mol_17li * li6_frac
        mols['3007']  = mol_17li * li7_frac

    elif breeder == 'pebble':
        vol_breeder = vol_total - BISO_VOLUME
        vol_li4sio4 = (0.1304 / (0.1304+0.3790)) * vol_breeder
        vol_be      = (0.3790 / (0.1304+0.3790)) * vol_breeder
        vol_he      = (0.423/(1-0.423)) * vol_breeder # 0.423 --ppark 2025-11-25

        amu_li4sio4 = 4*amu_li_avg + AMU_SI + 4*AMU_O
        mol_li4sio4 = vol_li4sio4 * DENSITY_LI4SIO4 / amu_li4sio4
        mol_be      = vol_be * DENSITY_Be / AMU_Be9
        mol_he      = vol_he * DENSITY_He / AMU_He4

        # Lithium isotopes
        mols['3006'] = 4*mol_li4sio4 * li6_frac
        mols['3007'] = 4*mol_li4sio4 * li7_frac

        # Silicon isotopes / notice '+=' operator!
        mols['14028'] += mol_li4sio4 * SI28_ABUNDANCE 
        mols['14029'] += mol_li4sio4 * SI29_ABUNDANCE 
        mols['14030'] += mol_li4sio4 * SI30_ABUNDANCE 

        # Oxygen isotopes  / notice '+=' operator!
        mols['8016'] += 4*mol_li4sio4 * O16_ABUNDANCE 
        mols['8017'] += 4*mol_li4sio4 * O17_ABUNDANCE 
        mols['8018'] += 4*mol_li4sio4 * O18_ABUNDANCE 

        # Beryllium (Be-9, 100% natural abundance)
        mols['4009'] = mol_be 

        # Helium (He-4, 100% natural abundance)
        mols['2004'] = mol_he

        density_breeder = (vol_li4sio4*DENSITY_LI4SIO4 + vol_be*DENSITY_Be + vol_he*DENSITY_He)/ (vol_breeder+vol_he)
        rho_total = (mass_uo2 + mass_sic + vol_li4sio4*DENSITY_LI4SIO4 + vol_be*DENSITY_Be + vol_he*DENSITY_He)/ (vol_total+vol_he)

        cube_vol = vol_total + vol_he 
        cube_half = cube_vol**(1/3) * 1/2
        print(f"*****Try cube volume of {cube_vol:.4e} [cm]")
        print(f"*****Try cube half-length of {cube_half:.6f} [cm]")
    
    # Carbon isotopes from SiC
    mols['6012'] = mol_sic * C12_ABUNDANCE 
    mols['6013'] = mol_sic * C13_ABUNDANCE 

    # Silicon isotopes from SiC / notice '+=' operator!
    mols['14028'] += mol_sic * SI28_ABUNDANCE 
    mols['14029'] += mol_sic * SI29_ABUNDANCE 
    mols['14030'] += mol_sic * SI30_ABUNDANCE 

    # Oxygen isotopes from UO2 / notice '+=' operator!
    mols['8016'] += mol_o * O16_ABUNDANCE 
    mols['8017'] += mol_o * O17_ABUNDANCE 
    mols['8018'] += mol_o * O18_ABUNDANCE 
    
    # Uranium isotopes
    mols['92234'] = mol_uo2 * U234_ABUNDANCE
    mols['92235'] = mol_uo2 * U235_ABUNDANCE
    mols['92238'] = mol_uo2 * U238_ABUNDANCE
    
    # Calculate total mols
    total_mols = sum(mols.values())
    
    # Calculate atomic fractions
    atomic_fractions = {}
    for zaid, num_mols in mols.items():
        atomic_fractions[zaid] = num_mols / total_mols

    # Print results
    print(f"\nInput parameters:")
    print(f"  BISO volume fraction: {biso_vol_frac:.4f}")
    print(f"  Li-6 enrichment: {li6_enrichment_at_percent:.2f} at%")
    print()
    print(f"{breeder.upper()} mass density [g/cm³]: {density_breeder:.6f}")
    print(f"Mixture mass density [g/cm³]: {rho_total:.6f}")
    
    return atomic_fractions

def print_pebble_het():
    """
    Print atomic fractions in MCNP format for Li4SiO4-Be WITHOUT any BISO
    """

    # Volume fractions of Li4SiO4 and Be based on EU DEMO HCPB
    vf_li4sio4 = (0.1304 / (0.1304+0.3790+0.3730)) 
    vf_be      = (0.3790 / (0.1304+0.3790+0.3730)) 
    vf_he      = (0.3730 / (0.1304+0.3790+0.3730))

    # Moles of Li4SiO4 and Be
    li6_frac = 0.60
    li7_frac = 1 - li6_frac
    amu_li_avg = li6_frac * AMU_LI6 + li7_frac * AMU_LI7
    amu_li4sio4 = 4*amu_li_avg + AMU_SI + 4*AMU_O
    mol_li4sio4 = vf_li4sio4 * DENSITY_LI4SIO4 / amu_li4sio4
    mol_be      = vf_be * DENSITY_Be / AMU_Be9
    mol_he      = vf_he * DENSITY_He / AMU_He4

    # Molar fractions of Li4SiO4 and Be
    mf_li4sio4 = mol_li4sio4 / (mol_li4sio4+mol_be+mol_he)
    mf_be      =      mol_be / (mol_li4sio4+mol_be+mol_he)
    mf_he      =      mol_he / (mol_li4sio4+mol_be+mol_he)
    mol_total  = 9*mol_li4sio4 + mol_be + mol_he

    ''' Multiply molar ratios to molar fractions '''
    mols = {}

    # Lithium isotopes
    mols['3006'] = 4*mol_li4sio4/mol_total * li6_frac
    mols['3007'] = 4*mol_li4sio4/mol_total * li7_frac

    # Silicon isotopes / notice '+=' operator!
    mols['14028'] = mol_li4sio4/mol_total * SI28_ABUNDANCE 
    mols['14029'] = mol_li4sio4/mol_total * SI29_ABUNDANCE 
    mols['14030'] = mol_li4sio4/mol_total * SI30_ABUNDANCE 

    # Oxygen isotopes  / notice '+=' operator!
    mols['8016'] = 4*mol_li4sio4/mol_total * O16_ABUNDANCE 
    mols['8017'] = 4*mol_li4sio4/mol_total * O17_ABUNDANCE 
    mols['8018'] = 4*mol_li4sio4/mol_total * O18_ABUNDANCE 

    # Beryllium (Be-9, 100% natural abundance)
    mols['4009'] = mol_be/mol_total

    # Helium (He-4, 100% natural abundance)
    mols['2004'] = mol_he/mol_total

    print(f"c       ----------------------------------------------------------------------------")
    print(f"c       Li4SiO4 ({vf_li4sio4*100:.2f} vol% = {mf_li4sio4*100:.2f} mol%, 60 at% Li-6)")
    print(f"c       Be      ({vf_be*100:.2f} vol% = {mf_be*100:.2f} mol%)                       ")
    print(f"c       He      ({vf_he*100:.2f} vol% = {mf_he*100:.2f} mol%)                       ")
    print(f"c       ----------------------------------------------------------------------------")
    print(f"c                                                                              ")
    print(f" m300    3006.02c  {mols['3006']:.6f}  $  {4*mol_li4sio4/mol_total:.4f} *  60.000 at% in Li")
    print(f"         3007.02c  {mols['3007']:.6f}  $  {4*mol_li4sio4/mol_total:.4f} *  10.000 at% in Li")
    print(f"        14028.02c  {mols['14028']:.6f}  $  {mol_li4sio4/mol_total:.4f} *  92.223 at% in Si")
    print(f"        14029.02c  {mols['14029']:.6f}  $  {mol_li4sio4/mol_total:.4f} *   4.685 at% in Si")
    print(f"        14030.02c  {mols['14030']:.6f}  $  {mol_li4sio4/mol_total:.4f} *   3.092 at% in Si")
    print(f"         8016.02c  {mols['8016']:.6f}  $  {4*mol_li4sio4/mol_total:.4f} *  99.757 at% in O")
    print(f"         8017.02c  {mols['8017']:.6f}  $  {4*mol_li4sio4/mol_total:.4f} *   0.038 at% in O")
    print(f"         8018.02c  {mols['8018']:.6f}  $  {4*mol_li4sio4/mol_total:.4f} *   0.205 at% in O")
    print(f"         4009.02c  {mols['4009']:.6f}  $  {mol_be/mol_total:.4f} * 100.000 at% in Be")
    print(f"         2004.02c  {mols['2004']:.6f}  $  {mol_he/mol_total:.4f} * 100.000 at% in He")

 




def print_atomic_fractions(atomic_fractions):
    """
    Print atomic fractions in MCNP format with isotopic abundances
    """
    # Define the order and format for output with isotopic abundances
    isotope_data = {
         '2004': (' 2004.02c', 'He', 100.000),
         '9019': (' 9019.02c', 'F',  100.000),   # F-19 is 100% of natural F
         '3006': (' 3006.02c', 'Li', None),      # Will be calculated
         '3007': (' 3007.02c', 'Li', None),      # Will be calculated
         '4009': (' 4009.02c', 'Be', 100.000),   # Be-9 is 100% of natural Be
         '6012': (' 6012.02c', 'C',   98.93),    # C-12 abundance
         '6013': (' 6013.02c', 'C',    1.07),    # C-13 abundance
         '8016': (' 8016.02c', 'O',   99.757),   # O-16 abundance
         '8017': (' 8017.02c', 'O',    0.038),   # O-17 abundance
         '8018': (' 8018.02c', 'O',    0.205),   # O-18 abundance
        '14028': ('14028.02c', 'Si',  92.223),   # Si-28 abundance
        '14029': ('14029.02c', 'Si',   4.685),   # Si-29 abundance
        '14030': ('14030.02c', 'Si',   3.092),   # Si-30 abundance
        '82204': ('82204.02c', 'Pb (83Pb-17Li by at%)',   1.40),    # Pb-204 abundance
        '82206': ('82206.02c', 'Pb (83Pb-17Li by at%)',  24.10),    # Pb-206 abundance
        '82207': ('82207.02c', 'Pb (83Pb-17Li by at%)',  22.10),    # Pb-207 abundance
        '82208': ('82208.02c', 'Pb (83Pb-17Li by at%)',  52.40),    # Pb-208 abundance
        '92234': ('92234.02c', 'U',    0.0054),  # U-234 at 0.0053 wt% ≈ 0.0054 at%
        '92235': ('92235.02c', 'U',    0.7204),  # U-235 at 0.7114 wt% ≈ 0.7204 at%
        '92238': ('92238.02c', 'U',   99.2742)   # U-238 abundance
    }
    
    # Calculate Li isotope percentages from the actual fractions
    li6_atoms = atomic_fractions.get('3006', 0)
    li7_atoms = atomic_fractions.get('3007', 0)
    li_total = li6_atoms + li7_atoms
    if li_total > 0:
        li6_percent = (li6_atoms / li_total) * 100
        li7_percent = (li7_atoms / li_total) * 100
        isotope_data['3006'] = (' 3006.02c', 'Li', li6_percent)
        isotope_data['3007'] = (' 3007.02c', 'Li', li7_percent)
    
    print("\nAtomic fractions in MCNP format:")
    print("-" * 70)
    
    # Print in specified order
    output_order = ['2004','82204', '82206', '82207', '82208', 
                    '9019', '3006', '3007','4009',
                    '14028', '14029', '14030',
                    '6012', '6013', 
                    '8016', '8017', '8018', 
                    '92234', '92235', '92238']
    
    for zaid_key in output_order:
        if zaid_key in atomic_fractions and atomic_fractions[zaid_key] > 0:
            zaid_mcnp, element, abundance = isotope_data[zaid_key]
            fraction = atomic_fractions[zaid_key]
            if abundance is not None:
                print(f"{zaid_mcnp:12s}  {fraction:.11f}    $ {abundance:7.3f} at% in {element}")
            else:
                print(f"{zaid_mcnp:12s}  {fraction:.11f}")


def calculate_weight_fraction(atom_fraction_li6):
    """
    Convert Li-6 atomic fraction to weight fraction
    
    Parameters:
    -----------
    atom_fraction_li6 : float
        Atomic fraction of Li-6 in lithium (0 to 1)
    
    Returns:
    --------
    float : Weight fraction of Li-6
    """
    # Atomic fractions
    f6_at = atom_fraction_li6
    f7_at = 1.0 - atom_fraction_li6
    
    # Calculate average atomic weight
    avg_atomic_weight = f6_at * AMU_LI6 + f7_at * AMU_LI7
    
    # Calculate weight fractions
    f6_wt = (f6_at * AMU_LI6) / avg_atomic_weight
    f7_wt = (f7_at * AMU_LI7) / avg_atomic_weight
    
    return f6_wt, f7_wt, avg_atomic_weight




if __name__ == "__main__":

    # Calculate atomic fractions
    biso_vol_frac  = float(0.01) # 0.01 # 0.50
    li6_enrichment = float(60) # 7.5 wt% = 8.640288 at% # 20 wt% = 22.576908 at%
    atomic_fractions = calculate_atomic_fractions('pebble', biso_vol_frac, li6_enrichment) # biso_vol_frac, li6_enrichment) #  
    print_atomic_fractions(atomic_fractions)

    print_pebble_het()

    # Print summary statistics
    # total_frac = sum(atomic_fractions.values())
    # print(f"\nTotal atomic fraction sum: {total_frac:.12f}")
    # print(f"Deviation from 1.0: {abs(1.0 - total_frac):.2e}")


    # Convert Li-6 weight fraction to atomic fraction
    # print("="*60)
    # print("Li-6 atomic fraction to weight fraction conversion")
    # for w6 in [0.075, 0.20]:
    #     a6 = w6/6.015/(w6/6.015+(1-w6)/7.016)
    #     w7, a7 = 1-w6, 1-a6
    #     print(f"Li-6: {w6:.6f} wt% → {a6:.6f} at% ({a6*100:.6f}%)")
    #     print(f"Li-7: {w7:.6f} wt% → {a7:.6f} at% ({a7*100:.6f}%)")
    #     print(f".")

    # print(AMU_UO2)
    bv = 4/3*np.pi*(0.05)**3
    
    # for p in [0.01, 0.50]:
    #     l, _ = cube_dims(p)
    #     breedv = l**3 - bv
    #     print(f"== biso vol frac {p} ==")
    #     print(f"cube length [cm]: {l:.6e}")
    #     print(f"cube half [cm]: {(l/2):.6f}")
    #     print(f"cube volume [cc]: {l**3:.4e}")
    #     print(f"biso volume [cc]: {bv:.4e}")
    #     print(f"UO2 volume  [cc]: {(KERNEL_VOLUME):.4e}")
    #     print(f"SiC volume  [cc]: {(BISO_VOLUME-KERNEL_VOLUME):.4e}")
    #     print(f"breeder vol [cc]: {breedv:.4e}")
    # print(KERNEL_VOLUME)
    # print(BISO_VOLUME-KERNEL_VOLUME)
    # print(l**3-BISO_VOLUME)
