import math, csv, random

"""
This script generates origins of 3 different kinds of spheres
all packed into the same volume for a given overall packing factor

--ppark 2025-11-15

"""

# Volume per sphere
def sph_vol(r): 
    return 4/3*math.pi*(r**3)


def main():

    """ Problem parameters """

    # Sphere radii and volumes - BE SURE TO ADD +0.00001 BUFFER TO RADII !!!
    r_A = 0.015 + 0.00001  # [cm] A = Li4SiO4 / OD = 0.4 mm
    r_B = 0.03 + 0.00001  # [cm] B = Be       / OD = 1.0 mm
    r_C = 0.05 + 0.00001  # [cm] C = BISO     / OD = 1.0 mm

    V_A0 = sph_vol(r_A)
    V_B0 = sph_vol(r_B)
    V_C0 = sph_vol(r_C)

    # Target volume fractions 
    frac_A_in_AB  = 0.256
    frac_B_in_AB  = 0.744 
    frac_C_in_ABC = 0.5   # 0.010  # BISO is 1 vol% of total

    # Scaling factor - for testing purposes 
    CUBE_SCALE = 1.00


    """ Derived parameters """

    # Target volumes
    target_vol_ABC = V_C0 / frac_C_in_ABC
    target_vol_AB  = target_vol_ABC - V_C0
    target_vol_A   = target_vol_AB * frac_A_in_AB
    target_vol_B   = target_vol_AB * frac_B_in_AB

    # Number of A and B spheres
    n_A = int(round(target_vol_A / V_A0))
    n_B = int(round(target_vol_B / V_B0))

    print(f"Target volume *of material*: {target_vol_ABC:.6f} [cm³]")
    print(f"Number of pebbles: Li4SiO4: {n_A}, Be: {n_B}")


    # Generate coordinates of B spheres
    B_coords = generate_bcc_coords(n_B, r_B, r_C)
    C_coords = (0.0,0.0,0.0)
    print(f"Number of generated B spheres: {n_B}")

    # Calculate unit cell cube size
    max_coord = 0.0
    for x, y, z in B_coords:
        max_coord = max(max_coord, abs(x), abs(y), abs(z))

    cube_half  = max_coord + r_B
    cube_side  = 2*cube_half
    cube_vol   = cube_side**3
    remain_vol = cube_vol - V_C0

    print(f"Cube half: {cube_half:.6f} [cm]")
    print(f"Cube side: {cube_side:.6f} [cm]")

    pf_B  = V_B0*n_B / remain_vol
    pf_AB = (V_A0*n_A + V_B0*n_B) / remain_vol

    vol_A = n_A*V_A0
    vol_B = n_B*V_B0
    vf_A_AB = vol_A / (vol_A + vol_B)
    vf_B_AB = vol_B / (vol_A + vol_B)

    print(f"Be packing fraction in breeder vol:    {pf_B:.6f}")
    print(f"Li4SiO4 + Be pack frac in breeder vol: {pf_AB:.6f}")

    print(f"Actual volume of Li4SiO4: {vol_A:.6f} [cm³]")
    print(f"Actual volume of Be     : {vol_B:.6f} [cm³]")
    print(f"Actual vol frac Li4SiO4 : {vf_A_AB:.6f} ")
    print(f"Actual vol frac Be      : {vf_B_AB:.6f} ")

    ''' TURN ON WHEN YOU NEED TO ACTUALLY GENERATE NEW SET OF POINTS
    """ Populate empty spaces between Be (B) and BISO (C) with Li4SiO4 pebbles (A) """
    GRID_STEP = 0.0005  # [cm]

    # Build initial grid for A centers
    print("\nBuilding grid of candidate A centers...")
    grid_points = grid_build(cube_half, r_A, GRID_STEP)
    print(f"Initial grid points: {len(grid_points):,}")

    # Existing spheres at this stage: central C + all B
    existing_spheres = [(C_coords[0], C_coords[1], C_coords[2], r_C)] + \
                       [(x, y, z, r_B) for (x, y, z) in B_coords]

    # Prune grid for C and B (remove points where an A would overlap them)
    print("Pruning grid points that overlap C or any B sphere...")
    grid_points = list(grid_pruned_iter(grid_points, existing_spheres, r_A))
    print(f"Grid points after pruning for C and B: {len(grid_points):,}")

    # Place A spheres by scanning the grid
    print(f"Placing {n_A} A spheres by scanning grid...")
    A_coords, grid_points = place_spheres(n_A, r_A, existing_spheres, grid_points)
    print("Finished placing A spheres.")



    """ Output to CSV """

    csv_name = "sphere_packing.csv"
    print(f"Writing coordinates to: {csv_name}")

    all_spheres = []
    all_spheres.append(("BISO", 0.0, 0.0, 0.0, r_C))
    for (x, y, z) in B_coords: all_spheres.append(("Be",      x, y, z, r_B))
    for (x, y, z) in A_coords: all_spheres.append(("Li4SiO4", x, y, z, r_A))

    with open("sphere_packing.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["type", "x_cm", "y_cm", "z_cm", "radius_cm"])
        for row in all_spheres:
            writer.writerow([row[0],
                             f"{row[1]:.8f}",
                             f"{row[2]:.8f}",
                             f"{row[3]:.8f}",
                             f"{row[4]:.8f}"])

    print("\nCentral C sphere center (cm):", C_coords)
    print("Number of A spheres placed:", len(A_coords))
    print("Number of B spheres placed:", len(B_coords))
    '''




def generate_bcc_coords(n_B, r_B, r_C):
    """ 
    Generate at least n_B BCC lattice positions for Be (B) pebbles,
    avoid overlap with BISO (C) sphere at origin.

    Args:
        n_B [int]
            Number of B spheres to draw
        r_B, r_C
            Radii of B and C spheres in cm

    Return :
        list of (x,y,z) in cm.
    """

    V_B0 = sph_vol(r_B)
    a = 4/math.sqrt(3) * r_B

    n = 2
    while True:
        coords = []

        # Corners
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    coords.append((i,j,k))

        # Body-centered sites
        for i in range(n-1):
            for j in range(n-1):
                for k in range(n-1):
                    coords.append((i+0.5, j+0.5, k+0.5))

        # Center lattice at origin
        shift = (n-1)/2
        coords_centered = [ ((x-shift)*a, (y-shift)*a, (z-shift)*a) for (x,y,z) in coords ]

        # Remove coords that overlap with C 
        min_dr = r_B + r_C
        coords_pruned = []
        for (x,y,z) in coords_centered:
            if x**2 + y**2 + z**2 >= min_dr**2: # 
                coords_pruned.append((x,y,z))

        if len(coords_pruned) >= n_B:
            return coords_pruned[:n_B]

        n += 1


def compute_pf_from_coords(coords, r):
    """
    Compute packing fraction of spheres of radius r
    occupying the minimal cube that encloses them.
    """
    V_sphere = sph_vol(r)
    n = len(coords)

    max_coord = 0.0
    for x, y, z in coords:
        max_coord = max(max_coord, abs(x), abs(y), abs(z))

    cube_half = max_coord + r
    cube_side = 2 * cube_half
    cube_vol  = cube_side**3

    return n * V_sphere / cube_vol


def grid_build(cube_half, r, dx):
    """
    Build a grid 
    """
    min_coord = -cube_half + r
    max_coord =  cube_half - r

    n_steps = int(math.floor((max_coord - min_coord)/dx)) + 1

    points = []
    count, next_print = 0, 0.05

    for i in range(n_steps):
        x = min_coord + i * dx
        for j in range(n_steps):
            y = min_coord + j * dx
            for k in range(n_steps):
                z = min_coord + k * dx
                points.append((x,y,z))
                count += 1
                progress = count / (n_steps**3)
                if progress >= next_print:
                    print(f"  Progress: {progress*100:.1f}% ({count:,}/{(n_steps**3):,})")
                    next_print += 0.05

    return points


def grid_pruned_iter(grid_points, sphere_coords, r, print_every=0.05):
    """
    Yields grid points that do NOT overlap any sphere in sphere_coords.
    
    Uses streaming (no second full list) and prints progress updates.
    print_every = fraction of total scanned at which to print a status line.
    """
    # Precompute squared exclusion radii
    min_d2_list = [(r + cr)**2 for (_,_,_,cr) in sphere_coords]

    total = len(grid_points)
    next_print = print_every
    count = 0

    for (x,y,z) in grid_points:
        count += 1
        ok = True

        # overlap check
        for (cx,cy,cz,cr), d2 in zip(sphere_coords, min_d2_list):
            dx = x - cx
            dy = y - cy
            dz = z - cz
            if dx*dx + dy*dy + dz*dz < d2:
                ok = False
                break

        if ok:
            yield (x,y,z)

        # progress print
        progress = count / total
        if progress >= next_print:
            print(f"  Grid prune: {progress*100:.1f}% ({count:,}/{total:,})")
            next_print += print_every



def place_spheres(n_spheres, r, sphere_coords, grid_points):
    """
    Greedy placement of n_spheres with radius r onto grid_points,
    avoiding overlap with sphere_coords and previously placed spheres.
    """
    points = []
    grid_candidates = list(grid_points)  # make a copy so we can modify
    random.shuffle(grid_candidates)      # randomize order of candidates to avoid spheres being unbalanced to one side
    existing = list(sphere_coords)       # start with C + B spheres

    for i in range(n_spheres):
        placed = False

        # Loop over candidate grid points
        for idx, (x, y, z) in enumerate(grid_candidates):
            overlap = False

            # Check against all existing spheres (C, B, and placed A)
            for (cx, cy, cz, cr) in existing:
                dx = x - cx
                dy = y - cy
                dz = z - cz
                if dx*dx + dy*dy + dz*dz < (r + cr)**2:
                    overlap = True
                    break

            if overlap:
                continue  # try next candidate

            # If we reach here, this candidate works
            points.append((x, y, z))
            existing.append((x, y, z, r))   # include this A sphere in collision checks
            placed = True

            # Keep ONLY points far enough from this new A sphere
            new_grid_points = []
            for (gx, gy, gz) in grid_candidates:
                dx = gx - x
                dy = gy - y
                dz = gz - z
                if dx*dx + dy*dy + dz*dz >= (2*r)**2:
                    new_grid_points.append((gx, gy, gz))
            grid_candidates = new_grid_points

            print(f"Placed A-sphere {i+1} / {n_spheres} "
                  f"(remaining candidates: {len(grid_candidates):,})")

            # Done placing this A-sphere, move on to the next one
            break

        if not placed:
            print(f"Could not place A-sphere {i+1}/{n_spheres} "
                  f"(no valid grid point left)")
            break

    return points, grid_candidates




if __name__ == '__main__':
    main()