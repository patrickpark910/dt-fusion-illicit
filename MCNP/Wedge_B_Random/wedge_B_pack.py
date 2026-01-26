import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

"""
SPHERE-PACKING IN RECTANGULAR PRISMS
My script to generate (x,y,z) coordinates for 1-mm wide BISO spheres
packed in a trapezoidal prism, which represents a borehole wedge through
the tokamak blankets.

"""

def pack(length, height, n_target, offset=0, r=0.05, name=None, spheres_existing=None, max_attempts=int(1e8)):
    """
    Packs spheres into a rectangular prism using Random Sequential Packing.
    --ppark 2025-12-21
    
    Args:
        length (float): Side length of the square base. Spheres placed in [-length/2, length/2] for x and y).
        height (float): Height of the packing region along the z-axis.
        offset (float): Vertical offset. Spheres are placed with z in [offset, height + offset].
        n_target (int): Target number of spheres to pack.
        r (float): Radius of each sphere. Default is 0.05 cm.
        name (str): Name of particle; should be 'biso', 'be', or 'li4sio4'.
        spheres_existing (pd.DataFrame or None): DataFrame with columns ['name', 'r', 'x', 'y', 'z'] of existing spheres
            to avoid overlap with. These are not included in the returned DataFrame. Default is None.
        max_attempts (int): Maximum placement attempts per sphere before giving up. Default is 10000.
        
    Returns:
        pd.DataFrame: DataFrame with columns ['name', 'r', 'x', 'y', 'z'] of successfully packed spheres.
    """
    spheres = []
    
    margin = r  # Safety margin from walls
    min_dist_sq = (2*r)**2 * 1.02  # Minimum distance squared (with buffer)
    
    # Convert existing spheres to numpy array for efficient distance calculation
    if spheres_existing is not None and len(spheres_existing) > 0:
        existing = spheres_existing[['x', 'y', 'z']].values.astype(float)
    else:
        existing = None
    
    print(f"Attempting to pack {n_target} {name} spheres...")
    start_time = time.time()

    for i in range(n_target):
        for _ in range(max_attempts):

            # Generate random point: 
            # x and y in (-length/2 + r, length/2 - r) 
            # z in (offset + r, height + offset - r)
            pt = np.array([np.random.uniform(-length/2 + margin, length/2 - margin),
                           np.random.uniform(-length/2 + margin, length/2 - margin),
                           np.random.uniform(offset + margin, height + offset - margin)])
            
            # Check overlap with existing spheres
            if existing is not None:
                dists_sq = np.sum((existing - pt)**2, axis=1)
                if np.any(dists_sq < min_dist_sq):
                    continue  # Overlap detected
            
            # Check overlap with newly placed spheres 
            if len(spheres) > 0:
                coords = np.array([[s['x'], s['y'], s['z']] for s in spheres], dtype=float)
                dists_sq = np.sum((coords - pt)**2, axis=1)
                if np.any(dists_sq < min_dist_sq):
                    continue  # Overlap detected
            
            spheres.append({'name': name, 'r': r, 'x': pt[0], 'y': pt[1], 'z': pt[2]})
            break  # Successfully placed

    print(f"Packed {len(spheres)}/{n_target} {name} spheres in {time.time()-start_time:.2f}s")

    return pd.DataFrame(spheres, columns=['name', 'r', 'x', 'y', 'z'])

