import numpy as np
import matplotlib.pyplot as plt
import time

def fertile_to_biso_density(fertile_kgm3, isotope='U238'):
    """Calculates and prints BISO density metrics."""
    vol_kernel = 4/3 * np.pi * 0.04**3
    vol_biso   = 4/3 * np.pi * 0.05**3
    biso_per_cc = fertile_kgm3 * (238.02+32)/238.05 / vol_kernel / 10.5 / 100**3 * 1000
    biso_vol_percent = biso_per_cc * vol_biso * 100
    
    print(f"{fertile_kgm3:>8} kg/m³ | {biso_per_cc:>10.4f} part/cc | {biso_vol_percent:>8.4f} Vol %")


def pack_trapezoidal_prism(a, b, h, n_target, r=0.05, max_attempts=10000):
    """Packs spheres into a trapezoidal prism using Random Sequential Adsorption."""
    spheres = []
    
    # Pre-calculate wall slope and safety margin
    slope = (a - b) / (2.0 * h)
    margin = r / np.sqrt(1 + slope**2) # Adjust for wall slant
    min_dist_sq = (diam * 1.05) ** 2   # Minimum distance squared (with buffer)
    
    print(f"\nAttempting to pack {n_target} spheres...")
    start_time = time.time()

    for i in range(n_target):
        for _ in range(max_attempts):
            # 1. Generate Random Point
            w_max = max(a, b)
            pt = np.array([np.random.uniform(-w_max/2, w_max/2),
                           np.random.uniform(-w_max/2, w_max/2),
                           np.random.uniform(0, h)])
            
            # 2. Check Boundary (Inside Prism?)
            current_half_width = (b + (a - b) * (pt[2] / h)) / 2.0
            limit = current_half_width - margin
            
            if (pt[2] < r or pt[2] > h - r or abs(pt[0]) > limit or abs(pt[1]) > limit):
                continue # Out of bounds

            # 3. Check Overlap (Vectorized)
            if len(spheres) > 0:
                # Compute distance to all existing spheres at once
                dists_sq = np.sum((np.array(spheres) - pt)**2, axis=1)
                if np.any(dists_sq < min_dist_sq):
                    continue # Overlap detected

            spheres.append(pt)
            break # Successfully placed

    print(f"Packed {len(spheres)}/{n_target} spheres in {time.time()-start_time:.2f}s")
    return np.array(spheres)

def visualize_and_save(spheres, a, b, h, filename="packing_viz.png"):
    """Visualizes the packing and saves the figure."""
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot Spheres (Scatter is faster than surface plots)
    if len(spheres) > 0:
        ax.scatter(spheres[:,0], spheres[:,1], spheres[:,2], s=20, alpha=0.6, c='b')

    # Draw Wireframe
    wb, wa = b/2, a/2
    corners = np.array([[-wb,-wb,0], [wb,-wb,0], [wb,wb,0], [-wb,wb,0], [-wb,-wb,0],
                        [-wa,-wa,h], [wa,-wa,h], [wa,wa,h], [-wa,wa,h], [-wa,-wa,h]])
    
    # Draw top/bottom squares and connecting lines
    ax.plot(corners[:5,0], corners[:5,1], corners[:5,2], 'k-')
    ax.plot(corners[5:,0], corners[5:,1], corners[5:,2], 'k-')
    for i in range(4):
        ax.plot([corners[i,0], corners[i+5,0]], [corners[i,1], corners[i+5,1]], [corners[i,2], corners[i+5,2]], 'k-')

    ax.set_title(f"Packed {len(spheres)} Spheres"); ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    max_range = max(a, b, h) / 2
    ax.set_xlim(-max_range, max_range); ax.set_ylim(-max_range, max_range); ax.set_zlim(0, h)
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"Visualization saved to {filename}")

if __name__ == "__main__":
    # --- Configuration ---
    a, b, h = 0.6370, 0.7155, 45.0000
    N_target = round(1647 * 75/80)
    densities = [0, 0.03, 0.3, 0.6, 1.5, 3, 7.5, 10, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 250, 500, 750, 1000]

    # --- Run Density Calcs ---
    print("--- Density Calculations ---")
    for f in densities:
        fertile_to_biso_density(f)

    # --- Run Packing ---
    spheres = pack_trapezoidal_prism(a, b, h, N_target, r=0.05)
    
    # --- Save Results ---
    np.savetxt('sphere_positions.csv', spheres, delimiter=',', header='x,y,z', comments='')
    visualize_and_save(spheres, a, b, h)