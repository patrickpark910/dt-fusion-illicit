import numpy as np
import matplotlib.pyplot as plt

from utilities import *

# Set parameters
R0 = 4.0      # Major radius in meters
a = 1.2       # Minor radius in meters
kappa = 1.84  # Elongation (>1 for vertical elongation)
delta = 0.5   # Triangularity (positive for D-shape)
d = 1     # Offset distance (1 cm = 0.01 m)

# Create parameter array
t = np.linspace(0, 2*np.pi, 1000)

# Calculate original shape
R_orig, Z_orig = miller_d_shape(t, R0, a, kappa, delta)

# Calculate offset shape
R_offset, Z_offset = miller_d_shape_offset(t, R0, a, kappa, delta, d)

# Create the plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot 1: Both contours
ax1.plot(R_orig, Z_orig, 'b-', linewidth=2, label='Original Miller D-shape')
ax1.plot(R_offset, Z_offset, 'r-', linewidth=2, label=f'Offset ({d*100:.1f} cm)')
ax1.set_xlabel('R (m)', fontsize=12)
ax1.set_ylabel('Z (m)', fontsize=12)
ax1.set_title('Miller D-shape with Offset Contour', fontsize=14)
ax1.grid(True, alpha=0.3)
ax1.axis('equal')
ax1.legend(loc='best', fontsize=10)

# Plot 2: Zoomed in view to see the offset clearly
# Find a region to zoom (e.g., the outer midplane)
zoom_center_R = R0 + a
zoom_center_Z = 0
zoom_window = 0.5

ax2.plot(R_orig, Z_orig, 'b-', linewidth=2, label='Original')
ax2.plot(R_offset, Z_offset, 'r-', linewidth=2, label=f'Offset ({d*100:.1f} cm)')

# Add some normal vectors for visualization
n_vectors = 20
t_vectors = np.linspace(0, 2*np.pi, n_vectors)
for t_val in t_vectors:
    R_pt, Z_pt = miller_d_shape(t_val, R0, a, kappa, delta)
    R_off_pt, Z_off_pt = miller_d_shape_offset(t_val, R0, a, kappa, delta, d)
    ax2.arrow(R_pt, Z_pt, R_off_pt - R_pt, Z_off_pt - Z_pt,
              head_width=0.02, head_length=0.01, fc='gray', ec='gray', alpha=0.5)

ax2.set_xlim([zoom_center_R - zoom_window, zoom_center_R + zoom_window])
ax2.set_ylim([zoom_center_Z - zoom_window, zoom_center_Z + zoom_window])
ax2.set_xlabel('R (m)', fontsize=12)
ax2.set_ylabel('Z (m)', fontsize=12)
ax2.set_title('Zoomed View with Normal Vectors', fontsize=14)
ax2.grid(True, alpha=0.3)
ax2.axis('equal')
ax2.legend(loc='best', fontsize=10)

plt.tight_layout()

# Print some information
print(f"Miller D-shape parameters:")
print(f"  Major radius R0 = {R0} m")
print(f"  Minor radius a = {a} m")
print(f"  Elongation κ = {kappa}")
print(f"  Triangularity δ = {delta}")
print(f"  Offset distance = {d*100} cm")
print(f"\nAspect ratio = {R0/a:.2f}")
print(f"Plasma height ≈ {2*kappa*a:.2f} m")
print(f"Plasma width ≈ {2*a:.2f} m")

# Calculate and print the perimeter (approximate)
dR = np.diff(R_orig)
dZ = np.diff(Z_orig)
perimeter_orig = np.sum(np.sqrt(dR**2 + dZ**2))
dR_off = np.diff(R_offset)
dZ_off = np.diff(Z_offset)
perimeter_offset = np.sum(np.sqrt(dR_off**2 + dZ_off**2))
print(f"\nApproximate perimeter:")
print(f"  Original: {perimeter_orig:.3f} m")
print(f"  Offset: {perimeter_offset:.3f} m")
print(f"  Difference: {perimeter_offset - perimeter_orig:.3f} m")

plt.show()