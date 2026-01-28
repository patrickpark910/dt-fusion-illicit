import openmc
import numpy as np
import matplotlib.pyplot as plt

# 1. Load the tracks
tracks = openmc.Tracks('tracks.h5')

# 2. Extract all coordinates to find the bounds
all_coords = []
for track in tracks:
    for p_track in track.particle_tracks:
        # p_track.states is a numpy record array containing 'xyz'
        all_coords.append(p_track.states['xyz'])

# Combine all points into one large array (N, 3)
coords = np.concatenate(all_coords)

# Calculate min/max for X, Y, and Z
min_xyz = coords.min(axis=0)
max_xyz = coords.max(axis=0)

# 3. Create the plot
ax = tracks.plot()

# 4. Apply a tight fit with a small margin (e.g., 5%)
margin = 0.05
for i, set_lim in enumerate([ax.set_xlim, ax.set_ylim, ax.set_zlim]):
    span = max_xyz[i] - min_xyz[i]
    # If the particle is stuck at a single point, span might be 0
    if span == 0: span = 1.0 
    
    set_lim(min_xyz[i] - margin * span, max_xyz[i] + margin * span)

plt.show()