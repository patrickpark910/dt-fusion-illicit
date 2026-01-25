import numpy as np

vf_fertile = 100 / (6.7 * 1000)  # 1 g/cm³ = 1000 kg/m³
fractions  = np.array([1, vf_fertile])
fractions /= np.array([1, vf_fertile]).sum()
vf_flibe_new, vf_fertile_new = fractions

print(vf_flibe_new, vf_fertile_new)