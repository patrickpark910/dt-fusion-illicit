import subprocess
import sys
import os

os.environ["OPENMC_CROSS_SECTIONS"] = "/Users/gretali/Desktop/research2025/FissileDependence/endfb-viii.0-hdf5/cross_sections.xml"

fertile_densities = [150, 120, 30, 15, 1.5, 0.5]
dir = "/Users/gretali/Desktop/research2025/FissileDependence/dt-fusion-illicit/SpatialSelfShielding/HCPB_het/U_"

command = ["openmc"]

for density in fertile_densities:
    try:
        # Run the command and wait for it to complete
        print(density)
        process = subprocess.Popen(command, cwd=dir + str(density), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(f"[{density}] {line}", end="")

        process.wait()
        if process.returncode == 0:
            print(f"Command for density {density} finished successfully.\n")
        else:
            print(f"Command for density {density} failed with return code {process.returncode}\n")

    except FileNotFoundError:
        print(f"Command not found: {' '.join(command)}")