import openmc
import os, sys
import numpy as np

# Import helper functions
from Python.device import *
from Python.parameters import *
from Python.utilities import *


def main():
    print("\n\n")
    print("="*42)
    print(f"Running FLiBe Thorium OpenMC model for Li-6 enrichment: {7.5} wt%")

    current_run = Reactor()

    print(f"Check if '{current_run.path}' exists: {os.path.isdir(current_run.path)}")

    if os.path.isdir(current_run.path):
        print(f"Warning. Directory {current_run.path} already exists, so running OpenMC will fail. Skipping...")
        return
    else:
        current_run.set_xs_path()
        current_run.run_openmc()

if __name__ == '__main__':
    main()