import openmc
import os, sys
import numpy as np
import pandas as pd

# Import helper functions
from Python.reactor import *
from Python.parameters import *
from Python.utilities import *


@timer
def main():

    for breeder in ['flibe']:
        for fertile_element in ['U']:
            for fd_kgm3 in FERTILE_BULK_DENSITY_KGM3[0:1]:
                current_run = Reactor(breeder=breeder,fertile_element=fertile_element,fertile_density_kgm3=fd_kgm3, run_openmc=True)

                print(f"Check if '{current_run.path}' exists: {os.path.isdir(current_run.path)}")

                if current_run.run:
                    if os.path.isdir(current_run.path):
                        print(f"{Colors.YELLOW}Warning. {Colors.END}Directory {current_run.path} already exists, so the OpenMC run will be skipped...")
                        continue
                    else:
                        current_run.run_openmc()

                os.mkdir(f"./Figures/Data/", exists_ok=True)



if __name__ == '__main__':
    main()