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

    os.makedirs(f"./Figures/Data/", exist_ok=True)
    os.makedirs(f"./OpenMC/", exist_ok=True)


    for breeder in ['flibe']:
        for fertile_element in ['U']:
            for fd_kgm3 in FERTILE_BULK_DENSITY_KGM3:
                
                current_run = Reactor(breeder=breeder,fertile_element=fertile_element,fertile_density_kgm3=fd_kgm3, calc_volumes=False, run_openmc=True)

                print(f"Check if '{current_run.path}' exists: {os.path.isdir(current_run.path)}")

                if current_run.calc_volumes:
                    if os.path.exists(f"{current_run.path}/volume_1.h5"):
                        print(f"{Colors.YELLOW}Warning. {Colors.END}File {current_run.path}/volume_1.h5 already exists, so this OpenMC volume calculation will be skipped...")
                    else:
                        current_run.volumes()

                if current_run.run_openmc:
                    if os.path.isdir(current_run.path):
                        print(f"{Colors.YELLOW}Warning. {Colors.END}Directory {current_run.path} already exists, so this OpenMC run will be skipped...")
                        continue
                    else:
                        current_run.openmc()

                



if __name__ == '__main__':
    main()