import openmc
import os, sys
import argparse
import numpy as np
import pandas as pd

# Import helper functions
from Python.reactor import *
# from Python.plot import *
from Python.parameters import *
from Python.utilities import *


@timer
def main():

    os.makedirs(f"./Figures/Data/", exist_ok=True)
    os.makedirs(f"./OpenMC/", exist_ok=True)

    parser = argparse.ArgumentParser(description="Choose run type with -r flag: ['tallies', 'volume', 'plot']")
    parser.add_argument("-r", "--run_type", 
                        type=str, default='tallies', 
                        help="Specify run type: ['tallies', 'volume', 'plot']" )
    args = parser.parse_args()
    run_type = args.run_type


    if run_type == 'plot':
        for breeder in ['arc','flibe','ll']:
            current_run = Reactor(breeder=breeder, run_type='plot')

            if os.path.exists(f"{current_run.path}/{current_run.breeder_name}_xz.ppm"):
                print(f"{Colors.YELLOW}Warning.{Colors.END} File {current_run.path}/{current_run.breeder_name}_xz.ppm already exists, so this OpenMC volume calculation will be skipped...")
            else:
                current_run.plot()


    elif run_type == 'volume':
        for breeder in ['arc','flibe','ll']:
            current_run = Reactor(breeder=breeder, run_type='volume')

            if os.path.exists(f"{current_run.path}/volume_1.h5"):
                print(f"{Colors.YELLOW}Warning.{Colors.END} File {current_run.path}/volume_1.h5 already exists, so this OpenMC volume calculation will be skipped...")
            else:
                current_run.volumes()

    elif run_type == 'tallies':

        for breeder in ['arc','flibe','ll']: # 'arc','flibe',
            for fertile_element in ['U']:
                for fbd_kgm3 in FERTILE_BULK_DENSITY_KGM3: # [FERTILE_BULK_DENSITY_KGM3[0]]: # 

                    cv, pp = False, False
                    #if fbd_kgm3 == 0:
                    #    cv, pp = True, True
                    
                    current_run = Reactor(breeder=breeder,
                                          fertile_element=fertile_element,
                                          fertile_bulk_density_kgm3=fbd_kgm3, 
                                          run_type='tallies')

                    print(f"Check if '{current_run.path}' exists: {os.path.isdir(current_run.path)}")

                    if has_statepoint(current_run.path):
                        print(f"{Colors.YELLOW}Warning.{Colors.END} File {current_run.path}/statepoint.h5 already exists, so this OpenMC run will be skipped...")
                    else:
                        current_run.openmc()

                



if __name__ == '__main__':
    main()