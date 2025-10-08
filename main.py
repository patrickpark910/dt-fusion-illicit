import openmc
import os, sys
import argparse
import numpy as np
import pandas as pd

# Import helper functions
from Python.reactor import *
from Python.arc     import *
from Python.flibe   import *
from Python.pbli    import *
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
        for breeder in ['ARC','FLiBe','LL']: # make this match class name
            
            current_run = build_reactor(breeder, breeder_name=breeder, run_type='plot')

            if os.path.exists(f"{current_run.path}/{current_run.breeder_name}_xz.ppm"):
                print(f"{Colors.YELLOW}Warning.{Colors.END} File {current_run.path}/{current_run.breeder_name}_xz.ppm already exists, so this OpenMC volume calculation will be skipped...")
            else:
                current_run.plot()


    elif run_type == 'volume':
        for breeder in ['ARC','FLiBe','LL']: # make this match class name

            current_run = build_reactor(breeder, breeder_name=breeder, run_type='volume')
            
            if os.path.exists(f"{current_run.path}/volume_1.h5"):
                print(f"{Colors.YELLOW}Warning.{Colors.END} File {current_run.path}/volume_1.h5 already exists, so this OpenMC volume calculation will be skipped...")
            else:
                current_run.volumes()


    elif run_type == 'tallies':

        for breeder in ['ARC','FLiBe','LL']: # make this match class name
            for fertile_element in ['U']:
                for fbd_kgm3 in FERTILE_BULK_DENSITY_KGM3: # [FERTILE_BULK_DENSITY_KGM3[0]]: # 
                    
                    current_run = build_reactor(breeder, breeder_name=breeder,
                                                fertile_element=fertile_element,
                                                fertile_bulk_density_kgm3=fbd_kgm3, 
                                                run_type='tallies')

                    print(f"Check if '{current_run.path}' exists: {os.path.isdir(current_run.path)}")

                    if has_statepoint(current_run.path):
                        print(f"{Colors.YELLOW}Warning.{Colors.END} File {current_run.path}/statepoint.h5 already exists, so this OpenMC run will be skipped...")
                    else:
                        current_run.openmc()

                    current_run.extract_tallies()


def build_reactor(breeder:str, **kwargs):
    """
    Create reactor model from mapping into abstract base class (ABC)
    """
    try:
        cls = eval(breeder) 
        reactor = cls(**kwargs)
        reactor.initialize()
    except KeyError:
        raise ValueError(f"Unknown breeder '{breeder}'.")
    return reactor



if __name__ == '__main__':

    main()