import openmc
import os, sys
import argparse
import numpy as np
import pandas as pd

# Import helper functions
from Python.reactor    import *
from Python.arc        import *
from Python.ball       import *
from Python.flibe      import *
from Python.dcll       import *
from Python.hcpb       import *
from Python.output     import *
from Python.parameters import *
from Python.utilities  import *


@timer
def main():

    os.makedirs(f"./Figures/Data/", exist_ok=True)
    os.makedirs(f"./OpenMC/", exist_ok=True)

    parser = argparse.ArgumentParser(description=f"Choose ONE run type with -r flag: ['tallies', 'volume', 'plot']. Choose blankets with -b flag, multiple separated by spaces: {BLANKETS}. Choose fertile isotopes with -i flag, multiple separated by spaces: {ISOTOPES}")
    
    parser.add_argument("-r", "--run_type", 
                        type=str, default='tallies', 
                        help="Specify run type: ['tallies', 'volume', 'plot']" )

    parser.add_argument("-b", "--blankets", 
                        type=str, nargs="+", default=BLANKETS,
                        help=f"Specify blankets, separated by space, among this list: {BLANKETS}" )

    parser.add_argument("-i", "--isotopes", 
                        type=str, nargs="+", default=ISOTOPES,
                        help=f"Specify fertile isotopes, separated by space, among this list: {ISOTOPES}" )
    
    parser.add_argument("-f", "--fertile",
                        type=float, nargs="+", default=FERTILE_KGM3,
                        help=f"Specify fertile loadings in kg/m3 (default: {FERTILE_KGM3})")
    
    parser.add_argument("-l", "--lithium", 
                        type=float, default=None,
                        help=f"Specify lithium enrichment. Only one value allowed" )
    
    parser.add_argument("-p", "--n_particles", 
                        type=lambda x: int(float(x)), default=N_PARTICLES,
                        help=f"Specify number of particles in integer scientific notation, e.g. 1e6" )
    
    parser.add_argument("-c", "--n_cycles", 
                        type=lambda x: int(float(x)), default=N_CYCLES,
                        help=f"Specify number of cycles in integers, e.g. 25" )
    
    parser.add_argument("--no_xml", 
                        dest="print_xml", action="store_false",
                        help="Disable printing model.xml files")
    
    parser.add_argument("--no_run", 
                        dest="run_openmc", action="store_false",
                        help="Disable running the OpenMC calculation")
    
    parser.add_argument("--no_debug", 
                        dest="run_debug", action="store_false",
                        help="Disable debugging statements")

    # Set defaults for the 'dest' targets above so they are True unless the flag is used
    parser.set_defaults(print_xml=True, run_openmc=True, run_debug=True)

    run_type   = parser.parse_args().run_type
    blankets   = parser.parse_args().blankets   # [b.upper() if isinstance(b, str) else b for b in args.blankets]  # too extra to make exception for FLiBe, just type it right yourself :D --ppark
    isotopes   = parser.parse_args().isotopes
    densities  = parser.parse_args().fertile
    lithium    = parser.parse_args().lithium
    print_xml  = parser.parse_args().print_xml
    run_openmc = parser.parse_args().run_openmc 
    run_debug  = parser.parse_args().run_debug
    n_particles = parser.parse_args().n_particles
    n_cycles    = parser.parse_args().n_cycles


    if run_type == 'plot':
        for blanket in blankets:  # make this match class name
            
            current_run = build_reactor(blanket, blanket_name=blanket, run_type='plot')

            if os.path.exists(f"{current_run.path}/{current_run.blanket_name}_xz.ppm"):
                print(f"{C.YELLOW}Warning.{C.END} File {current_run.path}/{current_run.blanket_name}_xz.ppm already exists, so this OpenMC volume calculation will be skipped...")
            elif current_run.run_openmc:
                current_run.plot()


    elif run_type == 'volume':
        for blanket in blankets:  # make this match class name

            current_run = build_reactor(blanket, blanket_name=blanket, run_type='volume')
            
            if os.path.exists(f"{current_run.path}/volume_1.h5"):
                print(f"{C.YELLOW}Warning.{C.END} File {current_run.path}/volume_1.h5 already exists, so this OpenMC volume calculation will be skipped...")
            elif current_run.run_openmc:
                current_run.volumes()


    elif run_type == 'tallies':

        for blanket in blankets:  # ['FLiBe', 'DCLL', 'HCPB', 'ARC', 'ARCB']  # make each blanket str match class name
            for fertile_isotope in isotopes:
                for fertile_kgm3 in densities: 
                    current_run = build_reactor(blanket,
                                                blanket_name=blanket,
                                                fertile_isotope=fertile_isotope,
                                                fertile_kgm3=fertile_kgm3, 
                                                lithium_enrich=lithium,
                                                run_type='tallies',
                                                n_particles=n_particles,
                                                n_cycles=n_cycles,
                                                print_xml=print_xml,
                                                run_openmc=run_openmc,
                                                run_debug=run_debug)

                    print(f"Check if '{current_run.path}' exists: {os.path.isdir(current_run.path)}")

                    current_run.compile()
                    current_run.openmc()
                    current_run.extract_tallies()
                    
                print(f"Collating tallies for {blanket} {fertile_isotope} at {current_run.temp_k}")
                collate_tallies(blanket, fertile_isotope, current_run.breeder_enrich, current_run.temp_k, current_run.breeder_volume)


def build_reactor(blanket:str, **kwargs):
    """
    Create reactor model from mapping into abstract base class (ABC)
    """
    try:
        cls = eval(blanket) 
        reactor = cls(**kwargs)
        reactor.initialize()
    except KeyError:
        raise ValueError(f"Unknown blanket '{blanket}'.")
    return reactor





if __name__ == '__main__':
    main()