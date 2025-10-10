import openmc
import os, sys
import argparse
import numpy as np
import pandas as pd

# Import helper functions
from Python.reactor import *
from Python.arc     import *
from Python.ball     import *
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
        for breeder in ['ARC','ARCBall','FLiBe','LL']: # make this match class name
            
            current_run = build_reactor(breeder, breeder_name=breeder, run_type='plot', run_openmc=True)

            if os.path.exists(f"{current_run.path}/{current_run.breeder_name}_xz.ppm"):
                print(f"{Colors.YELLOW}Warning.{Colors.END} File {current_run.path}/{current_run.breeder_name}_xz.ppm already exists, so this OpenMC volume calculation will be skipped...")
            elif current_run.run_openmc:
                current_run.plot()


    elif run_type == 'volume':
        for breeder in ['ARC','ARCBall','FLiBe','LL']: # make this match class name

            current_run = build_reactor(breeder, breeder_name=breeder, run_type='volume', run_openmc=True)
            
            if os.path.exists(f"{current_run.path}/volume_1.h5"):
                print(f"{Colors.YELLOW}Warning.{Colors.END} File {current_run.path}/volume_1.h5 already exists, so this OpenMC volume calculation will be skipped...")
            elif current_run.run_openmc:
                current_run.volumes()


    elif run_type == 'tallies':

        for breeder in ['ARC','ARCBall','FLiBe','LL']: # make this match class name
            for fertile_element in ['U']:
                for fbd_kgm3 in FERTILE_BULK_DENSITY_KGM3: # [FERTILE_BULK_DENSITY_KGM3[0]]: # 
                    
                    current_run = build_reactor(breeder, breeder_name=breeder,
                                                fertile_element=fertile_element,
                                                fertile_bulk_density_kgm3=fbd_kgm3, 
                                                run_type='tallies',
                                                run_openmc=True)

                    print(f"Check if '{current_run.path}' exists: {os.path.isdir(current_run.path)}")

                    if has_statepoint(current_run.path):
                        print(f"{Colors.YELLOW}Warning.{Colors.END} File {current_run.path}/statepoint.h5 already exists, so this OpenMC run will be skipped...")
                    elif current_run.run_openmc:
                        current_run.openmc()

                    if os.path.exists(f"{current_run.path}/tallies_summary.csv"): 
                        print(f"{Colors.YELLOW}Warning.{Colors.END} File {current_run.path}/tallies_summary.csv already exists, so tally extraction will be skipped...")
                    elif current_run.run_openmc:
                        current_run.extract_tallies()

            collate_tallies(breeder, current_run.temp_k, current_run.breeder_volume)




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


def collate_tallies(breeder,temp_k,vol_m3):

    df_all = pd.DataFrame(columns=['filename','fertile_kg/m3', 'fertile_mt', 'tbr','Pu239_kg/yr'])
    
    tally_folders = [x for x in os.listdir("./OpenMC/") if x.startswith(f"tallies_{breeder}_{temp_k}K")]

    for folder in tally_folders:

        # Extract the fertile loading
        for part in folder.split("_"):
            if part.startswith("U") or part.startswith("Th"):
                fertile = float(part.replace("kgm3", "").lstrip("UTh"))
                mt = fertile*vol_m3/1e3 # metric tons of fertile isotope

        tally_summary = f"./OpenMC/{folder}/tallies_summary.csv"
        df = pd.read_csv(tally_summary)

        tbr = df[ df['cell']=='total' ]['tbr'].values[0]
        pu  = df[ df['cell']=='total' ]['Pu239_kg/yr'].values[0]

        df_all.loc[len(df_all)] = {'filename': folder,
                              'fertile_kg/m3': fertile,
                                 'fertile_mt': mt,
                                        'tbr': tbr,
                                'Pu239_kg/yr': pu }

    dst = f"./Figures/Data/{breeder}_total_rxns_{temp_k}K.csv"
    df_all.to_csv(dst, index=False)
    print(f"{Colors.GREEN}Comment.{Colors.END} Collated tallies for {breeder} at {temp_k} K to: {dst}")
                          



if __name__ == '__main__':

    main()