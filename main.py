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
            for fertile_element in ['U','Th']: # ,'Th']:
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

                print(f"Collating tallies for {breeder} {fertile_element} at {current_run.temp_k} and breeder vol {current_run.breeder_volume} m3")
                collate_tallies(breeder, fertile_element, current_run.temp_k, current_run.breeder_volume)




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


def collate_tallies(breeder,fertile_element,temp_k,vol_m3):

    if fertile_element == 'U':
        df_all   = pd.DataFrame(columns=['filename','fertile_kg/m3', 'fertile_mt', 'Li6(n,t)', 'Li7(n,Xt)','U238(n,g)','tbr','Pu239_kg/yr'])
        fertile_isotope = 'U238'
        fissile_isotope = 'Pu239'

    elif fertile_element == 'Th':
        df_all = pd.DataFrame(columns=['filename','fertile_kg/m3', 'fertile_mt', 'Li6(n,t)', 'Li7(n,Xt)','Th232(n,g)','tbr','Pu239_kg/yr'])
        fertile_isotope = 'Th232'
        fissile_isotope =  'U233'

    df_ebins = pd.DataFrame(columns=['filename','fertile_kg/m3', 'fertile_mt', 'energy mid [eV]', 'mean'])
    tally_folders = [x for x in os.listdir("./OpenMC/") if (x.startswith(f"tallies_{breeder}_{temp_k}K")) and x.split("_")[-1].startswith(fertile_element)]

    for folder in tally_folders:

        # Extract the fertile loading
        part = folder.split("_")[-1]
        fertile = float(part.replace("kgm3", "").lstrip(fertile_element))
        mt = fertile*vol_m3/1e3 # metric tons of fertile isotope

        tally_summary = f"./OpenMC/{folder}/tallies_summary.csv"
        tally_Ebins   = f"./OpenMC/{folder}/{fertile_isotope}_n-gamma_Ebins.csv"
        
        try:
            df = pd.read_csv(tally_summary)
        except:
            print(f"{Colors.YELLOW}Warning.{Colors.END} File 'tallies_summary.csv' not found in {folder}, skipping...")
            continue

        try:
            df_bins = pd.read_csv(tally_Ebins)
        except:
            print(f"{Colors.YELLOW}Warning.{Colors.END} File '{fertile_isotope}_n-gamma_Ebins.csv' not found in {folder}, skipping...")
            continue

        li6  = df[ df['cell']=='total' ]['Li6(n,t)'].values[0]
        li7  = df[ df['cell']=='total' ]['Li7(n,t)'].values[0]
        u238 = df[ df['cell']=='total' ][f'{fertile_isotope}(n,g)'].values[0]
        tbr  = df[ df['cell']=='total' ]['tbr'].values[0]
        pu   = df[ df['cell']=='total' ][f'{fissile_isotope}_kg/yr'].values[0]

        Ebins = df_bins.groupby('energy mid [eV]')['mean'].sum().reset_index()
        Ebins['filename']      = folder
        Ebins['fertile_kg/m3'] = fertile
        Ebins['fertile_mt']    = mt
        df_ebins = pd.concat([df_ebins, Ebins])
        df_ebins.to_csv(f"./Figures/Data/{breeder}_{fertile_element}_n-gamma_{temp_k}K.csv",index=False)

        df_all.loc[len(df_all)] = {'filename': folder,
                              'fertile_kg/m3': fertile,
                                 'fertile_mt': mt,
                                   'Li6(n,t)': li6,
                                  'Li7(n,Xt)': li7,
                    f'{fertile_isotope}(n,g)': u238,
                                        'tbr': tbr,
                   f'{fissile_isotope}_kg/yr': pu }

    dst = f"./Figures/Data/{breeder}_{fertile_element}_rxns_{temp_k}K.csv"
    # df_all.to_csv(dst, index=False)
    # print(f"{Colors.GREEN}Comment.{Colors.END} Collated tallies for {breeder} at {temp_k} K to: {dst}")
                          



if __name__ == '__main__':

    main()