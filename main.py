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

    args = parser.parse_args()
    run_type = args.run_type
    blankets = args.blankets # [b.upper() if isinstance(b, str) else b for b in args.blankets]  # too extra to make exception for FLiBe, just type it right yourself :D --ppark
    isotopes = args.isotopes


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
                for fertile_kgm3 in FERTILE_KGM3: 
                    
                    current_run = build_reactor(blanket,
                                                blanket_name=blanket,
                                                fertile_isotope=fertile_isotope,
                                                fertile_kgm3=fertile_kgm3, 
                                                run_type='tallies',
                                                run_openmc=True,
                                                run_debug=True)

                    print(f"Check if '{current_run.path}' exists: {os.path.isdir(current_run.path)}")

                    if current_run.run_openmc:
                        if has_statepoint(current_run.path):
                            print(f"{C.YELLOW}Warning.{C.END} File {current_run.path}/statepoint.h5 already exists, so this OpenMC run will be skipped...")
                        else:
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


def collate_tallies(blanket, fertile_isotope, breeder_enrich,temp_k, vol_m3):
    """
    Collates all the tallies for given [blanket, fertile isotope, temperature] across multiple fertile kg/m³

    Args:
        blanket (str): one of ['FLiBe', 'DCLL', 'HCPB', 'ARC', 'ARCB']
        fertile_isotope (str): one of ['U238', 'Th232']
        breeder_enrich (float): enrichment of lithium-6 in breeder
        temp_k (float): temperature of the system for which you want to collate data 
        breeder_volume (float): [m³] volume of breeder
    """

    if fertile_isotope == 'U238':        
        fissile_isotope = 'Pu239'

    elif fertile_isotope == 'Th232':
        fissile_isotope =  'U233'

    df_all = pd.DataFrame(columns=['filename','fertile_kg/m3', 'fertile_mt', 
                                   'Li6(n,t)', 'Li7(n,Xt)', f'{fertile_isotope}(n,g)',
                                   'tbr', f'{fissile_isotope}_kg/yr'])

    tally_folders = [x for x in os.listdir("./OpenMC/") if (x.startswith(f"tallies_{blanket}_{temp_k}K_Li{breeder_enrich:04.1f}")) and x.split("_")[-2].startswith(fertile_isotope)]

    flux_list, ng_list = [], []  # Use lists instead of empty DataFrames

    for folder in tally_folders:

        # Extract the fertile loading
        part = folder.split("_")[-1]
        fertile = float(part.replace("kgm3", ""))
        mt = fertile*vol_m3/1e3 # metric tons of fertile isotope

        tally_summary = f"./OpenMC/{folder}/tallies_summary.csv"
        tally_ng      = f"./OpenMC/{folder}/{fertile_isotope}_n-gamma_Ebins.csv"
        tally_flux    = f"./OpenMC/{folder}/{fertile_isotope}_flux_Ebins.csv"


        try:
            df = pd.read_csv(tally_summary)
        except FileNotFoundError:
            print(f"{C.YELLOW}Warning.{C.END} File 'tallies_summary.csv' not found in {folder}, skipping...")
            continue

        try:
            df_ng = pd.read_csv(tally_ng)
        except FileNotFoundError:
            print(f"{C.YELLOW}Warning.{C.END} File '{fertile_isotope}_n-gamma_Ebins.csv' not found in {folder}, skipping...")
            continue

        try:
            df_flux = pd.read_csv(tally_flux)
        except FileNotFoundError:
            print(f"{C.YELLOW}Warning.{C.END} File '{fertile_isotope}_flux_Ebins.csv' not found in {folder}, skipping...")
            continue


        li6  = df[ df['cell']=='total' ]['Li6(n,t)'].values[0]
        li7  = df[ df['cell']=='total' ]['Li7(n,t)'].values[0]
        u238 = df[ df['cell']=='total' ][f'{fertile_isotope}(n,g)'].values[0]
        tbr  = df[ df['cell']=='total' ]['tbr'].values[0]
        pu   = df[ df['cell']=='total' ][f'{fissile_isotope}_kg/yr'].values[0]

        li6_stdev  = df[ df['cell']=='total' ]['Li6(n,t)_stdev'].values[0]
        li7_stdev  = df[ df['cell']=='total' ]['Li7(n,t)_stdev'].values[0]
        u238_stdev = df[ df['cell']=='total' ][f'{fertile_isotope}(n,g)_stdev'].values[0]
        tbr_stdev  = df[ df['cell']=='total' ]['tbr_stdev'].values[0]
        pu_stdev   = df[ df['cell']=='total' ][f'{fissile_isotope}_kg/yr_stdev'].values[0]

        df_all.loc[len(df_all)] = {'filename': folder,
                              'fertile_kg/m3': fertile,
                                 'fertile_mt': mt,
                                   'Li6(n,t)': li6,
                             'Li6(n,t)_stdev': li6_stdev,
                                  'Li7(n,Xt)': li7,
                            'Li7(n,Xt)_stdev': li7_stdev,
                    f'{fertile_isotope}(n,g)': u238,
              f'{fertile_isotope}(n,g)_stdev': u238_stdev,
                                        'tbr': tbr,
                                  'tbr_stdev': tbr_stdev,
                   f'{fissile_isotope}_kg/yr': pu,
             f'{fissile_isotope}_kg/yr_stdev': pu_stdev, }

        dst = f"./Figures/Data/{blanket}_{temp_k}K_Li{breeder_enrich:04.1f}_{fertile_isotope}"
        df_all.to_csv(f"{dst}_rxns.csv", index=False)


        cols = ["energy low [eV]", "energy high [eV]", "energy mid [eV]", "mean"]     

        # Group by energy midpoint and sum mean values across all cells while preserving energy bin boundaries  
        # groupby("energy mid [eV]", as_index=False) - Groups all rows that have the same "energy mid [eV]" value together
        # .agg(**{...}) - "output_column_name": ("input_column_name", "aggregation_function")
        # "energy low [eV]": ("energy low [eV]", "first") - Takes the first "energy low [eV]" value in the group
        # "energy high [eV]": ("energy high [eV]", "first") - Takes the first "energy high [eV]" value in the group
        # "mean": ("mean", "sum") - Sums all "mean" values in the group

        ng =  (df_ng[cols].groupby("energy mid [eV]", as_index=False)
                          .agg(**{"energy low [eV]" : ("energy low [eV]", "first"),
                                  "energy high [eV]": ("energy high [eV]", "first"),
                                              "mean": ("mean", "sum"),} ) )

        flux = (df_flux[cols].groupby("energy mid [eV]", as_index=False)
                             .agg(**{"energy low [eV]" : ("energy low [eV]", "first"),
                                    "energy high [eV]": ("energy high [eV]", "first"),
                                                "mean": ("mean", "sum"),} ) )

        ng['filename'],      flux['filename']      = folder, folder
        ng['fertile_kg/m3'], flux['fertile_kg/m3'] = fertile, fertile
        ng['fertile_mt'],    flux['fertile_mt']    = mt, mt
        ng['br_vol_m3'],     flux['br_vol_m3']     = vol_m3, vol_m3

        # avoid deprecation warning smh this change by pandas feels unnecessary
        ng_list.append(ng)
        flux_list.append(flux)

        # Concatenate once at the end
        df_ng_collated = pd.concat(ng_list, ignore_index=True) if ng_list else pd.DataFrame()
        df_flux_collated = pd.concat(flux_list, ignore_index=True) if flux_list else pd.DataFrame()

        # Reorder columns
        df_ng_collated = df_ng_collated[['filename', 'fertile_kg/m3', 'fertile_mt', 'br_vol_m3', 
                                         'energy mid [eV]', 'mean', 'energy low [eV]', 'energy high [eV]']]

        df_flux_collated = df_flux_collated[['filename', 'fertile_kg/m3', 'fertile_mt', 'br_vol_m3', 
                                             'energy mid [eV]', 'mean', 'energy low [eV]', 'energy high [eV]']]

        df_ng_collated.to_csv(f"{dst}_n-gamma.csv",index=False)
        df_flux_collated.to_csv(f"{dst}_flux.csv",index=False)


    print(f"{C.GREEN}Comment.{C.END} Collated tallies for {blanket} at {temp_k} K to: {dst}")


if __name__ == '__main__':
    main()