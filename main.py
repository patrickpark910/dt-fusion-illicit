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


def collate_tallies(blanket, fertile_isotope, breeder_enrich, temp_k, vol_m3):
    """
    Collates all the tallies for given [blanket, fertile isotope, temperature]
    across multiple fertile kg/m³ into combined CSVs in ./Figures/Data/.

    Args:
        blanket (str): one of ['FLiBe', 'DCLL', 'HCPB']
        fertile_isotope (str): one of ['U238', 'Th232']
        breeder_enrich (float): enrichment of lithium-6 in breeder
        temp_k (float): temperature of the system
        vol_m3 (float): [m³] volume of breeder
    """

    rows_all = []
    rxnsE_list, fluxE_list, leakE_list = [], [], []

    tally_folders = [x for x in os.listdir("./OpenMC/")
                     if x.startswith(f"tallies_{blanket}_{temp_k}K_Li{breeder_enrich:04.1f}")
                     and x.split("_")[-3].startswith(fertile_isotope)]

    dst = f"./Figures/Data/{blanket}_{temp_k}K_Li{breeder_enrich:04.1f}_{fertile_isotope}"

    for folder in tally_folders:

        # Extract fertile loading from folder name
        part = folder.split("_")[-2]
        fertile = float(part.replace("kgm3", ""))
        mt = fertile * vol_m3 / 1e3

        # File paths (must match output.py exports)
        tally_summary = f"./OpenMC/{folder}/tallies_summary.csv"
        tally_leak    = f"./OpenMC/{folder}/tallies_leakage.csv"
        tally_rxnsE   = f"./OpenMC/{folder}/Ebins_rxns.csv"
        tally_fluxE   = f"./OpenMC/{folder}/Ebins_flux.csv"
        tally_leakE   = f"./OpenMC/{folder}/Ebins_leakage.csv"

        # Collate tallies_summary.csv
        try:
            df = pd.read_csv(tally_summary)
        except FileNotFoundError:
            print(f"{C.YELLOW}Warning.{C.END} File 'tallies_summary.csv' not found in {folder}, skipping...")
            continue

        try:
            df_leak = pd.read_csv(tally_leak)
        except FileNotFoundError:
            print(f"{C.YELLOW}Warning.{C.END} File 'tallies_leakage.csv' not found in {folder}, skipping...")
            continue

        tot  = df[df['cell'] == 'total']
        leak = df_leak[df_leak['void_cell'] == 'total']

        rows_all.append({
            'filename':        folder,
            'fertile_kg/m3':   fertile,
            'fertile_mt':      mt,
            'Li6(n,t)':        tot['Li6(n,t)'].values[0],
            'Li6(n,t)_sd':     tot['Li6(n,t)_stdev'].values[0],
            'Li7(n,Xt)':       tot['Li7(n,t)'].values[0],
            'Li7(n,Xt)_sd':    tot['Li7(n,t)_stdev'].values[0],
            'Be9(n,2n)':       tot['Be9(n,2n)'].values[0],
            'Be9(n,2n)_sd':    tot['Be9(n,2n)_stdev'].values[0],
            'Pb(n,2n)':        tot['Pb(n,2n)'].values[0],
            'Pb(n,2n)_sd':     tot['Pb(n,2n)_stdev'].values[0],
            'U238(n,g)':       tot['U238(n,g)'].values[0],
            'U238(n,g)_sd':    tot['U238(n,g)_stdev'].values[0],
            'Th232(n,g)':      tot['Th232(n,g)'].values[0],
            'Th232(n,g)_sd':   tot['Th232(n,g)_stdev'].values[0],
            'tot(n,fis)':      tot['tot(n,fis)'].values[0],
            'tot(n,fis)_sd':   tot['tot(n,fis)_stdev'].values[0],
            'tbr':             tot['tbr'].values[0],
            'tbr_sd':          tot['tbr_stdev'].values[0],
            'Pu239_kg/yr':     tot['Pu239_kg/yr'].values[0],
            'Pu239_kg/yr_sd':  tot['Pu239_kg/yr_stdev'].values[0],
            'U233_kg/yr':      tot['U233_kg/yr'].values[0],
            'U233_kg/yr_sd':   tot['U233_kg/yr_stdev'].values[0],
            'leak [n/src-n]':  leak['leakage'].values[0],
            'leak_sd':         leak['leakage_stdev'].values[0],
            'heat [MW]':       tot['heating [MW]'].values[0],
            'heat_sd':         tot['heating_stdev'].values[0],
            'fisq [MW]':       tot['fisq [MW]'].values[0],
            'fisq_sd':         tot['fisq_stdev'].values[0],
        })

        # Collate Ebins_rxns.csv (already summed over cells)
        try:
            df_rxnsE = pd.read_csv(tally_rxnsE)
            df_rxnsE['filename']   = folder
            df_rxnsE['fertile_mt'] = mt
            df_rxnsE['br_vol_m3']  = vol_m3
            rxnsE_list.append(df_rxnsE)
        except FileNotFoundError:
            print(f"{C.YELLOW}Warning.{C.END} File 'Ebins_rxns.csv' not found in {folder}, skipping...")

        # Collate Ebins_flux.csv (still per-cell, needs groupby) 
        try:
            df_fluxE = pd.read_csv(tally_fluxE)
            cols = ['energy low [eV]', 'energy high [eV]', 'energy mid [eV]', 'mean']
            fluxE = (df_fluxE[cols]
                     .groupby('energy mid [eV]', as_index=False)
                     .agg(**{'energy low [eV]':  ('energy low [eV]', 'first'),
                             'energy high [eV]': ('energy high [eV]', 'first'),
                             'mean':             ('mean', 'sum')}))
            fluxE['filename']      = folder
            fluxE['fertile_kg/m3'] = fertile
            fluxE['fertile_mt']    = mt
            fluxE['br_vol_m3']     = vol_m3
            fluxE_list.append(fluxE)
        except FileNotFoundError:
            print(f"{C.YELLOW}Warning.{C.END} File 'Ebins_flux.csv' not found in {folder}, skipping...")

        # Collate Ebins_leakage.csv (already summed over cells)
        try:
            df_leakE = pd.read_csv(tally_leakE)
            df_leakE['filename']      = folder
            df_leakE['fertile_kg/m3'] = fertile
            df_leakE['fertile_mt']    = mt
            df_leakE['br_vol_m3']     = vol_m3
            leakE_list.append(df_leakE)
        except FileNotFoundError:
            print(f"{C.YELLOW}Warning.{C.END} File 'Ebins_leakage.csv' not found in {folder}, skipping...")

    # Export collated files
    if rows_all:
        df_all = pd.DataFrame(rows_all).sort_values(by='fertile_kg/m3', ascending=True)
        df_all.to_csv(f"{dst}_rxns.csv", index=False)

    if rxnsE_list:
        pd.concat(rxnsE_list, ignore_index=True).to_csv(f"{dst}_Ebins_rxns.csv", index=False)

    if fluxE_list:
        df_fluxE_collated = pd.concat(fluxE_list, ignore_index=True)
        df_fluxE_collated = df_fluxE_collated[['filename', 'fertile_kg/m3', 'fertile_mt', 'br_vol_m3',
                                               'energy mid [eV]', 'mean', 'energy low [eV]', 'energy high [eV]']]
        df_fluxE_collated.to_csv(f"{dst}_flux.csv", index=False)

    if leakE_list:
        pd.concat(leakE_list, ignore_index=True).to_csv(f"{dst}_leak.csv", index=False)

    print(f"{C.GREEN}Comment.{C.END} Collated tallies for {blanket} at {temp_k} K to: {dst}")


if __name__ == '__main__':
    main()