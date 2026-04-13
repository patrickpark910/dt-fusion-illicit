import os
import openmc
import pandas as pd
import numpy as np

def extract_rates(breeder, breeder_enrich, isotope):
    """
    Extracts the TBR and (n,gamma) reaction rates along with their 
    standard deviations from OpenMC statepoint files.
    """
    base_dir = "./OpenMC"
    results = []
    breeder = breeder.lower()
    
    if not os.path.exists(base_dir):
        print(f"Directory {base_dir} not found. Run from Jupyter/Wedge")
        return

    print(f"{breeder}_Li{breeder_enrich}_wedgeA")

    folders = [f for f in os.listdir(base_dir) if f.startswith((f"{breeder}_Li{breeder_enrich}_wedgeA", f"{breeder}_Li{breeder_enrich}_wedgeC")) and isotope in f]
    
    for folder in sorted(folders):
        
        path = os.path.join(base_dir, folder)
        if not os.path.isdir(path):
            continue
  
        parts = folder.split('_')

        case = parts[2][-1] 
        folder_iso = parts[3]
        
        if folder_iso != isotope:
            continue
            
        loading_str = parts[4].replace('kgm3', '')
        try:
            loading = float(loading_str)
        except ValueError:
            continue
            
        sp_files = [f for f in os.listdir(path) if f.startswith('statepoint.') and f.endswith('.h5')]
        if not sp_files:
            print(f"No statepoint found in {folder}")
            continue
        
        sp_path = os.path.join(path, max(sp_files, key=lambda x: int(x.split('.')[1])))
            
        try:
            sp = openmc.StatePoint(sp_path)
            
            # Extract (n,Xt) / (n,t) — tritium production
            # Try new tally names first, fall back to old names
            try:
                li_tally = sp.get_tally(name="Li rxn rates total")
                li_df = li_tally.get_pandas_dataframe()
                li_xt = li_df[li_df["score"] == "(n,Xt)"]
                tbr_mean = li_xt["mean"].sum()
                tbr_std  = np.sqrt((li_xt["std. dev."]**2).sum())
            except LookupError:
                li_tally = sp.get_tally(name="Total (n,t) rxn rate")
                li_df = li_tally.get_pandas_dataframe()
                tbr_mean = li_df["mean"].values[0]
                tbr_std  = li_df["std. dev."].values[0]
            
            # Extract (n,gamma) for the fertile isotope
            try:
                fertile_tally = sp.get_tally(name="Fertile rxn rates total")
            except LookupError:
                fertile_tally = sp.get_tally(name="Total fertile rxn rate")
            
            fertile_slice = fertile_tally.get_slice(nuclides=[isotope], scores=['(n,gamma)'])
            f_df = fertile_slice.get_pandas_dataframe()
            ng_mean = f_df["mean"].values[0]
            ng_std  = f_df["std. dev."].values[0]
            
            # Extract fission rate for the fertile isotope
            fission_slice = fertile_tally.get_slice(nuclides=[isotope], scores=['fission'])
            fission_df = fission_slice.get_pandas_dataframe()
            fission_mean = fission_df["mean"].values[0]
            fission_std  = fission_df["std. dev."].values[0]
            
            results.append({
                "folder": folder,
                "case": case,
                "loading_kg_m3": loading,
                "tbr": tbr_mean,
                "tbr_std": tbr_std,
                f"{isotope.lower()}_gamma": ng_mean,
                f"{isotope.lower()}_gamma_std": ng_std,
                f"{isotope.lower()}_fis": fission_mean,
                f"{isotope.lower()}_fis_std": fission_std
            })
            
            print(f"Extracted {folder}")
                
        except Exception as e:
            print(f"Error processing {folder}: {e}")
            
    if results:
        df = pd.DataFrame(results)
        out_path = f"./Figures/extracted_rates_{breeder}_{isotope}.csv"
        df.to_csv(out_path, index=False)
        print(f"\nSaved {len(results)} records to {os.path.abspath(out_path)}")
    else:
        print("No data extracted.")

if __name__ == "__main__":
    for (b, e) in [('dcll', 90.0), ('hcpb', 60.0)]:
        for isotope in ['U238', 'Th232']:
            extract_rates(b, e, isotope)