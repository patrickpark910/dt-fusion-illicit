import os
import openmc
import pandas as pd
import numpy as np

def extract_rates(breeder, isotope):
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

    folders = [f for f in os.listdir(base_dir) if f.startswith((f"{breeder}_wedgeA", f"{breeder}_wedgeC")) and isotope in f]
    
    for folder in sorted(folders):
        path = os.path.join(base_dir, folder)
        if not os.path.isdir(path):
            continue
  
        parts = folder.split('_')
        if len(parts) < 3:
            continue
        case = parts[1][-1] 
        folder_iso = parts[2]
        
        if folder_iso != isotope:
            continue
            
        loading_str = parts[3].replace('kgm3', '')
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
            
            # --- Extract TBR (n,t) ---
            li_tally = sp.get_tally(name="Total (n,t) rxn rate")
            li_df = li_tally.get_pandas_dataframe()
            
            tbr_mean = float(li_df["mean"].sum())
            # Propagate error: sqrt(sum of squares of std. dev.)
            tbr_std = float(np.sqrt((li_df["std. dev."]**2).sum()))
            
            # --- Extract (n,gamma) for the specific isotope ---
            fertile_tally = sp.get_tally(name="Total fertile rxn rate")
            f_df = fertile_tally.get_pandas_dataframe()
            
            ng_mask = (f_df["nuclide"] == isotope) & (f_df["score"] == "(n,gamma)")
            
            ng_mean = float(f_df.loc[ng_mask, "mean"].sum())
            # Propagate error for the filtered subset
            ng_std = float(np.sqrt((f_df.loc[ng_mask, "std. dev."]**2).sum()))
            
            results.append({
                "folder": folder,
                "case": case,
                "loading_kg_m3": loading,
                "tbr": tbr_mean,
                "tbr_std": tbr_std,
                f"{isotope.lower()}_n_gamma": ng_mean,
                f"{isotope.lower()}_n_gamma_std": ng_std
            })
            
            print(f"Extracted {folder}")
            
        except Exception as e:
            print(f"Error processing {folder}: {e}")
            
    if results:
        df = pd.DataFrame(results)
        out_path = f"extracted_rates_{breeder}_{isotope}.csv"
        df.to_csv(out_path, index=False)
        print(f"\nSaved {len(results)} records to {os.path.abspath(out_path)}")
    else:
        print("No data extracted.")

if __name__ == "__main__":
    for b in ['hcpb', 'dcll']:
        for isotope in ['U238', 'Th232']:
            extract_rates(b, isotope)