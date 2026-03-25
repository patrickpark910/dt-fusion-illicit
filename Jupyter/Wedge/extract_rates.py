import os
import openmc
import pandas as pd

def extract_rates(breeder, isotope):
    """
    Extracts the TBR (Total (n,t) rxn rate) and (n,gamma) reaction rate 
    from OpenMC statepoint files for wedge configurations.
    
    This function scans the OpenMC output directories for specified isotope loadings, 
    dynamically determines the most recent statepoint file in each matching folder,
    extracts the specified reaction rates using OpenMC pandas dataframes, 
    and saves the collated results into an isotope-specific CSV file.
    """
    base_dir = "./OpenMC"
    results = []
    breeder = breeder.lower()
    
    # Check if directory exists
    if not os.path.exists(base_dir):
        print(f"Directory {base_dir} not found. Run from Jupyter/Wedge")
        return

    # Find matching folders that include both the breeder and the isotope
    folders = [f for f in os.listdir(base_dir) if f.startswith((f"{breeder}_wedgeA", f"{breeder}_wedgeC")) and isotope in f]
    
    for folder in sorted(folders):
        path = os.path.join(base_dir, folder)
        if not os.path.isdir(path):
            continue
  
        # Parse case and loading from folder name
        # Format: dcll_wedgeC_U238_150.00kgm3
        parts = folder.split('_')
        if len(parts) < 3:
            continue
        case = parts[1][-1] # gets 'A' or 'C'
        folder_iso = parts[2]
        
        # Only process folders matching the current isotope argument
        if folder_iso != isotope:
            continue
            
        loading_str = parts[3].replace('kgm3', '')
        try:
            loading = float(loading_str)
        except ValueError:
            continue
            
        # Look for statepoint files
        sp_files = [f for f in os.listdir(path) if f.startswith('statepoint.') and f.endswith('.h5')]
        if not sp_files:
            print(f"No statepoint found in {folder}")
            continue
        # Get the latest statepoint by batch number
        sp_path = os.path.join(path, max(sp_files, key=lambda x: int(x.split('.')[1])))
            
        try:
            sp = openmc.StatePoint(sp_path)
            
            # Extract TBR
            li_tally = sp.get_tally(name="Total (n,t) rxn rate")
            li_df = li_tally.get_pandas_dataframe()
            # print(li_df)
            tbr = float(li_df["mean"].sum())
            
            # Extract (n,gamma) for the specific isotope
            fertile_tally = sp.get_tally(name="Total fertile rxn rate")
            f_df = fertile_tally.get_pandas_dataframe()
            # print(f_df)
            ng_mask = (f_df["nuclide"] == isotope) & (f_df["score"] == "(n,gamma)")
            ng_rate = float(f_df.loc[ng_mask, "mean"].sum())
            
            results.append({
                "folder": folder,
                "case": case,
                "loading_kg_m3": loading,
                "tbr": tbr,
                f"{isotope.lower()}_n_gamma": ng_rate
            })
            
            print(f"Extracted {folder}")
            
        except Exception as e:
            print(f"Error processing {folder}: {e}")
            
    # Save to CSV
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
