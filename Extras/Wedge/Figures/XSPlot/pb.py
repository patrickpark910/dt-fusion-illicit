import numpy as np

def main():
    # Define the input files and their corresponding isotopic abundances
    isotopes = {
        'Pb204n2n.txt': 0.014,  # 1.4%
        'Pb206n2n.txt': 0.241,  # 24.1%
        'Pb207n2n.txt': 0.221,  # 22.1%
        'Pb208n2n.txt': 0.524   # 52.4%
    }

    data = {}
    all_energies = []

    # 1. Read the data from each file
    for filename in isotopes.keys():
        try:
            # skiprows=2 ignores the title and header rows
            en, xs = np.loadtxt(filename, skiprows=2, unpack=True)
            data[filename] = (en, xs)
            all_energies.extend(en)
        except FileNotFoundError:
            print(f"Error: {filename} not found. Please ensure it is in the same directory.")
            return

    # 2. Create a common, sorted energy grid containing all unique energy points
    common_energy = np.unique(np.sort(all_energies))
    
    # 3. Calculate the weighted elemental cross section
    elemental_xs = np.zeros_like(common_energy)

    for filename, weight in isotopes.items():
        en, xs = data[filename]
        
        # Interpolate the isotopic XS onto the common grid.
        # left=0.0: The (n,2n) reaction has a threshold energy; below the minimum data point, XS is 0.
        # right=xs[-1]: If the common grid goes slightly past an isotope's max energy, hold the last value.
        interp_xs = np.interp(common_energy, en, xs, left=0.0, right=xs[-1])
        
        # Add the weighted contribution
        elemental_xs += weight * interp_xs

    # 4. Write the results to the output file
    output_filename = 'Pbn2n.txt'
    with open(output_filename, 'w') as f:
        f.write("Natural Pb(n,2n) ENDFB-8.0 (Weighted Average)\n")
        f.write("Energy(eV) XS(b)\n")
        for e, x in zip(common_energy, elemental_xs):
            f.write(f"{e:g}\t{x:g}\n")

    print(f"Successfully generated {output_filename} with {len(common_energy)} energy points.")

if __name__ == "__main__":
    main()