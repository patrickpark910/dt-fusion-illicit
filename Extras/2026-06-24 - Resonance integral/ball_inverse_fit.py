import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define the inverted reciprocal function
def production_rate(rho, c1, c2, c3):
    return rho / (c1 * rho**2 + c2 * rho + c3)

# List of files and corresponding labels for the plot
datasets = [
    {'file': 'FLiBe_900K_Li07.5_U238_summary.csv', 'label': 'FLiBe (Li6 7.5%)', 'color': '#66b420'},
    {'file': 'HCPB_900K_Li60.0_U238_summary.csv',  'label': 'HCPB (Li6 60.0%)', 'color': '#b41f24'},
    {'file': 'DCLL_900K_Li90.0_U238_summary.csv',  'label': 'DCLL (Li6 90.0%)', 'color': '#0047ba'}
]

# Constants for extracting A, B, and C
M_SQ = 8.0 
V = 1.0  

plt.figure(figsize=(7, 4.5))

for data in datasets:
    # 1. Load the data
    df = pd.read_csv(data['file'])
    
    # Extract x and y
    rho_fert = df['fertile_kg/m3'].values
    R_fiss = df['Pu239_kg/yr'].values
    
    # 2. Fit the data
    # Adjust p0 if the algorithm struggles to converge for a specific dataset
    try:
        popt, pcov = curve_fit(production_rate, rho_fert, R_fiss, p0=[1e-6, 1e-3, 1.0])
        c1, c2, c3 = popt
        
        # 3. Print Results
        print(f"--- FITTING RESULTS FOR {data['label']} ---")
        print(f"c1 = {c1:.4e} | c2 = {c2:.4e} | c3 = {c3:.4e}")
        
        # Extract original ARC parameters
        B = -(c1 * M_SQ) / V
        C = c2 * M_SQ
        A = c3 * M_SQ * V
        print(f"A = {A:.4e} | B = {B:.4e} | C = {C:.4e}\n")
        
        # 4. Plotting
        # Scatter plot for raw OpenMC data
        plt.scatter(rho_fert, R_fiss, color=data['color'], marker='o', alpha=0.6, label=f"{data['label']} Data")
        
        # Smooth line for the fitted curve
        rho_smooth = np.linspace(min(rho_fert), max(rho_fert), 200)
        plt.plot(rho_smooth, production_rate(rho_smooth, *popt), 
                 color=data['color'], linestyle='-', linewidth=2, label=f"{data['label']} Fit")
                 
    except Exception as e:
        print(f"Could not fit data for {data['label']}. Error: {e}\n")

# Finalize the plot formatting
plt.xlabel('Fertile Isotope Density, $\\rho_{fert}$ (kg/m$^3$)')
plt.ylabel('Pu-239 Production Rate, $R_{fiss}$ (kg/yr)')
plt.title('Fissile Production Rate vs. Fertile Density (Blanket Comparison)')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('fit.pdf')
plt.savefig('fit.png')
plt.show()