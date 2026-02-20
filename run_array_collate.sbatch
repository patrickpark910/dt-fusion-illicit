#!/bin/bash
#SBATCH -J dt_fusion_collate          # Job name
#SBATCH --account=paceship-tritium   # Charge account
#SBATCH -N 1                          # One node
#SBATCH --ntasks-per-node=1           # Single task for collation
#SBATCH -p cpu-small                 # Partition
#SBATCH --mem=16G                     # Memory allocation (less needed)
#SBATCH -t 05:00:00                   # Duration (1 hour should be plenty)
#SBATCH -q inferno                    # Queue (inferno = resource-guaranteed, uses compute credits)
#SBATCH --mail-type=BEGIN,END,FAIL    # Email notifications
#SBATCH --mail-user=ezoccoli3@gatech.edu
#SBATCH -o report-collate-%j.out      # Combined stdout/stderr
#SBATCH --dependency=afterok:JOBID    # Replace JOBID with your array job ID

# Always run from the folder you submitted in
cd $SLURM_SUBMIT_DIR

# Load and activate OpenMC environment
module load anaconda3
conda activate openmc-env

# Move to your project directory (scratch)
cd /storage/scratch1/7/ezoccoli3/dt-fusion-illicit

# Collate all tallies after array jobs complete (FLiBe, DCLL, HCPB)
echo "Collating tallies for FLiBe, DCLL, HCPB..."

python << 'PYTHON_SCRIPT'
import sys
sys.path.insert(0, '.')

from Python.parameters import *
from Python.utilities import *
import main
from main import build_reactor

for blanket in ['FLiBe', 'DCLL', 'HCPB']:
    for fertile_isotope in ['U238', 'Th232']:
        print(f"Collating tallies for {blanket} {fertile_isotope}...")
        sample_run = build_reactor(blanket,
                                  blanket_name=blanket,
                                  fertile_isotope=fertile_isotope,
                                  fertile_kgm3=0.0,
                                  run_type='tallies',
                                  run_openmc=False)
        sample_run.initialize()
        vol_m3 = sample_run.blanket_volume
        main.collate_tallies(blanket, fertile_isotope, sample_run.breeder_enrich, sample_run.temp_k, vol_m3)

# Run rob.py for analysis
print("Running rob.py for analysis...")
import rob
rob.run_all()

# Redo all plots (plot.py) for all results in Figures/Data
# Note: plot.py expects FLiBe, DCLL, HCPB CSVs; after FLiBe-only array run only FLiBe files exist.
print("Running plot.py for all results...")
try:
    import plot
    plot.plot_all()
    print("plot.py completed.")
except Exception as e:
    print(f"plot.py skipped or failed (e.g. missing DCLL/HCPB data): {e}")

print("Collation and plots complete!")
PYTHON_SCRIPT

echo "Collation job complete. Outputs are in $SLURM_SUBMIT_DIR"