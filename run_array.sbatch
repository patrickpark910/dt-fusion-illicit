#!/bin/bash
#SBATCH -J dt_fusion_array            # Job name
#SBATCH --account=paceship-tritium    # Charge account
#SBATCH --array=0-39                  # 40 jobs: 0-39 (2 isotopes × 20 loadings)
#SBATCH -N 1                          # One node per job
#SBATCH --ntasks-per-node=1           # 1 task, 24 CPUs for OpenMP
#SBATCH --cpus-per-task=24            # Use all 24 cores for OpenMP threads
#SBATCH -p cpu-small                  # Partition
#SBATCH --mem=128G                    # Memory allocation
#SBATCH -t 10:00:00                   # Duration per job (8 hours; inferno allows longer)
#SBATCH -q inferno                    # Queue (inferno = resource-guaranteed, uses compute credits)
#SBATCH --mail-type=BEGIN,END,FAIL    # Email notifications
#SBATCH --mail-user=ezoccoli3@gatech.edu
#SBATCH -o report-array-%A-%a.out     # Combined stdout/stderr (%A=jobID, %a=array index)

# Always run from the folder you submitted in
cd $SLURM_SUBMIT_DIR

# Load and activate OpenMC environment
module load anaconda3
conda activate openmc-env

# Move to your project directory (scratch)
cd /storage/scratch1/7/ezoccoli3/dt-fusion-illicit

# Calculate which configuration to run based on array index
# Array index maps to: (isotope_index * 20) + loading_index
# Isotopes: 0=U238, 1=Th232
# Loadings: 0-19 (20 values)

ISOTOPES=("U238" "Th232")
FERTILE_KGM3=(0 0.1 0.5 1.5 3 7.5 15 30 45 60 75 90 105 120 135 150 250 500 750 999.99)

ISOTOPE_IDX=$((SLURM_ARRAY_TASK_ID / 20))
LOADING_IDX=$((SLURM_ARRAY_TASK_ID % 20))

ISOTOPE=${ISOTOPES[$ISOTOPE_IDX]}
LOADING=${FERTILE_KGM3[$LOADING_IDX]}

echo "Array job $SLURM_ARRAY_TASK_ID: Running FLiBe with $ISOTOPE at ${LOADING} kg/m³"

# Export variables for Python script
export ISOTOPE LOADING

# Run single configuration
python << 'PYTHON_SCRIPT'
import sys
sys.path.insert(0, '.')

from Python.reactor import *
from Python.flibe import *
from Python.parameters import *
from Python.utilities import *
from main import build_reactor  # ADD THIS LINE
import os

# Override to run only this specific configuration
ISOTOPE = os.environ['ISOTOPE']
LOADING = float(os.environ['LOADING'])
BLANKET = "FLiBe"

print(f"Running: {BLANKET} {ISOTOPE} {LOADING} kg/m³")

current_run = build_reactor(BLANKET,
                            blanket_name=BLANKET,
                            fertile_isotope=ISOTOPE,
                            fertile_kgm3=LOADING,
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

print(f"{C.GREEN}Completed: {BLANKET} {ISOTOPE} {LOADING} kg/m³{C.END}")
PYTHON_SCRIPT

echo "Array job $SLURM_ARRAY_TASK_ID complete: FLiBe $ISOTOPE ${LOADING} kg/m³"