#!/bin/bash
#SBATCH -J dt_fusion_array            # Job name
#SBATCH --account=paceship-tritium    # Charge account
#SBATCH --array=0-119                 # 120 jobs: FLiBe(0-39) + DCLL(40-79) + HCPB(80-119); each block 2 isotopes × 20 loadings
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
# 0-39 FLiBe, 40-79 DCLL, 80-119 HCPB; within each block: 2 isotopes × 20 loadings

BLANKETS=("FLiBe" "DCLL" "HCPB")
ISOTOPES=("U238" "Th232")
FERTILE_KGM3=(0 0.1 0.5 1.5 3 7.5 15 30 45 60 75 90 105 120 135 150 250 500 750 999.99)

BLANKET_IDX=$((SLURM_ARRAY_TASK_ID / 40))
SUB_IDX=$((SLURM_ARRAY_TASK_ID % 40))
ISOTOPE_IDX=$((SUB_IDX / 20))
LOADING_IDX=$((SUB_IDX % 20))

BLANKET=${BLANKETS[$BLANKET_IDX]}
ISOTOPE=${ISOTOPES[$ISOTOPE_IDX]}
LOADING=${FERTILE_KGM3[$LOADING_IDX]}

echo "Array job $SLURM_ARRAY_TASK_ID: Running $BLANKET with $ISOTOPE at ${LOADING} kg/m³"

# Export variables for Python script
export BLANKET ISOTOPE LOADING

# Run single configuration
python << 'PYTHON_SCRIPT'
import sys
sys.path.insert(0, '.')

from Python.reactor import *
from Python.parameters import *
from Python.utilities import *
import main  # pulls in FLiBe, DCLL, HCPB via main's imports
from main import build_reactor
import os

BLANKET = os.environ['BLANKET']
ISOTOPE = os.environ['ISOTOPE']
LOADING = float(os.environ['LOADING'])

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

echo "Array job $SLURM_ARRAY_TASK_ID complete: $BLANKET $ISOTOPE ${LOADING} kg/m³"