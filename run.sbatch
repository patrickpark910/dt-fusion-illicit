#!/bin/bash
#SBATCH -J dt_fusion_flibe             # Job name
#SBATCH --account=paceship-tritium     # Charge account
#SBATCH -N 1                           # One node
#SBATCH --ntasks-per-node=1            # One task (OpenMP uses all 24 CPUs)
#SBATCH --cpus-per-task=24             # Use all 24 cores for OpenMP threads
#SBATCH -p cpu-small                   # Partition (192GB node, 24 cores)
#SBATCH --mem=128G                     # Memory allocation
#SBATCH -t 24:00:00                    # Duration (24 hours; inferno allows longer than embers)
#SBATCH -q inferno                     # Queue (inferno = resource-guaranteed, uses compute credits)
#SBATCH --mail-type=BEGIN,END,FAIL     # Email notifications
#SBATCH --mail-user=ezoccoli3@gatech.edu
#SBATCH -o report-%j.out               # Combined stdout/stderr

# Always run from the folder you submitted in
cd $SLURM_SUBMIT_DIR

# Load and activate OpenMC environment
module load anaconda3
conda activate openmc-env

# Move to your project directory (scratch)
cd /storage/scratch1/7/ezoccoli3/dt-fusion-illicit

# Run FLiBe tallies only (1e7 particles, 10 batches)
python main.py -r tallies -b FLiBe

# Run analysis and generate plots
python rob.py

echo "Job complete. Outputs are in $SLURM_SUBMIT_DIR"