#!/bin/bash
#SBATCH -J dt_fusion_volume_7          # Job name
#SBATCH --account=paceship-tritium     # Charge account
#SBATCH -N 1                           # Number of nodes
#SBATCH --ntasks-per-node=24           # Number of cores per node
#SBATCH -p cpu-small                   # Memory per core
#SBATCH -t 99:99:00                    # Duration (30 min)
#SBATCH -q inferno                     # Queue (embers = free)
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ezoccoli3@gatech.edu

# Always run from the folder you submitted in
cd $SLURM_SUBMIT_DIR

# Load and activate OpenMC environment
module load anaconda3
conda activate openmc-env

cd /storage/scratch1/7/ezoccoli3/dt-fusion-illicit

# Run OpenMC via your Python driver
python main.py -r tallies
#python plot.py 
echo "Job complete. Outputs are in $SLURM_SUBMIT_DIR"
