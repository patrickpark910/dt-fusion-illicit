#!/bin/bash
#SBATCH --job-name=u238_hcpb
#SBATCH --nodes=1                   
#SBATCH --ntasks=32                 # total number of tasks across all nodes
#SBATCH --cpus-per-task=1           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=4G                    # total memory
#SBATCH --time=48:00:00             # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin,end,fail  # receive email notifications
#SBATCH --mail-user=pp0191@princeton.edu

module purge
module load anaconda3/2025.6
conda activate openmc-env
python prism_hcpb.py