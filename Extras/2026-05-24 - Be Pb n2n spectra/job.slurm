#!/bin/bash
#SBATCH --job-name=spectra             # create a short name for your job
#SBATCH --nodes=1                      # node count
#SBATCH --ntasks=1                     # total number of tasks across all nodes
#SBATCH --cpus-per-task=64             # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=4G                       # memory total
#SBATCH --time=24:00:00                # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin,end,fail     # receive email notifications
#SBATCH --mail-user=pp0191@princeton.edu

module purge
module load anaconda3/2025.6
conda activate openmc-env
python -u n2n.py
