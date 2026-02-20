# Clean run (FLiBe, DCLL, HCPB – 1e6 particles)

1. **On cluster, remove old OpenMC output** (so array jobs run fresh):
   ```bash
   cd /storage/scratch1/7/ezoccoli3/dt-fusion-illicit
   rm -rf OpenMC/tallies_* OpenMC/plot_* OpenMC/volume_*
   ```

2. **Pull latest** (if you pushed from elsewhere):
   ```bash
   git pull origin emmaflibe
   ```

3. **Submit array** (120 jobs: FLiBe + DCLL + HCPB, each 2 isotopes × 20 loadings; 1e6 particles × 10 batches):
   ```bash
   sbatch run_array.sbatch
   ```

4. **After array completes, submit collate** (replace JOBID with the array job ID from step 3):
   ```bash
   sbatch --dependency=afterok:JOBID run_array_collate.sbatch
   ```
   Collate rebuilds rxns CSVs for all three blankets, runs rob.py, then plot.py.

5. **Push results to git** (from cluster; no OpenMC – it’s gitignored):
   ```bash
   git add main.py Python/ plot.py rob.py Figures/Data/
   git add -f run.sbatch run_array.sbatch run_array_collate.sbatch
   git commit -m "FLiBe DCLL HCPB run 1e6 particles, rob and plot"
   git push origin emmaflibe
   ```
