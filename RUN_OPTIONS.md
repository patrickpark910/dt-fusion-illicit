# Running Options for FLiBe Tallies

This document explains the different ways to run your OpenMC simulations.

## Option 1: Single Job with Parallelization (Recommended for Full Run)

**File:** `run.sbatch`

- Uses all 24 cores per OpenMC run (OpenMP threads)
- Runs all configurations sequentially
- Queue: **inferno** (resource-guaranteed; uses compute credits)
- Time limit: 24 hours
- **Expected runtime:** ~20-30 hours with 24 cores (vs ~64 hours sequential)

**Usage:**
```bash
sbatch run.sbatch
```

**What it does:**
- Runs FLiBe tallies for U238 and Th232
- 20 different fertile loadings each (40 total runs)
- Each OpenMC run uses 24 threads for parallelization
- Runs `rob.py` after completion

---

## Option 2: Test Run with Reduced Loadings

**File:** `run_test.sbatch`

- Uses 24 cores per OpenMC run
- Only runs 9 loadings instead of 20: `[0, 0.5, 3, 15, 30, 60, 150, 500, 999.99]`
- Time limit: 2 hours
- **Expected runtime:** ~2-3 hours

**Usage:**
```bash
sbatch run_test.sbatch
```

**Use this to:**
- Test that everything works before running the full set
- Verify parallelization is working
- Check that outputs are correct

---

## Option 3: Job Array (Fastest - Runs in Parallel)

**Files:** `run_array.sbatch` and `run_array_collate.sbatch`

- Runs 40 jobs in parallel (one per configuration)
- Each job uses 24 cores
- Queue: **inferno** (resource-guaranteed; uses compute credits)
- Time limit per job: 8 hours
- **Expected runtime:** ~4-8 hours total (all jobs run simultaneously)

**Usage:**
```bash
# Step 1: Submit the array job (use --parsable for reliable job ID)
JOBID=$(sbatch --parsable run_array.sbatch)
echo "Array job ID: $JOBID"

# Step 2: Submit collation job that waits for array to finish
sbatch --dependency=afterok:$JOBID run_array_collate.sbatch
```

**Or manually:**
```bash
# Submit array
sbatch run_array.sbatch
# Note the job ID (e.g., 12345678)

# Edit run_array_collate.sbatch and replace JOBID with your array job ID
# Then submit:
sbatch run_array_collate.sbatch
```

**What it does:**
- Array job runs 40 independent jobs (one per configuration)
- Each job runs a single (blanket, isotope, loading) combination
- Collation job runs after all array jobs complete
- Collates results and runs `rob.py`

**Advantages:**
- Fastest option (all runs in parallel)
- If one configuration fails, others continue
- Can monitor progress per configuration

**Disadvantages:**
- Uses more cluster resources (40 nodes simultaneously)
- May have longer queue wait time

---

## Comparison

| Option | Runtime | Cores Used | Best For |
|--------|---------|------------|----------|
| Option 1 (Single) | ~20-30 hrs | 1 node × 24 cores | Full production run, moderate speed |
| Option 2 (Test) | ~2-3 hrs | 1 node × 24 cores | Testing before full run |
| Option 3 (Array) | ~4-6 hrs | 40 nodes × 24 cores | Fastest, if cluster resources available |

---

## Phoenix queues: inferno vs embers

Scripts use the **inferno** queue (resource-guaranteed):
- **inferno:** Uses compute credits ($68/month free tier ≈ 10,000 CPU-hours on a 192GB node); longer time limits, no preemption. Unused credits roll over (max 4 months).
- **embers:** Free, preemptable backfill; shorter time limits (e.g. 4 h), jobs can be preempted.

To use the free embers queue instead, change `-q inferno` to `-q embers` in the sbatch files and reduce `-t` to the embers limit (e.g. 4 hours).

---

## Notes

- All scripts use `--cpus-per-task=24` and `--ntasks-per-node=1` for OpenMP on 24 cores
- The `reactor.py` file uses `SLURM_CPUS_PER_TASK` for thread count
- Convert line endings on the cluster: `sed -i 's/\r$//' run*.sbatch`

---

## Troubleshooting

**If jobs are preempted:**
- You are likely on **embers**. Switch to **inferno** (edit `-q inferno` in the sbatch files) to use compute credits for non-preemptible runs.

**QOSMaxWallDurationPerJobLimit:**
- The requested time exceeds the queue limit. On **inferno** you can request longer (e.g. 8–24 h). On **embers**, use 4 h or whatever the QOS allows.

**If time limit exceeded:**
- Increase `-t` time limit
- Use Option 3 (array) to run faster
- Reduce number of loadings for testing

**If parallelization not working:**
- Check that OpenMC was compiled with OpenMP support
- Verify `SLURM_CPUS_PER_TASK` is set: `echo $SLURM_CPUS_PER_TASK`
- Check OpenMC output for "threads" or "OpenMP" messages
