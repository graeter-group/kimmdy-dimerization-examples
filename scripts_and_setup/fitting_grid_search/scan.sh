#!/bin/bash
#SBATCH --job-name=grid_test
#SBATCH --output=logs/test_job_%A_%a.out
#SBATCH --error=logs/test_job_%A_%a.err
#SBATCH --array=0-99
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

source /hits/fast/mbm/schuepbs/miniconda3/etc/profile.d/conda.bash
conda activate ./myenv

# Parameters
total_k1=1000
total_k2=1000
total_jobs=100
combinations_per_job=$(( (total_k1 * total_k2) / total_jobs ))

task_id=${SLURM_ARRAY_TASK_ID}
start_idx=$(( task_id * combinations_per_job ))
end_idx=$(( (task_id + 1) * combinations_per_job - 1 ))

start_k1_idx=$(( start_idx / total_k2 ))
end_k1_idx=$(( end_idx / total_k2 ))

spacing_k1=$(( end_k1_idx - start_k1_idx + 1 ))
spacing_k2=$total_k2

# Map indices to actual values
k1_min=0.0
k1_max=5.0
k2_min=0.0
k2_max=1.0

step_k1=$(python3 -c "print(($k1_max - $k1_min)/($total_k1 - 1))")
start_k1=$(python3 -c "print($k1_min + $start_k1_idx * $step_k1)")
end_k1=$(python3 -c "print($k1_min + $end_k1_idx * $step_k1)")

start_k2=$k2_min
end_k2=$k2_max

# Run your Python script
echo "Starting python script"
# Change XXX to system TdT or dT20
python scan_standalone_XXX.py "$start_k1" "$end_k1" "$spacing_k1" "$start_k2" "$end_k2" "$spacing_k2" "$task_id"
