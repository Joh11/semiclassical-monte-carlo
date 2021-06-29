#!/bin/bash -l
#SBATCH --job-name="zoomed 8x8"
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cscsci
#SBATCH --constraint=mc
#SBATCH --account=mr27
#SBATCH --cpus-per-task=36
#SBATCH --hint=nomultithread
#SBATCH --output=skl_dyn_zoomed.out
#SBATCH --error=skl_dyn_zoomed.err

module load daint-mc
module load Julia
module load JuliaExtensions

export JULIA_NUM_THREADS=8

echo running julia ...
srun julia --project=.. -O3 --check-bounds=no $1 # skl_dyn_factor.jl
