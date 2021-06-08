#!/bin/bash -l
#SBATCH --job-name="scmc"
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --account=mr27
#SBATCH --cpus-per-task=36
#SBATCH --hint=nomultithread
#SBATCH --output=skl_dyn.out
#SBATCH --error=skl_dyn.err

module load daint-mc
module load Julia
module load JuliaExtensions

export JULIA_NUM_THREADS=8

echo running julia ...
srun julia --project=.. -O3 --check-bounds=no skl_dyn_factor.jl
