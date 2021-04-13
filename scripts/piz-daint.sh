#!/bin/bash -l
#SBATCH --job-name="scmc"
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=debug
#SBATCH --constraint=mc
#SBATCH --account=mr27
#SBATCH --cpus-per-task=36
#SBATCH --hint=nomultithread
#SBATCH --output=fig4.out
#SBATCH --error=fig4.err

module load daint-mc
module load Julia
module load JuliaExtensions

export JULIA_NUM_THREADS=72
srun julia --project=.. -O3 --check-bounds=no fig4.jl
