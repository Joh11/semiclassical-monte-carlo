#!/bin/bash -l
#SBATCH --job-name="scmc"
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --account=mr27
#SBATCH --cpus-per-task=72
#SBATCH --hint=nomultithread

module load daint-gpu
module load Julia
module load JuliaExtensions

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun julia -O3 --check-bounds=no fig4.jl
