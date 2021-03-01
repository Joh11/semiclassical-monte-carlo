#!/bin/bash -l
#SBATCH --job-name="scmc"
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --account s1008

module load daint-gpu
module load Julia
module load JuliaExtensions

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun julia -O3 --check-bounds=no scmc.jl
