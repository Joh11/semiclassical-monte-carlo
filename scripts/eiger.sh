#!/bin/bash -l
#SBATCH --job-name="UUD finite size scaling"
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --account=mr27
#SBATCH --cpus-per-task=36
#SBATCH --hint=nomultithread
##SBATCH --output=uud_qsl_scan.out
##SBATCH --error=uud_qsl_scan.err

module load cpeGNU
module load Julia
module load JuliaExtensions

export JULIA_NUM_THREADS=16

echo running julia ...
srun julia --project=.. -O3 --check-bounds=no $@
