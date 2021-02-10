#!/bin/bash
#SBATCH --job-name=wacommplusplus
#SBATCH --output=wacommplusplus.out
#SBATCH --error=wacommplusplus.err
#SBATCH --partition=xhicpu
#SBATCH --nodes=2
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=32

module load cuda/10.1
module load mvapich2-2.3.5-gcc-8.3.1

if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi

export OMP_NUM_THREADS=$omp_threads

mpirun ./wacommplusplus