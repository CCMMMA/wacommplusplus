#!/bin/bash
#SBATCH --job-name=wacommplusplus
#SBATCH --output=wacommplusplus.out
#SBATCH --error=wacommplusplus.err
#SBATCH --partition=xhicpu
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=16

# Set the scratch directory base
SCRATCH_DIR_BASE=$HOME/wacommplusplus_work

# Here the full path to the WaComM++ executable
WACOMMPLUSPLUS_EXEC=""

# The path to the config file
CONFIG_FILE=""

# The path to the restart file
RESTART_FILE=""

# Load the CUDA module
module load cuda/10.1

# Load the MVAPICH2 module
module load mvapich2-2.3.5-gcc-8.3.1

# If the thread number is not specified, set it to 1
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi

# Set the path for the scratch directory
SCRATCH_DIR=$SCRATCH_DIR_BASE/$SLURM_JOB_ID

# Create the scratch directory
mkdir -p $SCRATCH_DIR

# Change the current directory to the scratch
cd $SCRATCH_DIR

# Create the output and the restarts directories
mkdir output
mkdir restarts

# Copy here the processed input files directory
cp -r $PROCESSED_DIR .

# Copy the restart file
cp $RESTART_FILE WACOMM_rst_20201130Z00.txt

# Link the executable
ln -sf $WACOMMPLUSPLUS_EXEC .

# Link the configuration file
ln -sf $CONFIG_FILE wacomm.json

# Export the number of threads
export OMP_NUM_THREADS=$omp_threads

# Execute WaComM++
mpirun ./wacommplusplus