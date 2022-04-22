#!/bin/bash
#SBATCH --partition=xgpu

# Set the scratch directory base
SCRATCH_DIR_BASE=$HOME/build_wacomm

# Here the full path to the WaComM++ executable
WACOMMPLUSPLUS_EXEC="$HOME/wacommplusplus/build/wacommplusplus"

# The path to the config file
CONFIG_FILE="$HOME/wacommplusplus/build/wacomm.json"

# The path to the sources file
SOURCE_FILE="$HOME/wacommplusplus/build/sources.json"

# The path to the processed input file directory
PROCESSED_DIR="$HOME/wacommplusplus/build/processed"

# The path to the restart file
#RESTART_FILE="$HOME/wacommplusplus/build/WACOMM_rst_20210701Z00.txt"

# Load the CUDA module
module load cuda/10.1

# Load the MVAPICH2 module
#module load mvapich2-2.3.5-gcc-8.3.1
module load ompi-4.1.0-gcc-8.3.1


# For PurpleJeans ONLY!
# Set Stack size:

ulimit -s 10240


# If the thread number is not specified, set it to 1
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi

CODE_NAME="N"${SLURM_NNODES}_"P"${SLURM_TASKS_PER_NODE}_"T"${SLURM_CPUS_PER_TASK}_"G0"
echo Creating $CODE_NAME

# Set the path for the scratch directory
SCRATCH_DIR=$SCRATCH_DIR_BASE/"GENERAL"

# Change the current directory to the scratch
cd $SCRATCH_DIR

# Link the executable
ln -sf $WACOMMPLUSPLUS_EXEC .

# Link the configuration file
ln -sf $CONFIG_FILE wacomm.json

# Link the configuration file
ln -sf $SOURCE_FILE sources.json

# Export the number of threads
export OMP_NUM_THREADS=$omp_threads
echo THREADS = $OMP_NUM_THREADS

mpirun --bind-to none --mca orte_base_help_aggregate 0 ./wacommplusplus
#mpirun --bind-to none ./wacommplusplus
