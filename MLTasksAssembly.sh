#!/bin/bash
#SBATCH --job-name=MLTasksAssembly
#SBATCH --account=ACF-UTHSC0001
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --output=MLTasksAssembly.o%j
#SBATCH --error=MLTasksAssembly.e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=01-00:00:00

###########################################
cd $SLURM_SUBMIT_DIR
echo "Current directory: $(pwd)"
###########################################
eval "$(conda shell.bash hook)"

source activate multiethnic

echo "The environment has been activated."

python MLTasksAssembly.py gender_race combination 2 5 BLACK

echo "The execution has been done."


