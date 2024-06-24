# MLTasksAssembly

This repository contains scripts for running machine learning tasks in a high-performance computing environment. The primary components of this repository are:

- `MLTasksAssembly.py`: A Python script containing the machine learning tasks.
- `MLTasksAssembly.sh`: A Bash script for submitting the Python script as a job on an HPC cluster using SLURM.

## Files

### MLTasksAssembly.py

This file contains the main code for performing various machine learning tasks. It is intended to be run in an environment with the necessary dependencies installed.

The `MLTasksAssembly.py` script requires the following input arguments:

1. **which_count**: Specifies the type of count to be performed. Possible values are `gender_race`, `ancestry`, `gender`.
2. **count_type**: Specifies the type of count. Possible values are `single` or `combination`.
3. **num_features**: Specifies the number of features to be combined. Possible values are integers such as `2`, `3`, `4`.
4. **min_cases**: Specifies the minimum number of cases. Normally, it is `5`.
5. **DDP_group**: Specifies the data disadvantaged population (DDP) as target ethnic group. Possible values are `BLACK`, `ASIAN`, `NAT_A` for the TCGA dataset.

### Example Command

    ```sh
    python MLTasksAssembly.py gender_race combination 2 5 BLACK
    ```

This command will run the `MLTasksAssembly.py` script with the specified input arguments.

### MLTasksAssembly.sh

This file is a SLURM job submission script. It sets up the job parameters, activates the required Python environment, and runs the `MLTasksAssembly.py` script.

## Environment Setup

The SLURM script activates a Conda environment named `multiethnic`. Ensure this environment is set up with the necessary dependencies before submitting the job.

    ```sh
    conda create --name multiethnic
    conda activate multiethnic
    # Install necessary packages
    conda install numpy pandas
    pip install scipy
    ```

## Acknowledgement

This work has been supported by NIH R01 grant.

## Contact

For any queries, please contact:

Prof. Yan Cui (ycui2@uthsc.edu)

Dr. Teena Sharma (tee.shar6@gmail.com)
