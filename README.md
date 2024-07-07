# MLTasksAssembly

This repository contains the script for obtaining the number of samples associated with each machine learning task that includes various input combinations of cancer type, omics feature type from the TCGA dataset, clinical outcome endpoint, and event time threshold (in years).

## Files

### MLTasksList.py

This file contains the main code for obtaining the number of samples associated with each machine learning task. It is intended to be run in an environment with the necessary dependencies installed.

The `MLTasksList.py` script requires the following input arguments:

1. **count_category**: Specifies the type of count to be performed. Possible values are `GenderRace`, `Race`, `Gender`.
2. **omicsConfiguration**: Specifies the configuration of omics data to be used. Possible values are `single` or `combination`.
3. **num_features**: Specifies the number of features to be combined. This is required if `omicsConfiguration` is `combination`. Possible values are integers such as `2`, `3`, `4`.
4. **min_cases**: Specifies the minimum number of cases. Normally, it is `5`.
5. **DDP_group**: Specifies the data disadvantaged population (DDP) as the target ethnic group. Possible values are `BLACK`, `ASIAN`, `NAT_A` for the TCGA dataset. This is required if `count_category` is `GenderRace` or `Race`.

### Example Commands

    python MLTasksList.py --count_category Gender --omicsConfiguration single --min_cases 5
    
    python MLTasksList.py --count_category Race --omicsConfiguration single --min_cases 5 --DDP_group BLACK
    
    python MLTasksList.py --count_category GenderRace --omicsConfiguration single --min_cases 5 --DDP_group BLACK
    
    python MLTasksList.py --count_category Gender --omicsConfiguration combination --num_features 2 --min_cases 5

It will run the `MLTasksList.py` script with the specified input arguments.


## Environment Setup

Create a virtual environment named `multiethnic`. Ensure this environment is set up with the necessary dependencies before submitting the job.

    **Python Version**: Python 3.6+

**Required Python Packages**:

`scipy`

`pandas`

`numpy`

`argparse`

`openpyxl` (for reading/writing Excel files)


## Acknowledgement

This work has been supported by NIH R01 grant.


## Contact

For any queries, please contact:

Prof. Yan Cui (ycui2@uthsc.edu)

Dr. Teena Sharma (tee.shar6@gmail.com)
