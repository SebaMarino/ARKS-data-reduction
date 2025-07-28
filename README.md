# ARKS-data-reduction
Repository of reduction scripts used for ARKS to be run in CASA 6.4.1.12

Each directory is dedicated to a particular system and must contain the following subdirectories:

- calibrated: where the calibrated data should be placed in subsubdirectories with a name to identify the antenna configuration (e.g. calibrated/12mLB/arc/projects/ARKS/data/reduction/HD10647/calibrated/12mLB/uid___A002_X1171dca_X5df0.ms.split.cal). Note that you will need to create this directory and its subdirectories.
- reduced: it contains a script reduce_*target*.py that needs to be run to produce the reduced data. Note that this script calls functions in the file functions_data_reduction.py, which is in the main directory.
- corrected: it contains a script correct_*target*.py that generates the corrected data using as input the file MCMC_results.json that is in this subsubdirectory.
