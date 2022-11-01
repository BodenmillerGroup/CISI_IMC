## Import libraries
import numpy as np
import sys
import os


## Import/Specify path to CISI code
# Specify CISI code directory path
CODE_DIR = snakemake.params['cisi_path']
# Add code directory to systems paths
sys.path.append(CODE_DIR)

# Import CISI training fnc.
from train_dictionary_and_compositions import train_U_and_A


## Output/error into log file
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f


    ## Read in snakemake parameters
    k = int(snakemake.wildcards['k'])
    d = int(snakemake.wildcards['d'])
    nmeasurements = int(snakemake.wildcards['m'])

    sce = snakemake.params['sce']
    outpath = snakemake.params['out_path']                       
    default_params = snakemake.params['default_params']


    ## Train CISI
    # Catch SMAF failing (U being empty) and create empty dummy file for snakemake not to fail
    try:
        (training_res, training_res_no_noise,
         U_best, Phi_best, X_test) = train_U_and_A(sce,
                                                   outpath,
                                                   split_by=default_params['split_by'],
                                                   k_cv=default_params['k_cv'],
                                                   test_set=tuple(default_params['test_names']),
                                                   lda1=k,
                                                   normalization=default_params['normalization'],
                                                   d=d,
                                                   maxItr=default_params['maxItr'],
                                                   nmeasurements=nmeasurements,
                                                   maxcomposition=default_params['maxcomposition'],
                                                   save=default_params['save'],
                                                   analysis_normalization=default_params['analysis_normalization'],
                                                   best_A_method=default_params['best_A_method'])
    except ValueError as ve:
        if 'dictionary' in str(ve):
            with open(os.path.join(outpath, 'no_noise_simulation_results.txt'), 'w') as f:
                pass
        else:
            raise ve.with_traceback(sys.exc_info()[2])
        


