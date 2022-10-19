## Import libraries
import numpy as np
import sys

## Import/Specify path to CISI code
# Specify CISI code directory path
CODE_DIR = snakemake.params['cisi_path']
# Add code directory to systems paths
sys.path.append(CODE_DIR)

# Import CISI training fnc.
from train_dictionary_and_compositions import train_U_and_A


## Read in snakemake parameters
k = int(snakemake.wildcards['k'])
d = int(snakemake.wildcards['d'])
nmeasurements = int(snakemake.wildcards['m'])

sce = snakemake.params['sce']
outpath = snakemake.params['out_path']                       
default_params = snakemake.params['default_params']


## Train CISI
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
                                           save=default_params['add_noise'],
                                           analysis_normalization=default_params['analysis_normalization'])
