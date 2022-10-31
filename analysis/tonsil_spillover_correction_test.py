# Import libraries
import anndata as ad
from pathlib import Path
import errno
import os
import numpy as np
import sys


# Helper fncs
import helpers.analysis_utils


## CISI
# Configure code directory as module

# Find code directory relative to our directory
THIS_DIR = os.path.dirname('__file__')
CODE_DIR = os.path.abspath(os.path.join(THIS_DIR, '..', 'code'))
# Add code directory to systems paths
sys.path.append(CODE_DIR)

# Import CISI training fnc.
from train_dictionary_and_compositions import train_U_and_A
from analyze_dictionary_and_compositions import analyze_U_and_A


# Specify input paths
data_path = Path('/mnt/bb_dqbm_volume')
tonsil_path = Path(os.path.join(data_path,
                                'data/Tonsil_th152/preprocessed_data/spe.h5ad'))

# Specify output path
out_path = Path(os.path.join(data_path, 'analysis/Tonsil_th152/training/full'))
# Create output directory if it doesn't exist
out_path.mkdir(parents=True, exist_ok=True)


# Check that input files/dictionary exist
if not helpers.analysis_utils.is_valid_file(tonsil_path, ['.h5ad']):
    # If file is not found, throw error
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            tonsil_path)


# Read in SingleCellExperiment converted to anndata by cellconverter in R
sce_tonsil = ad.read_h5ad(tonsil_path)

# Remove uninteresting proteins/channels
sce_tonsil = sce_tonsil[:, ~sce_tonsil.var.index.str.contains('Histone|Ir[0-9]|E-Cad',
                                                              regex=True, case=False)]


# Define test rois
test_names_tonsil = ('20220520_TsH_th152_cisi1_002',)


# Specify k-sparsity of dictionary used in training and name of final folder
# for each dataset, where results will be saved
k = 1
folder_name_cor = 'sm_corrected/no_norm'
folder_name_uncor = 'sm_uncorrected/no_norm'
folder_name_combi = 'sm_combi/no_norm'
normalization = 'none'

# Train CISI
(training_res_cor, training_res_comp_cor,
 U_best_cor, Phi_best_cor,
 X_test_cor) = train_U_and_A(sce_tonsil,
                             os.path.join(out_path, folder_name_cor),
                             split_by='percentage', k_cv=4,
                             test_set=test_names_tonsil,
                             lda1=k, normalization=normalization,
                             layer='compcounts')
(training_res_uncor, training_res_comp_uncor,
 U_best_uncor, Phi_best_uncor,
 X_test_uncor) = train_U_and_A(sce_tonsil,
                               os.path.join(out_path, folder_name_uncor),
                               split_by='percentage', k_cv=4,
                               test_set=test_names_tonsil,
                               lda1=k, normalization=normalization)

# Test spillover corrected results not spillover corrected data for simulation and analysis
training_res_combi, training_res_noisy_combi = analyze_U_and_A(sce_tonsil[sce_tonsil.obs.index.isin(X_test_cor), ],
                                                               U_best_cor,
                                                               [Phi_best_cor], ['none'],
                                                               os.path.join(out_path, folder_name_combi),
                                                               None, norm=normalization)

