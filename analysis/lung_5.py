# Import libraries
import anndata as ad
from pathlib import Path
import errno
import os
import numpy as np


# Helper fncs
import helpers.analysis_utils


## CISI
# Import system libraries to configure code directory as module
from os.path import dirname, abspath, join
import sys

# Find code directory relative to our directory
THIS_DIR = dirname('__file__')
CODE_DIR = abspath(join(THIS_DIR, '..', 'code'))
# Add code directory to systems paths
sys.path.append(CODE_DIR)

# Import CISI training fnc.
from train_dictionary_and_compositions import train_U_and_A
from analyze_dictionary_and_compositions import analyze_U_and_A


# Specify input paths
data_path = Path('/mnt/bb_dqbm_volume')
lung_path = Path(os.path.join(data_path,
                                'data/Immucan_lung/Lung_sce_original.h5ad'))

# Specify output path
out_path = Path(os.path.join(data_path, 'analysis/Immucan_lung/training/subset'))
# Create output directory if it doesn't exist
out_path.mkdir(parents=True, exist_ok=True)


# Check that input files/dictionary exist
if not helpers.analysis_utils.is_valid_file(lung_path, ['.h5ad']):
    # If file is not found, throw error
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            lung_path)


# Read in SingleCellExperiment converted to anndata by cellconverter in R
sce_lung = ad.read_h5ad(lung_path)

# Remove uninteresting proteins/channels
sce_lung = sce_lung[:, ~sce_lung.var.index.str.contains('Histone|Ir[0-9]|DNA',
                                                              regex=True, case=False)]


# Define test rois
# Get all roi names for indexing
roi_names_lung = sce_lung.obs['sample_id'].unique().to_list()
test_names_lung = tuple(np.random.choice(roi_names_lung,
                                           int(0.95*len(roi_names_lung)),
                                           replace=False))


# Specify k-sparsity of dictionary used in training and name of final folder
# for each dataset, where results will be saved
k = 1
folder_name_norm = 'normalized'
folder_name_unnorm = 'unnormalized'

# Train CISI
(training_res_norm, training_res_comp_norm,
 U_best_norm, Phi_best_norm,
 X_test_norm) = train_U_and_A(sce_lung,
                             os.path.join(out_path, folder_name_norm),
                             split_by='percentage', k_cv=4,
                             test_set=test_names_lung,
                             lda1=k, normalization='paper_norm')
(training_res_unnorm, training_res_comp_unnorm,
 U_best_unnorm, Phi_best_unnorm,
 X_test_unnorm) = train_U_and_A(sce_lung,
                               os.path.join(out_path, folder_name_unnorm),
                               split_by='percentage', k_cv=4,
                               test_set=test_names_lung,
                               lda1=k, normalization='none')
