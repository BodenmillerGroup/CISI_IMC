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
CODE_DIR = os.path.abspath(os.path.join(THIS_DIR, '../..', 'code'))
# Add code directory to systems paths
sys.path.append(CODE_DIR)

# Import CISI training fnc.
from train_dictionary_and_compositions import train_U_and_A


# Specify input paths
data_path = Path('/mnt/bb_dqbm_volume')
tiffs_path = Path(os.path.join(data_path, 'data/Tonsil_th152/steinbock/img'))
panel_path = Path(os.path.join(data_path, 'data/Tonsil_th152/steinbock/panel.csv'))


# Specify output path
out_path = Path(os.path.join(data_path, 'analysis/Tonsil_th152/training_pixel/full'))
# Create output directory if it doesn't exist
out_path.mkdir(parents=True, exist_ok=True)


# Check that input files/dictionary exist
if not helpers.analysis_utils.is_valid_file(panel_path, ['.csv']):
    # If file is not found, throw error
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            panel_path)
if not helpers.analysis_utils.is_valid_directory(tiffs_path):
    # If directory is not found or doesn't contain TIFF files, throw an error
    raise Exception('Input TIFFS path {0} directory does not exist or does\
    not contain any valid TIFF files.'.format(tiffs_path)) from FileNotFoundError

    
# Read in images and convert to anndata 
images = helpers.analysis_utils.anndata_from_tiff(tiffs_path, panel_path)
# Remove uninteresting proteins/channels
images = images[:, ~images.var.index.str.contains('Histone|Ir[0-9]|E-Cad', regex=True, case=False)]


# Define test rois
test_names = ('20220520_TsH_th152_cisi1_002', '20220520_TsH_th152_cisi1_001', '20220520_TsH_th152_cisi1_003',
             '20220520_TsH_th152_cisi1_005',)


# Specify k-sparsity of dictionary used in training and name of final folder
# for each dataset, where results will be saved
k = 3
folder_name = 'k_3'
normalization = 'none'
num_blocks=50

# Train CISI
(training_res, training_res_comp,
 U_best, Phi_best,
 X_test) = train_U_and_A(images,
                         os.path.join(out_path, folder_name),
                         split_by='percentage', k_cv=4,
                         test_set=test_names,
                         lda1=k, normalization=normalization,
                         num_blocks=num_blocks)

