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
tonsil_path = Path(os.path.join(data_path,
                                'data/Tonsil_th152/preprocessed_data/spe.h5ad'))

# Specify output path
out_path = Path(os.path.join(data_path, 'analysis/Tonsil_th152/training/subset'))
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
immune_channels = ['CD15', 'CD20', 'CD3', 'CD38', 'CD4', 'CD68', 'CD8a', 
                   'ICOS', 'Ki-67', 'MPO', 'panCK', 'SMA', 'CD303', 'FOXP3', 
                   'GranzymeB']

# immune_channels = ["CD8a", "CD20", "CD4", "CD11c", "CD3", "Ki-67", "SMA", "panCK"]
sce_tonsil = sce_tonsil[:, sce_tonsil.var.index.isin(immune_channels)]


# Define test rois
test_names_tonsil = ('20220520_TsH_th152_cisi1_002',)


# Specify k-sparsity of dictionary used in training, name of final folder
# for each dataset, where results will be saved and if (test-)analysis is done
# on normalized data or not
k = 1
folder_name = 'pilot_immune_channels_normalized'
#folder_name = 'test_proof_of_concept'
normalization = 'none'
analysis_normalization = False

# Train CISI
(training_res, training_res_comp,
 U_best, Phi_best,
 X_test) = train_U_and_A(sce_tonsil,
                         os.path.join(out_path, folder_name),
                         split_by='percentage', k_cv=4,
                         test_set=test_names_tonsil,
                         lda1=k, normalization=normalization,
                         analysis_normalization=analysis_normalization)
'''
                         d=80,
                         maxItr=10,
                         nmeasurements=4,
                         maxcomposition=2)
'''