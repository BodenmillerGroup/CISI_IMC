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
import sys

# Find code directory relative to our directory
THIS_DIR = os.path.dirname('__file__')
CODE_DIR = os.path.abspath(os.path.join(THIS_DIR, '..', 'code'))
# Add code directory to systems paths
sys.path.append(CODE_DIR)

# Import CISI training fnc.
from analyze_dictionary_and_compositions import analyze_U_and_A


## Specify input paths
data_path = Path('/mnt/bb_dqbm_volume')
tonsil_path = Path(os.path.join(data_path,
                                'data/Tonsil_th152/preprocessed_data/spe.h5ad'))
####
## TODO: Change!!!
U_path = Path(os.path.join(data_path,
                                'analysis/Tonsil_th152/training/subset/....'))
A_path = Path(os.path.join(data_path,
                                'analysis/Tonsil_th152/training/subset/....'))

## Specify output path
out_path = Path(os.path.join(data_path, 'analysis/Tonsil_th152/training/subset...'))
# Create output directory if it doesn't exist
out_path.mkdir(parents=True, exist_ok=True)
####

# Check that input files/dictionary exist
if not helpers.analysis_utils.is_valid_file(tonsil_path, ['.h5ad']):
    # If file is not found, throw error
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            tonsil_path)
if not helpers.analysis_utils.is_valid_file(U_path, ['.csv']):
    # If file is not found, throw error
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            U_path)
if not helpers.analysis_utils.is_valid_file(A_path, ['.txt']):
    # If file is not found, throw error
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            A_path)
    
## Specify parameters
### TODO: change!! 
dictionary_size = 10

# Define test rois
test_names_tonsil = ('20220520_TsH_th152_cisi1_002',)

# Specify k-sparsity of dictionary used in training and name of final folder
# for each dataset, where results will be saved
normalization = 'none'


## Read in data
# Read in SingleCellExperiment converted to anndata by cellconverter in R
sce_tonsil = ad.read_h5ad(tonsil_path)
# Channels of interest
channels_of_interest = ['CD8a', 'CD20', 'CD4', 'CD11c', 'CD3', 'Ki-67', 'SMA', 'panCK']
# Remove uninteresting proteins/channels
sce_tonsil = sce_tonsil[:, sce_tonsil.var.index.isin(channels_of_interest)]

# Read U
U_tonsil = np.genfromtxt(U_path, delimiter=',', skip_header=True,
                         usecols=list(range(1, (dictionary_size)+1)))
U_names = np.genfromtxt(U_path, delimiter=',', usecols=0, skip_header=1,
                        dtype='S20').tolist()

# Read A
A_tonsil = np.loadtxt(A_path, skiprows=1, usecols=list(range(2, len(channels_of_interest)+2))
A_names = np.loadtxt(A_path, max_rows=1, dtype='S20').tolist()
  
                      
## Throw error if A and U don't have the same proteins or they are in a different order
if not U_names!=A_names:
    # If file is not found, throw error
    raise ValueError(('A and U do not have the same proteins or they are in a ' +
                      'different order. Please check both files, to prevent CISI ' +
                      'from computing wrong results.\n' +
                      'Proteins in A: {0}\n'.format(A_names) +
                      'Proteins in U: {0}'.format(U_names))
                      

## Predict performance of A given A and SCE
# Test spillover corrected results not spillover corrected data for simulation and analysis
predicted_res, predicted_res_noisy = analyze_U_and_A(sce_tonsil[sce_tonsil.obs.index.isin(test_names_tonsil), ],
                                                               U_tonsil, [A_tonsil], ['none'], out_path,
                                                               None, norm=normalization)

