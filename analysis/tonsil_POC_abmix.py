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
from analyze_dictionary_and_compositions import analyze_U_and_A
from decompress import decompress
from utils import compare_results


## Specify input paths
data_path = Path('/mnt/bb_dqbm_volume')
training_path = Path(os.path.join(data_path,
                                'data/Tonsil_th152/preprocessed_data/spe.h5ad'))
## TODO: Change!!!
experiment_path = Path(os.path.join(data_path,
                                '......'))

U_path = Path(os.path.join(data_path,
                                'analysis/.../experiment/proof_of_concept_fixed_A/gene_modules.csv'))
A_path = Path(os.path.join(data_path,
                                'analysis/.../experiment/proof_of_concept_fixed_A/version_1.txt'))

## Specify output path
out_path_training = Path(os.path.join(data_path, 'analysis/.../experiment/proof_of_concept_fixed_A/training'))
out_path_experiment = Path(os.path.join(data_path, 'analysis/.../experiment/proof_of_concept_fixed_A/experiment'))
out_path_experiment_simulation = Path(os.path.join(data_path, 'analysis/.../experiment/proof_of_concept_fixed_A/experiment-simulation'))
# Create output directory if it doesn't exist
out_path_training.mkdir(parents=True, exist_ok=True)
out_path_experiment.mkdir(parents=True, exist_ok=True)
out_path_experiment_simulation.mkdir(parents=True, exist_ok=True)


# Check that input files/dictionary exist
if not helpers.analysis_utils.is_valid_file(training_path, ['.h5ad']):
    # If file is not found, throw error
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            training_path)
if not helpers.analysis_utils.is_valid_file(experiment_path, ['.h5ad']):
    # If file is not found, throw error
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            experiment_path)    
if not helpers.analysis_utils.is_valid_file(U_path, ['.csv']):
    # If file is not found, throw error
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            U_path)
if not helpers.analysis_utils.is_valid_file(A_path, ['.txt']):
    # If file is not found, throw error
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            A_path)
    
## Specify parameters for CISI
dictionary_size = 20
normalization = 'none'

# Define test rois
test_names_training = ('20220520_TsH_th152_cisi1_002',)


## Read in data
# Read in SingleCellExperiment converted to anndata by cellconverter in R
sce_training = ad.read_h5ad(training_path)
sce_experiment = ad.read_h5ad(experiment_path)
# Channels of interest
proteins_of_interest = ['CD8a', 'CD20', 'CD4', 'CD11c', 'CD3', 'Ki-67', 'SMA', 'panCK']
channels_of_interest = ['channel 141', 'channel 145', 'channel 148', 'channel 152',
                       'channel 154', 'channel 165']
# Remove uninteresting proteins/channels
sce_training = sce_training[:, sce_training.var.index.isin(proteins_of_interest)]

# Read U
U = np.genfromtxt(U_path, delimiter=',', skip_header=True,
                         usecols=list(range(1, (dictionary_size)+1)))
U_names = [x.decode() for x in np.genfromtxt(U_path, delimiter=',', usecols=0, skip_header=1,
                                             dtype='S20')]

# Read A
A = np.loadtxt(A_path, skiprows=1, usecols=list(range(2, len(proteins_of_interest)+2)))
A_names = [x.decode() for x in np.loadtxt(A_path, max_rows=1, dtype='S20')]
A_channels = [x.decode() for x in np.loadtxt(A_path, usecols=1, dtype='S20')]
  
                      
## Throw error if A and U don't have the same proteins or A has different channels than the specified channels
## of interest
if (all(e in U_names for e in A_names) & all(e in A_names for e in U_names) &
   all(e in A_channels for e in channels_of_interest) & all(e in channels_of_interest for e in A_channels)):
    U_index = [U_names.index(ind) for ind in proteins_of_interest]
    A_index = [A_names.index(ind) for ind in proteins_of_interest]
    A_channel_index = [A_channels.index(ind) for ind in channels_of_interest]
    
    U = U[U_index, :]
    A = A[A_index, :]
    A = A[:, A_channel_index]
    
else:
    # Throw error
    raise ValueError(('A and U do not have the same proteins or A does ' +
                      'not have the channels of interest.' +
                      'Please check all files, to prevent CISI ' +
                      'from computing wrong results.\n' +
                      'Proteins in U: {0}'.format(U_names) +
                      'Proteins in A: {0}\n'.format(A_names) +
                      'Channels in A: {0}\n'.format(A_channels))
                      

## Predict performance of A given A,U and SCE
# Test for training data and real experiment data
(predicted_res_training, 
 predicted_res_noisy_training) = analyze_U_and_A(sce_training[sce_training.obs.index.isin(test_names_training), ],
                                                 U, [A], ['none'], out_path_training,
                                                 None, norm=normalization)
(predicted_res_training, 
 predicted_res_noisy_training) = analyze_U_and_A(sce_experiment[:, sce_experiment.var.index.isin(proteins_of_interest)],
                                                 U, [A], ['none'], out_path_experiment_simulation,
                                                 None, norm=normalization)

## Decompression
# Decompress composite channels                     
decompressed_x = decompress((sce_experiment[:, sce_experiment.var.index.isin(channels_of_interest)].X).T, 
                            U, [A])
# Remove infinit values                   
decompressed_x[np.isnan(decompressed_x)] = 0
# Compute statistics
decompression_results = compare_results((sce_experiment[:, sce_experiment.var.index.isin(proteins_of_interest)].X).T, 
                                        decompressed_x)
d_gene = np.percentile(1 - distance.pdist(A.dot(U).T, 'correlation'), 90)

# Open file to save decompression results                 
decompression_file = open(os.path.join(out_path_experiment, 'simulation_results.txt'), 'w')
colnames = ['version', 'Overall pearson', 'Overall spearman', 'Gene average',
            'Sample average', 'Sample dist pearson', 'Sample dist spearman',
            'Gene dist pearson', 'Gene dist spearman',
            'Matrix coherence (90th ptile)']
decompression_file.write('\t'.join(colnames) + '\n')                     
decompression_file.write('\t'.join([str(x) for x in ['']+decompression_results+[d_gene]]) + '\n')
decompression_file.close()
                     
# Write decomposed X to anndata
decompressed_anndata = sce_experiment[:, sce_experiment.var.index.isin(proteins_of_interest)].copy()
for k in list(decompressed_anndata.layers.keys()):
    del decompressed_anndata.layers[k]
    decompressed_anndata.X = decompressed_x.T
    decompressed_anndata.write(os.path.join(out_path_experiment, 'X_decomposed.h5ad'))

# Write original X subseted to individual protein expression levels to the same place as the decomposed X
(sce_experiment[:, sce_experiment.var.index.isin(proteins_of_interest)]).write(os.path.join(out_path_experiment, 'X_test.h5ad'))
                                    
                     
                     
                     
                     
                     
                     
