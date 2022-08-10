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
lung_path = Path(os.path.join(data_path,
                              'data/Immucan_lung/Lung_sce.h5ad'))

# Specify output path
out_path = Path(os.path.join(data_path, 'analysis'))
# Create output directory if it doesn't exist
out_path.mkdir(parents=True, exist_ok=True)


# Check that input files/dictionary exist
if not helpers.analysis_utils.is_valid_file(tonsil_path, ['.h5ad']):
    # If file is not found, throw error
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            tonsil_path)
if not helpers.analysis_utils.is_valid_file(lung_path, ['.h5ad']):
    # If file is not found, throw error
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT),
                            lung_path)


# Read in SingleCellExperiment converted to anndata by cellconverter in R
sce_tonsil = ad.read_h5ad(tonsil_path)
sce_lung = ad.read_h5ad(lung_path)

# Remove uninteresting proteins/channels
sce_tonsil = sce_tonsil[:, ~sce_tonsil.var.index.str.contains('Histone|Ir[0-9]|E-Cad',
                                                              regex=True, case=False)]
sce_lung = sce_lung[:, ~sce_lung.var.index.str.contains('Histone|Ir[0-9]|E-Cad',
                                                        regex=True, case=False)]

# Get intersection of common proteins/channels and randomly select the same amount
# of cells in lung as in tonsil
intersection_proteins = (sce_tonsil.var.index).intersection(sce_lung.var.index)
sce_tonsil = sce_tonsil[:, intersection_proteins]
sce_lung = sce_lung[(np.random.choice(range(sce_lung.n_obs), sce_tonsil.n_obs,
                                      replace=False)), intersection_proteins]

# Define test rois
test_names_tonsil = ('20220520_TsH_th152_cisi1_002',)
# Get all roi names for indexing
roi_names_lung = sce_lung.obs['sample_id'].unique().to_list()
test_names_lung = tuple(np.random.choice(roi_names_lung,
                                           int(0.2*len(roi_names_lung)),
                                           replace=False))


# Specify k-sparsity of dictionary used in training and name of final folder
# for each dataset, where results will be saved
k = 1
folder_name = "k_1"

# Train CISI
(training_res_tonsil, training_res_comp_tonsil,
 U_best_tonsil, Phi_best_tonsil,
 X_test_tonsil) = train_U_and_A(sce_tonsil,
                                          os.path.join(out_path,
                                                       'Tonsil_th152/training/full', folder_name),
                                          split_by='percentage', k_cv=4,
                                          test_set=test_names_tonsil,
                                          lda1=k)
(training_res_lung, training_res_comp_lung,
 U_best_lung, Phi_best_lung,
 X_test_lung) = train_U_and_A(sce_lung,
                                        os.path.join(out_path,
                                                     'Immucan_lung/training/subset', folder_name),
                                        split_by='percentage', k_cv=4,
                                        test_set=test_names_lung,
                                        lda1=k)


# Test training results on the opposite dataset
training_res_lung, training_res_comp_lung = analyze_U_and_A(sce_lung[sce_lung.obs.index.isin(X_test_lung), ],
                                                            U_best_tonsil,
                                                            [Phi_best_tonsil], ["none"],
                                                            os.path.join(out_path,
                                                                         'tonsil_vs_lung/test_lung', folder_name),
                                                            "none")
training_res_tonsil, training_res_comp_tonsil = analyze_U_and_A(sce_tonsil[sce_tonsil.obs.index.isin(X_test_tonsil), ],
                                                                U_best_lung,
                                                                [Phi_best_lung], ["none"],
                                                                os.path.join(out_path,
                                                                             'tonsil_vs_lung/test_tonsil',
                                                                             folder_name), "none")
