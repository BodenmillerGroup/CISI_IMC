# Import libraries (from original script)
import numpy as np
from analyze_predictions import *

# Import libraries (additionaly added)
from utils import sparse_decode, compare_results, get_observations, get_observations_no_noise
from decompress import decompress
from pathlib import Path
import sys
import os
import pandas as pd


'''
Analyze performance of computed U and A
(Nine different metrics, correlations and distances gene- and sample-wise)

For given X, U and A, simulate composite results y and analyze performance:
inputs:
    X: anndata object containing numpy array X (cells/pixels x proteins)
       with proteins names as X.var_names
    U: a dictionary of gene modules (proteins x modules)
    Phi: list containing (best) composition matrices (composite channels (measurements) x proteins,
         binary)
    lasso_sparsity: error tolerance applied to l1 sparsity constraint (default: 0.2)
    THREADS: # Number of threads used (default=20)
    outpath: If specified, the best 50 U versions are saved as .txt files in
             compositions_A folder
    k: Which crossvalidation led to best results
    versions: list containing versions of phi's in Phi
    layer: which layer in anndata object to use (default: X)
    norm: which normalization used before simulating decomposition measurements
    save: Which decomposed X is saved.
          Either the X decomposed from noisy simulated data or simulated data
          without noise (default: no_noise)
    snr: Signal to noise ratio used to simulate noisy composite data
    correct_comeasured: Correct comeasured genes that are not coexpressed (default: False)
    train_corr: Correlations between genes in training data X (genes that are coexpressed)
                used in correct_comeasured (default: None)

outputs:
    results_df: All performance analyis in one pandas df
    simulation_results.txt: Performance analysis
'''

def analyze_U_and_A(X_input, U, Phi, versions, outpath, k, lasso_sparsity=0.2,
                    THREADS=20, layer=None, norm='none', save='no_noise', snr=5,
                    correct_comeasured=False, train_corr=None):
    # Select layer of anndata object that should be used and transpose it
    # to proteins x cells/channels
    if layer is not None:
        X_test = (X_input.layers[layer]).T
    else:
        X_test = (X_input.X).T

    # Extract data matrix X from anndata object and apply selected normalization
    match norm:
        case 'paper_norm':
            # Normalize X? (otherwise spams.nmf underneath will estimate U to only contain
            # 0s and smaf won't work)
            # Normalizaiton: Every row in X_input is divided by the corresponding vector
            # element in the rowwise norm (proteins are normalized accross all cells/pixels)
            X_test = (X_test.T / np.linalg.norm(X_test, axis=1)).T
        case 'min_max_norm':
            X_test = (X_test-X_test.min(axis=1, keepdims=True)) / (
                X_test.max(axis=1, keepdims=True)-X_test.min(axis=1, keepdims=True)
                )
        case _:
            # In case no valid normalization is given, an error is thrown
            raise ValueError(('The normalization {0} used by smaf is not valid.'.format(norm) +
                              'Please use one of the following: paper_norm, min_max_norm, ' +
                              'or none.'))

    # Empirical observation: using a sparsity constraint that is softer than
    # that used during training slightly improves results
    sparsity = lasso_sparsity / 10

    # Correlation between genes?
    d_gene = np.array([np.percentile(1 - distance.pdist(phi.dot(U).T, 'correlation'), 90)
                    for phi in Phi])

    xs = np.argsort(d_gene)

    # Create output directories if necessary
    path = Path(outpath)
    path.mkdir(parents=True, exist_ok=True)

    # Write x_test to file for analysis
    X_save = X_input.copy()
    X_save.X = X_test.T
    for k in list(X_save.layers.keys()):
        del X_save.layers[k]
    X_save.write(os.path.join(path, 'X_test.h5ad'))

    # Write output to file
    f2 = open(os.path.join(path, 'simulation_results.txt'), 'w')
    f3 = open(os.path.join(path, 'composite_simulation_results.txt'), 'w')
    colnames = ['version', 'Overall pearson', 'Overall spearman', 'Gene average',
                'Sample average', 'Sample dist pearson', 'Sample dist spearman',
                'Gene dist pearson', 'Gene dist spearman',
                'Matrix coherence (90th ptile)', 'Best crossvalidation fold']
    f2.write('\t'.join(colnames) + '\n')
    f3.write('\t'.join(colnames) + '\n')

    # Create empty pandas df to return all analysis
    results_df = pd.DataFrame(columns=colnames)
    results_noNoise_df = pd.DataFrame(columns=colnames)

    for i in xs:
        phi = Phi[i]

        y = get_observations(X_test, phi, snr, normalization=norm)
        x2 = decompress(y, U, phi, correct_comeasured, train_corr)
        x2[np.isnan(x2)] = 0
        results = compare_results(X_test, x2)
        f2.write('\t'.join([str(x) for x in [versions[i]]+results+[d_gene[i]]+[k]]) + '\n')

        y_noNoise = get_observations_no_noise(X_test, phi, normalization=norm)
        x2_noNoise = decompress(y_noNoise, U, phi, correct_comeasured, train_corr)
        x2_noNoise[np.isnan(x2_noNoise)] = 0
        results_noNoise = compare_results(X_test, x2_noNoise)
        f3.write('\t'.join([str(x) for x in [versions[i]]+results_noNoise+[d_gene[i]]+[k]]) + '\n')

        # Either safe the decomposed results X computed from noisy simulated data
        # or simulated data without noise
        if save=='no_noise':
            # Write x_noNoise to anndata
            x2_noNoise_anndata = X_input.copy()
            for k in list(x2_noNoise_anndata.layers.keys()):
                del x2_noNoise_anndata.layers[k]
            x2_noNoise_anndata.X = x2_noNoise.T
            x2_noNoise_anndata.write(os.path.join(path, 'X_simulated_'+str(i)+'.h5ad'))
        elif save=='noise':
            # Write x_noNoise to anndata
            x2_anndata = X_input.copy()
            for k in list(x2_anndata.layers.keys()):
                del x2_anndata.layers[k]
            x2_anndata.X = x2.T
            x2_anndata.write(os.path.join(path, 'X_simulated_'+str(i)+'.h5ad'))
        else:
            # In case no valid type of simulated data to save is given, an error is thrown
            raise ValueError(('The simulation type of data used to decompose the test '
                              'data {0} is not valid.'.format(save) +
                              'Please use one of the following: no_noise or noise'))


        # Add results to pandas df
        results_df.loc[len(results_df)] = [versions[i]]+results+[d_gene[i]] + [k]
        results_noNoise_df.loc[len(results_noNoise_df)] = [versions[i]]+results_noNoise+[d_gene[i]] + [k]

    f2.close()
    f3.close()

    return results_df, results_noNoise_df
