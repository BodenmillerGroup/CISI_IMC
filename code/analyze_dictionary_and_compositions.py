# Import libraries (from original script)
import numpy as np
from analyze_predictions import *

# Import libraries (additionaly added)
from utils import sparse_decode, compare_results, get_observations, simulate_composite_measurements
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

outputs:
    results_df: All performance analyis in one pandas df
    simulation_results.txt: Performance analysis
'''

def analyze_U_and_A(X_input, U, Phi, versions, outpath, k, lasso_sparsity=0.2,
                    THREADS=20, layer=None):
    # Select layer of anndata object that should be used and transpose it
    # to proteins x cells/channels
    if layer is not None:
        X_test = (X_input.layers[layer]).T
    else:
        X_test = (X_input.X).T

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
    X_input.write(os.path.join(path, 'X_test.h5ad'))

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
    results_comp_df = pd.DataFrame(columns=colnames)

    for i in xs:
        phi = Phi[i]
        y = get_observations(X_test, phi, snr=5)
        w = sparse_decode(y, phi.dot(U), sparsity, method='lasso', numThreads=THREADS)
        x2 = U.dot(w)
        x2_anndata = X_input
        x2_anndata.X = x2.T
        x2_anndata.write(os.path.join(path, 'X_simulated_'+str(i)+'.h5ad'))
        results = compare_results(X_test, x2)
        f2.write('\t'.join([str(x) for x in [versions[i]]+results+[d_gene[i]]+[k]]) + '\n')

        y_comp = simulate_composite_measurements(X_test, phi)
        x2_comp = decompress(y_comp, U, phi)
        results_comp = compare_results(X_test, x2_comp)
        f3.write('\t'.join([str(x) for x in [versions[i]]+results_comp+[d_gene[i]]+[k]]) + '\n')

        # Add results to pandas df
        results_df.loc[len(results_df)] = [versions[i]]+results+[d_gene[i]] + [k]
        results_comp_df.loc[len(results_comp_df)] = [versions[i]]+results_comp+[d_gene[i]] + [k]

    f2.close()
    f3.close()

    return results_df, results_comp_df
