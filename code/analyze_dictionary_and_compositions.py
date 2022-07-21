# Import libraries (from original script)
import numpy as np
from analyze_predictions import *

# Import libraries (additionaly added)
from utils import sparse_decode, compare_results, get_observations
from pathlib import Path
import sys
import os
import pandas as pd


'''
Analyze performance of computed U and A
(Nine different metrics, correlations and distances gene- and sample-wise)

For given X, U and A, simulate composite results y and analyze performance:
inputs:
    X: anndata object containing numpy array X (proteins x cells/pixels)
       with proteins names as X.var_names
    U: a dictionary of gene modules (proteins x modules)
    Phi: list containing (best) composition matrices (composite channels (measurements) x proteins,
         binary)
    lasso_sparsity: error tolerance applied to l1 sparsity constraint (default: 0.2)
    THREADS: # Number of threads used (default=20)
    outpath: If specified, the best 50 U versions are saved as .txt files in
             compositions_A folder

outputs:
    results_df: All performance analyis in one pandas df
    simulation_results.txt: Performance analysis
'''

def analyze_U_and_A(X_input, U, Phi, outpath, lasso_sparsity=0.2, THREADS=20,
                    layer=None):
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
    # Write output to file
    f2 = open(os.path.join(path, 'simulation_results.txt'), 'w')
    colnames = ['version', 'Overall pearson', 'Overall spearman', 'Gene average',
                'Sample average', 'Sample dist pearson', 'Sample dist spearman',
                'Gene dist pearson', 'Gene dist spearman',
                'Matrix coherence (90th ptile)']
    f2.write('\t'.join(colnames) + '\n')

    # Create empty pandas df to return all analysis
    results_df = pd.DataFrame(columns=colnames)

    for i in xs:
        phi = Phi[i]
        y = get_observations(X_test, phi, snr=5)
        w = sparse_decode(y, phi.dot(U), sparsity, method='lasso', numThreads=THREADS)
        x2 = U.dot(w)
        results = compare_results(X_test, x2)
        f2.write('\t'.join([str(x) for x in [i]+results+[d_gene[i]]]) + '\n')
        print(d_gene[i], compare_results(X_test, x2))

        # Add results to pandas df
        results_df.loc[len(results_df)] = [i]+results+[d_gene[i]]
    f2.close()

    return results_df
