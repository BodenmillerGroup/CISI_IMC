# Import libraries
import numpy as np
from utils import sparse_decode_blocks


'''
Analyze performance of computed U and A
(Nine different metrics, correlations and distances gene- and sample-wise)

For given X, U and A, simulate composite results y and analyze performance:
inputs:
    y: Compressed measurements (composite channels (measurements) x samples)
    U: a dictionary of gene modules (proteins x modules)
    phi: list containing (best) composition matrices (composite channels (measurements) x proteins,
         binary)
    sparsity: an error threshold - when optimizing over U we will search for the
              sparsest fit while tolerating at most this error (default: 0.1)
    method: 'lasso' or 'omp' used to decompress W
    numThreads: # Number of threads used (default=20)
    outpath: If specified, the best 50 U versions are saved as .txt files in
             compositions_A folder

outputs:
    Estimated X (proteins x cells/pixels)
'''

def decompress(y, U, phi, sparsity=0.1, method='lasso', numThreads=20,
               worstFit=1., mink=0, nonneg=True, num_blocks=20):
    # Call fnc sparse_decode to compute W
    w = sparse_decode_blocks(y, phi.dot(U), sparsity, numThreads, method,
                             worstFit, mink, nonneg, num_blocks)
    # Calculate estimated X
    x2 = U.dot(w)
    #x2[np.isnan(x2)] = 0

    return x2
