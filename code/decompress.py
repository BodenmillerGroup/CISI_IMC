# Import libraries
import numpy as np
from utils import sparse_decode_blocks, select_and_correct_comeasured


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
    correct_comeasured: Correct comeasured genes that are not coexpressed (default: False)
    train_corr: Correlations between genes in training data X (genes that are coexpressed)
                used in correct_comeasured (default: None)

outputs:
    Estimated X (proteins x cells/pixels)
'''

def decompress(y, U, phi, sparsity=0.1, method='lasso', numThreads=20,
               worstFit=1., mink=0, nonneg=True, num_blocks=20,
               correct_comeasured=False, train_corr=None):
    # Call fnc sparse_decode to compute W
    w = sparse_decode_blocks(y, phi.dot(U), sparsity, numThreads, method,
                             worstFit, mink, nonneg, num_blocks)
    # Calculate estimated X
    x2 = U.dot(w)

    # Set nan and values smaller than 0 to 0 in x2 (not possible values)
    x2[np.isnan(x2)] = 0
    x2[x2 < 0] = 0

    # Correct with phi_corr and training_corr
    if correct_comeasured:
        # Correlations between genes in Phi/A (genes that are comeasured) used
        # in correct_comeasured
        phi_corr = (np.einsum('ij,ik->jk', phi, phi)/phi.sum(0)).T - np.eye(phi.shape[1])
        # Correct X
        x2, cp = select_and_correct_comeasured(x2, y, phi, phi_corr, train_corr)

    return x2
