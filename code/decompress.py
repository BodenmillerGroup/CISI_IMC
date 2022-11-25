# Import libraries
import numpy as np
from utils import sparse_decode_blocks, select_and_correct_comeasured


'''
Decompress Composite Data
(using specified U and A)

For given y, U and A, decompose composite measurements y using the given U and A:
inputs:
    y: Compressed measurements (composite channels (measurements) x samples)
    U: a dictionary of gene modules (proteins x modules)
    phi: composite matrix used for composite measurements y (composite channels
         (measurements) x proteins)
    sparsity: an error threshold when computing W using method='lasso' and the sparsity
              k when using method='lasso'
    method: 'lasso' or 'omp' used to decompress W (default: lasso)
    numThreads: # Number of threads used (default=20)
    correct_comeasured: Correct comeasured genes that are not coexpressed (default: False)
    train_corr: Correlations between genes in training data X (genes that are coexpressed)
                used in correct_comeasured (default: None)
    num_blocks: number of blocks used to calculate W (should be bigger for pixel-wise?)
                (default: 20)

outputs:
    x2: Numpy array containing decompressed X (proteins x cells/pixels)
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
