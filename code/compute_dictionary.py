# Import libraries (from original script)
import numpy as np
import spams
from scipy.stats import entropy
from scipy.spatial import distance

# Import libraries (additionaly added)
from pathlib import Path
from utils import sparse_decode
import os
# from scipy.stats import zscore
import pandas as pd


'''
Sparse Module Activity Factorization (SMAF)
(Computes a dictionary (U) from training data)

Find X = UW with special constraints:
inputs:
    X_input: anndata object containing numpy array X (cells/pixels x proteins)
    d: the number of features (columns) in the dictionary U
    lda1: in mode 1 (recommended) the number of nonzeros per column in W
    lda2: an error threshold - when optimizing over U we will search for the
          sparsest fit while tolerating at most this error
    mode: Mode in which spams lasso fnc. regularizes U and W (constraint of min.)
          (default: mode=1, Note: mode=2 seems to produce empty matrices in some cases,
          maybe constraint are to stringthend?)
          Find columns alpha for each column x in X:
          Mode=1
            min_{alpha} + ||alpha||_1 s.t. ||x-W*alpha||_2^2 <= lambda1
          Mode=2
            min_{alpha} + 0.5||x-W*alpha||_2^2 + lambda1||alpha||_1 +0.5 lambda2||alpha||_2^2
    THREADS: # Number of threads used (default=4)
    doprint: Print correlations between predicted X and real X and some additional info
    normalization: How data is normalized before running smaf (default: paper_norm)
                   Options: paper_norm (normalization used in paper, protein-wise),
                   min_max_norm (protein-wise), zscore_norm (protein-wise) or none
    layer: Name of layer in anndata object to be used as X (default: None, =anndata.X)
    others: All other parameters belong to the fnc. call of lasso from spams and
            further information is available on the respective website

outputs:
    U: a dictionary of gene modules (proteins x modules)
    W: the module activity levels in each cell (/pixel?) of training data
       (modules x pixels)

Example:
    d = 500
    k = 15
    UW = (np.random.random((X.shape[0],d)),np.random.random((d,X.shape[1])))
    U,W = smaf(X,d,k,0.1,maxItr=10,activity_lower=0.,module_lower=min(20,X.shape[0]/20),donorm=True,mode=1,use_chol=True,mink=5,doprint=True)
'''
def smaf(X_input, d, lda1, lda2, maxItr=10, UW=None, posW=False, posU=True,
         use_chol=False, module_lower=1, activity_lower=1, donorm=False,
         mode=1, mink=0, U0=[], U0_delta=0.1, doprint=False, THREADS=4,
         outpath=None, normalization='paper_norm', layer=None):

    # Select layer of anndata object that should be used in SMAF and transpose it
    # to proteins x cells/channels
    if layer is not None:
        X_mat = (X_input.layers[layer]).T
    else:
        X_mat = (X_input.X).T

    # Extract data matrix X from anndata object and apply selected normalization
    match normalization:
        case 'paper_norm':
            # Normalize X? (otherwise spams.nmf underneath will estimate U to only contain
            # 0s and smaf won't work)
            # Normalizaiton: Every row in X_input is divided by the corresponding vector
            # element in the rowwise norm (proteins are normalized accross all cells/pixels)
            X = (X_mat.T / np.linalg.norm(X_mat, axis=1)).T
        case 'min_max_norm':
            X = (X_mat-X_mat.min(axis=1, keepdims=True)) / (
                X_mat.max(axis=1, keepdims=True)-X_mat.min(axis=1, keepdims=True)
                )
        case 'zscore_norm':
            X = ((X_mat.T-np.mean(X_mat, axis=1)) / np.std(X_mat, axis=1)).T # zscore(X_mat, axis=1)
        case 'none':
            X = X_mat
        case _:
            # In case no valid normalization is given, an error is thrown
            raise ValueError(('The normalization {0} used by smaf is not valid.'.format(normalization) +
                              'Please use one of the following: paper_norm, min_max_norm, ' +
                              'min_max_norm_channelwise, zscore_norm or none.'))

    # use Cholesky when we expect a very sparse result
    # this tends to happen more on the full vs subsampled matrices
    if UW == None:
        # Initialze U, W randomly
        U, W = spams.nmf(np.asfortranarray(X), return_lasso=True, K = d, numThreads=THREADS)
        W = np.asarray(W.todense())
    else:
        U, W = UW
    Xhat = U.dot(W)
    Xnorm = np.linalg.norm(X)**2 / X.shape[1]

    # For maxItr iterations approximate U, W by:
    # a. Update the module dictionary as U=LassoNegative()
    # b. Normalize each module with L2-norm=1
    # c. Update the activity levels as W=OMP()
    for itr in range(maxItr):
        # Default mode using lasso
        if mode == 1:
            # In this mode the ldas correspond to an approximate desired fit
            # Higher lda will be a worse fit, but will result in a sparser sol'n
            U = spams.lasso(np.asfortranarray(X.T), D=np.asfortranarray(W.T),
                            lambda1=lda2*Xnorm, mode=1, numThreads=THREADS,
                            cholesky=use_chol, pos=posU)
            U = np.asarray(U.todense()).T
        elif mode == 2:
            if len(U0) > 0:
                U = projected_grad_desc(W.T, X.T, U.T, U0.T, lda2, U0_delta, maxItr=400)
                U = U.T
            else:
                U = spams.lasso(np.asfortranarray(X.T), D=np.asfortranarray(W.T),
                                lambda1=lda2, lambda2=0.0, mode=2, numThreads=THREADS,
                                cholesky=use_chol, pos=posU)
                U = np.asarray(U.todense()).T
        # Default: False
        if donorm:
            U = U / np.linalg.norm(U, axis=0)
            U[np.isnan(U)] = 0
        # Default mode=1: estimate W
        if mode == 1:
            wf = (1 - lda2)
            W = sparse_decode(X, U, lda1, worstFit=wf, mink=mink, numThreads=THREADS)
        elif mode == 2:
            if len(U0) > 0:
                W = projected_grad_desc(U, X, W, [], lda1, 0., nonneg=posW, maxItr=400)
            else:
                W = spams.lasso(np.asfortranarray(X), D=np.asfortranarray(U),
                                lambda1=lda1, lambda2=1.0, mode=2, numThreads=THREADS,
                                cholesky=use_chol, pos=posW)
                W = np.asarray(W.todense())

        Xhat = U.dot(W)
        module_size = np.average([np.exp(entropy(u)) for u in U.T if u.sum()>0])
        activity_size = np.average([np.exp(entropy(abs(w))) for w in W.T])
        # Print correlation between X and estimated X, module size, activity size,
        # lda1 and lda2
        if doprint:
            print(distance.correlation(X.flatten(), Xhat.flatten()), module_size,
                  activity_size, lda1, lda2)
        if module_size < module_lower:
            lda2 /= 2.
        if activity_size < activity_lower:
            lda2 /= 2.

    # Remove empty columns (proteins that are never chosen?)
    U = U[:, (U.sum(0) > 0)]

    # If an outpath is given
    if outpath!=None:
        path = Path(outpath)
        path.mkdir(parents=True, exist_ok=True)
        # Save U
        np.save(os.path.join(path, 'gene_modules.npy'), U)
        pd.DataFrame(U, columns=list(range(1, U.shape[1]+1)),
                     index=X_input.var_names).to_csv(os.path.join(path, 'gene_modules.csv'))

    return U, W
