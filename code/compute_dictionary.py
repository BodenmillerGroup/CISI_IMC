# Import libraries (from original script)
import numpy as np
import spams
from scipy.stats import entropy
from scipy.spatial import distance

# Import libraries (additionaly added)
from pathlib import Path
from utils import sparse_decode
import os


'''
Sparse Module Activity Factorization (SMAF)
(Computes a dictionary (U) from training data)

Find X = UW with special constraints:
inputs:
    X: proteins x cells (/pixels?) -> numpy array
    d: the number of features (columns) in the dictionary U
    lda1: in mode 1 (recommended) the number of nonzeros per column in W
    lda2: an error threshold - when optimizing over U we will search for the
          sparsest fit while tolerating at most this error
    THREADS: # Number of threads used (default=4)
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
         outpath=None):

    # Normalize X? (otherwise spams.nmf underneath will estimate U to only contain
    # 0s and smaf won't work)
    X = (X_input.T / np.linalg.norm(X_input, axis=1)).T

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
            U = U / np.linalg.norm(U,axis=0)
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
        np.savetxt(os.path.join(path, 'gene_modules.csv'), U, delimiter=',')

    return U, W
