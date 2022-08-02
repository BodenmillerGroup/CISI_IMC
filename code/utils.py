# Import libraries
import numpy as np
import spams
import os
from pathlib import Path
import pandas as pd
from re import search
from analyze_predictions import *
import math
import random


'''
Helper file
(Contains general helper files used in more than one file)
'''


'''
New functions
'''

# Function that checks that file_path is an existing file and has a certain extension
def check_file(file_path, extension):
    file = Path(file_path)
    if file.exists() and file.is_file() and file.suffix in extension:
        return True
    else:
        return False


# Simulate composite measurements using composite matrix A assuming no noise
def simulate_composite_measurements(X, phi):
    # Simulate composite measurements
    return phi.dot(X)


'''
Old functions (from original publication)
'''

# Fixing U, approximating W
def sparse_decode(Y, D, k, numThreads, method='omp', worstFit=1., mink=0, nonneg=False):
    if method == 'omp':
        while k > mink:
            W = spams.omp(np.asfortranarray(Y), np.asfortranarray(D), L=k,
                 numThreads=numThreads)
            W = np.asarray(W.todense())
            fit = 1 - np.linalg.norm(Y - D.dot(W))**2 / np.linalg.norm(Y)**2
            if fit < worstFit:
                break
            else:
                k -= 1
    elif method == 'lasso':
        Ynorm = np.linalg.norm(Y)**2 / Y.shape[1]
        W = spams.lasso(np.asfortranarray(Y), np.asfortranarray(D), lambda1=k*Ynorm,
                  mode=1, numThreads=numThreads, pos=nonneg)
        W = np.asarray(W.todense())

    return W


# Fixing U, approximating W (in blocks)
def sparse_decode_blocks(Y, D, lda=0.1, numThreads=20, method='omp', worstFit=1.,
                         mink=0, nonneg=False, num_blocks=20):
	W = np.zeros((D.shape[1], Y.shape[1]))
	ynorm = np.linalg.norm(Y, axis=0)
	xs = np.argsort(ynorm)
	block_size = int(len(xs) / num_blocks)
	for i in range(0, len(xs), block_size):
		idx = xs[i:i+block_size]
		w = sparse_decode(Y[:, idx], D, lda, numThreads, method, worstFit,
                    mink, nonneg)
		W[:, idx] = w
	return W


# Simulate random composite observations Y from X
def get_observations(X0, Phi, snr=5, return_noise=False):
    noise = np.array([np.random.randn(X0.shape[1]) for _ in range(X0.shape[0])])
    noise *= np.linalg.norm(X0)/np.linalg.norm(noise)/snr
    if return_noise:
        return Phi.dot(X0 + noise), noise
    else:
        return Phi.dot(X0 + noise)


# Compare X to predicted X by correlations and distances between them
def compare_results(A, B):
    results = list(correlations(A, B, 0))[:-1]
    results += list(compare_distances(A, B))
    results += list(compare_distances(A.T, B.T))
    return results
