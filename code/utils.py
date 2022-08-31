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
def get_observations_no_noise(X_input, phi, normalization='none'):
    # Extract data matrix X from anndata object and apply selected normalization
    match normalization:
        case 'paper_norm':
            # Normalize X? (otherwise spams.nmf underneath will estimate U to only contain
            # 0s and smaf won't work)
            # Normalizaiton: Every row in X_input is divided by the corresponding vector
            # element in the rowwise norm (proteins are normalized accross all cells/pixels)
            X = (X_input.T / np.linalg.norm(X_input, axis=1)).T
        case 'min_max_norm':
            X = (X_input-X_input.min(axis=1, keepdims=True)) / (
                X_input.max(axis=1, keepdims=True)-X_input.min(axis=1, keepdims=True)
                )
        # case 'zscore_norm':
        #     X = ((X_mat.T-np.mean(X_mat, axis=1)) / np.std(X_mat, axis=1)).T # zscore(X_mat, axis=1)
        case 'none':
            X = X_input
        case _:
            # In case no valid normalization is given, an error is thrown
            raise ValueError(('The normalization {0} used by smaf is not valid.'.format(normalization) +
                              'Please use one of the following: paper_norm, min_max_norm, ' +
                              'or none.'))
    # Simulate composite measurements
    return phi.dot(X)


'''
Old functions (from original publication)
'''

# Fixing U, approximating W
def sparse_decode(Y, D, k, numThreads, method='omp', worstFit=1., mink=0, nonneg=True):
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
                         mink=0, nonneg=True, num_blocks=20):
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
def get_observations(X0, Phi, snr=5, return_noise=False, normalization='none'):
    # Extract data matrix X from anndata object and apply selected normalization
    match normalization:
        case 'paper_norm':
            # Normalize X? (otherwise spams.nmf underneath will estimate U to only contain
            # 0s and smaf won't work)
            # Normalizaiton: Every row in X_input is divided by the corresponding vector
            # element in the rowwise norm (proteins are normalized accross all cells/pixels)
            X = (X0.T / np.linalg.norm(X0, axis=1)).T
        case 'min_max_norm':
            X = (X0-X0.min(axis=1, keepdims=True)) / (
                X0.max(axis=1, keepdims=True)-X0.min(axis=1, keepdims=True)
                )
        # case 'zscore_norm':
        #     X = ((X_mat.T-np.mean(X_mat, axis=1)) / np.std(X_mat, axis=1)).T # zscore(X_mat, axis=1)
        case 'none':
            X = X0
        case _:
            # In case no valid normalization is given, an error is thrown
            raise ValueError(('The normalization {0} used by smaf is not valid.'.format(normalization) +
                              'Please use one of the following: paper_norm, min_max_norm, ' +
                              'or none.'))

    noise = np.array([np.random.randn(X.shape[1]) for _ in range(X.shape[0])])
    noise *= np.linalg.norm(X)/np.linalg.norm(noise)/snr
    if return_noise:
        return Phi.dot(X + noise), noise
    else:
        return Phi.dot(X + noise)


# Compare X to predicted X by correlations and distances between them
def compare_results(A, B):
    results = list(correlations(A, B, 0))[:-1]
    results += list(compare_distances(A, B))
    results += list(compare_distances(A.T, B.T))
    return results


# Find genes that are comeasured but are not coexpressed (correlation<threshold)
# and correct their expressions
def select_and_correct_comeasured(x, xc, phi, phi_corr, training_corr,
                                  phi_thresh=0.6, train_thresh=0.1):
	# find comeasured genes that are not coexpressed
	comeasured = []
	for i in range(phi_corr.shape[0]):
		xs = np.argsort(-phi_corr[i])
		for j in xs:
			if phi_corr[i, j] < phi_thresh:
				break
			comeasured.append((phi_corr[i, j], i, j))
	corrected_pairs = []
	for c in sorted(comeasured, reverse=True):
		if training_corr[c[1], c[2]] < train_thresh:
			x, both_nz = correct_coexpression(x, xc, phi, c[1], c[2])
			corrected_pairs.append((c[1], c[2], both_nz))
	return x, corrected_pairs


# Set closest not coexpressed genes to 0?
def correct_coexpression(x, xc, phi, i, j):
	# pick the gene with nearest expression pattern in scRNA
	thresh_i = np.percentile(x[i], 99.9) / 100
	thresh_j = np.percentile(x[j], 99.9) / 100
	both_nz = (x[i] > thresh_i)*(x[j] > thresh_j)
	dist = distance.cdist([phi[:, i], phi[:, j]], xc[:, both_nz].T, 'correlation')
	i_closer = np.where(both_nz)[0][dist[0] < dist[1]]
	j_closer = np.where(both_nz)[0][dist[0] > dist[1]]
	x[i, j_closer] = 0
	x[j, i_closer] = 0
	return x, both_nz
