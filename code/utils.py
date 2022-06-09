# Import libraries
import numpy as np
import spams
import os
from pathlib import Path
import pandas as pd
from re import search
from analyze_predictions import *


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


# Function that reads in a cell.csv type file and prepares it for decomposition
# training (selects only meanIntensity, ImageNumber and ObjectNumber columns and
# determines protein names)
def read_input(X_input, panel_input):
    # Check if input is a valid csv file
    if check_file(X_input, ['.csv']) and check_file(panel_input, ['.csv']):
        # Read intensities form .csv file
        X = np.genfromtxt(X_input, delimiter=',', skip_header=1)
        # Read in header names
        X_header = pd.read_csv(X_input, nrows=0).columns.tolist()
        # Get list of mean intensity columns and Image/Object number columns
        selected_columns = list(map(lambda col:
            bool(search('(Intensity_MeanIntensity_FullStackFiltered_c|ImageNumber|ObjectNumber)', col)),
            X_header))
        # Filter for above selected columns
        X = X[:, selected_columns]
        # Remove everything before channel numbers in X_header
        X_header = [x.replace("Intensity_MeanIntensity_FullStackFiltered_c", "")
                    for x in X_header]

        # Read in panel data
        panel = pd.read_csv(panel_input)
        proteins = [panel.loc[panel['channel'] == int(i), 'Target'].values[0]
                    for i in X_header if i not in ['ImageNumber', 'ObjectNumber']]

    else:
        # If it is none of the above input types an error is
        raise ValueError('''No valid input for X was given.
              The input for X needs to be a .csv, .npy or numpy array with
              dimensions: proteins x cells (/pixels?).''')

    return X, X_header, proteins


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
