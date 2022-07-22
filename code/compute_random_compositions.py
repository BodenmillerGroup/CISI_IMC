# Import libraries (from original script)
import numpy as np
import spams
from scipy.stats import entropy
from scipy.spatial import distance
import os

# Import libraries (additionaly added)
from utils import sparse_decode, compare_results, get_observations
from pathlib import Path
from alive_progress import alive_bar


'''
Estimating best random composition matrix
(Randomly creates possible A matrices satisfying constraints,
 estimating W and analising the fit)

Find A given U, X with special constraints:
inputs:
    X: anndata object containing numpy array X (proteins x cells/pixels)
       with proteins names as X.var_names
    U: a dictionary of gene modules (proteins x modules)
    mode: 'M', max genes per composition is constrained
          'G', maximum times each gene is represented is constrained (default)
    lasso_sparsity: error tolerance applied to l1 sparsity constraint (default: 0.2)
    nmeasurements: number of channels (# rows of A)
    maxcomposition: maximum times each gene is represented (mode G),
                    or max genes per composition (mode M)
    THREADS: # Number of threads used (default=20)
    outpath: If specified, the best 50 U versions are saved as .txt files in
             compositions_A folder
    num_phi: the number of best phi's that are saved (default: 1, at most: 50)
    binary: boolean variable that specifies if Phi is binary or not (default: False)
            If True, then Phi will be binarized directly after computation and the
            best Phi is chosen according to the results of its binarized version

outputs:
    Phi/A: composition matrix (composite channels (measurements) x proteins,
           binary)
    compositions_A/version_<i>.txt: If outpath is specified, the best 50
                                    composition matrices are saved
'''


def compute_A(X_input, U, nmeasurements, maxcomposition, mode='G', lasso_sparsity=0.2,
              outpath=None, THREADS=20, layer=None, num_phi=1, binary=False):
    # Raise error if unsupported mode is given
    if (mode != 'M') and (mode != 'G'):
        raise AssertionError("Unsupported mode 'mode' given!", mode)

    # Select layer of anndata object that should be used and transpose it
    # to proteins x cells/channels
    if layer is not None:
        X = (X_input.layers[layer]).T
    else:
        X = (X_input.X).T

    # Assert that number of phi's to be returned is between 1 and 50
    assert (num_phi>=1 and num_phi<=50), ('Number of Phis to be returned in ' +
                                          'compute_random_compositions is not ' +
                                          'between 1 and 50.')

    # Empirical observation: using a sparsity constraint that is softer than
    # that used during training slightly improves results
    sparsity = lasso_sparsity / 10

    # Find a good composition matrix by generating a bunch of random matrices
    # and testing each for preservation of distances
    print('%d measurements' % nmeasurements)
    best = np.zeros(50)
    Phi = [None for _ in best]

    # Progress bar depicting the progress of
    with alive_bar(2000, title="Iterating over random A's...") as bar:
        for _ in range(2000):
            # Initialze random A with constraints
            while True:
                if mode == 'M':
                    phi = random_phi_subsets_m(nmeasurements, X.shape[0],
                                   (2, maxcomposition), d_thresh=0.8)
                elif mode == 'G':
                    phi = random_phi_subsets_g(nmeasurements, X.shape[0],
                                   (1, maxcomposition), d_thresh=0.8)
                if check_balance(phi):
                    break
            if _%100 == 0:
                print(_, best)

            # If binary= True then binarize phi (A)
            if binary:
                phi = np.where(phi > 0, 1, 0)

            # Initialze composite oservations Y using X and noise
            y = get_observations(X, phi, snr=5)
            # Compute W given Y, A and U
            w = sparse_decode(y, phi.dot(U), sparsity, method='lasso', numThreads=THREADS)
            x2 = U.dot(w)
            # Compare X to predicted X
            r = compare_results(X, x2)
            # Update best A
            if r[2] > best.min():
                i = np.argmin(best)
                best[i] = r[2]
                Phi[i] = phi

            # Mark end on one iteration for progressbar
            bar()

    xs = np.argsort(best)
    best = best[xs[::-1]]

    # Sort Phi, such that best composite matrix is at the beginning
    Phi = [Phi[i] for i in xs]

    # If outpath is specified, then the best num_phi composite matrices are saved
    # into files
    if outpath!=None:
        path = Path(os.path.join(outpath, 'compositions_A'))
        path.mkdir(parents=True, exist_ok=True)
        for i in xs[:num_phi]:
            f1 = open(os.path.join(path, 'version_%d.txt' % i), 'w')
            phi = Phi[i]

            # Add protein names
            f1.write('\t'.join(X_input.var_names) + '\n')

            for j in range(nmeasurements):
                genes = ['channel %d' % j]
                for k in range(X.shape[0]):
                    genes.append(str(phi[j, k]))
                f1.write('\t'.join(genes) + '\n')
            f1.close()

    return Phi[0], best[0]

# For mode 'M': Create random A matrix with at most n[1] selected proteins
# per channel (at most n[1] 1s per row)
def random_phi_subsets_m(m, g, n, d_thresh=0.4):
    Phi = np.zeros((m, g))
    Phi[0, np.random.choice(g, np.random.randint(n[0], n[1]), replace=False)] = 1
    Phi[0] /= Phi[0].sum()
    # Add random rows with at most n[1] 1's which have a min distance from last rows
    # of 1 - dist > d_thresh (such that rows/compositions are not to similar) until
    # Phi has m measurements
    for i in range(1, m):
        dmax = 1
        while dmax > d_thresh:
            p = np.zeros(g)
            p[np.random.choice(g, np.random.randint(n[0], n[1]), replace=False)] = 1
            p /= p.sum()
            dmax = 1 - distance.cdist(Phi[:i], [p], 'cosine').min()
        Phi[i] = p
    return Phi


# For mode 'G': Create random A matrix with proteins being selected at most n[1]
# times over all channels (at most n[1] 1s per column)
def random_phi_subsets_g(m, g, n, d_thresh=0.4):
    Phi = np.zeros((m, g))
    Phi[np.random.choice(m,np.random.randint(n[0], n[1]), replace=False), 0] = 1
    for i in range(1, g):
        dmax = 1
        # Add random columns with at most n[1] 1's which have a min distance from
        # last column of 1 - dist > d_thresh (such that columns are not to similar)
        # until Phi has g columns
        while dmax > d_thresh:
            p = np.zeros(m)
            p[np.random.choice(m,np.random.randint(n[0], n[1]), replace=False)] = 1
            dmax = 1 - distance.cdist(Phi[:,:i].T, [p], 'correlation').min()
        Phi[:,i] = p
    Phi = (Phi.T / Phi.sum(1)).T
    return Phi


# Check that no channel is unused or that difference between the number of times
# each protein is used is not to big
def check_balance(Phi, thresh=4):
    x = Phi.sum(0) + 1e-7
    if (x.max() / x.min() > thresh) or (Phi.sum(1).min() == 0):
        return False
    else:
        return True
