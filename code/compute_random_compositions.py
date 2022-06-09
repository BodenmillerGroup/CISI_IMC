# Import libraries (from original script)
import numpy as np
import spams
from scipy.stats import entropy
from scipy.spatial import distance
import os

# Import libraries (additionaly added)
from utils import sparse_decode, compare_results, get_observations
from pathlib import Path


'''
Estimating best random composition matrix
(Randomly creates possible A matrices satisfying constraints,
 estimating W and analising the fit)

Find A given U, X with special constraints:
inputs:
    X: proteins x cells (/pixels?) -> .csv / .npy / numpy array
    U: a dictionary of gene modules (proteins x modules)
    mode: 'M', max genes per composition is constrained
          'G', maximum times each gene is represented is constrained (default)
    lasso_sparsity: error tolerance applied to l1 sparsity constraint (default: 0.2)
    nmeasurements: number of channels (# rows of A)
    maxcomposition: maximum times each gene is represented (mode G),
                    or max genes per composition (mode M)
    THREADS: # Number of threads used (default=20)

outputs:
    A: composition matrix (composite channels (measurements) x proteins)


    W: the module activity levels in each cell (/pixel?) of training data
       (modules x pixels)
'''


def compute_A(X, U, nmeasurements, maxcomposition, mode='G', lasso_sparsity=0.2,
              outpath=None, proteins=[], THREADS=20):
    # Raise error if unsupported mode is given
    if (mode != 'M') and (mode != 'G'):
        raise AssertionError("Unsupported mode 'mode' given!", mode)

    # Empirical observation: using a sparsity constraint that is softer than
    # that used during training slightly improves results
    sparsity = lasso_sparsity / 10

    # Find a good composition matrix by generating a bunch of random matrices
    # and testing each for preservation of distances
    print('%d measurements' % nmeasurements)
    best = np.zeros(50)
    Phi = [None for _ in best]
    for _ in range(100):#2000):
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

    xs = np.argsort(best)
    best = best[xs[::-1]]
    # Sort Phi, such that best composite matrix is at the beginning
    Phi = [Phi[i] for i in xs]
    # Correlation between genes?
    d_gene = np.array([np.percentile(1 - distance.pdist(phi.dot(U).T, 'correlation'), 90)
                    for phi in Phi])
    xs = np.argsort(d_gene)

    # Binarize phi (A)
    Phi_binary = [np.where(Phi[i] > 0, 1, 0) for i in xs]

    # If outpath is specified, then the best 50 composite matrices are saved
    # into files
    # If proteins is specified, protein names are added as a header
    if outpath!=None:
        path = Path(os.path.join(outpath, 'compositions_A'))
        path.mkdir(parents=True, exist_ok=True)
        for i in xs:
            f1 = open(os.path.join(path, 'version_%d.txt' % i), 'w')
            phi = Phi_binary[i]

            # Add protein names if available
            if proteins:
                genes = ['channels'] + proteins
                f1.write('\t'.join(genes) + '\n')

            for j in range(nmeasurements):
                genes = ['channel %d' % j]
                for k in range(X.shape[0]):
                    genes.append(str(phi[j, k]))
                f1.write('\t'.join(genes) + '\n')
            f1.close()

    return Phi_binary

# For mode 'M': Create random A matrix with at most n[1] selected proteins
# per channel (at most n[1] 1s per row)
def random_phi_subsets_m(m, g, n, d_thresh=0.4):
    Phi = np.zeros((m, g))
    Phi[0, np.random.choice(g, np.random.randint(n[0], n[1]), replace=False)] = 1
    Phi[0] /= Phi[0].sum()
    for i in range(1, m):
        dmax = 1
        # TODO: Check what happens if dmax < d_thresh, does it just add the
        # same p until nmeasurements?
        while dmax > d_thresh:
            p = np.zeros(g)
            p[np.random.choice(g,np.random.randint(n[0], n[1]), replace=False)] = 1
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
        # TODO: Check what happens if dmax < d_thresh, does it just add the
        # same p until nmeasurements?
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
