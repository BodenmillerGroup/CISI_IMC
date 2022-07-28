# Import libraries (from original script)
import numpy as np
from scipy.spatial import distance

# Import libraries (additionaly added)
import os


'''
Compute pairwise correlations between proteins/channels

Save into correlations.csv file and return correlations:
inputs:
    X_input: anndata object containing numpy array X (cells/pixels x proteins)
             Will be divided into: training, validate and test set
    outpath: Specify where output files should be saved to
'''

def make_cor(X_input, outpath):

    X = (X_input.X).T
    corr = 1 - distance.pdist(X, 'correlation')
    np.savetxt(os.path.join(outpath, 'correlations.csv'), corr, delimiter=',')

    return corr
