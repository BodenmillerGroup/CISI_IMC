# Import libraries
import numpy as np

# Import libraries (additionaly added)
import os


'''
Compute pairwise conditional probabilities between proteins/channels

Matrix (proteins x proteins) contains if protein i is True (bigger than threshold)
given that protein j is True (bigger than threshold) for all possible combinations

inputs:
    X_input: anndata object containing numpy array X (cells/pixels x proteins)
             Will be divided into: training, validate and test set
    outpath: Specify where output files should be saved to
    threshold:
'''

def make_cond_prob(X_input, outpath, threshold=10.0):

        X = (X_input.X).T
        cond_prob = gene_conditional_probs(X, threshold)
        np.savetxt(os.path.join(outpath, 'conditional_probability.csv'), cond_prob,
                delimiter=',')

        return cond_prob


# Function that actually does the calculations
def gene_conditional_probs(X, thresh):
        CP = np.eye(X.shape[0])
        for i in range(X.shape[0]):
                a = (X[i] > thresh)
                for j in range(X.shape[0]):
                        b = (X[j] > thresh)
                        CP[i,j] = np.average(a*b) / np.average(a)
        return CP
