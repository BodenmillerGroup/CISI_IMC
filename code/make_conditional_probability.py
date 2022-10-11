# Import libraries
import numpy as np

# Import libraries (additionaly added)
import os


'''
Compute pairwise conditional probabilities between proteins/channels

Matrix (proteins x proteins) contains conditional probability of protein i given
protein j (calculated as a cell being positiv for a specified protein, if the
protein has an expression of at least > threshold)

inputs:
    X_input: anndata object containing numpy array X (cells/pixels x proteins)
    outpath: Specify where output files should be saved to
    threshold: threshold to determine if a cell is 'positive' for a given protein
'''

def make_cond_prob(X_input, outpath, threshold=10.0):

        # Extract expression matrix from anndata object and transpose
        X = (X_input.X).T
        # Calculate conditional probabilities
        cond_prob = gene_conditional_probs(X, threshold)
        # Save probabilities
        np.savetxt(os.path.join(outpath, 'conditional_probability.csv'), cond_prob,
                delimiter=',')

        return cond_prob


# Function that actually does the calculations
def gene_conditional_probs(X, thresh):
        # Create empty matrix of the shape protein x protein
        CP = np.eye(X.shape[0])
        # Loop through all protein combinations
        for i in range(X.shape[0]):
                a = (X[i] > thresh)
                for j in range(X.shape[0]):
                        b = (X[j] > thresh)
                        # Compute CP by multiplying the average number of times
                        # cells were positive for both given cells by the average
                        # number of times protein a was positive
                        CP[i,j] = np.average(a*b) / np.average(a)
        return CP
