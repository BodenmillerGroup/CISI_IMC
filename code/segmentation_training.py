# Import libraries
from utils import read_input, split_X
from compute_dictionary import smaf
from compute_random_compositions import compute_A
from analyze_dictionary_and_compositions import analyze_U_and_A


## Main
if __name__ == '__main__':

    # Define input paths
    X_input = '/Users/lelo/Desktop/Munchkin/ETH/Master Thesis/CISI_IMC/tests/data/cell.csv'
    panel_input = '/Users/lelo/Desktop/Munchkin/ETH/Master Thesis/CISI_IMC/tests/data/TH152_panel.csv'

    # Define output path
    outpath = None

    # Set parameters
    d = 80
    lda1 = 3
    lda2 = 0.2
    nmeasurements = 10
    maxcomposition = 3

    ## Start training
    # Read in inputs
    X, X_header, proteins = read_input(X_input, panel_input)

    # Split into training, validate and test set
    X_training, X_validate, X_test = split_X(X, X_header, set_sizes=[[3], [5], [2]])

    # Compute U
    U, W = smaf(X_training, d, lda1, lda2, maxItr=10, UW=None, posW=False, posU=True,
                use_chol=False, module_lower=1, activity_lower=1, donorm=False,
                mode=1, mink=0, U0=[], U0_delta=0.1, doprint=False, THREADS=4,
                outpath=outpath)

    # Compute A
    Phi = compute_A(X_validate, U, nmeasurements, maxcomposition, mode='G',
                    lasso_sparsity=0.2, outpath=outpath, proteins=proteins,
                    THREADS=20)

    # Analize training
    analyze_U_and_A(X_test, U, Phi, outpath='data', lasso_sparsity=0.2, THREADS=20)
