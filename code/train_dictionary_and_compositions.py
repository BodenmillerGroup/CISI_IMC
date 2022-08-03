# Import libraries
from compute_dictionary import smaf
from compute_random_compositions import compute_A
from analyze_dictionary_and_compositions import analyze_U_and_A
from alive_progress import alive_bar
import numpy as np
import os
from make_correlations import make_cor
from make_conditional_probability import make_cond_prob


'''
Compute dictionary U and random composite matrix A/Phi from training data and
save results from crossvalidation

Find U, A/phi given X and return goodness of fit:
inputs:
    X_input: anndata object containing numpy array X (cells/pixels x proteins)
             Will be divided into: training, validate and test set
    outpath: Specify where output files should be saved to (used in all fnc)
    layer: Name of layer in anndata object to be used as X (default: None, =anndata.X)
           (used in all fnc)
    lasso_sparsity: error tolerance applied to l1 sparsity constraint (default: 0.2)
                    (used in compute_A and analyze_U_and_A)

    For fnc smaf():
        d: the number of features (columns) in the dictionary U
        lda1: in mode 1 (recommended) the number of nonzeros per column in W
        lda2: an error threshold - when optimizing over U we will search for the
              sparsest fit while tolerating at most this error
        mode_smaf: Mode in which spams lasso fnc. regularizes U and W (constraint of min.)
                   (default: mode=1, Note: mode=2 seems to produce empty matrices in some cases,
                   maybe constraint are to stringthend?)
                   Find columns alpha for each column x in X:
                   Mode=1
                    min_{alpha} + ||alpha||_1 s.t. ||x-W*alpha||_2^2 <= lambda1
                   Mode=2
                    min_{alpha} + 0.5||x-W*alpha||_2^2 + lambda1||alpha||_1 + 0.5 lambda2||alpha||_2^2
        THREADS_smaf: # Number of threads used (default=4)
        doprint: Print correlations between predicted X and real X and some additional info
        normalization: How data is normalized before running smaf (default: paper_norm)
                       Options: paper_norm (normalization used in paper, protein-wise),
                       min_max_norm (protein-wise), zscore_norm (protein-wise) or none
        others: All other parameters belong to the fnc. call of lasso from spams and
                further information is available on the respective website
                (maxItr, UW, posW, posU, use_chol, module_lower, activity_lower,
                donorm, mode, mink, U0, U0_delta)

    For fnc compute_A():
        nmeasurements: number of channels (# rows of A) (default: 10)
        maxcomposition: maximum times each gene is represented (mode G),
                        or max genes per composition (mode M)
                        (default: 3)
        mode_phi: 'M', max genes per composition is constrained
                  'G', maximum times each gene is represented is constrained (default)
        THREADS_A: # Number of threads used (default=20)
        num_phi: the number of best phi's that are returned/saved (default: 1, at most: 50)
        binary: boolean variable that specifies if Phi is binary or not (default: False)
                If True, then Phi will be binarized directly after computation and the
                best Phi is chosen according to the results of its binarized version

    For fnc analyze_U_and_A():
        THREADS_A_and_U: # Number of threads used (default=20)

    For fnc train_validate_test_split():
        split_by: either split by 'roi' or 'percentage' (default: 'roi')
        k_cv: number k cross-validations (default: 4)
              if split_by='roi, k_cv needs to be smaller and a multiple of the
              number of rois (not including the test roi)
        Only one of the next two parameters needs to be set!
        test_set: tuple of rois used as test sets with the same names as in
                  the anndata object column sample_id
                  (if split_by='roi', then test_set must be set and test_size
                  can't be used)
        test_size: size of test set (eg. 0.2)

    For fnc. make_cond_prob():
        threshold_cond_prob: Threshold to stratify if particular protein in particular
                             cell is counted as positive (default: 10.0)


outputs:
    training_res: results of training evaluated on test set (pandas dataframe)
    cor: pairwise protein correlations (numpy array)
    cond_prob: pairwise conditional probabilities between proteins/channels
               (numpy array)

'''


def train_U_and_A(X, outpath, layer=None, d = 80, lda1 = 3, lda2 = 0.2, maxItr=10,
                  UW=None, posW=False, posU=True,use_chol=False, module_lower=1,
                  activity_lower=1, donorm=False, mode_smaf=1, mink=0, U0=[],
                  U0_delta=0.1, doprint=False, normalization='paper_norm',
                  THREADS_smaf=4, nmeasurements = 10, maxcomposition = 3, mode_phi='G',
                  lasso_sparsity=0.2, THREADS_A=20, num_phi=1, binary=False,
                  THREADS_A_and_U=20, split_by='roi', k_cv=4, test_set=(), test_size=None,
                  threshold_cond_prob=10.0):

    # Add a seed to use by numpy for reproducibility
    np.random.seed(11)

    # Create test, validate and training sets for k_cv times crossvalidation
    X_test, remaining = train_validate_test_split(X, split_by, k_cv, test_set, test_size)

    ## Start training
    # Do crossvalidation to determine best U, Phi/A
    U_best = None
    Phi_best = None
    k_best = None
    version_best = None
    all_pearson_cor = [0] * k_cv
    for k in range(k_cv):
        # Initialize progress bar
        with alive_bar(4, title="Fold {0}".format(str(k)), force_tty=True) as bar:
            # Split into training and validating sets
            X_validate = X[X.obs.index.isin(remaining[k]), ]
            X_training = X[X.obs.index.isin(np.concatenate(remaining[:k]+remaining[(k+1):])), ]
            # Mark step for progressbar
            bar()

            # Compute U
            U, W = smaf(X_training, d, lda1, lda2, maxItr, UW, posW, posU, use_chol,
                        module_lower, activity_lower, donorm, mode_smaf, mink, U0, U0_delta,
                        doprint, THREADS_smaf, os.path.join(outpath, 'crossvalidation_'+str(k)),
                        normalization, layer)
            # Mark step for progressbar
            bar()

            # Compute A
            Phi, pearson_cor, versions = compute_A(X_validate, U, nmeasurements, maxcomposition, 
                                                  mode_phi, lasso_sparsity, 
                                                  os.path.join(outpath, 'crossvalidation_'+str(k)),
                                                  THREADS_A, layer, num_phi, binary)
            # Mark step for progressbar
            bar()

            # Keep best U, Phi and all pearson correlations of training data
            if pearson_cor > max(cor for cor in all_pearson_cor if cor is not None):
                U_best = U
                Phi_best = Phi
                k_best = k
                version_best
            all_pearson_cor[k] = pearson_cor

            # Mark step for progressbar
            bar()

    # Print training correlations
    print(('Cross-fold validations done.\n' +
          'The pearson correlations are: {0}'.format(all_pearson_cor)))

    # Analize training
    training_res, training_res_comp = analyze_U_and_A(X_test, U_best, [Phi_best],
                                                      [version_best], outpath, k_best, 
                                                      lasso_sparsity,
                                                      THREADS_A_and_U)

    # Calculate pairwise protein correlations
    cor = make_cor(X, outpath)
    # Compute pairwise conditional probabilities between proteins/channels
    cond_prob = make_cond_prob(X, outpath, threshold_cond_prob)

    return training_res, training_res_comp, U_best, Phi_best, X_test.obs.index


# Function splitting X into training, validate and test either by 'roi' or 'percentage'
# If one or multiple specific rois should be used as test cases, then their names
# should be specified in test_set as a tuple with the same names as in the anndata
# object column sample_id
# k_cv is the number of crossvalidations (if split_by='roi': k_cv needs to be smaller
# and a multiple of the number of rois not including the test roi)
# Only one of the two test_set and test_size (eg. 0.7) needs to be set
# (if split_by='roi': only option test_set can be used)
def train_validate_test_split(X, split_by='roi', k_cv=4, test_set=(), test_size=None):
    # Split into training, validate and test set
    if split_by=='roi':

        # Assert that test_set is not empty
        if test_set == ():
            raise ValueError(('split_by is set to "roi", but no test_set is given. ' +
                              'Please specify a tuple of rois used as test sets '  +
                              'in test_set.'))

        # Get all roi names for indexing
        roi_names = X.obs['sample_id'].unique().to_list()

        # Get all roi names not used for testing and assert that k_cv is smaller
        # and a multiple of the number of remaining rois
        i_remaining = list(set(roi_names) - set(test_set))
        if (len(i_remaining)<k_cv) or (k_cv%len(i_remaining)!=0):
            raise ValueError(('The number of crossvalidations k_cv ({0}) '.format(k_cv) +
                              'is not smaller or a multiple of the number of ' +
                              'remaining rois!'))

        X_test = X[X.obs['sample_id'].isin(test_set), ]
        remaining_rois = np.array_split(np.random.permutation(np.array(i_remaining)), k_cv)

        # Convert sets of rois to sets of row names of according sets
        remaining = [None] * len(remaining_rois)
        for i in range(len(remaining_rois)):
            remaining[i] = X.obs.index[X.obs["sample_id"].isin(remaining_rois[i])].to_numpy()

    elif split_by=='percentage':
        # If test_set is set, then the specified test roi(s) will be set aside for
        # testing (test_size will not be used)
        if test_set:
            # Get all roi names for indexing
            roi_names = X.obs['sample_id'].unique().to_list()

            # Index test set and remaining sets
            i_remaining = list(set(roi_names) - set(test_set))
            X_test = X[X.obs['sample_id'].isin(test_set), ]
            X_remaining = X[~X.obs["sample_id"].isin(test_set), ]

            # Split remainder according to k_cv
            remaining = np.array_split(np.random.permutation(X_remaining.obs.index.to_numpy()), k_cv)

        else:
            # Get all obs names for indexing
            obs_names = X.obs.index.to_numpy()
            # Index test set according to test_size
            test_names = np.random.choice(obs_names, int(test_size*len(obs_names)),
                                          replace=False)
            X_test = X[X.obs.index.isin(test_names), ]

            # Split remainder according to k_cv
            i_remaining = list(set(obs_names) - set(test_names))
            remaining = np.array_split(np.random.permutation(i_remaining), k_cv)

    else:
        # If it is none of the above splitting methods, return ValueError
        raise ValueError(('No valid splitting method for X was given. ' +
                         '<split_by> needs to be either "roi" or "percentage".'))

    return X_test, remaining
