# Import libraries
import numpy as np
import spams


'''
Helper file
(Contains general helper files used in more than one file)
'''


# Fixing U, approximating W
def sparse_decode(Y, D, k, numThreads, worstFit=1., mink=0):
    while k > mink:
        W = spams.omp(np.asfortranarray(Y), np.asfortranarray(D), L=k, numThreads=numThreads)
        W = np.asarray(W.todense())
        fit = 1 - np.linalg.norm(Y - D.dot(W))**2 / np.linalg.norm(Y)**2
        if fit < worstFit:
            break
        else:
            k -= 1
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
