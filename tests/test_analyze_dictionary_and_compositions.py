# Import system libraries to configure code directory as module
from os.path import dirname, abspath, join
import sys

# Find code directory relative to our directory
THIS_DIR = dirname(__file__)
CODE_DIR = abspath(join(THIS_DIR, '..', 'code'))
# Add code directory to systems paths
sys.path.append(CODE_DIR)

# Import libraries
import unittest
import numpy as np
from analyze_dictionary_and_compositions import analyze_U_and_A
from utils import check_file
import os
import anndata as ad


# Test inputs to compute_dictionary
class TestPerformanceUandA(unittest.TestCase):

    def test_analyze_dictionary_and_composition(self):
        P = [np.loadtxt(os.path.join('data/compositions_A/version_16.txt'),
                        skiprows=1, usecols=list(range(2, 45)))]
        results = analyze_U_and_A(X_input=ad.read_h5ad('data/test.h5ad'),
                                  U=np.array(np.load('data/gene_modules.npy')),
                                  Phi=P, outpath='data')
        print(results)
        self.assertTrue(check_file('data/simulation_results.txt', [".txt"]))


## Main
if __name__ == '__main__':
    unittest.main()
