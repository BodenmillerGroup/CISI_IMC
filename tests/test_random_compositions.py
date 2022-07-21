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
import compute_random_compositions
from utils import check_file
import anndata as ad


# Test SMAF
class TestCompositeA(unittest.TestCase):

    def test_A(self):
        A = compute_random_compositions.compute_A(
            X_input=ad.read_h5ad('data/test.h5ad'),
            U=np.array(np.load('data/gene_modules.npy')),
            nmeasurements=10, maxcomposition=3, outpath='data')
        self.assertFalse(np.all((A == 0)))

    def test_A_binary(self):
        A = compute_random_compositions.compute_A(
            X_input=ad.read_h5ad('data/test.h5ad'),
            U=np.array(np.load('data/gene_modules.npy')),
            nmeasurements=10, maxcomposition=3, outpath='data', binary=True)
        self.assertFalse(np.all((A == 0)))


## Main
if __name__ == '__main__':
    unittest.main()
