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
from decompress import decompress
from utils import check_file, simulate_composite_measurements
import os
import anndata as ad


# Test inputs to analyze_U_and_A
class TestDecompression(unittest.TestCase):

    def test_decompression(self):
        phi = np.loadtxt(os.path.join('data/compositions_A/version_16.txt'),
                         skiprows=1, usecols=list(range(2, 45)))
        X = ad.read_h5ad('data/test.h5ad')[:50, ]
        y = simulate_composite_measurements(X, phi)

        results = decompress(y, np.array(np.load('data/gene_modules.npy')), phi)
        self.assertFalse(np.all((results < 0)))


## Main
if __name__ == '__main__':
    unittest.main()
