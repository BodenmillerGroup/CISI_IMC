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
import compute_dictionary
from utils import check_file
import anndata as ad

# Test SMAF
class TestSMAF(unittest.TestCase):

    test_data = ad.read_h5ad('data/test.h5ad')

    def test_smaf(self):
        U, W = compute_dictionary.smaf(self.test_data,
                                       80, 3, 0.2)
        self.assertFalse(np.all((U == 0)))

    def test_smaf_with_output(self):
        U, W = compute_dictionary.smaf(self.test_data,
                                       80, 3, 0.2, outpath='data')
        self.assertTrue(check_file('data/gene_modules.csv', ['.csv']) and
                        check_file('data/gene_modules.npy', ['.npy']))

    def test_smaf_with_layer(self):
        U, W = compute_dictionary.smaf(self.test_data,
                                       80, 3, 0.2, layer='exprs')
        self.assertFalse(np.all((U == 0)))

    def test_smaf_with_zscore(self):
        U, W = compute_dictionary.smaf(self.test_data,
                                       80, 3, 0.2, normalization='zscore_norm')
        self.assertFalse(np.all((U == 0)))

    def test_smaf_with_min_max(self):
        U, W = compute_dictionary.smaf(self.test_data,
                                       80, 3, 0.2, normalization='min_max_norm')
        self.assertFalse(np.all((U == 0)))

    def test_smaf_with_min_max_channelwise(self):
        U, W = compute_dictionary.smaf(self.test_data,
                                       80, 3, 0.2, normalization='min_max_norm_channelwise')
        self.assertFalse(np.all((U == 0)))

## Main
if __name__ == '__main__':
    unittest.main()
