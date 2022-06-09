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


# Test SMAF
class TestSMAF(unittest.TestCase):

    def test_smaf_csv(self):
        U, W = compute_dictionary.smaf(np.array(np.load('data/training_data_sce.npy')),
                                       80, 3, 0.2)
        self.assertFalse(np.all((U == 0)))

    def test_smaf_csv_with_output(self):
        U, W = compute_dictionary.smaf(np.array(np.load('data/training_data_sce.npy')),
                                       80, 3, 0.2, outpath='data')
        self.assertTrue(check_file('data/gene_modules.csv', ['.csv']) and
                        check_file('data/gene_modules.npy', ['.npy']))


## Main
if __name__ == '__main__':
    unittest.main()
