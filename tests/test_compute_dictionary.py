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


# Test inputs to compute_dictionary
class TestInputX(unittest.TestCase):

    def test_csv_validity(self):
        t = compute_dictionary.check_input('data/cell.csv', [".csv"])
        self.assertTrue(t)

    def test_npy_validity(self):
        t = compute_dictionary.check_input('data/training_data_sce.npy', [".npy"])
        self.assertTrue(t)

    def test_csv(self):
        t = compute_dictionary.read_input('data/cell.csv')[:5, 0]
        self.assertIsNone(np.testing.assert_array_equal(t, np.array([0.000158776,
                                                                     0.000232407,
                                                                     0.000186816,
                                                                     0.000181502,
                                                                     0.000392305])))

    def test_input_error(self):
        with self.assertRaises(ValueError):
            compute_dictionary.read_input('dat')


# Test SMAF
class TestSMAF(unittest.TestCase):

    def test_smaf_csv(self):
        U, W = compute_dictionary.smaf('data/cell.csv', 80, 3, 0.2)
        self.assertFalse(np.all((U == 0)))








## Main
if __name__ == '__main__':
    unittest.main()
