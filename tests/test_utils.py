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
import utils


# Test inputs to compute_dictionary
class TestInputX(unittest.TestCase):

    def test_csv_validity(self):
        t = utils.check_file('data/cell.csv', [".csv"])
        self.assertTrue(t)

    def test_npy_validity(self):
        t = utils.check_file('data/training_data_sce.npy', [".npy"])
        self.assertTrue(t)

    def test_X_input(self):
        X, X_header, proteins = utils.read_input('data/cell.csv', 'data/TH152_panel.csv')
        self.assertIsNone(np.testing.assert_array_equal(X[2, :5],
                                                        np.array([0.000158776,
                                                                  0.000232407,
                                                                  0.000186816,
                                                                  0.000181502,
                                                                  0.000392305])))

    def test_X_header(self):
        t = ['ImageNumber', 'ObjectNumber', '1', '10', '11', '12', '13', '14', '15',
             '16', '17', '18', '19', '2', '20', '21', '22', '23', '24', '25', '26',
             '27', '28', '29', '3', '30', '31', '32', '33', '34', '35', '36', '37',
             '38', '39', '4', '40', '41', '42', '43', '5', '6', '7', '8', '9']
        X, X_header, proteins = utils.read_input('data/cell.csv', 'data/TH152_panel.csv')
        self.assertListEqual(X_header, t)

    def test_proteins(self):
        t = ['Histone H3', 'CD45RA', 'CD20', 'CD11b', 'CD56 ', 'GranzymeB', 'CD3',
             'DC-LAMP', 'CD11c', 'PD-1', 'GITR', 'SMA', 'MMP9', 'CD14', 'PD-L1',
             'TCF1/TCF7', 'CD45RO', 'FOXP3', 'ICOS', 'Ki-67', 'CD8a', 'Tim-3',
             'E-Cad/P-Cad', 'IRF4', 'VISTA', 'CD1c', 'CD4', 'CD31', 'CXCL12',
             'CCL21', 'panCK', 'CXCL13', 'Ir191', 'CCR7', 'Ir193', 'Vimentin',
             'CD15', 'MPO', 'CD38', 'HLA-DR', 'CD27', 'CD303', 'CD68']
        X, X_header, proteins = utils.read_input('data/cell.csv', 'data/TH152_panel.csv')
        self.assertListEqual(proteins, t)

    def test_input_error(self):
        with self.assertRaises(ValueError):
            utils.read_input('dat', 'data/TH152_panel.csv')


# Test inputs to compute_dictionary
class TestSplitX(unittest.TestCase):

    def test_split_roi(self):
        X_test = np.array([[0, 1, 2, 3, 2],
                           [0, 0, 0, 0, 1],
                           [1, 2, 3, 4, 5],
                           [6, 7, 8, 9, 10]])
        X_header_test = ['ImageNumber', 'ObjectNumber', '1', '10']
        X_training, X_validate, X_test = utils.split_X(X_test, X_header_test,
                                                       by='roi',
                                                       set_sizes=[[2, 3], [0], [1]])
        self.assertTrue(X_training.shape==(2, 3) and
                          X_validate.shape==(2, 1) and
                          X_test.shape==(2, 1))

    def test_split_percent(self):
        X_test = np.array([[0, 1, 2, 3, 2],
                           [0, 0, 0, 0, 1],
                           [1, 2, 3, 4, 5],
                           [6, 7, 8, 9, 10]])
        X_header_test = ['ImageNumber', 'ObjectNumber', '1', '10']
        X_training, X_validate, X_test = utils.split_X(X_test, X_header_test,
                                                       by='percentage',
                                                       set_sizes=[0.5, 0.25, 0.25])
        self.assertTrue(X_training.shape==(2, 3) and
                          X_validate.shape==(2, 1) and
                          X_test.shape==(2, 1))


## Main
if __name__ == '__main__':
    unittest.main()
