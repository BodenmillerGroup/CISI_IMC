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


# Test SMAF
class TestCompositeA(unittest.TestCase):

    '''
    def test_A(self):
        A = compute_random_compositions.compute_A(
            X=np.array(np.load('data/training_data_sce.npy')),
            U=np.array(np.load('data/gene_modules.npy')),
            nmeasurements=10, maxcomposition=3)
        self.assertFalse(np.all((A[1] == 0)))
    '''

    def test_A_with_output(self):
        p = ['Histone H3', 'CD45RA', 'CD20', 'CD11b', 'CD56 ', 'GranzymeB', 'CD3',
             'DC-LAMP', 'CD11c', 'PD-1', 'GITR', 'SMA', 'MMP9', 'CD14', 'PD-L1',
             'TCF1/TCF7', 'CD45RO', 'FOXP3', 'ICOS', 'Ki-67', 'CD8a', 'Tim-3',
             'E-Cad/P-Cad', 'IRF4', 'VISTA', 'CD1c', 'CD4', 'CD31', 'CXCL12',
             'CCL21', 'panCK', 'CXCL13', 'Ir191', 'CCR7', 'Ir193', 'Vimentin',
             'CD15', 'MPO', 'CD38', 'HLA-DR', 'CD27', 'CD303', 'CD68']
        A = compute_random_compositions.compute_A(
            X=np.array(np.load('data/training_data_sce.npy')),
            U=np.array(np.load('data/gene_modules.npy')),
            nmeasurements=10, maxcomposition=3,
            outpath='data', proteins=p)
        self.assertTrue(check_file('data/compositions_A/version_49.txt', ['.txt']))


## Main
if __name__ == '__main__':
    unittest.main()
