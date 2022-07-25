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
import train_dictionary_and_compositions
from utils import check_file
import anndata as ad


# Test inputs to train_validate_test_split
class TestSplitX(unittest.TestCase):

    test_data = ad.read_h5ad('data/test.h5ad')

    def test_split_by_roi(self):
        test, remaining = train_dictionary_and_compositions.train_validate_test_split(
            self.test_data, split_by='roi', k_cv=4, test_set=('20220520_TsH_th152_cisi1_001',))
        self.assertTrue(np.all(test.obs['sample_id']=='20220520_TsH_th152_cisi1_001') and
                        not any('20220520_TsH_th152_cisi1_001' in row for row in remaining))

    def test_split_by_percentage(self):
        test, remaining = train_dictionary_and_compositions.train_validate_test_split(
            self.test_data, split_by='percentage', k_cv=2, test_set=('20220520_TsH_th152_cisi1_002',))
        self.assertTrue(np.all(test.obs['sample_id']=='20220520_TsH_th152_cisi1_002') and
                        not any('20220520_TsH_th152_cisi1_002' in row for row in remaining) and
                        len(remaining)==2)

    def test_split_by_percentage_test_size(self):
        test, remaining = train_dictionary_and_compositions.train_validate_test_split(
            self.test_data, split_by='percentage', k_cv=4, test_size=0.2)
        self.assertTrue(test.X.shape[0]==int(1000*0.2) and
                        len(remaining)==4)

# # Test train_U_and_A
# class TestTraining(unittest.TestCase):
#
#     test_data = ad.read_h5ad('data/test.h5ad')
#
#     def test_train_U_and_A(self):
#         train_dictionary_and_compositions.train_U_and_A(self.test_data, 'data/training',
#                                                         test_set=('20220520_TsH_th152_cisi1_001',))
#         self.assertTrue(check_file('data/training/simulation_results.txt', ['.txt']))



## Main
if __name__ == '__main__':
    unittest.main()
