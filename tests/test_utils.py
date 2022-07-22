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
class TestCheckFile(unittest.TestCase):

    def test_check_file(self):
        t = utils.check_file('data/test.h5ad', [".h5ad"])
        self.assertTrue(t)


## Main
if __name__ == '__main__':
    unittest.main()
