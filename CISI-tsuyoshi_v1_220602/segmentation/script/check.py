#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 15:07:38 2022

@author: tsuyoshi
"""
import numpy as np


z = np.load("/Users/tsuyoshi/Documents/prj_CISI/CISI-test_melanoma_220126/preprocessing/FilteredCompositeImages/tissue1/segmented/cell_masks.npz")
#x = np.load("/Users/tsuyoshi/Documents/prj_CISI/CISI-test_melanoma_220126/training/data/relative_abundance.npy")
#y = np.load("/Users/tsuyoshi/Documents/prj_CISI/CISI-test_melanoma_220126/training/data/average_on_expression_level.npy")

x = z.files
y = z['data']
