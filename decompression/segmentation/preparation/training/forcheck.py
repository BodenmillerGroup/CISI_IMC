#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 19:05:08 2022

@author: tsuyoshi
"""

#%%

import numpy as np
import argparse
from collections import defaultdict

f = open("/Users/tsuyoshi/Documents/prj_CISI/CISI-test_melanoma_220126/training/data/selected_genes.txt")
Genes = [line.strip() for line in f]
f.close()
f = open("/Users/tsuyoshi/Documents/prj_CISI/CISI-test_melanoma_220126/training/result/10_measurements/3_max/compositions_v48_2.csv")
header = f.readline()
Phi_dict = defaultdict(list)
for line in f:
   		ls = line.strip().split(',') 
   		Phi_dict[int(ls[0])].append(ls[1])
Phi = np.zeros((len(Phi_dict),len(Genes)))
for k,v in Phi_dict.items():
   		gi = [Genes.index(g) for g in v]
   		Phi[k,gi] = 1

#%%

f = open("/Users/tsuyoshi/Documents/prj_CISI/CISI-test_melanoma_220126/training/result/10_measurements/3_max/compositions_v48_3.csv")

lines = [line.strip().split(',') for line in f]


gt=[Genes.index(g) for g in ["SMA","gp100"]]
gt=[Genes.index(g) for g in ['HLA-DR', 'S100', 'Glut1', 'CD8a', 'caveoli', 'MPO']]

Pt = np.zeros((len(Phi_dict),len(Genes)))
Pt[0,gt] = 1
x=Phi_dict[1]
x = range(1,3)
ls[1:]

#%%
import numpy as np
import argparse
from collections import defaultdict

f = open("/Users/tsuyoshi/Documents/prj_CISI/CISI-test_melanoma_220126/training/data/selected_genes.txt")
Genes = [line.strip() for line in f]
f.close()
f = open("/Users/tsuyoshi/Documents/prj_CISI/CISI-test_melanoma_220126/training/result/10_measurements/3_max/compositions_v48_2.csv")
header = f.readline()
Phi_dict = defaultdict(list)
for line in f:
   		ls = line.strip().split(',')
   		ls = ' '.join(ls).split()
   		Phi_dict[int(ls[0])].extend(ls[1:])
Phi = np.zeros((len(Phi_dict),len(Genes)))
for k,v in Phi_dict.items():
   		gi = [Genes.index(g) for g in v]
   		Phi[k,gi] = 1