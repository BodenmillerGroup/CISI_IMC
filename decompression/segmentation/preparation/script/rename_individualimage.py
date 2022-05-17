#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 13:38:48 2022

@author: tsuyoshi
"""

##import imageio
import matplotlib.pyplot as plt
import PIL
from PIL import Image
import os
from numpy import asarray
import pandas as pd
import numpy as np
#%%
# input zone
tissue_post = list(range(1,6))
tissue_pre = ["r1_1","r1_2","r2_2","r2_3","r2_4"]
panel = "TH108_panel.csv"

#%%
cwd = os.getcwd()
imd = cwd.replace("script","individual_image")
trd = cwd.replace("script","training")
tags = pd.read_csv("%s/data/%s" % (trd,panel))
z = tags[["Metal Tag","Target"]]
isotope_target_dict = dict(z.values)

f = open('%s/data/selected_genes.txt' % trd)
genes = np.array([line.strip() for line in f])

#%%
for i in tissue_post:
    os.mkdir(os.path.join(imd,"tissue%d"%i))
    os.mkdir(os.path.join(imd,"tissue%d"%i,"original"))
    opd = os.path.join(imd,"tissue%d"%i,"original")
    cimd = "%s/%s" % (imd,tissue_pre[i-1])
    images = os.listdir(cimd) 
    for ci in images:
        target = isotope_target_dict.get(ci.split('_')[-1].split('.')[0])
        newname = "%s.tiff" % target
        if target in genes:
            os.rename(os.path.join(cimd,ci),os.path.join(opd,newname))
        














