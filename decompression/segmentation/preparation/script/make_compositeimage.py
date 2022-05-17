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
tissue = "r2_4"
panel = "TH108_panel.csv"
nms = 10
nmax = 2

#%%
cwd = os.getcwd()
trd = cwd.replace("script","training")
dtd = os.path.join(trd,"data")
imd = cwd.replace("script","individual_image")
opd = cwd.replace("script","composite_image")
cimd = "%s/%s" % (imd,tissue)

f = open('%s/selected_genes.txt' % dtd)
genes = np.array([line.strip() for line in f])

tags = pd.read_csv("%s/%s" % (dtd,panel))

phid = os.path.join(trd,"result","%d_measurements","%d_max") % (nms,nmax)
phi = np.load('%s/phi.npy' % phid)

images = os.listdir(cimd)


#%%
for m in range(phi.shape[0]):
    x = phi[m].astype(bool)
    targets = genes[x]

    z = tags[["Target","Metal Tag"]]
    w = dict(z.values)
    for k,v in w.items():
        if k in targets:
            isotopes = [v for k, v in w.items() if k in targets]

    cimages = []
    for j in isotopes:
        y = [i for i in images if j  in i]
        cimages =  cimages + y

    cmpia = 0
    for l in range(len(cimages)):
        ci = Image.open("%s/%s" % (cimd,cimages[l]))
        cia = asarray(ci)
        cmpia = cmpia + cia
    cmpi = Image.fromarray(cmpia).save("%s/%s_composite_%d.tiff" % (opd,tissue,m))


