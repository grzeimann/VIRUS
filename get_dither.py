# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 20:07:35 2016

@author: gregz
"""

import warnings
import numpy as np
import os.path as op
from utils import biweight_location

ucam = ["004", "008", "012", "013",  "016", "017", "020", "024", "025", "027",
        "032", "037", "038", "041", "047", "051"]
ucam = ["013","016","017","027","041"]        
fields = ["HF10", "HF11", "HF12", "HF15", "HF16", "HF17", "HF23", "HF47",
          "HF49", "HF50", "HF55", "HF56"]

mnx1 = 0
mnx2 = 0
mnx3 = 0
mny1 = 0.
mny2 = 0
mny3 = 0

thresh = 3
dither1 = []
dither2 = []
dither3 = []
for field in fields:
    matches = []        
    for specid in ucam:
        folder = op.join('/work/00115/gebhardt/maverick/sci', field, 'c'+specid)
        D1 = op.join(folder, 'd1_cont.dat')
        D2 = op.join(folder, 'd2_cont.dat')
        D3 = op.join(folder, 'd3_cont.dat')
        if not op.exists(D1):
            continue
        if not op.exists(D2):
            continue
        if not op.exists(D3):
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cat1 = np.loadtxt(D1)
            cat2 = np.loadtxt(D2)
            cat3 = np.loadtxt(D3)
        if not cat1.size:
            continue
        if not cat2.size:
            continue
        if not cat3.size:
            continue
        if cat1.ndim < 2:
            cat1 = cat1[np.newaxis,:]
        if cat2.ndim < 2:
            cat2 = cat2[np.newaxis,:]
        if cat3.ndim < 2:
            cat3 = cat3[np.newaxis,:]
        for i, icx in enumerate(cat1[:,1]):
            for j, jcx in enumerate(cat2[:,1]):
                for k, kcx in enumerate(cat3[:,1]):
                    if (np.abs(icx-jcx+mnx1)<thresh and np.abs(icx-kcx+mnx2)<thresh 
                                 and np.abs(jcx-kcx+mnx3)<thresh):
                         if (np.abs(cat1[i,2]-cat2[j,2]+mny1)<thresh 
                                 and np.abs(cat1[i,2]-cat3[k,2]+mny2)<thresh 
                                 and np.abs(cat2[j,2]-cat3[k,2]+mny3)<thresh):
                             matches.append(np.array([icx, jcx, kcx, cat1[i,2],
                                                      cat2[j,2], cat3[k,2]]))
                             #print(icx, jcx, kcx, cat1[i,2],
                             #                         cat2[j,2], cat3[k,2])
                         else:
                             continue
                    else:
                        continue
    print(matches)
    if not matches:
        continue
    matches = np.vstack(matches)
    mnx = biweight_location(matches[:,:3],axis=(1,))
    mny = biweight_location(matches[:,3:],axis=(1,))
    matches[:,:3] = matches[:,:3] - mnx[:,np.newaxis]
    matches[:,3:] = matches[:,3:] - mny[:,np.newaxis]
    print(matches[:,0:1] - matches[:,:3])
    print(matches[:,3:4] - matches[:,3:])
    dx = biweight_location(matches[:,:3],axis=(0,))
    dy = biweight_location(matches[:,3:],axis=(0,))
    dither1.append(np.array([dx[0]-dx[0],dy[0]]-dy[0]))
    dither2.append(np.array([dx[0]-dx[1],dy[0]]-dy[1]))
    dither3.append(np.array([dx[0]-dx[2],dy[0]]-dy[2]))
    print('''Field: %s
                   D1: (%0.2f, %0.2f)
                   D2: (%0.2f, %0.2f)
                   D3: (%0.2f, %0.2f)
          ''' % (field, dx[0]-dx[0], dy[0]-dy[0], dx[0] - dx[1], dy[0]-dy[1], 
                 dx[0]-dx[2], dy[0]-dy[2]))
    
    
                        