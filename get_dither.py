# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 20:07:35 2016

@author: gregz
"""

import numpy as np
import os.path as op
from utils import biweight_location

ucam = ["004", "008", "012", "013",  "016", "017", "020", "024", "025", "027",
        "032", "037", "038", "041", "047", "051"]
        
fields = ["HF10", "HF11", "HF12", "HF15", "HF16", "HF17", "HF23", "HF47",
          "HF49", "HF50", "HF55", "HF56"]

thresh = 3.
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
        cat1 = np.loadtxt(D1)
        cat2 = np.loadtxt(D1)
        cat3 = np.loadtxt(D1)
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
                    if (np.abs(icx-jcx)<thresh and np.abs(icx-kcx)<thresh 
                                 and np.abs(jcx-kcx)<thresh):
                         if (np.abs(cat1[i,2]-cat2[j,2])<thresh 
                                 and np.abs(cat1[i,2]-cat3[k,2])<thresh 
                                 and np.abs(cat2[j,2]-cat3[k,2])<thresh):
                             matches.append(np.array([icx, jcx, kcx, cat1[i,2],
                                                      cat2[j,2], cat3[k,2]]))
                         else:
                             continue
                    else:
                        continue
    matches = np.vstack(matches)
    mnx = biweight_location(matches[:,:3],axis=(1,))
    mny = biweight_location(matches[:,3:],axis=(1,))
    matches[:,:3] -= mnx
    matches[:,3:] -= mny
    dx = biweight_location(matches[:,:3],axis=(0,))
    dy = biweight_location(matches[:,3:],axis=(0,))
    dither1.append(np.array([dx[0],dy[0]]))
    dither2.append(np.array([dx[1],dy[1]]))
    dither3.append(np.array([dx[2],dy[2]]))
    print('''Field: %s
                   D1: (%0.2f, %0.2f)
                   D2: (%0.2f, %0.2f)
                   D3: (%0.2f, %0.2f)
          ''' % (field, dx[0], dy[0], dx[1], dy[1], dx[2], dy[2]))
    
    
                        