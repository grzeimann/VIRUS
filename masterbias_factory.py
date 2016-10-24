# -*- coding: utf-8 -*-
"""

Science script for "quick" reduction.

The script requires an output directory.  It also requires key information
for selecting the science and calibration directories.  
It is required to input the date and obsid for the science folders

@author: gregz
"""

from __future__ import print_function
import argparse as ap
from astropy.io import fits
import numpy as np
import pandas as pd
from pandas import DataFrame as DF
import matplotlib.pyplot as plt
import glob
import os
import os.path as op
import textwrap
import re
from utils import biweight_location

_dataframe = DF()    
AMPS = ["LL", "LU", "RL", "RU"]
SPECID = ["004", "008", "012", "013", "016", "017", "020", "024", "025", "027",
          "032", "037", "038", "041", "047", "051"]


def parse_args(argv=None):
    """Parse the command line arguments

    Parameters
    ----------
    argv : list of string
        arguments to parse; if ``None``, ``sys.argv`` is used

    Returns
    -------
    Namespace
        parsed arguments
    """
    description = textwrap.dedent('''Make MasterBias Script - 
    
                     This scripts create, examines, and compares the following:
                         masterbias frames for evolution in the bias

                     The script places the files in the given output directory.
                     
                     ''')
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.RawTextHelpFormatter)
                        
    parser.add_argument("--specid", nargs='?', type=str, 
                        help='''List of SPECID's for processing. 
                        Ex: "008,012,027"''', default = None)

    parser.add_argument("--amps", nargs='?', type=str, 
                        help='''List of SPECID's for processing. 
                        Ex: "LL,RL"
                        Default: "LL,LU,RL,RU"''', default = None)

    parser.add_argument("--output", nargs='?', type=str, 
                        help='''Output Directory [REQUIRED]''', 
                        default=None)

    parser.add_argument("--cal_dirs", nargs='?', type=str, 
                        help='''Calibration Directory [REQUIRED]
                        Needs to be a string in quotes if wildcard used.''', 
                        default=None)
                                                                         
    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.specid:
        args.specid = args.specid.replace(" ", "").split(',')
    else:
        args.specid = SPECID

    # Check that the arguments are filled
    if args.cal_dirs:
        nargs = args.cal_dirs.replace(" ", "").split('*')
        if len(nargs)>1:
            args.cal_dirs = args.cal_dirs.replace(" ", "")
            args.cal_dirs = glob.glob(args.cal_dirs)
        else:
            args.cal_dirs = args.cal_dirs.replace(" ", "").split(',')  
    else:
        msg = 'No calibration directory was provided'
        parser.error(msg)
         
    if args.output is None:
        msg = 'No output directory was provided'
        parser.error(msg) 
    else:
        if not op.exists(args.output):
            print("Making directory: {:s}".format(args.output))
            os.mkdir(args.output)
    
    if args.amps:
        args.amps = args.amps.replace(" ", "").split(',')
    else:
        args.amps = AMPS   
            
    return args
    

def build_dataframe(_dataframe, date, fn, args):
    F = fits.open(fn)
    specid = F[0].header['SPECID']
    if specid in args.specid:
        blank, txlow, txhigh, tylow, tyhigh, blank = re.split('[: \[ \] ,]', 
                                                            F[0].header['TRIMSEC'])
        txlow                = int(txlow)-1
        txhigh               = int(txhigh)
        tylow                = int(tylow)-1
        tyhigh               = int(tyhigh)
        blank, bxlow, bxhigh, bylow, byhigh, blank = re.split('[: \[ \] ,]', 
                                                            F[0].header['BIASSEC'])
        bxlow                = int(bxlow)-1
        bxhigh               = int(bxhigh)
        bylow                = int(bylow)-1
        byhigh               = int(byhigh)    
        over = biweight_location(F[0].data[bylow:byhigh,bxlow:bxhigh])
        amp = (F[0].header['CCDPOS'].replace(" ", "") 
               + F[0].header['CCDHALF'].replace(" ", ""))
        A = {'filename' : pd.Series(fn, index=[date+'T'+F[0].header['UT']]),
             'TRIM_XL' : pd.Series(txlow, index=[date+'T'+F[0].header['UT']]),
             'TRIM_XH' : pd.Series(txhigh, index=[date+'T'+F[0].header['UT']]),
             'TRIM_YL' : pd.Series(tylow, index=[date+'T'+F[0].header['UT']]),
             'TRIM_YH' : pd.Series(tyhigh, index=[date+'T'+F[0].header['UT']]),
             'BIAS_XL' : pd.Series(bxlow, index=[date+'T'+F[0].header['UT']]),
             'BIAS_XH' : pd.Series(bxhigh, index=[date+'T'+F[0].header['UT']]),
             'BIAS_YL' : pd.Series(bylow, index=[date+'T'+F[0].header['UT']]),
             'BIAS_YH' : pd.Series(byhigh, index=[date+'T'+F[0].header['UT']]),
             'overscan' : pd.Series(over, index=[date+'T'+F[0].header['UT']]),
             'SPECID' : pd.Series(specid, index=[date+'T'+F[0].header['UT']]),
             'AMP' : pd.Series(amp, index=[date+'T'+F[0].header['UT']])} 
        data = DF(A)
        print(data)
        _dataframe.append(data)
    
    
def main():
    args = parse_args()
    print(args.specid)
    lower_folder_struct = op.join('virus','virus*','exp*','virus','2*zro.fits')
    for date in args.cal_dirs:
        files = glob.glob(op.join(date,lower_folder_struct))
        for fn in files:
            build_dataframe(_dataframe, op.basename(date), fn, args)
    print(_dataframe)
            

   
if __name__ == '__main__':
    main() 