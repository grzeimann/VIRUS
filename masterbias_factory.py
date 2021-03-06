# -*- coding: utf-8 -*-
"""

Science script for "quick" reduction.

The script requires an output directory.  It also requires key information
for selecting the science and calibration directories.  
It is required to input the date and obsid for the science folders

@author: gregz
"""

from __future__ import print_function
import matplotlib
matplotlib.use('agg')
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
import warnings
from utils import biweight_location, is_outlier
from progressbar import ProgressBar

plt.ioff()   

AMPS = ["LL", "LU", "RL", "RU"]

SPECID = ["004", "008", "012", "013", "016", "017", "020", "024", "025", "027",
          "032", "037", "038", "041", "047", "051"]

CAM_IFUSLOT_DICT = {'004':'093',
                    '037':'074',
                    '027':'075',                
                    '047':'076',
                    '024':'073',
                    '013':'084',
                    '016':'085',
                    '041':'086',
                    '051':'083',
                    '008':'094',
                    '025':'095',
                    '038':'096',
                    '020':'103',
                    '032':'104',
                    '012':'105',
                    '017':'106',}

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
                        
    parser.add_argument("--zero", 
                        help='''Work on zero frames.''',
                        action="count", default=0)

    parser.add_argument("--dark", 
                        help='''Work on dark frames.''',
                        action="count", default=0)
                                                                         
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


def get_basename_depth(filename, depth=1):
    cnt=0
    while cnt < depth:
        filename = op.dirname(filename)
        cnt+=1
    return op.basename(filename)


def get_region_values(F, dx=0.05):
    a,b = F[0].data.shape
    avg = np.zeros((9,))
    cnt=0
    for j in [0.9, 0.5, 0.1]:
        for i in [0.1, 0.5, 0.9]:
            xl = int(i*b - dx*b/2.)
            xh = int(i*b + dx*b/2.)
            yl = int(j*a - dx*a/2.)
            yh = int(j*a + dx*a/2.)
            avg[cnt] = biweight_location(F[0].data[yl:yh,xl:xh]) 
            cnt+=1 
    return avg

def build_dataframe(_dataframe, date, fn, spec):
    F = fits.open(fn)
    specid = "%03d" %F[0].header['SPECID']
    if specid == spec:
        exp = get_basename_depth(fn, depth=2)
        obs = get_basename_depth(fn, depth=3)
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
        temp = F[0].header['DETTEMP']
        mjd = F[0].header['MJD']
        avg = get_region_values(F) 
        if np.all(np.isfinite(avg)):
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
                 'temp' : pd.Series(temp, index=[date+'T'+F[0].header['UT']]),
                 'mjd' : pd.Series(mjd, index=[date+'T'+F[0].header['UT']]),
                 'SPECID' : pd.Series(specid, index=[date+'T'+F[0].header['UT']]),
                 'AMP' : pd.Series(amp, index=[date+'T'+F[0].header['UT']]),
                 'OBS' : pd.Series(obs, index=[date+'T'+F[0].header['UT']]),
                 'EXP' : pd.Series(exp, index=[date+'T'+F[0].header['UT']])}
            for i,val in enumerate(avg):
                strv = 'VAL' + str(i)
                A[strv] = pd.Series(val, index=[date+'T'+F[0].header['UT']])
            data = DF(A)
            _dataframe = _dataframe.append(data)
    return(_dataframe)
    
    
def main():
    args = parse_args()
    _dataframe = DF()
    for spec in args.specid:
        if args.dark:
            if op.exists(op.join(args.outfolder, 'bias_DF_%s.csv' %spec)):
                _bias_dataframe = pd.read_csv(op.join(args.outfolder, 
                                          'bias_DF_%s.csv' %spec))
            else:
                print("Cannot find the dataframe for the biases.")
                
        ifuslot = CAM_IFUSLOT_DICT[spec]
        lower_folder_struct = op.join('virus', 'virus*', 'exp*', 'virus',
                                      '2*_%s*zro.fits' %ifuslot)
        progress = ProgressBar(len(args.cal_dirs), spec, fmt=ProgressBar.FULL)
        for date in args.cal_dirs:
            files = glob.glob(op.join(date,lower_folder_struct))
            for fn in files:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    _dataframe = build_dataframe(_dataframe, op.basename(date), fn, spec)
            progress.current+=1
            progress()
        progress.done()
    norm = plt.Normalize()
    colors = plt.cm.viridis_r(norm(np.arange(9+2)))
    if args.zero:
        _dataframe.to_csv('bias_DF_%s.csv' %spec)
    for amp in AMPS:
        fig1 = plt.figure(figsize=(8,6))
        fig2 = plt.figure(figsize=(8,6))
        fig3 = plt.figure(figsize=(8,6))
        df = _dataframe.query('AMP=="%s"'%amp)
        for i in xrange(9):
            strv = 'VAL' + str(i)
            df = df[(is_outlier(df['overscan'])<1)*(is_outlier(df[strv])<1)*
                    (is_outlier(df['temp'])<1)]
            plt.figure(fig1.number)
            plt.scatter(df['overscan'],df[strv]-df['overscan'], edgecolor='none',
                        s=25, color=colors[i,0:3], alpha=0.3)
            plt.xlabel('Overscan [ADU]')
            plt.ylabel('BIAS [ADU]')
            plt.figure(fig2.number)
            plt.scatter(df['temp'],df[strv]-df['overscan'], edgecolor='none',
                        s=25, color=colors[i,0:3], alpha=0.3)
            plt.xlabel('TEMP')
            plt.ylabel('BIAS [ADU]')
            plt.figure(fig3.number)
            plt.scatter(df['mjd'],df[strv]-df['overscan'], edgecolor='none',
                        s=25, color=colors[i,0:3], alpha=0.3)
            plt.xlabel('MJD')
            plt.ylabel('BIAS [ADU]')
        plt.figure(fig1.number)
        plt.savefig(op.join(args.output,'bias_struct_%s_%s_overscan.pdf' %(spec, amp)),dpi=150)
        plt.figure(fig2.number)
        plt.savefig(op.join(args.output,'bias_struct_%s_%s_temp.pdf' %(spec, amp)),dpi=150)
        plt.figure(fig3.number)
        plt.savefig(op.join(args.output,'bias_struct_%s_%s_mjd.pdf' %(spec, amp)),dpi=150)
        plt.close(fig1)
        plt.close(fig2)
        plt.close(fig3)

   
if __name__ == '__main__':
    main() 