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
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import os.path as op
import textwrap
from astropy.stats import biweight_location, biweight_midvariance

data_dir = '/work/03946/hetdex/maverick'


SPEC = ["LL","LU","RL","RU"]
SPECBIG = ["L","R"]  

# Default set of spectrographs for reduction
SPECID = ["004","008","012","013","016","017","020","024","025","027",
          "032","037","038","041","047","051"]

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
    description = textwrap.dedent('''Calibration Inspection Script - 
    
                     This scripts examines and compares the following:
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

    parser.add_argument("--side", nargs='?', type=str, 
                        help='''List of SPECID's for processing. 
                        Ex: "L"
                        Default: "L,R"''', default = None)

    parser.add_argument("--output", nargs='?', type=str, 
                        help='''Output Directory [REQUIRED]''', 
                        default=None)
                        
    parser.add_argument("-d","--inspect_darks", 
                        help='''Inspect darks.''',
                        action="count", default=0) 

    parser.add_argument("--dates", nargs='?', type=str, 
                        help='''Dates to examine.''', 
                        default=None)
                                                                         
    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.specid:
        args.specid = args.specid.replace(" ", "").split(',')
    else:
        args.specid = SPECID

    # Check that the arguments are filled
    if args.dates:
        nargs = args.dates.replace(" ", "").split('*')
        if len(nargs)>1:
            args.dates = args.dates.replace(" ", "")
            args.files = glob.glob(op.join(data_dir, args.dates, 'virus',
                                           'exp*','virus','2*.fits'))
        else:
            args.files = []
            args.dates = args.dates.replace(" ", "").split(',')  
            for date in args.date:
                for f in glob.glob(op.join(data_dir, date, 'virus',
                                           'exp*','virus','2*.fits')):
                   args.files.append(f)
                                               
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
        args.amps = SPEC

    if args.side:
        args.side = args.side.replace(" ", "").split(',')
    else:
        args.side = SPECBIG      
            
    return args
    
def make_file_dict(args, amp_exam=False, side_exam=False):
    file_dict = {}
    file_dict.setdefault("DATA", [])
    file_dict.setdefault("DATE-OBS", [])
    file_dict.setdefault("UT", [])
    file_dict.setdefault("AMBTEMP", [])
    file_dict.setdefault("DETTEMP", [])
    file_dict.setdefault("SPECID", [])
    if amp_exam:
        file_dict.setdefault("AMP", [])
    if side_exam:
        file_dict.setdefault("SIDE", [])

    for fn in args.files:
        p = pyfits.open(args.filename)
        file_dict['UT'].append(p[0].header['UT'])
        file_dict['DATE-OBS'].append(p[0].header['DATE-OBS'])
        file_dict['AMBTEMP'].append(p[0].header['AMBTEMP'])
        file_dict['DETTEMP'].append(p[0].header['DETTEMP'])
        file_dict['SPECID'].append(p[0].header['SPECID'])
        file_dict['AMP'].append(p[0].header['CCDPOS']+p[0].header['CCDHALF'])
        file_dict['DATA'].append(p[0].data)         
    return file_dict
    
def main():
    args = parse_args()
    if args.inspect_darks:
        file_dict = make_file_dict(args, amp_exam=True)
        name = [file_dict['DATE-OBS'][i]+file_dict['UT'][i] 
                for i in xrange(len(file_dict['UT']))]
        udate = list(np.unique(np.array(name)))
        cmap = plt.get_cmap('plasma_r')
        colors = cmap(np.arange(256))
        ncols = len(udate)
        nint  = 256 / ncols
        colors = colors[nint/2:(256-nint/2+1):nint]
        for cam in args.specid:
            plt.figure(figsize=(16,12))
            for ampind, amp in enumerate(args.amps):
                matches = [i for i in xrange(len(file_dict['SPECID']))
                          if file_dict['SPECID'][i]==cam 
                          and file_dict['AMP'][i]==amp]
                mny = []
                mxy = []
                for ind in matches:
                    data = file_dict["DATA"][ind]
                    N,D = data.shape
                    dettemp = file_dict["DETTEMP"][ind]
                    DATE_OBS = file_dict["DATE-OBS"][ind]
                    UT = file_dict["UT"][ind]
                    dind = udate.index(DATE_OBS+UT)
                    x = np.arange(300,D-300)
                    y = np.median(data[(N/2-250):(N/2+250),:],axis=0)
                    plt.subplot(2,2,ampind+1)
                    #plt.scatter(x, y[x], 
                    #            c = 'none', edgecolor='k', 
                    #            alpha=1, s=25)
                    plt.scatter(x, y[x], 
                                c = colors[dind], edgecolor='none', 
                                alpha=0.2, s=15)
                    plt.xlim([300,D-300])
                    plt.annotate("Date Obs:{:s}, UT:{:s}, Det Temp:{:02.1f}".format(
                                 DATE_OBS, UT, float(dettemp)), xy=(0.05, 
                                 (0.95-dind*0.02)), xycoords='axes fraction',
                                 color=colors[dind],fontsize=8)
                    plt.xlabel('X',fontsize=14)
                    plt.ylabel('ADU',fontsize=14)
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    rx = np.arange(200,D-200)
                    mny.append(np.percentile(y[rx],2))
                    mxy.append(np.percentile(y[rx],98))
                plt.annotate("Amp: {:s}".format(amp), xy=(0.45, 0.05), 
                             xycoords='axes fraction', color='k', fontsize=14)   
                ymin = np.min(np.array(mny))
                ymax = np.max(np.array(mxy))
                yran = ymax - ymin
                ylow = ymin - 0.1*yran
                yhigh = ymax + 0.7*yran
                #plt.ylim([-2, 50])
            plt.suptitle("Camera: {:s}".format(cam),fontsize=20)
            plt.savefig(op.join(args.output, 
                                "bias_inspection_cam{:s}.png".format(cam)),
                                dpi=150)     
                
    
if __name__ == '__main__':
    main() 