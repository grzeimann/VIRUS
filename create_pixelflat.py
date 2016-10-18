# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 16:06:44 2016

@author: gregz
"""

import argparse as ap
import glob
import os
from astropy.io import fits
import os.path as op
import shutil
import re
from scipy.stats import sigmaclip
from scipy.signal import medfilt2d
from utils import imbox, biweight_location, biweight_midvariance
import numpy as np
import datetime
    
SPEC    = [ "LL", "LU", "RL", "RU" ]
SPECBIG = [ "L", "R" ]


class VirusFrame:
    def __init__(self, filename=None):
        '''
        Initializing a VirusFrame for a given filename.
        This includes reading the header and setting reduction keywords.
        From this, the file with have attributes that can be tracked.
        '''
        
        if filename is None:
            print ( "No filename was given for VirusFrame to be initialized." )
            return None
        else:
            ######## OPEN AND KEEP FILES IN SOME SMART WAY ########
            self.filename               = filename
            self.origloc                = op.dirname(self.filename)
            try:
                self.basename, temp1, temp2 = op.basename(self.filename).split('_')
            except ValueError:
                print(self.filename)
                raise ValueError
            self.type                   = temp2.split('.fits', 1)[0]
            self.ifuslot                = temp1.split('LL')[0]
            self.time                   = self.basename.split('T')[1]
            self.hr                     = self.basename.split('T')[1][0:2]
            self.min                    = self.basename.split('T')[1][2:4]
            self.sec                    = self.basename.split('T')[1][4:8]

            ###### READ IN THE HEADER AND NECESSARY KEYWORDS ######
            self.trimsec = {}
            self.biassec = {}
            self.currentbase = {}
            self.fits = {}
            self.indname = {}

            for amp in SPEC:
                self.currentbase[amp] = self.basename
                rootname = op.join(self.origloc, self.basename + '_' 
                                                + self.ifuslot + amp + '_'
                                                + self.type + '.fits')
                self.indname[amp] = (self.basename + '_' + self.ifuslot + amp 
                                     + '_' + self.type + '.fits')
                hdulist = fits.open(rootname)
                trim = re.split('[\[ \] \: \,]', hdulist[0].header['TRIMSEC'])[1:-1]
                self.trimsec[amp] = [int(t)-((i+1)%2) for i,t in enumerate(trim)]
                bias = re.split('[\[ \] \: \,]', hdulist[0].header['BIASSEC'])[1:-1]
                self.biassec[amp] = [int(b)-((i+1)%2) for i,b in enumerate(bias)]                
                self.fits[amp] = hdulist[0].data   

    def addbase(self, action, amp):
        self.currentbase[amp] = action + self.currentbase[amp] 

        
def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")    

        
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
    description = "Pixel Flat Script"
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--flt", nargs='?', type=str, help='''Location of the flt files (either
                        pixel flats or fiber flats).  For example: 
                        "/work/03946/hetdex/maverick/virus_lab/20160309/camra0000005"''')
                        
    parser.add_argument("--drk", nargs='?', type=str, help='''Location of the drk files.  
                        For example: 
                        "/work/03946/hetdex/maverick/virus_lab/20160309/camra0000006"''')
                        
    parser.add_argument("--camra", nargs='?', type=int, help='''Number of the camra you plan to examine.
                        For example:                     
                        12 or 012."''')
                        
    parser.add_argument("--binning", nargs='?', type=str, help='''Whether the pixel flats should be binned or not?"''',
                        default = False)  

    parser.add_argument("--join", nargs='?', type=str, help='''Whether the pixel flats should be binned or not?"''',
                        default = False)  
                        
    parser.add_argument("--newlab", nargs='?', type=str, help='''Was the data taken on or after March 07, 2016?  
                        Should be a True or False input.  If none is given, the default is True."''',
                        default = True)

    parser.add_argument("--runclean", nargs='?', type=str, 
                        help='''Copy files from flt and drk folders, 
                        if False, then run from existing folders''',
                        default = True) 
                        
    parser.add_argument("--smoothing", type=str, help="Smoothing size (y,x)",
                        default = '(101,101)')

    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.flt is None:
        msg = 'The flt folder location was not provided'
        parser.error(msg)
    if args.drk is None:
        msg = 'The drk folder location was not provided'
        parser.error(msg)
    if args.camra is None:
        msg = 'The camra number was not provided'
        parser.error(msg)
    if args.binning is not False:
        args.binning = str2bool(args.binning)        
    if args.newlab is not True:
        args.newlab  = str2bool(args.newlab)    
    vals = re.split('[\( \, \)]',args.smoothing)
    if len(vals) is not 4:
        msg = 'The smooth parameter must be of the form: "(y,x)"'
        parser.error(msg)
    else:
        args.smoothy=int(vals[1])
        args.smoothx=int(vals[2])
        
    return args
    

def addkeywordHeader(fits, keyword, value, comment=None):
    if comment:
        fits.header[keyword] = (value, comment)
    else:
        fits.header[keyword] = value 
        

def subtractfits(frame, amp, image=None, value=None):
    
    if value is not None:
        frame.fits[amp] -= value    
    if image is not None:
        if frame.fits[amp].shape == image.data.shape:
            frame.fits[amp] -= image.data
        else:
            print("[ERROR] Tried to subtract an image of size: %s\n"
                  "[ERROR] From an image of size: %s\n" %(image.shape, 
                                                    frame.fits[amp].shape,))
    
def overscansub(frames, amp):
    for f in frames:
        overscan_data = f.fits[amp][f.biassec[amp][2]:f.biassec[amp][3],
                                    f.biassec[amp][0]:f.biassec[amp][1]]
        overscan_value = np.mean(sigmaclip(overscan_data.flatten(), low=3.5, 
                                           high=3.5)[0])
        print("Subtracting %0.3f from %s" %(overscan_value, 
                                            op.basename(f.indname[amp])))
        subtractfits(f, amp, value=overscan_value)
        
        
def meanfits(frames, amp):
    bigarray = np.array([f.fits[amp][f.trimsec[amp][2]:f.trimsec[amp][3],
                                     f.trimsec[amp][0]:f.trimsec[amp][1]] 
                                     for f in frames])
    return biweight_location(bigarray, axis=(0,))   
    
def biweight_filter_2d(image, width, adaptive=False, tol=0.001):
    a,b = image.shape
    if len(width) > 1:
        xw = width[0]
        yw = width[1]
    else:
        xw = width[0]
        yw = width[0]
    smooth = np.zeros(image.shape)
    var = np.zeros(image.shape)

    for i in xrange(b):
        for j in xrange(a):
            xl = np.max([0,i-xw])
            yl = np.max([0,j-yw])
            xh = np.min([b-1,i+xw])
            yh = np.min([a-1,j+yw])
            loc = (j-yl) * (xh+1-xl) + (i-xl)
            smooth[j,i] = biweight_location(np.delete(image[yl:yh+1,xl:xh+1],(loc)))
            var[j,i] = biweight_midvariance(np.delete(image[yl:yh+1,xl:xh+1],(loc)))
    return smooth, var
    
    
def writefits(fitsfile, args, outfile):
    p = fits.PrimaryHDU(data=fitsfile)
    p.header['SPECID'] = ("%i" %(args.camra), 'Spectrograph ID')
    p.header['HISTORY'] = ("Created on %s" %(datetime.datetime.now().strftime(
                                                      "%y/%m/%d at %H:%M:%S")))
    p.header['HISTORY'] = ("From pixel flat data taken on: %s" 
                                          %(op.basename(op.dirname(args.flt))))                                                  
    p.header['HISTORY'] = ("Found in folder: %s" 
                                          %(op.basename(args.flt)))                
    p.header['HISTORY'] = ("Corrected for scattered light using "
                           "data taken on: %s" 
                                          %(op.basename(op.dirname(args.drk))))                                                  
    p.header['HISTORY'] = ("Found in folder: %s" 
                                          %(op.basename(args.drk)))
    p.header['HISTORY'] = "Created by pixelflat2.py."
    p.header['HISTORY'] = "Steps include:"
    p.header['HISTORY'] = "1) Subtract overscan from all frames" 
    p.header['HISTORY'] = "   (mean value with 3.5 sigma clipping)"
    p.header['HISTORY'] = "2) Create median flat frame"
    p.header['HISTORY'] = "3) Create median dark frame"
    p.header['HISTORY'] = "4) Subtract average dark from average flat"
    p.header['HISTORY'] = "5) Smooth resultant using a 2d median filter"
    p.header['HISTORY'] = "   The width of the filter is (25,3) in the (y,x) direction"
    p.header['HISTORY'] = "6) Divide the smoothed flat from the unsmoothed"
    p.header['HISTORY'] = "7) Bin 2x1 in the x-direction"
    p.header['HISTORY'] = "8) Join Amplifiers to form a single CCD"

    p.writeto(outfile,clobber=True)
    
def main():
    args = parse_args()
    # March 7th, 2016 and beyond had a change in the naming convention going from "virus" to "camra" for the lab
    if args.newlab:
        subfolder = "camra"
        print("Looking at data on or after March 7, 2016.")
    else:
        subfolder = "virus"
        print("Looking at data before March 7, 2016.")
    
    newdir = "PixelFlats_cam%03d" %(args.camra)
    # If the directory does not exist, then create it.
    if not op.exists(newdir):
        os.mkdir(newdir)
    if not op.exists(op.join(newdir, 'flats')):
        os.mkdir(op.join(newdir, 'flats'))            
    if not op.exists(op.join(newdir, 'darks')):
        os.mkdir(op.join(newdir, 'darks'))
    if args.runclean:
        for sp in SPEC:
            # Get all of the pixel flat and dark images and copy to temporary directory.
            flat_names = glob.glob(op.join(args.flt, 'exp*', subfolder, 
                                           '*' + sp + '*.fits'))
            for f in flat_names:
                shutil.copy(f, op.join(newdir, 'flats'))
            dark_names = glob.glob(op.join(args.drk, 'exp*', subfolder, 
                                           '*' + sp + '*.fits'))
            for f in dark_names:
                shutil.copy(f, op.join(newdir, 'darks'))
    sp = SPEC[0]
    fframe = [] # list of flat frames
    dframe = [] # list of dark frames
    for f in glob.glob(op.join(newdir, 'flats', '2*' + sp + '*.fits')):
        fframe.append(VirusFrame(f))
    for f in glob.glob(op.join(newdir, 'darks', '2*' + sp + '*.fits')):
        dframe.append(VirusFrame(f))
    
    pixelflat_amp = {}
    for sp in SPEC:
        overscansub(fframe, sp)
        overscansub(dframe, sp)
        meanflat = meanfits(fframe, sp)
        meandark = meanfits(dframe, sp)
        #[subtractfits(frame, sp, image=meandark) for frame in fframe]
        flat_minus_dark = meanflat - meandark
        smooth = medfilt2d(flat_minus_dark,(1,args.smoothx))
        smooth = medfilt2d(smooth,(args.smoothy,1))
        pixelflat = flat_minus_dark / smooth
        if args.binning:
            pixelflat_amp[sp] = pixelflat[:,0::2]/2. + pixelflat[:,1::2]/2.
            writefits(smooth, args, op.join(newdir,'smooth_%s.fits' %(sp)))
            writefits(flat_minus_dark, args, 
                      op.join(newdir,'amp_%s.fits' %(sp)))
        else:
            pixelflat_amp[sp] = pixelflat
    
    a,b = pixelflat_amp[SPEC[0]].shape
    pixelflat_L = np.zeros((2*a,b))
    pixelflat_R = np.zeros((2*a,b))
    for sp in SPEC:
        name = "pixelflat_cam%03d_%s.fits" %(args.camra, sp)
        writefits(np.array(pixelflat_amp[sp]), args, op.join(newdir, name))
    if args.join:
        if args.newlab:
            pixelflat_L[:a,:] = pixelflat_amp["LL"]
            pixelflat_L[a:,:] = pixelflat_amp["LU"][::-1,::-1] # Flip x and y
            pixelflat_R[:a,:] = pixelflat_amp["RU"]
            pixelflat_R[a:,:] = pixelflat_amp["RL"][::-1,::-1] # Flip x and y
        else:
            pixelflat_L[:a,:] = pixelflat_amp["LL"][:,::-1] # Flip x
            pixelflat_L[a:,:] = pixelflat_amp["LU"][:,::-1] # Flip x
            pixelflat_R[:a,:] = pixelflat_amp["RU"][::-1,::-1] # Flip x and y
            pixelflat_R[a:,:] = pixelflat_amp["RL"][::-1,::-1] # Flip x and y       
        name_L = "pixelflat_cam%03d_L.fits" %(args.camra)
        writefits(pixelflat_L, args, op.join(newdir, name_L))
        name_R = "pixelflat_cam%03d_R.fits" %(args.camra)
        writefits(pixelflat_R, args, op.join(newdir, name_R))


if __name__ == '__main__':
    main()  