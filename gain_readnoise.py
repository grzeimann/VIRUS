# -*- coding: utf-8 -*-
"""
Created on Mon April 4 13:52:02 2016

Used to measure the gain and read noise in both lab pixel flat, lab fiber
flat, and HET fiber flat data.  

@author: gregz
"""

import numpy as np
import pyfits
import argparse as ap
import glob
import re
import os.path as op
import scipy.stats as stats
import time
import copy
import scipy
import matplotlib.pyplot as plt

SPEC = ["LL","LU","RL","RU"]

IFUSLOT_DICT = {'093':['004','023'],
                '074':['037','024'],
                '075':['027','001'],
                '076':['047','016'],
                '073':['024','051'],
                '084':['013','019'],
                '085':['016','026'],
                '086':['041','015'],
                '083':['051','033'],
                '094':['008','054'],
                '095':['025','020'],
                '096':['038','014'],
                '103':['020','004'],
                '104':['032','028'],
                #'105':['012','055'],
                '106':['017','022'],
                '300016':['012','999']} 

class Fibers:
    def __init__(self, N, D, order = 3):
        self.flag = -99999
        self.x = self.flag * np.ones((D,),dtype = np.int)
        self.y = np.zeros((D,))
        self.trace = np.zeros((D,))
        self.norm = np.zeros((D,))
        self.order = order
        self.polyvals = np.zeros((order,))
        self.deadfiber = False
        
    def fitpoly(self):
        sel = self.x != self.flag
        self.polyvals = np.polyfit(self.x[sel] / 1032.,self.y[sel],self.order)
        
    def evalpoly(self):
        self.trace = np.polyval(self.polyvals, np.arange(len(self.x)) / 1032.)         

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
    description = "Gain and Readnoise Measurement"
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--flt", nargs='?', type=str, help='''Location of the flt files (either
                        pixel flats or fiber flats).  For example: 
                        "/work/03946/hetdex/maverick/virus_lab/20160321/camra0000020"''')
                        
    parser.add_argument("--zro", nargs='?', type=str, help='''Location of the zro files.  
                        For example: "/work/03946/hetdex/maverick/virus_lab/20160322/camra0000021"''')
                        
    parser.add_argument("--lab", help="Is this lab data?", action="count", 
                        default=0)
        
    parser.add_argument("--fiber_trace", help='''Use pixels near trace.''',
                        action="count", default=0)

    parser.add_argument("--debug", help='''Debug.''',
                        action="count", default=0)
                        
    parser.add_argument("--specid", type=str, help="Camera or SPECID", 
                        default=None)

    parser.add_argument("--bin_bias", action="count", default=0, 
                        help="Bin the Bias Frames?")


    parser.add_argument("--dist_from_trace", type=float, 
                        help='''Pixel distance from trace''',
                        default=2.)

    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.flt is None:
        msg = 'The flt folder location was not provided'
        parser.error(msg)
    if args.zro is None:
        msg = 'The zro folder location was not provided'
        parser.error(msg)
    if args.specid is None:
        if not args.lab:
            msg = 'The specid was not provided and it is not lab data.'
            parser.error(msg)
    else:
        try:
            args.specid = "%03d" %(int(args.specid))
        except ValueError:
            msg = 'The specid was not proper format.  Try "008" or "8".'
            parser.error(msg)            
    return args
    
class Trace():
    def __init__(self, fits, debug=False):
        self.flat = fits
        self.debug = debug
        self.a, self.b = self.flat.shape
        self.y, self.x = np.indices(self.flat.shape)
        self.y = self.y.ravel()
        self.x = self.x.ravel()
        
    def gauss(self, x, a, b, c):
        return a / np.sqrt(2*np.pi*c**2) * np.exp(-1./2.*(x-b)**2/c**2)
    
    def calculate_trace(self, window=8.7, step=8.7, order=3, interactive=False):
        N, D     = self.flat.shape
        cols     = np.arange(30 , self.b-30, self.b/50)
        init_col = self.b/2
        fdist    = 1.5
        if self.debug:
            t1=time.time()
        if interactive:
            vl, gs, pk, dd = self.find_fibers(init_col, window, step, 
                                              interactive=interactive)
            allfibers = []
            for i in xrange(len(pk)):
                f = Fibers(N, D, order)
                f.x[init_col] = init_col
                f.y[init_col] = pk[i]
                f.norm[init_col] = gs[i]
                f.deadfiber = dd[i]>0
                allfibers.append(copy.deepcopy(f))
        else:
            vl, gs, pk = self.find_fibers(init_col, window, step, 
                                              interactive=interactive)
            allfibers = []
            for i in xrange(len(pk)):
                f = Fibers(N, D, order)
                f.x[init_col] = init_col
                f.y[init_col] = pk[i]
                f.norm[init_col] = gs[i]
                f.deadfiber = False
                allfibers.append(copy.deepcopy(f))            
        brcol = np.argmin(np.abs(cols-init_col))
        cols1 = cols[brcol::-1]
        cols2 = cols[(brcol+1)::1]
        for c in cols1:
            vl, gs, pk = self.find_fibers(c, window, step, interactive=False)
            xloc  = np.argmin(np.abs(allfibers[0].x - c))
            for i in xrange(len(pk)):
                yval = np.hstack([fi.y[xloc] for fi in allfibers])
                floc = np.argmin(np.abs(pk[i] - yval))
                if np.abs(pk[i] - yval[floc]) < fdist:
                    allfibers[floc].x[c] = c
                    allfibers[floc].y[c] = pk[i]
                    allfibers[floc].norm[c] = gs[i]
        for c in cols2:
            vl, gs, pk = self.find_fibers(c, window, step, interactive=False)
            xloc  = np.argmin(np.abs(allfibers[0].x - c))
            for i in xrange(len(pk)):
                yval = np.hstack([fi.y[xloc] for fi in allfibers])
                floc = np.argmin(np.abs(pk[i] - yval))
                if np.abs(pk[i] - yval[floc]) < fdist:
                    allfibers[floc].x[c] = c
                    allfibers[floc].y[c] = pk[i]
                    allfibers[floc].norm[c] = gs[i]

        cnt = 0     
        for fiber in allfibers:
            if fiber.deadfiber:
                if cnt>0:
                    fiber.x = allfibers[cnt-1].x
                    fiber.y = allfibers[cnt-1].y + step
                    fiber.fitpoly()
                    fiber.evalpoly()
                else:
                    fiber.x = allfibers[cnt+1].x
                    fiber.y = allfibers[cnt+1].y - step
                    fiber.fitpoly()
                    fiber.evalpoly()
            else:
                fiber.fitpoly()
                fiber.evalpoly()
            cnt += 1
        if self.debug:
            t2 = time.time()
            print("[Trace] Time Taken per Fiber: %0.2f ms" 
                  %((t2-t1)/len(pk)*1e3))
            print("[Trace] Time Taken over range: %0.2f s" %((t2-t1))) 
        
        return allfibers
        
        
    def find_fibers(self, col, window, step, interactive=False):
        y = self.flat[:,col]
        x = self.y.reshape((self.a,self.b))[:,col]
        mx = np.sort(y)[int(len(y)*0.99)]
        high = np.where(y>0.1*mx)[0]
        start = high[np.searchsorted(high,8,side='right')] + 3
        flag = True
        while flag:
            try:
                guessheight = (np.max(y[int(start-window/2.):
                                       int(start+window/2.)])
                               *np.sqrt(2*np.pi*2.2**2))
                popt = scipy.optimize.curve_fit(self.gauss,
                                                x[int(start-window/2.):
                                                  int(start+window/2.)],
                                                y[int(start-window/2.):
                                                  int(start+window/2.)],
                                                p0=(guessheight,start,2))
                flag = False
            except RuntimeError:
                start += 4
        peaks = []
        vals  = []    
        guess = []
        dead  = []
        peaks.append(popt[0][1])
        guess.append(guessheight)
        vals.append(popt[0][0])
        dead.append(0)
        loc = popt[0][1] + step
        pwindow = window * 5
        while loc < (len(y)-window):
            guessheight = (np.max(y[int(loc-window/2.):int(loc+window/2.)])
                           *np.sqrt(2*np.pi*2.2**2))
            if np.sum(y[int(loc-window/4.):int(loc+window/4.)]>0.10*mx)>2:
                try:
                    popt = scipy.optimize.curve_fit(self.gauss,
                                                    x[int(loc-window/2.):
                                                      int(loc+window/2.)],
                                                    y[int(loc-window/2.):
                                                      int(loc+window/2.)],
                                                    p0=(guessheight,loc,2), 
                                                    maxfev=2000)
                    if np.abs(popt[0][1]-loc) < 2:
                        loc = popt[0][1]
                        vals.append(popt[0][0])
                        guess.append(guessheight)
                        peaks.append(popt[0][1])
                        dead.append(0)
                        
                    loc += step
                except RuntimeError:
                    print("Scipy curvefit didn't find solution for single "
                          "trace postion at {:0.0f} (y "
                          "pos) in column {:0.0f} (x pos). ".format(loc,col))
                    loc += step
            else:
                if interactive:
                    plt.plot(x,y)
                    ax = plt.gca()
                    sel = (x > (loc-pwindow)) * (x< (loc+pwindow))
                    ax.set_xlim([loc-pwindow,loc+pwindow])
                    mn = y[sel].min()
                    mx = y[sel].max()
                    rn = mx - mn
                    ax.set_ylim([-.1*rn + mn, 1.1*rn+mn])
                    plt.show()
                    answer = raw_input("Is there a fiber at {:0.0f} (y pos) in"
                                       " column {:0.0f} (x pos)? ".format(loc,
                                                                          col))
                    plt.close()

                    if answer.lower() in ("yes", "y", "true", "t", "1"):
                        vals.append(0)
                        guess.append(0)
                        peaks.append(loc)
                        dead.append(1)
                loc += step
        if interactive:
            return vals, guess, peaks, dead
        else:
            return vals, guess, peaks

def main():
    args = parse_args()
    
    # If it is lab data, the sub folder is different than if it is HET data
    if args.lab:
        subfolder = "camra"
        print("Looking at LAB data.")
    else:
 	subfolder = "virus"
        print("Looking at HET data.")
    
    if args.specid is not None:
        for key in IFUSLOT_DICT:
            if IFUSLOT_DICT[key][0] == args.specid:
                ifuslot = key 
        print(ifuslot)
        flt_names = glob.glob(op.join(args.flt, 'exp*', subfolder, '*' + ifuslot + SPEC[0] + '*.fits'))
        zro_names = glob.glob(op.join(args.zro, 'exp*', subfolder, '*' + ifuslot + SPEC[0] + '*.fits'))
    else:
        flt_names = glob.glob(op.join(args.flt, 'exp*', subfolder, '*' + SPEC[0] + '*.fits'))
        zro_names = glob.glob(op.join(args.zro, 'exp*', subfolder, '*' + SPEC[0] + '*.fits'))
    # Sort the flt names, not sure why they aren't sorted from glob.glob, but they aren't!
    flt_names = sorted(flt_names)
    
    npairs   = len(flt_names)  / 2 
    nbiases  = len(zro_names)
    if npairs < 1:
        print("Must have at least two flt images for gain/readnoise measurement.")
        print("None found in: %s" %(op.join(args.flt, 'exp*', subfolder, '*' + SPEC[0] + '*.fits')))
        return None
        
    print ( "Examining %d pairs of flt files for gain and readnoise." %( npairs ) )
    print ( "Looking at flt files in %s" %( args.flt ) )
    print ( "Looking at zro files in %s" %( args.zro ) )

        
    flow    = 500
    fhigh   = 45000
    lthresh = 1500
    hthresh = 32000
                   
    gain = np.zeros((4,npairs))
    read = np.zeros((4,npairs))
    
    gainunit = {}
    readunit = {}
    readnoiseavg = {}
    gainhead = {}
    rdnoisehead = {}
    spcount = 0

    for sp in SPEC:
        beginning, ending    = zro_names[0].split(SPEC[0]) # using first zero frame
        filenameb1           = beginning + sp + ending
        p = pyfits.open(filenameb1)
        blank, xlow, xhigh, ylow, yhigh, blank = re.split('[: \[ \] ,]',p[0].header['BIASSEC'])
        xlow                 = int(xlow)
        xhigh                = int(xhigh)
        ylow                 = int(ylow)
        yhigh                = int(yhigh)
        blank, txlow, txhigh, tylow, tyhigh, blank = re.split('[: \[ \] ,]',p[0].header['TRIMSEC'])
        txlow                = int(txlow)-1
        txhigh               = int(txhigh)
        tylow                = int(tylow)-1
        tyhigh               = int(tyhigh)
        
        bias1                = np.array(p[0].data.copy(),dtype=np.float)
        overscan             = np.mean(bias1[ylow:yhigh,xlow:xhigh])
        bias1                -= overscan
        bias1                = bias1[tylow:tyhigh,txlow:txhigh] 
        if args.bin_bias:
            bias1 = (bias1[0::2,0::2] + bias1[0::2,1::2] + bias1[1::2,0::2] 
                        + bias1[1::2,1::2])

        beginning, ending    = flt_names[0].split(SPEC[0]) # using first zero frame
        filenameb1           = beginning + sp + ending
        p = pyfits.open(filenameb1)
        blank, fxlow, fxhigh, fylow, fyhigh, blank = re.split('[: \[ \] ,]',p[0].header['BIASSEC'])
        fxlow                 = int(fxlow)
        fxhigh                = int(fxhigh)
        fylow                 = int(fylow)
        fyhigh                = int(fyhigh)
        blank, ftxlow, ftxhigh, ftylow, ftyhigh, blank = re.split('[: \[ \] ,]',p[0].header['TRIMSEC'])
        ftxlow                = int(ftxlow)-1
        ftxhigh               = int(ftxhigh)
        ftylow                = int(ftylow)-1
        ftyhigh               = int(ftyhigh)
        
        mf1 = np.zeros((npairs,))
        mf2 = np.zeros((npairs,))
        if args.fiber_trace:
            beginning, ending    = flt_names[0].split(SPEC[0]) # looping through the first of the frames
            filenamef1           = beginning + sp + ending
            print("Calculating trace from: %s" %(filenamef1))
            p                    = pyfits.open(filenamef1)
            gainhead[sp]         = p[0].header['GAIN']
            rdnoisehead[sp]      = p[0].header['RDNOISE']
            flat1                = np.array(p[0].data.copy(),dtype=np.float)
            overscan             = np.mean(flat1[fylow:fyhigh,fxlow:fxhigh])
            flat1                -= overscan
            flat1                = flat1[ftylow:ftyhigh,ftxlow:ftxhigh]
            T = Trace(flat1, debug=args.debug)
            allfibers = T.calculate_trace(interactive=False)            
            y = np.array([F.trace for F in allfibers])
            a,b = flat1.shape
            yarr,xarr = np.indices(flat1.shape)
            mask = np.zeros((a,b),dtype=bool)
            if args.debug:
                t1 = time.time()
            for i in xrange(a):
                for j in xrange(b):
                    mask[i,j] = np.any(np.abs(y[:,j]-yarr[i,j])<args.dist_from_trace)
            if args.debug:
                t2 = time.time()
                print("Time taken creating mask: %0.2f s" %(t2-t1))
        else:
            mask = np.ones(bias1.shape,dtype=bool)

        bigbias = np.zeros(bias1.shape + (nbiases,))
        for i in xrange(nbiases):       
            beginning, ending   = zro_names[i].split(SPEC[0]) # using second zero frame
            filenameb           = beginning + sp + ending
            p = pyfits.open(filenameb)
            bias                = np.array(p[0].data.copy(), dtype=np.float)
            overscan            = np.mean (bias[ylow:yhigh,xlow:xhigh])
            bias               -= overscan
            bias                = bias[tylow:tyhigh,txlow:txhigh]
            if args.bin_bias:
                bias = (bias[0::2,0::2] + bias[0::2,1::2] + bias[1::2,0::2] 
                        + bias[1::2,1::2])
            bigbias[:,:,i] = bias
        rdnoiseimage     = np.std(bigbias, axis=2 )
        avgbiasimage     = np.mean(bigbias, axis=2 )
        rdnoiseclipped   = stats.sigmaclip(rdnoiseimage)[0]
        readnoiseavg[sp] = np.mean(rdnoiseclipped)
        
        for i in xrange(npairs):    
    
           beginning, ending    = flt_names[2*i].split(SPEC[0]) # looping through the first of the frames
           filenamef1           = beginning + sp + ending
           p                    = pyfits.open(filenamef1)
           exptime1 = p[0].header['EXPTIME']
           readtime1 = p[0].header['READTIME']
           flat1                = np.array(p[0].data.copy(),dtype=np.float)
           overscan             = np.mean(flat1[fylow:fyhigh,fxlow:fxhigh])
           flat1                -= overscan
           flat1                = flat1[ftylow:ftyhigh,ftxlow:ftxhigh]
           
           beginning, ending    = flt_names[2*i+1].split(SPEC[0]) # looking at consecutive pairs
           filenamef2           = beginning + sp + ending
           p = pyfits.open(filenamef2)
           exptime2 = p[0].header['EXPTIME']
           readtime2 = p[0].header['READTIME']
           flat2                = np.array(p[0].data.copy(),dtype=np.float)
           overscan             = np.mean(flat2[fylow:fyhigh,fxlow:fxhigh])
           flat2                -= overscan
           flat2                = flat2[ftylow:ftyhigh,ftxlow:ftxhigh]

           x, y = np.where((flat1 > flow) * (flat1 < fhigh) * (mask))
           if len(x)>10:
               mf1[i]  = np.mean(flat1[x,y])
               mf2[i]  = np.mean(flat2[x,y])
               mb  = np.mean(avgbiasimage[x,y])
               df   = flat1[x,y] - flat2[x,y]*mf1[i]/mf2[i]
               cdf  = stats.sigmaclip(df)[0]
               sdf  = np.std(cdf)
               mn   = (mf1[i] + mf2[i] - 2*mb) / 2.
               vr   = (sdf**2 - 2.*readnoiseavg[sp]**2) / 2.
               gain[spcount,i] = mn / vr
               read[spcount,i] = gain[spcount,i] * readnoiseavg[sp]
               print("%s | Gain: %01.3f | RDNOISE: %01.3f | F1: %5d | F2: %5d | Var: %05.1f | E1: %3.2f | E2: %3.2f | R1: %3.2f | R2: %3.2f " %(sp, gain[spcount,i],read[spcount,i], mf1[i], mf2[i], vr, exptime1, exptime2, readtime1, readtime2))
        gainunit[sp] = np.median(gain[spcount,( mf1 > lthresh ) * (mf2 > lthresh ) * ( mf1 < hthresh ) * (mf2 < hthresh ) * (mf1/vr<1.0)]) # Only include pixel flats above lthresh
        readunit[sp] = np.median(read[spcount,( mf1 > lthresh ) * (mf2 > lthresh ) * ( mf1 < hthresh ) * (mf2 < hthresh ) * (mf1/vr<1.0)])
        readnoiseavg[sp] *= gainunit[sp]
        spcount = spcount+1

    if args.fiber_trace:
	print('SPECID_AMP, GAIN, RDNOISE, GAIN_HEADER, RDNOISE_HEADER:')
	for sp in SPEC:
	    	print("%s_%s: %0.3f, %0.3f, %0.3f, %0.3f" %(args.specid, sp, gainunit[sp], readnoiseavg[sp], gainhead[sp], rdnoisehead[sp]))
    else:         
    	print('GAIN, RDNOISE:')
    	for sp in SPEC: 
        	print ("%s : %0.3f, %0.3f" %(sp, gainunit[sp], readnoiseavg[sp]))



        
if __name__ == '__main__':
    main()          

