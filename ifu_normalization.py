"""

Script for quick look at wavelenght solution

@author: gregz

"""
from __future__ import print_function

import sys
import os
import time

import matplotlib.pyplot as plt
import numpy as np
import argparse as ap
import os.path as op

from astropy.io import fits
from os import environ
from utils import biweight_location
from bspline import Bspline
from pyhetdex.cure.fibermodel import FiberModel
from scipy.optimize import lsq_linear

virus_config = "/work/03946/hetdex/maverick/virus_config"

SIDE = ["L","R"]
ucam = ["004", "008", "012", "013",  "016", "017", "020", "024", "025", "027",
        "032", "037", "038", "041", "047", "051"]

CUREBIN = None
if not CUREBIN:
    CUREBIN = environ.get('CUREBIN')
if not CUREBIN:
    print("Please set CUREBIN as  environment variable or in the script")
    sys.exit(1)
        
def run_cure_command(command, suppress_output=0, debug=1):
    '''
       Run a cure command
       
       Parameter suppress_output is used to quiet down. (default)
       Return the value of os.system(command)
    '''
    # store info in log?

    if debug:
        print('Running: \'%s\'' % command)
    if not suppress_output:
        return os.system(op.join(CUREBIN, command) +' 1>>output.log  2>> error.log')
    else:
        return os.system(op.join(CUREBIN, command))
        
def fiberextract(filename, distmodel, fibermodel, outname, opts):

    command = ('fiberextract %s -d %s -f %s -o "%s" %s' 
                %(opts, distmodel, fibermodel, outname, filename))
        
    run_cure_command(command, 0)
    
    return command
    
def throughput_fiberextract(Felist, args):
    nifu = len(Felist)
    nw = len(Felist[0][0].data[0,:])
    xp = np.linspace(0, 1, num=nw)
    nbspline = 12
    a = np.linspace(0, 1, nbspline)
    knots = np.hstack([0,0,np.vstack([a,a]).T.ravel(),1,1])
    b = Bspline(knots, 3)
    basis = np.array([b(xi) for xi in xp])
    B = np.zeros((nifu,nw))
    for i in xrange(nifu):
        spec = biweight_location(Felist[i][0].data,axis=(0,))
        mask = np.where((~np.isnan(spec))*(~np.isinf(spec))*(spec!=0))[0]
        sol = np.linalg.lstsq(basis[mask,:], spec[mask])[0]
        B[i,:] = np.dot(basis,sol)
#        if args.plot:  
#            pltfile = op.join(args.outfolder, 'spectrum_%i.pdf' %i)
#            fig = plt.figure(figsize=(8, 6))
#            plt.plot(xp, spec)
#            plt.plot(xp, B[i,:],'r--')
#            plt.xticks([])
#            plt.xlabel('Wavelength')
#            plt.ylabel('Arbitrary Units')
#            plt.xlim([0, 1])
#            fig.savefig(pltfile, dpi=150)
#            plt.close()
    if args.plot:
        norm = plt.Normalize()
        colors = plt.cm.viridis(norm(np.arange(nifu)+1))
        pltfile = op.join(args.outfolder, 'IFU_average_spectra.pdf')
        fig = plt.figure(figsize=(8, 6))
        avgB = biweight_location(B,axis=(0,))
        for i in xrange(nifu):
            with np.errstate(divide='ignore'):
                plt.plot(xp, B[i,:]/avgB, color=colors[i,0:3], alpha=0.9)
        plt.xticks([])
        plt.xlabel('Wavelength')
        plt.ylabel('Normalized Units')
        plt.xlim([0, 1])
        plt.ylim([0.5,1.5])
        fig.savefig(pltfile, dpi=150)
        plt.close()
    return B, avgB


def edit_fibermodel(Felist, Fiblist, header_tup, B, avgB, A, args):
    nifu = len(Fiblist)
    wave = np.arange(header_tup[0]) * header_tup[2] + header_tup[1]
    for i in xrange(nifu):
        F = FiberModel(Fiblist[i])
        nfib, nw = Felist[i][0].data.shape
        through = A[i] * B[i,:] / avgB
        for j in xrange(nfib):
            if args.debug:
                t1 = time.time()
                print("Working on Fibermodel, fiber: %i, %i" %(i,j))            
            mask = np.where(np.isfinite(through[j,:]))[0]
            basis = np.vstack([F.amplitudes[j].get_basis(F._scal_w(w)) 
                              for w in wave])
            x = through[j,mask]
            x[-1] = x[-2]                
            y = basis[mask,:]
            ax,ay = y.shape
            lb = -1.*np.ones((ay,))
            lb[-1]=-10.
            ub = 1.*np.ones((ay,))
            ub[-1]=10.
            sol = lsq_linear(np.array(basis[mask,:]), np.array(x), bounds=(lb,ub))
            F.amplitudes[j].A = sol['x'][:-1]
            F.amplitudes[j].mean = sol['x'][-1]
        filename = Fiblist[i][:-6]+'adjpy_'+Fiblist[i][-6]+'.fmod'
        if args.debug:
            print("Writing out %s" % filename)
        F.writeto(op.join(args.outfolder, filename))
        if args.debug:
            t2=time.time()
            print("Time Taken to reset amplitudes for IFU %i: %0.2f" %(i,t2-t1))

def normalize_fiberextract(Felist, Fe_e_list, Fenames, B, avgB, A, args):
    nifu = len(Felist)
    for i in xrange(nifu):
        outfile = op.join(op.dirname(Fenames[i]),'n'+op.basename(Fenames[i]))
        outfile_e = op.join(op.dirname(Fenames[i]),'e.n'+op.basename(Fenames[i]))
        outfile_t = op.join(op.dirname(Fenames[i]),'t'+op.basename(Fenames[i]))
        with np.errstate(divide='ignore'):
            norm = np.where(B[i,:]!=0., Felist[i][0].data / B[i,:] * avgB, 0.)
            norm_e = np.where(B[i,:]!=0., Fe_e_list[i][0].data / B[i,:] * avgB, 0.)
            through = np.where(B[i,:]!=0., A[i] * B[i,:] / avgB, 0.)
        Felist[i][0].data = norm
        Fe_e_list[i][0].data = norm_e
        Felist[i][0].header['HISTORY'] = 'Divided by Smoothed Average IFU Spectrum'
        Fe_e_list[i][0].header['HISTORY'] = 'Divided by Smoothed Average IFU Spectrum'
        Felist[i].writeto(outfile, clobber=True)
        Fe_e_list[i].writeto(outfile_e, clobber=True)
        Felist[i][0].data = through
        Felist[i].writeto(outfile_t, clobber=True)

def get_fiber_amps(Fiblist, Felist, header_tup, args):
    nifu = len(Fiblist)
    wave = np.arange(header_tup[0]) * header_tup[2] + header_tup[1]
    A = []
    for i in xrange(nifu):
        if args.debug:
            t1 = time.time()
            print("Working on IFU %i" %i)
        F = FiberModel(Fiblist[i])
        nfib, nw = Felist[i][0].data.shape
        a = np.zeros((nfib, nw))
        for j, w in enumerate(wave):
            for k in xrange(1, nfib+1):
                a[k-1, j] = F.get_wf_amplitude(w, k)
        A.append(a)
        if args.debug:
            t2=time.time()
            print("Time Taken to get amplitudes for IFU %i: %0.2f" %(i,t2-t1))
    return A

def plot_fiberextract(fibextract, psize, fsize, outfile):

    fig = plt.figure(figsize=(12, 9))
    p = fits.open(fibextract)[0].data
    crval1 = fits.open(fibextract)[0].header['CRVAL1']
    cdelt1 = fits.open(fibextract)[0].header['CDELT1']
    plots = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    
    solar = np.loadtxt(op.join(virus_config, 'solar_spec', 'sun.spec'))

    pos = 0
    xs = p.shape[1]
    nfib = p.shape[0]
    wave = cdelt1 * np.arange(xs) + crval1

    for i in [1, 0.5, 0]:
        for j in [0, 0.5, 1]:

            minx = int(j*xs - j*psize)
            maxx = minx + psize
            minf = int(i*nfib - i*fsize)
            maxf = minf + fsize
            minw = wave[minx]
            maxw = wave[maxx-1]
            miny = np.max([0,np.min(p[minf:maxf, minx:maxx])])
            maxy = np.max(p[minf:maxf, minx:maxx])
            solminx = np.searchsorted(solar[:,0],minw)
            solmaxx = np.searchsorted(solar[:,0],maxw,side='right')
            po = np.polyfit(wave[minx:maxx], 
                       np.median(p[minf:maxf,minx:maxx],axis=0),3)
            ps = np.polyfit(solar[solminx:solmaxx,0], 
                       solar[solminx:solmaxx,1],3)
            sub = fig.add_subplot(3, 3, plots[pos])
            sub.set_xticks(np.arange(int(minw/50)*50, int(maxw/50 + 1)*50, 50))
            sub.set_yticks(np.linspace(miny,maxy,5))
            mult = np.polyval(po, solar[solminx:solmaxx,0])
            div = np.polyval(ps, solar[solminx:solmaxx,0])
            sub.plot(solar[solminx:solmaxx,0], 
                     solar[solminx:solmaxx,1] * mult / div,
                     color=[0.35,0.23,0.35],linewidth=2)
            for k in xrange(minf,maxf):
                sub.plot(wave[minx:maxx], p[k,minx:maxx], 
                         color=[1.0,0.23,0.35], alpha=0.5)
                xran = maxw-minw
                yran = maxy-miny
                sub.text(0.1*xran+minw,0.85*yran+miny,
                         'Fibers: %03d - %03d' %(minf+1,maxf), color='red')
            sub.set_xlim([minw, maxw])
            sub.set_ylim([miny, maxy])
            pos += 1

    plt.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)

    
def parse_arg(args):
    """
    Command line parser
    Parameters
    ----------
    argv: list of strings
        list to be parsed
    output
    ------
    namespace:
        Parsed arguments
    """

    p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)

    p.add_argument('-t','--tracefile', dest='tracefile', 
                   help="""Traceframe for plotting""", default=None)
    p.add_argument('-d', '--dist', dest='dist',  default=None,
                   help="""Distortion model file""")
    p.add_argument('-f', '--fib', dest='fib', default=None,
                   help="""Fiber model file""")
    p.add_argument('-F', '--folder', dest='folder', default=None,
                   help="""Folder for twighlights and cals""")
    p.add_argument('-O', '--outfolder', dest='outfolder', default=None,
                   help="""Folder for twighlights and cals""")
    p.add_argument('-s', '--size', dest='psize', default=150,
                   help="""Size of the postage stamps to plot""")
    p.add_argument('-fs', '--nfibers', dest='fsize', default=20,
                   help="""Number of fibers used in each postage plot""")
    p.add_argument('-S', '--step', dest='step', default=10,
                   help="""Step size for the distortion overplot""")
    p.add_argument('-o', '--options', dest='opts', 
                   default="-n 1032 -W 3500,5500 -P",
                   help="""Fiberextract options.""")
    p.add_argument('-of', '--outfile', dest='outfile', default='FePlot.pdf',
                   help="""Size of the postage stamps to plot""")
    p.add_argument('-w', '--overwrite', dest='overwrite', default=0,
                   action="count", help="""Overwrite the existing Fe* file?""")
    p.add_argument('-l', '--loop', dest='loop', 
                   action="count", default=0,
                   help="""Loop through all cameras with default locations.""") 
    p.add_argument('-p', '--plot', dest='plot', 
                   action="count", default=0,
                   help="""Plot Fiber extracted twighlights vs. solar spectrum.""") 
    p.add_argument('-D', '--debug', dest='debug', 
                   action="count", default=0,
                   help="""Debug.""") 

                   
    args = p.parse_args(args=args)
               
    if not args.loop:
        if not args.tracefile:
            msg = ('If you are not running through a loop, '
                  'you need to provide a tracefile')
            p.error(msg) 
        if not args.dist:
            msg = ('If you are not running through a loop, '
                  'you need to provide a distortion model file')
            p.error(msg) 
        if not args.fib:
            msg = ('If you are not running through a loop, '
                  'you need to provide a fiber model file')
            p.error(msg)
    else:
        if not args.folder:
            msg = ('If you are running through a loop, '
                  'you need to provide a folder for the data/cals')
            p.error(msg)             
        if not args.outfolder:
            msg = ('If you are running through a loop, '
                  'you need to provide an out folder for the Fe files.')
            p.error(msg) 

            
    return args

def make_cal_filenames(folder, specid, side):
    mastertrace = op.join(folder, 'mastertrace_twi_%s_%s.fits' %(specid, side))    
    distfile = op.join(folder, 'mastertrace_twi_%s_%s.dist' %(specid, side))    
    fibfile = op.join(folder, 'mastertrace_twi_%s_%s.fmod' %(specid, side)) 
    return mastertrace, distfile, fibfile

def main():
    # parse the command line
    args = parse_arg(sys.argv[1:])
    if args.loop:
        Felist = []
        Fe_e_list = []
        Fiblist = []
        Fenames = []
        if args.debug:
            t1 = time.time()
        for uca in ucam:
            for sp in SIDE:
                mastertrace, distfile, fibfile = make_cal_filenames(
                                                          args.folder, uca, sp)
                outfile = op.join(args.outfolder, 'FePlot_cam%s_%s.pdf' 
                                                                    %(uca, sp))
                if (op.exists(mastertrace) and op.exists(distfile) 
                        and op.exists(fibfile)):
                    FeFile = op.join(args.outfolder,'Fe%s' %(op.basename(mastertrace)))
                    FeFile_e = op.join(args.outfolder,'e.Fe%s' %(op.basename(mastertrace)))
                    Fiblist.append(fibfile)
                    Fenames.append(FeFile)
                    if not op.exists(FeFile) or args.overwrite:
                        fiberextract(mastertrace, distfile, fibfile, 
                                  op.join(args.outfolder, "%s" %(op.basename(mastertrace))), 
                                  args.opts)
                    Felist.append(fits.open(FeFile))
                    Fe_e_list.append(fits.open(FeFile_e))
                    #if args.plot:
                    #    plot_fiberextract(FeFile, args.psize, args.fsize, 
                    #                      outfile)
        if args.debug:
            t2 = time.time()
            print("Time Taken Extracting Fibers or Gathering Data: %0.2f" %(t2-t1))
        outfile = op.join(args.outfolder, 'Avg_IFU_spectrum.fits')
        if args.debug:
            t1 = time.time()
        B, avgB = throughput_fiberextract(Felist, args)
        if args.debug:
            t2 = time.time()
            print("Time Taken fitting splines to Fe files: %0.2f" %(t2-t1))
        header_tup = (Felist[0][0].header['NAXIS1'], 
                      Felist[0][0].header['CRVAL1'], 
                      Felist[0][0].header['CDELT1'])
        if args.debug:
            t1 = time.time()      
        A = get_fiber_amps(Fiblist, Felist, header_tup, args)
        if args.debug:
            t2 = time.time()
            print("Time Taken getting fiber amps: %0.2f" %(t2-t1)) 
        hdu = fits.PrimaryHDU(np.array(avgB))
        hdu.header['CRVAL1'] = Felist[0][0].header['CRVAL1']
        hdu.header['CDELT1'] = Felist[0][0].header['CDELT1']
        hdu.writeto(outfile, clobber=True)
        if args.debug:
            t1 = time.time()  
        normalize_fiberextract(Felist, Fe_e_list, Fenames, B, avgB, A, args)
        if args.debug:
            t2 = time.time()
            print("Time Taken normalizing Fe files: %0.2f" %(t2-t1)) 
        if args.debug:
            t1 = time.time()             
        edit_fibermodel(Felist, Fiblist, header_tup, B, avgB, A, args)
        if args.debug:
            t2 = time.time()
            print("Time Taken rewriting Fib models: %0.2f" %(t2-t1))             
    else:
        FeFile = op.join(op.dirname(args.tracefile),
                         'Fe' + op.basename(args.tracefile))
        if not op.exists(FeFile) or args.overwrite:
            fiberextract(args.tracefile, args.dist, args.fib, 
                         args.tracefile, args.opts)
        if args.plot:
            plot_fiberextract(FeFile, args.psize, args.fsize, args.outfile)
        
if __name__ == "__main__":
    main()