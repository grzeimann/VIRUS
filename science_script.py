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
import shutil
import glob
import os
import os.path as op
import re
import call_cure as CC
import textwrap
import numpy as np
import warnings
from virus_cosmics import remove_cosmics


'''Global variables to be set for a given reduction'''
usemapping           = True # manual map IFUSLOT to SPECID 
CLEAN_AFTER_DONE     = True # clean the previous file after new one is created
configdir    = "/work/03946/hetdex/maverick/virus_config"
darkcurrentdir = "/work/00115/gebhardt/maverick/cal/lib_dark"
biasdir = "/work/00115/gebhardt/maverick/cal/lib_bias"


# Critical naming scheme for folders and object types
# Note the type doesn't have to initial match the folder because
# the code will change the type based on the folder it is sent.
sci_dir  = "sci"

# Dictionary of the mapping between IFUSLOT and SPECID, IFUID
IFUSLOT_DICT = {'073':['024','033'],
                '074':['037','024'],
                '075':['027','001'],
                '076':['047','016'],
                '083':['051','023'],
                '084':['013','019'],
                '085':['016','026'],
                '086':['041','015'],
                '093':['004','051'],
                '094':['008','054'],
                '095':['025','020'],
                '096':['038','014'],
                '103':['020','004'],
                '104':['032','028'],
                '105':['012','055'],
                '106':['017','022'],}

# Dictionary of the mapping between SPECID and IFUID

CAM_IFU_DICT = {'004':'051',
                '037':'024',
                '027':'001',                
                '047':'016',
                '024':'033',
                '013':'019',
                '016':'026',
                '041':'015',
                '051':'023',
                '008':'054',
                '025':'020',
                '038':'014',
                '020':'004',
                '032':'028',
                '012':'055',
                '017':'022',}

# Default set of spectrographs for reduction
SPECID = ["004","008","012","013","016","017","020","024","025","027","032",
          "037","038","041","047","051"]

SPEC = ["LL","LU","RL","RU"]
SPECBIG = ["L","R"]  

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
    description = textwrap.dedent('''Science Reduction Script - 
    
                     This script produces the following output:
                         S*.fits for sky subtracted frames
                         Fe*.fits for fiber extracted spectra
                         Cu*.fits for data cubes
                     
                     The script places the files in the given output directory.
                     
                     ''')
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.RawTextHelpFormatter)

    parser.add_argument("-r","--run_insitu", help='''Run cure commands now.''',
                        action="count", default=0)
                        
    parser.add_argument("--specid", nargs='?', type=str, 
                        help='''List of SPECID's for processing. 
                        Ex: "020,008".''', default = None)

    parser.add_argument("--output", nargs='?', type=str, 
                        help='''Output Directory [REQUIRED]''', 
                        default=None)
                        
    parser.add_argument("--skipbasic", help='''Skip basic science reduction.
                        Enter the current letter status of the reduction.
                        Ex: \"pses\"''',
                        type=str, default=None) 
                        
    parser.add_argument("-s","--subsky", help='''Subtract sky (no arg nec.)''',
                        action="count", default=0)                        

    parser.add_argument("-rc","--remove_cosmics", help='''Make new error frame with cosmics set to "-1".''',
                        action="count", default=0)  

    parser.add_argument("-f","--fiberextract", 
                        help='''Fiberextract (no arg nec.)''',
                        action="count", default=0) 

    parser.add_argument("-m","--makecube", 
                        help='''Make 3d Cube (no arg nec.)''',
                        action="count", default=0) 
                        
    parser.add_argument("-d","--detect", 
                        help='''Run detect''',
                        action="count", default=0) 

    parser.add_argument("-t","--use_twi", help='''Use twi distortion soln.''',
                        action="count", default=0) 

    parser.add_argument("-ud","--use_darks", help='''Subtract dark masters.''',
                        action="count", default=0) 

    parser.add_argument("--dark1_mult_val", type=float, 
                        help='''Multiplication factor applied to dark masters.
                        If this is not set, then it is calculated from the
                        science frame.''',
                        default=1.) 

    parser.add_argument("--dark2_mult_val", type=float, 
                        help='''Multiplication factor applied to dark masters.
                        If this is not set, then it is calculated from the
                        science frame.''',
                        default=1.) 

    parser.add_argument("--dark3_mult_val", type=float, 
                        help='''Multiplication factor applied to dark masters.
                        If this is not set, then it is calculated from the
                        science frame.''',
                        default=1.) 
                        
    parser.add_argument("--dark4_mult_val", type=float, 
                        help='''Multiplication factor applied to dark masters.
                        If this is not set, then it is calculated from the
                        science frame.''',
                        default=1.)                         

    parser.add_argument("--bias1_mult_val", type=float, 
                        help='''Multiplication factor applied to dark masters.
                        If this is not set, then it is calculated from the
                        science frame.''',
                        default=1.) 

    parser.add_argument("--cal_dir", nargs='?', type=str, 
                        help='''Calibration Directory [REQUIRED]''', 
                        default=None)

    parser.add_argument("-dd","--use_deformer_default", 
                        help='''Use default deformer solutions in:
                        "/home/03946/hetdex/virus_config/DeformerDefaults"''',
                        action="count", default=0) 

    parser.add_argument("-rd","--rerun_deformer", 
                        help="Rerun Deformer with science frame",
                        action="count", default=0) 

    parser.add_argument("--instr", nargs='?', type=str, 
                        help='''Instrument to process. 
                        Default: "virus"
                        Ex: "camra" for lab data.''', default = "virus")
                        
    parser.add_argument("--rootdir", nargs='?', type=str, 
                        help='''Root Directory
                        Default: \"/work/03946/hetdex/maverick\"''', 
                        default="/work/03946/hetdex/maverick")
                        
    parser.add_argument("-sd","--scidir_date", nargs='?', type=str,
                        help='''Science Directory Date.     [REQUIRED]
                        Ex: \"20160412\"''', default=None)

    parser.add_argument("-so","--scidir_obsid", nargs='?', type=str,
                        help='''Science Directory ObsID.    [REQUIRED]
                        Ex: \"3\" or \"102\"''', default=None)
                        
    parser.add_argument("-se","--scidir_expnum", nargs='?', type=str,
                        help='''Science Directory exposure number.
                        Ex: \"1\" or \"05\"''', default=None) 
                        
    parser.add_argument("--suboverscan_options", nargs="?", type=str, 
                        help='''Set subtract overscan options.
                        Default: \"-s -a -k 2.8 -t -z\"''', 
                        default="-s -a -k 2.8 -t -z")        
 
    parser.add_argument("--skysubtract_options", nargs="?", type=str, 
                        help='''Set sky subtraction options.
                        Default: \"-J -w 250 --output-both\".''', 
                        default="-J -w 250 --output-both")

    parser.add_argument("--fiberextract_options", nargs="?", type=str, 
                        help='''Set fiber extract options.
                        Default: \"-n 1032 -W 3500,5500 -P\".''', 
                        default="-n 1032 -W 3500,5500 -P")

    parser.add_argument("--makecube_options", nargs="?", type=str, 
                        help='''Set makecube options.
                        Default: \"\".''', 
                        default="")
                          
    args = parser.parse_args(args=argv)

    # Check that the arguments are filled
    if args.specid:
        args.specid = args.specid.replace(" ", "").split(',')
    else:
        args.specid = SPECID
        
    if args.run_insitu:
        args.run_insitu = True
    else:
        args.run_insitu = False
        
    if args.output is None:
        msg = 'No output directory was provided'
        parser.error(msg)    

    if args.cal_dir is None:
        msg = 'No calibration directory was provided'
        parser.error(msg)    
        
    if args.scidir_date is None:
        msg = 'No science directory date was provided'
        parser.error(msg) 
    else:
        args.scidir_date = args.scidir_date.replace(" ", "").split(',')

    if args.scidir_obsid is None:
        msg = 'No science directory ObsID was provided'
        parser.error(msg) 
    else:
        args.scidir_obsid = args.scidir_obsid.replace(" ", "").split(',')

    if args.scidir_expnum is not None:
        args.scidir_expnum = args.zrodir_expnum.replace(" ", "").split(',')
    
    args.sci_file_loc = []
    for date in args.scidir_date:
        for obsid in args.scidir_obsid:
            if args.scidir_expnum is not None:   
                for expnum in args.scidir_expnum:
                    args.sci_file_loc.append(op.join(args.rootdir, date, 
                                    args.instr, 
                                    "{:s}{:07d}".format(args.instr,int(obsid)), 
                                    "exp{:02d}".format(int(expnum))))
            else:
                args.sci_file_loc.append(op.join(args.rootdir, date, 
                                   args.instr,
                                   "{:s}{:07d}".format(args.instr,int(obsid))))
                
    return args  
    
class commandinfo(object):
    @classmethod
    def writecommand(cls, f, s):
        """Write the header to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        f.write('\n'.join(s) + "\n")

'''VirusFrame class.  This is critical to call cure commands from 
"call_cure.py", however this is a local definition due to the ever changing
nature of data to be read in.  Some have header key words set, others do not.
When finished, this will be a global class for all reductions, but for now,
I will just use it here, locally.'''
class VirusFrame:
    def __init__(self, initial_base = None, filename = None, dir_dict = None):
        '''
        Initializing a VirusFrame for a given filename.
        This includes reading the header and setting reduction keywords.
        From this, the file with have attributes that can be tracked.
        '''
        
        if filename is None:
            print("No filename was given for VirusFrame to be initialized.")
            return None
        if filename is None:
            print("No file type was given for VirusFrame to be initialized.")
            return None
        else:
            ######## OPEN AND KEEP FILES IN SOME SMART WAY ########
            self.filename             = filename
            self.origloc              = op.dirname(self.filename)
            self.basename, tmp1, tmp2 = op.basename(self.filename ).split('_')
            self.type                 = dir_dict
            self.ifuslot              = tmp1[0:3]
            self.time                 = self.basename.split('T')[1]
            self.hr                   = self.basename.split('T')[1][0:2]
            self.min                  = self.basename.split('T')[1][2:4]
            self.sec                  = self.basename.split('T')[1][4:8]
            self.clean                = CLEAN_AFTER_DONE

            
            ###### READ IN THE HEADER AND NECESSARY KEYWORDS ######
            self.trimsec    = {}
            self.biassec    = {}
            self.actionbase = {}
            for amp in SPEC:    
                rootname          = op.join(self.origloc, self.basename 
                                                          + '_' 
                                                          + self.ifuslot 
                                                          + amp 
                                                          + '_' 
                                                          + self.type 
                                                          + '.fits' )
                hdulist = fits.open(rootname)
                trimstr = re.split('[\[ \] ]',hdulist[0].header['TRIMSEC'])[1]
                biasstr = re.split('[\[ \] ]',hdulist[0].header['BIASSEC'])[1]
                self.trimsec[amp] = "\"" + trimstr + "\""
                self.biassec[amp] = "\"" + biasstr + "\""
                self.actionbase[amp] = ''     

            self.actionbase["L"] = initial_base  
            self.actionbase["R"] = initial_base  
            
            # Use mapping because if headers don't include SPECID
            if usemapping:
                self.specid = IFUSLOT_DICT[self.ifuslot][0]
                self.ifuid = IFUSLOT_DICT[self.ifuslot][1]

            else:
                self.specid = hdulist[0].header['SPECID']
                self.ifuid = hdulist[0].header['IFUID']

                
            self.object      = hdulist[0].header['OBJECT']
            self.orggain     = hdulist[0].header['GAIN']
            self.orgrdnoise  = hdulist[0].header['RDNOISE']
            self.exptime     = hdulist[0].header['EXPTIME']
                    
    def addbase(self, action, amp = None, side = None):
        s = []
        if self.clean:
            if amp is not None:
                filename   = op.join(self.origloc, self.actionbase[amp]
                                                   + self.basename 
                                                   + '_' 
                                                   + self.ifuslot 
                                                   + amp 
                                                   + '_' 
                                                   + self.type 
                                                   + '.fits')
                                                    
                filename_e = op.join(self.origloc, 'e.' 
                                                    + self.actionbase[amp] 
                                                    + self.basename 
                                                    + '_' 
                                                    + self.ifuslot 
                                                    + amp 
                                                    + '_' 
                                                    + self.type 
                                                    + '.fits') 
                                                    
                self.actionbase[amp] = action + self.actionbase[amp]
                
                if op.exists(filename):
                    os.remove(filename)
                s.append('{:s} '.format(filename))
                if op.exists(filename_e):
                    os.remove(filename_e)
                s.append('{:s} '.format(filename_e))

            if side is not None:
                filename = op.join(self.origloc, self.actionbase[side] 
                                                 + self.basename 
                                                 + '_' 
                                                 + self.ifuslot 
                                                 + '_' 
                                                 + self.type 
                                                 + '_' 
                                                 + side 
                                                 + '.fits')
                                                   
                filename_e = op.join(self.origloc, 'e.' 
                                                   + self.actionbase[side] 
                                                   + self.basename 
                                                   + '_' 
                                                   + self.ifuslot 
                                                   + '_' 
                                                   + self.type 
                                                   + '_' 
                                                   + side 
                                                   + '.fits')  
                                                   
                self.actionbase[side] = action + self.actionbase[side]
                
                if op.exists(filename):
                    os.remove(filename)
                s.append('{:s} '.format(filename))
                if op.exists(filename_e):
                    os.remove(filename_e) 
                s.append('{:s} '.format(filename_e))
        else:
            if amp is not None:
                self.actionbase[amp] = action + self.actionbase[amp]
                
            if side is not None:
                self.actionbase[side] = action + self.actionbase[side]
        
        return ''.join(s)
        
class ditherinfo(object):
    @classmethod
    def writeHeader(cls, f):
        """Write the header to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = []
        s.append("# basename is the base file name of the data files.")
        s.append("#          _{L,R}.fits is added for the left and right spectrographs")
        s.append("# modelbase is the base file name of the dist, fmod, pmode files corresponding to the data files")
        s.append("# $Id:$")
        s.append("#")
        s.append("# basename modelbase ditherx dithery seeing norm airmass")
        s.append("#")
        f.write('\n'.join(s) + "\n")

    @classmethod
    def writeDither(cls, f, filename, basename, ditherx, dithery, seeing, norm, airmass):
        """Write something to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = ("%s %s %4.2f %4.2f %4.2f %4.2f %4.2f" %
             (filename, basename, ditherx, dithery, seeing, norm, airmass))
        f.write(s)
        f.write("\n")
        f.flush()

def initial_setup(initial_base=None, file_loc_dir=None, redux_dir=None, 
                  DIR_DICT=None, uifuslot=None):
    '''
    Running the initial setup which includes:
    1) Building a standard reduction folder structure
    2) Copying files from specified location into the structure
    3) Creating class variable, VirusFrame, for each file which records
       keywords for other functions as well as time of the observation and 
       basename
    '''
    vframes = [] # list for VirusFrame objects
    if redux_dir is None:
        print("Please provide a reduction directory \"redux_dir\".")
        return None
    else:
        if not op.exists(redux_dir):
            os.mkdir(redux_dir)

    if initial_base is None:        
        print("WHOOPS! You forgot an initial base, even if it is just ''")
        return None

    if file_loc_dir is None:        
        print("Please provide a file location directory \"file_loc_dir\".")
        print("This is a list of the location directories for each file type.")
        return None
        
    if DIR_DICT is None:        
        print("Please provide a directory dictionary \"DIR_DICT\".")
        return None

    if uifuslot is None:        
        print("Please provide a IFUSLOT IDs for reduction (nested complaint).")
        return None
    
    '''Loop through the input directories and link files to the new structure.
       Create a VirusFrame class for each frame that can maintain info for 
       each original frame. The VirusFrame is critical for the rest of the 
       reduction pipeline.  Only gather VirusFrame objects for LL frames 
       as a single VirusFrame serves for all four amps.'''     
    for i in xrange(len(file_loc_dir)):
        if file_loc_dir[i] is not None:
            if not op.exists(op.join(redux_dir, DIR_DICT[i])):
                os.mkdir(op.join(redux_dir, DIR_DICT[i]))
            for file_loc in file_loc_dir[i]:
                # Example path: date/instr/instrOBSID/expnum/instr/files"
                # file_loc hold the date, instr, and instrOBSID section
                basename = op.basename(file_loc)
                if basename[0:3]=="exp":
                    filenames = glob.glob(op.join(file_loc, "*/*.fits")) 
                    if not filenames:
                        print("No files found for wildcard search: %s" 
                                        %(op.join(file_loc, "*/*.fits")))
                        return None
                else:
                    filenames = glob.glob(op.join(file_loc, "*/*/*.fits"))                  
                    if not filenames:
                        print("No files found for wildcard search: %s" 
                                        %(op.join(file_loc, "*/*/*.fits")))
                        return None
                for f in filenames:
                    temp, temp1, temp2 = op.basename(f).split('_')
                    newname = temp+'_'+temp1+'_'+DIR_DICT[i]+'.fits'
                    test = op.exists(op.join(redux_dir, DIR_DICT[i], 
                                             newname))
                    if not test:
                        if temp1[0:3] in uifuslot:
                            os.symlink(f, op.join(redux_dir, DIR_DICT[i], 
                                                  newname))
                for f in filenames:
                    temp, temp1, temp2 = op.basename(f).split('_')
                    amp = temp1[-2:]
                    newname = temp+'_'+temp1+'_'+DIR_DICT[i]+'.fits'
                    if amp == "LL":
                        if temp1[0:3] in uifuslot:
                            vframes.append(VirusFrame(initial_base,
                                                      op.join(redux_dir, 
                                                              DIR_DICT[i],
                                                              newname),
                                                              DIR_DICT[i]))
    
    return vframes          
                
def clean_folder(frames, cmd, side=None, amp=None):
    s = ['rm ']
    for f in frames:
        if f.clean:
            if side:
                for side in SPECBIG:
                    filename   = op.join(f.origloc, f.actionbase[side] 
                                                    + f.basename 
                                                    + '_' 
                                                    + f.ifuslot 
                                                    + '_' 
                                                    + f.type 
                                                    + '_' 
                                                    + side 
                                                    + '.fits')
                    filename_e = op.join(f.origloc, 'e.' 
                                                    + f.actionbase[side] 
                                                    + f.basename 
                                                    + '_' 
                                                    + f.ifuslot 
                                                    + '_' 
                                                    + f.type 
                                                    + '_' 
                                                    + side 
                                                    + '.fits')
                    if op.exists(filename):
                        os.remove(filename)
                    if op.exists(filename_e):
                        os.remove(filename_e)
                    s.append("{:s} {:s} ".format(filename,filename_e))                    
            if amp:
                for sp in SPEC:
                    filename   = op.join(f.origloc, f.actionbase[sp] 
                                                    + f.basename 
                                                    + '_' 
                                                    + f.ifuslot 
                                                    + sp
                                                    + '_' 
                                                    + f.type 
                                                    + '.fits')
                    filename_e = op.join(f.origloc, 'e.' 
                                                    + f.actionbase[sp] 
                                                    + f.basename 
                                                    + '_' 
                                                    + f.ifuslot 
                                                    + sp
                                                    + '_' 
                                                    + f.type 
                                                    + '.fits')
                    if op.exists(filename):
                        os.remove(filename)
                    if op.exists(filename_e):
                        os.remove(filename_e)
                    s.append("{:s} {:s} ".format(filename,filename_e))
    cmd.append(''.join(s))                    
    return cmd
    
def build_name(frame, side=None, amp=None, error=None, red_dir=None):
    if side:
        if not error:
            filename = op.join(frame.origloc, frame.actionbase[side] 
                                              + frame.basename 
                                              + '_' 
                                              + frame.ifuslot 
                                              + '_' 
                                              + frame.type 
                                              + '_' 
                                              + side 
                                              + '.fits')
        else:
            filename = op.join(frame.origloc, 'e.'
                                               + frame.actionbase[side] 
                                               + frame.basename 
                                               + '_' 
                                               + frame.ifuslot 
                                               + '_' 
                                               + frame.type 
                                               + '_' 
                                               + side 
                                               + '.fits')
    if red_dir:
        return op.join(red_dir, filename)
    else:
        return filename
        
def build_mastername(mastername, side, uca, red_dir, lamp=None):
    if lamp:
        filename = op.join(red_dir, mastername 
                                    + '_' 
                                    + lamp 
                                    + '_' 
                                    + uca 
                                    + '_' 
                                    + side 
                                    + '.fits')
    else:
        filename = op.join(red_dir, mastername 
                                    + '_' 
                                    + uca 
                                    + '_' 
                                    + side 
                                    + '.fits')  
    return filename
    
def flatten(x, ans):
    for i in x:
        if isinstance(i,list):
            ans = flatten(i, ans)
        else:
            ans.append(i)
    return ans
    
def main():
    '''
    Can run the basic reduction which includes:
    1) Make error frames from Readnoise/Gain (in ADUs)
    2) Overscan subtract
    3) Trim Overscan
    4) Subtract master bias from sci frames
    5) Ccdcombine frames which puts units in e-
    6) Add photon noise to combined frames

    Also, can run:
    1) Sky subtraction
    2) Fiber Extraction
    '''
    # Read in the arguments
    args = parse_args()
    # Reduction Directory
    redux_dir = args.output
    # file directories
    file_loc_dir = [args.sci_file_loc] # Keep same order with dir_dict
    DIR_DICT = {0:sci_dir} # Dictionary for looping through directory names
    # Set up VirusFrame objects for each file
    # Unique camera id's (aka SPECID)
    uifuslot = []
    ucam = args.specid
    for ifuslot, values in IFUSLOT_DICT.items():
        specid = values[0]
        if specid in ucam:    
            uifuslot.append(ifuslot)
    if args.skipbasic:
        vframes = initial_setup(args.skipbasic, file_loc_dir, redux_dir, 
                                DIR_DICT,uifuslot)
    else:
        vframes = initial_setup('', file_loc_dir, redux_dir, DIR_DICT,uifuslot)

    # Pooling the different frames for different purposes
    vframes = [v for v in vframes for uca in ucam if v.specid == uca]
    # Get the trimsec and biassec for each frame
    # Grab from the first frame and assume it is the same for all
    trimsec = vframes[0].trimsec["LL"] 
    biassec = vframes[0].biassec["LL"] 
    
    # loop over SPECID
    commandfile = []
    for ind,uca in enumerate(ucam):
        vframesselect = [v for v in vframes if v.specid == uca]                  

        commandfile_fn = op.join(redux_dir, "commands_{:s}".format(uca))
        print("Creating {:s}".format(commandfile_fn))
        commandfile.append(open(commandfile_fn, 'w'))
        cmd = []
        
        if not args.skipbasic:
            for amp in SPEC:
                # Make Error Frame
                cmd = CC.mkerrorframe(vframesselect, amp, cmd, 
                                      run=args.run_insitu) 
    
                # Subtract Overscan
                cmd = CC.subtractoverscan(biassec, args.suboverscan_options,
                                          vframesselect, amp, cmd, 
                                          run=args.run_insitu)  
                
                # Trim amplifier
                cmd = CC.extractfits(trimsec, vframesselect, amp, cmd,
                                     run=args.run_insitu)
    
                # Subtract Master Bias
                masterbiasname = op.join(biasdir, 'masterbias' 
                                                    + '_' 
                                                    + uca 
                                                    + '_' 
                                                    + amp 
                                                    + '.fits')
                F = fits.open(masterbiasname)
                F[0].data *= args.bias1_mult_val
                nmasterbiasname = op.join(redux_dir, 
                                          op.basename(masterbiasname))
                F.writeto(nmasterbiasname,clobber=True)
                cmd = CC.subtractbias(vframesselect, nmasterbiasname, amp, cmd, 
                                      run=args.run_insitu)
                if args.use_darks:
                    masterbasename = ('masterdark' + '_' + uca + '_' + amp 
                                         + '.fits') 
                    masterdarkname = op.join(darkcurrentdir, masterbasename) 
                    masterdarknamenew = op.join(redux_dir, masterbasename)

                    hdu = fits.open(masterdarkname)
                    if amp[1]=='LL':
                        hdu[0].data[:] = args.dark1_mult_val*hdu[0].data
                    if amp[1]=='LU':
                        hdu[0].data[:] = args.dark2_mult_val*hdu[0].data
                    if amp[1]=='RL':
                        hdu[0].data[:] = args.dark3_mult_val*hdu[0].data
                    if amp[1]=='RU':
                        hdu[0].data[:] = args.dark4_mult_val*hdu[0].data 
                       
                    hdu.writeto(masterdarknamenew, clobber=True)
                    cmd = CC.subtractbias(vframesselect, masterdarknamenew, 
                                          amp, cmd, run=args.run_insitu)
            # CCDCombine (multiplies by gain and joins amplifiers)
            cmd = CC.ccdcombine(vframesselect, cmd, run=args.run_insitu)
            '''Change from amps to side convention'''
            for v in vframesselect:
                v.actionbase["L"] = v.actionbase["LL"]
                v.actionbase["R"] = v.actionbase["RL"]
            # Clean the "sci" files that aren't joined
            if args.run_insitu:
                cmd = clean_folder(vframesselect, cmd, amp="LL")
                
            # Add photon noise
            for side in SPECBIG:    
                cmd = CC.addphotonnoise(vframesselect, side, cmd, 
                                        run=args.run_insitu)
                # Subtract Dark using PYFITS and not CURE
                # Change to CURE when implimented more seriously

#                if args.use_darks:
#                    if side in SPECBIG:
#                        masterbasename = ('masterdark' + '_' + uca + '_' + side 
#                                         + '.fits') 
#                        masterdarkname = op.join(darkcurrentdir, masterbasename) 
#                        masterdarknamenew = op.join(redux_dir, masterbasename) 
#                        if not op.exists(masterdarknamenew):
#                            shutil.copy(masterdarkname,masterdarknamenew)
#                        hdu = fits.open(masterdarknamenew)
#                        if args.dark_mult_val:
#                            for v in vframesselect:
#                                hdu1 = fits.open(build_name(v, side=side))
#                                hdu1[0].data = (hdu1[0].data 
#                                                - hdu[0].data * args.dark_mult_val)
#                                hdu1[0].header['HISTORY'] = 'subtracted:'
#                                hdu1[0].header['HISTORY'] = ('file: %s' 
#                                                            ' multiplied by %0.3f'
#                                                            % (masterbasename,
#                                                               args.dark_mult_val))
#                                with warnings.catch_warnings():
#                                    warnings.simplefilter("ignore")
#                                    hdu1.writeto(build_name(v, side=side),
#                                                 clobber=True)
#                                hdu1.close()
#                        else:
#                            # CRITICAL LINES FOR DARK SUBTRACTION
#                            mdd = np.median(hdu[0].data[1028:1036,416:616])
#                            mdv = []
#                            for v in vframesselect:
#                                hdu1 = fits.open(build_name(v, side=side))
#                                # CRITICAL LINES FOR DARK SUBTRACTION
#                                mdv = np.median(hdu1[0].data[1028:1036,416:616])
#                                hdu1[0].data = hdu1[0].data - hdu[0].data*mdv/mdd
#                                hdu1[0].header['HISTORY'] = 'subtracted:'
#                                hdu1[0].header['HISTORY'] = ('file: %s' 
#                                                            ' multiplied by %0.3f'
#                                                            % (masterbasename,
#                                                               mdv/mdd))
#                                with warnings.catch_warnings():
#                                    warnings.simplefilter("ignore")
#                                    hdu1.writeto(build_name(v, side=side),
#                                                 clobber=True)
#                                hdu1.close()
                
                opts = '-f ' + op.join(configdir, 'PixelFlats', 
                                       'pixelflat_cam%03d_%s.fits' %(int(uca), 
                                                                     side))
                cmd = CC.dividepixelflat(vframesselect, opts, side, cmd)               
        
        if args.use_twi:
            basetrace = 'mastertrace_twi'
        else:
            basetrace = 'mastertrace'
        
        if args.use_deformer_default:
            ifolder = op.join(configdir,"DeformerDefaults")
            basetrace = 'mastertrace_twi'
        else:
            ifolder = args.cal_dir

        if args.rerun_deformer:
            solar = op.join(configdir,'solar_spec','sun.spec')
            opts =  ("--debug --template-spectrum %s --corr-dw 100.0 "
                    "--corr-max-shift 0.5"  %(solar))
	    opts = ("--debug --max-trace-drift 3.0")
            for side in SPECBIG:
                distmodel = op.join(ifolder, basetrace 
                                                  + '_' 
                                                  + uca 
                                                  + '_' 
                                                  + side 
                                                  + '.dist')
                fibermodel = op.join(ifolder, basetrace 
                                                    + '_' 
                                                    + uca 
                                                    + '_' 
                                                    + side 
                                                    + '.fmod')
    
                CC.deformer(vframesselect, args.output, opts, cmd, rerun=True, 
                         dist=distmodel, fmod=fibermodel, run=args.run_insitu,
                         side=side)

        # Run sky subtraction            
        if args.subsky:  
            for side in SPECBIG:
                if args.rerun_deformer:
                    ifolder = args.output
                    add = '_adj'
                else:
                    add = ''
                distmodel = op.join(ifolder, basetrace 
                                                  + '_' 
                                                  + uca 
                                                  + '_' 
                                                  + side 
                                                  + add
                                                  + '.dist')
                fibermodel = op.join(ifolder, basetrace 
                                                    + '_' 
                                                    + uca 
                                                    + '_' 
                                                    + side 
                                                    + '.fmod')
                cmd = CC.subtractsky(vframesselect, side, distmodel, 
                                     fibermodel, args.skysubtract_options,
                                     cmd, run=args.run_insitu)  
                                     
        if args.remove_cosmics:
            for side in SPECBIG:
                for v in vframesselect:
                    filename = build_name(v, side=side) 
                    remove_cosmics(filename)
                    error_frame = op.join(op.dirname(filename),
                                          'e.c' + op.basename(filename))
                    new_error_frame = op.join(op.dirname(filename),
                                          'e.' + op.basename(filename))
                    shutil.copy(error_frame, new_error_frame)  
                          
        # Run fiberextract
        if args.fiberextract:
            for side in SPECBIG:                
                distmodel = op.join(ifolder, basetrace 
                                                  + '_' 
                                                  + uca 
                                                  + '_' 
                                                  + side 
                                                  + '.dist')
                if op.exists(distmodel):
                    fibermodel = op.join(ifolder, basetrace 
                                                       + '_' 
                                                       + uca 
                                                       + '_' 
                                                       + side 
                                                       + '.fmod')
                    cmd = CC.fiberextract(vframesselect, side, distmodel, 
                                          fibermodel, 
                                          args.fiberextract_options, cmd,
                                          run=args.run_insitu, 
                                          new_deformer=args.rerun_deformer,
                                          cal_folder=args.output)
        
        # Run mkcube
        if args.makecube:
            IFUcen_dir = op.join(configdir, 'IFUcen_files')
            IFUcen_file = 'IFUcen_VIFU{:s}.txt'.format(CAM_IFU_DICT[uca])
            IFUcen_file_fn = op.join(IFUcen_dir, IFUcen_file)
            shutil.copy(op.join(IFUcen_dir, IFUcen_file), 
                        op.join(redux_dir, IFUcen_file))
            cmd.append('cp {:s} {:s}'.format(op.join(IFUcen_dir, IFUcen_file), 
                                             op.join(redux_dir, IFUcen_file)))            
            side = SPECBIG[0]
            for v in vframesselect:
                ditherfile_fn = op.join(op.abspath(redux_dir), 
                                       "dither_{:s}_{:s}.txt".format(v.ifuslot, 
                                                                   v.basename))
                print("Creating {:s}".format(ditherfile_fn))
                ditherfile = open(ditherfile_fn, 'w')
                ditherinfo.writeHeader(ditherfile)        
                dither_fn = op.join(op.abspath(redux_dir), sci_dir, 
                                    "{:s}{:s}_{:s}_{:s}".format(
                                                        v.actionbase[side], 
                                                        v.basename, v.ifuslot, 
                                                        v.type))
                modelbase = op.join(op.abspath(args.cal_dir), 
                                    "{:s}_{:s}".format(basetrace, v.specid))
                if op.exists(modelbase+"_L.dist"): 
                    if op.exists(modelbase+"_R.dist"):
                        ditherinfo.writeDither(ditherfile, dither_fn, 
                                               modelbase, 0.00, 0.00, 1.50,
                                               1.00, 1.22) 
                        cmd = CC.mkcube(IFUcen_file_fn, ditherfile_fn, 
                                        args.makecube_options, cmd, 
                                        run=args.run_insitu)  

        # Run detect
        if args.detect:
            IFUcen_dir = op.join(configdir, 'IFUcen_files')
            IFUcen_file = 'IFUcen_VIFU{:s}.txt'.format(CAM_IFU_DICT[uca])
            IFUcen_file_fn = op.join(IFUcen_dir, IFUcen_file)
            shutil.copy(op.join(IFUcen_dir, IFUcen_file), 
                        op.join(redux_dir, IFUcen_file))
            cmd.append('cp {:s} {:s}'.format(op.join(IFUcen_dir, IFUcen_file), 
                                             op.join(redux_dir, IFUcen_file)))            
            side = SPECBIG[0]
            for v in vframesselect:
                output_fn = op.join(op.abspath(redux_dir), 
                                       "detect_{:s}_{:s}".format(v.ifuslot, 
                                                                   v.basename))
                ditherfile_fn = op.join(op.abspath(redux_dir), 
                                       "dither_{:s}_{:s}.txt".format(v.ifuslot, 
                                                                   v.basename))
                print("Creating {:s}".format(ditherfile_fn))
                ditherfile = open(ditherfile_fn, 'w')
                ditherinfo.writeHeader(ditherfile)        
                dither_fn = op.join(op.abspath(redux_dir), sci_dir, 
                                    "{:s}{:s}_{:s}_{:s}".format(
                                                        v.actionbase[side], 
                                                        v.basename, v.ifuslot, 
                                                        v.type))
                modelbase = op.join(op.abspath(args.cal_dir), 
                                    "{:s}_{:s}".format(basetrace, v.specid))
                if op.exists(modelbase+"_L.dist") and op.exists(modelbase+"_R.dist"):
                    ditherinfo.writeDither(ditherfile, dither_fn, modelbase,
                                           0.00, 0.00, 1.50, 1.00, 1.22) 
                    cmd = CC.detect(IFUcen_file_fn, ditherfile_fn, output_fn,
                                    args.detect_options, cmd, 
                                    run=args.run_insitu)  
            
        cmd = flatten(cmd,[])                             
        commandinfo.writecommand(commandfile[-1],cmd)
        commandfile[-1].close()
        
if __name__ == '__main__':
    main()  
