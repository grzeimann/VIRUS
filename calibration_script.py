# -*- coding: utf-8 -*-
"""

Calibration script for "quick" reduction.

The script requires an output directory.  It also requires key information
for selecting the zero, comp, flat, dark, and or twighlight directories.  
It is required to input the date and obsid for the zero, comp, and flat 
folders.
@author: gregz
"""
from __future__ import print_function
import argparse as ap
import numpy as np
import pyfits
import shutil
import glob
import os
import os.path as op
import re
import call_cure as CC
import textwrap

'''Global variables to be set for a given reduction'''
LAMP_DICT            = {0:'Hg',1:'Cd'}
usemapping           = True # manual map IFUSLOT to SPECID 
CLEAN_AFTER_DONE     = True # clean the previous file after new one is created
initial_base         = ''   # used to start reductions at a given place

pixflatdir   = "/Users/gregz/cure/virus_early/PixelFlats"
configdir    = "/home/03946/hetdex/virus_config"

# Critical naming scheme for folders and object types
# Note the type doesn't have to initial match the folder because
# the code will change the type based on the folder it is sent.
zro_dir  = "zro"
flt_dir  = "flt"
cmp_dir  = "cmp"
drk_dir  = "drk"
twi_dir  = "twi"

# Dictionary of the mapping between IFUSLOT and SPECID, IFUID
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
                '105':['012','055'],
                '106':['017','022'],
                '999':['031','999']}

# Dictionary of the mapping between SPECID and IFUID

CAM_IFU_DICT = {'004':'023',
                '037':'024',
                '027':'001',                
                '047':'016',
                '024':'051',
                '013':'019',
                '016':'026',
                '041':'015',
                '051':'033',
                '008':'054',
                '025':'020',
                '038':'014',
                '020':'004',
                '032':'028',
                '012':'055',
                '017':'022',
                '031':'999'}

# Default set of spectrographs for reduction
SPECID = ["004","008","012","013","016","017","020","024","025","027","032","037","038","041","047","051"]


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
    description = textwrap.dedent('''Calibration Reduction Script - 
    
                     This script produces the following output:
                         masterbias_{SPECID}_{AMP}.fits
                         masterarc_{LAMP}_{SPECID}_{SIDE}.fits
                         masterarc_{SPECID}_{SIDE}.fits
                         mastertrace_{SPECID}_{SIDE}.fits
                         mastertrace_{SPECID}_{SIDE}.dist
                         mastertrace_{SPECID}_{SIDE}.fmod
                         mastertrace_{SPECID}_{SIDE}.pmod
                     
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

    parser.add_argument("--instr", nargs='?', type=str, 
                        help='''Instrument to process. 
                        Default: "virus"
                        Ex: "camra" for lab data.''', default = "virus")
                        
    parser.add_argument("--rootdir", nargs='?', type=str, 
                        help='''Root Directory
                        Default: \"/work/03946/hetdex/maverick\"''', 
                        default="/work/03946/hetdex/maverick")
                        
    parser.add_argument("-mb","--make_masterbias", 
                        help='''Make Masterbias''',
                        action="count", default=0) 

    parser.add_argument("-md","--make_masterdark", 
                        help='''Make Masterdark''',
                        action="count", default=0) 

    parser.add_argument("-ma","--make_masterarc", 
                        help='''Make Masterarc''',
                        action="count", default=0) 

    parser.add_argument("-mt","--make_mastertrace", 
                        help='''Make Mastertrace''',
                        action="count", default=0) 
 
    parser.add_argument("-rd","--run_deformer", 
                        help='''Run Deformer''',
                        action="count", default=0)
                        
    parser.add_argument("--use_dist_file", nargs='?', type=str, 
                        help='''Dist file for using twighlight in Deformer.
                        Ex: "calib/20160504/mastertrace_twi_L.dist"''', 
                        default=None)

    parser.add_argument("--temp_spec", nargs='?', type=str, 
                        help='''Template solar spectrum for deformer.
                        Ex: "/home/03946/hetdex/virus_config/solar_spec/sun.spec""''', 
                        default="/home/03946/hetdex/virus_config/solar_spec/sun.spec")
                       
    parser.add_argument("--use_other_masterbias", nargs='?', type=str, 
                        help='''Master Bias Directory.
                        Ex: "calib/20160504"''', 
                        default=None)
                        
    parser.add_argument("-zd","--zrodir_date", nargs='?', type=str,
                        help='''Zero Directory Date.        [REQUIRED]
                        Ex: \"20160412\"''', default = None)

    parser.add_argument("-zo","--zrodir_obsid", nargs='?', type=str,
                        help='''Zero Directory ObsID.       [REQUIRED]
                        Ex: \"3\" or \"102\"''', default = None)
                        
    parser.add_argument("-ze","--zrodir_expnum", nargs='?', type=str,
                        help='''Zero Directory exposure number.
                        Ex: \"1\" or \"05\"''', default = None) 
 
    parser.add_argument("-dd","--drkdir_date", nargs='?', type=str,
                        help='''Dark Directory Date.       
                        Ex: \"20160412\"''', default = None)

    parser.add_argument("-do","--drkdir_obsid", nargs='?', type=str,
                        help='''Dark Directory ObsID.       
                        Ex: \"3\" or \"102\"''', default = None)
                        
    parser.add_argument("-de","--drkdir_expnum", nargs='?', type=str,
                        help='''Dark Directory exposure number.
                        Ex: \"1\" or \"05\"''', default = None) 

    parser.add_argument("-td","--twidir_date", nargs='?', type=str,
                        help='''Twighlight Directory Date.       
                        Ex: \"20160412\"''', default = None)

    parser.add_argument("-to","--twidir_obsid", nargs='?', type=str,
                        help='''Twighlight Directory ObsID.       
                        Ex: \"3\" or \"102\"''', default = None)
                        
    parser.add_argument("-te","--twidir_expnum", nargs='?', type=str,
                        help='''Twighlight Directory exposure number.
                        Ex: \"1\" or \"05\"''', default = None) 
                        
    parser.add_argument("-cd","--cmpdir_date", nargs='?', type=str,
                        help='''Comp Directory Date.        [REQUIRED]
                        Ex: \"20160412\"''', default = None)

    parser.add_argument("-co","--cmpdir_obsid", nargs='?', type=str,
                        help='''Comp Directory ObsID.       [REQUIRED]
                        Ex: \"3\" or \"102\"''', default = None)
                        
    parser.add_argument("-ce","--cmpdir_expnum", nargs='?', type=str,
                        help='''Comp Directory exposure number.
                        Ex: \"1\" or \"05\"''', default = None) 

    parser.add_argument("-fd","--fltdir_date", nargs='?', type=str,
                        help='''Flat Directory Date.        [REQUIRED]
                        Ex: \"20160412\"''', default = None)

    parser.add_argument("-fo","--fltdir_obsid", nargs='?', type=str,
                        help='''Flat Directory ObsID.       [REQUIRED]
                        Ex: \"3\" or \"102\"''', default = None)
                        
    parser.add_argument("-fe","--fltdir_expnum", nargs='?', type=str,
                        help='''Flat Directory exposure number.
                        Ex: \"1\" or \"05\"''', default = None) 

    parser.add_argument("-t","--twi_deformer", 
                        help='''Run deformer using master twighlight.''',
                        action="count", default=0)
   
    parser.add_argument("--deformer_options", nargs="?", type=str, 
                        help='''Set deformer options.
                        Default: \"-s 6 --debug --dump_psf_data\".''', 
                        default="-s 6 --debug --dump_psf_data")

    parser.add_argument("--twideformer_options", nargs="?", type=str, 
                        help='''Set twilight deformer options.
                        Default: \"-s 6 --debug --dump_psf_data\".''', 
                        default="-s 6 --debug --dump_psf_data")

    parser.add_argument("--headfits_options", nargs="?", type=str, 
                        help='''Set headfits options for normalization step.
                        Default: \"-m -k EXPTIME -v 1 -V \"normalized\"\"''', 
                        default="-m -k EXPTIME -v 1 -V \"normalized\"")

    parser.add_argument("--suboverscan_options", nargs="?", type=str, 
                        help='''Set subtract overscan options.
                        Default: \"-s -a -k 2.8 -t -z\"''', 
                        default="-s -a -k 2.8 -t -z")        

    parser.add_argument("--masterbias_options", nargs="?", type=str, 
                        help='''Set meanfits options for masterbias.
                        Default: \"-m -k 3.0\"''', 
                        default="-m -k 3.0")   

    parser.add_argument("--masterdark_options", nargs="?", type=str, 
                        help='''Set meanfits options for masterdark.
                        Default: \"-m -k 3.0\"''', 
                        default="-m -k 3.0")        

    parser.add_argument("--masterarc_options", nargs="?", type=str, 
                        help='''Set meanfits options for masterarc.
                        Default: \"--maxmem 1024 -s -t -m -k 2.8\"''', 
                        default="--maxmem 1024 -s -t -m -k 2.8")  

    parser.add_argument("--mastertrace_options", nargs="?", type=str, 
                        help='''Set meanfits options for mastertrace.
                        Default: \"--maxmem 1024 -s -t -m -k 2.8\"''',                         
                        default="--maxmem 1024 -s -t -m -k 2.8")
                          
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
        
    if not args.make_masterbias:
        if not args.make_masterarc:
            if not args.make_mastertrace:
                if not args.run_deformer:
                    print("Running masterbias, masterdark, masterarc, mastertrace, and deformer")
                    args.make_masterbias = 1
                    args.make_masterarc = 1
                    args.make_mastertrace = 1
                    args.run_deformer = 1

    if args.make_masterbias:    
        if args.zrodir_date is None:
            msg = 'No zero directory date was provided'
            parser.error(msg) 
        else:
            args.zrodir_date = args.zrodir_date.replace(" ", "").split(',')
    
        if args.zrodir_obsid is None:
            msg = 'No zero directory ObsID was provided'
            parser.error(msg) 
        else:
            args.zrodir_obsid = args.zrodir_obsid.replace(" ", "").split(',')
    
        if args.zrodir_expnum is not None:
            args.zrodir_expnum = args.zrodir_expnum.replace(" ", "").split(',')
        
        args.zro_file_loc = []
        for date in args.zrodir_date:
            for obsid in args.zrodir_obsid:
                if args.zrodir_expnum is not None:   
                    for expnum in args.zrodir_expnum:
                        args.zro_file_loc.append(op.join(args.rootdir, date, 
                                        args.instr, 
                                        "{:s}{:07d}".format(args.instr,int(obsid)), 
                                        "exp{:02d}".format(int(expnum))))
                else:
                    args.zro_file_loc.append(op.join(args.rootdir, date, 
                                       args.instr,
                                       "{:s}{:07d}".format(args.instr,int(obsid))))
    else:
        args.zro_file_loc = None

    if args.make_masterdark:
        if args.drkdir_date is None:
            msg = 'No dark directory date was provided'
            parser.error(msg)             
        else:
            args.drkdir_date = args.drkdir_date.replace(" ", "").split(',')
    
        if args.drkdir_obsid is None:
            msg = 'No dark directory ObsID was provided'
            parser.error(msg) 
        else:
            args.drkdir_obsid = args.drkdir_obsid.replace(" ", "").split(',')
    
        if args.drkdir_expnum is not None:
            args.drkdir_expnum = args.drkdir_expnum.replace(" ", "").split(',')
            
        args.drk_file_loc = []
        for date in args.drkdir_date:
            for obsid in args.drkdir_obsid:
                if args.drkdir_expnum is not None:   
                    for expnum in args.drkdir_expnum:
                        args.drk_file_loc.append(op.join(args.rootdir, date, 
                                        args.instr, 
                                        "{:s}{:07d}".format(args.instr,int(obsid)), 
                                        "exp{:02d}".format(int(expnum))))
                else:
                    args.drk_file_loc.append(op.join(args.rootdir, date, 
                                       args.instr,
                                       "{:s}{:07d}".format(args.instr,int(obsid))))
    else:
        args.drk_file_loc = None    

    if args.make_mastertrace:
        if args.twidir_date is None:
            print("No twighlight date provided so skipping twighlights.")
            args.twi_file_loc = None    
        else:
            args.twidir_date = args.twidir_date.replace(" ", "").split(',')

            if args.twidir_obsid is None:
                msg = 'No dark directory ObsID was provided'
                parser.error(msg) 
            else:
                args.twidir_obsid = args.twidir_obsid.replace(" ", "").split(',')
    
            if args.twidir_expnum is not None:
                args.twidir_expnum = args.twidir_expnum.replace(" ", "").split(',')
                
            args.twi_file_loc = []    
            for date in args.twidir_date:
                for obsid in args.twidir_obsid:
                    if args.twidir_expnum is not None:   
                        for expnum in args.twidir_expnum:
                            args.twi_file_loc.append(op.join(args.rootdir, date, 
                                        args.instr, 
                                        "{:s}{:07d}".format(args.instr,int(obsid)), 
                                        "exp{:02d}".format(int(expnum))))
                    else:
                        args.twi_file_loc.append(op.join(args.rootdir, date, 
                                       args.instr,
                                       "{:s}{:07d}".format(args.instr,int(obsid))))
    else:
        args.twi_file_loc = None    
       

    if args.make_masterarc:
        if args.cmpdir_date is None:
            msg = 'No cmp directory date was provided'
            parser.error(msg) 
        else:
            args.cmpdir_date = args.cmpdir_date.replace(" ", "").split(',')
    
        if args.cmpdir_obsid is None:
            msg = 'No cmp directory ObsID was provided'
            parser.error(msg) 
        else:
            args.cmpdir_obsid = args.cmpdir_obsid.replace(" ", "").split(',')
    
        if args.cmpdir_expnum is not None:
            args.cmpdir_expnum = args.cmpdir_expnum.replace(" ", "").split(',')
        
        args.cmp_file_loc = []
        for date in args.cmpdir_date:
            for obsid in args.cmpdir_obsid:
                if args.cmpdir_expnum is not None:   
                    for expnum in args.cmpdir_expnum:
                        args.cmp_file_loc.append(op.join(args.rootdir, date, 
                                        args.instr, 
                                        "{:s}{:07d}".format(args.instr,int(obsid)), 
                                        "exp{:02d}".format(int(expnum))))
                else:
                    args.cmp_file_loc.append(op.join(args.rootdir, date, 
                                       args.instr,
                                       "{:s}{:07d}".format(args.instr,int(obsid))))
    else:
        args.cmp_file_loc = None                       
    
    if args.make_mastertrace:
        if args.fltdir_date is None:
            msg = 'No flt directory date was provided'
            parser.error(msg) 
        else:
            args.fltdir_date = args.fltdir_date.replace(" ", "").split(',')
    
        if args.fltdir_obsid is None:
            msg = 'No flt directory ObsID was provided'
            parser.error(msg) 
        else:
            args.fltdir_obsid = args.fltdir_obsid.replace(" ", "").split(',')
    
        if args.fltdir_expnum is not None:
            args.fltdir_expnum = args.fltdir_expnum.replace(" ", "").split(',')
        
        args.flt_file_loc = []
        for date in args.fltdir_date:
            for obsid in args.fltdir_obsid:
                if args.fltdir_expnum is not None:   
                    for expnum in args.fltdir_expnum:
                        args.flt_file_loc.append(op.join(args.rootdir, date, 
                                        args.instr, 
                                        "{:s}{:07d}".format(args.instr,int(obsid)), 
                                        "exp{:02d}".format(int(expnum))))
                else:
                    args.flt_file_loc.append(op.join(args.rootdir, date, 
                                       args.instr,
                                       "{:s}{:07d}".format(args.instr,int(obsid))))
    else:
        args.flt_file_loc = None
                
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
    def __init__ ( self, filename = None, dir_dict = None):
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
                hdulist = pyfits.open(rootname)
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
            else:
                self.specid = hdulist[0].header['SPECID']
                
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


               

def initial_setup(file_loc_dir=None, redux_dir=None, DIR_DICT=None, 
                  uifuslot=None):
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
                            vframes.append(VirusFrame(op.join(redux_dir, 
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
    Running the basic calibration reduction which includes:
    1) Make error frames from Readnoise/Gain (in ADUs)
    2) Overscan subtract
    3) Trim Overscan
    4) Create master bias frame
    5) Subtract master bias from cmp/flt/drk frames
    7) Ccdcombine frames which puts units in e-
    8) Add photon noise to combined frames
    9) Create master dark frames (e-)
    10) Normalize cmps and flts 
    11) Combine cmps and flts into masterarc and mastertrace
    12) Run deformer
    '''
    # Read in the arguments
    args = parse_args()
    # Reduction Directory
    redux_dir = args.output
    # file directories
    file_loc_dir = [args.zro_file_loc, 
                args.drk_file_loc, 
                args.cmp_file_loc, 
                args.flt_file_loc,
                args.twi_file_loc] # Keep same order with directory dictionary
    DIR_DICT = {0:zro_dir,    
                1:drk_dir,    
                2:cmp_dir,    
                3:flt_dir,
                4:twi_dir} # Dictionary for looping through directory names

    # Unique camera id's (aka SPECID)
    ucam = args.specid
    uifuslot = []
    for ifuslot, values in IFUSLOT_DICT.items():
        specid = values[0]
        if specid in ucam:    
            uifuslot.append(ifuslot)
    if args.make_masterbias or args.make_masterdark or args.make_masterarc or args.make_mastertrace:
        # Set up VirusFrame objects for each file
        
        vframes = initial_setup(file_loc_dir, redux_dir, DIR_DICT, uifuslot)
    
        # Pooling the different frames for different purposes
        vframes = [v for v in vframes for uca in ucam if v.specid == uca]
        oframes = [v for v in vframes if v.type != "zro"] 
        dframes = [v for v in vframes if v.type == "drk"] 
        cframes = [v for v in vframes if v.type == "flt" or v.type == "cmp"] 
        lframes = [v for v in vframes if v.type == "cmp"] 
        fframes = [v for v in vframes if v.type == "flt"] 
        tframes = [v for v in vframes if v.type == "twi"] 
    
    
        # Get the trimsec and biassec for each frame
        # Grab from the first frame and assume it is the same for all
        trimsec = vframes[0].trimsec["LL"] 
        biassec = vframes[0].biassec["LL"] 
        
        # loop over SPECID
        commandfile = []
        for uca in ucam:
            vframesselect = [v for v in vframes if v.specid == uca]                 
            zframesselect = [v for v in vframes if v.type == "zro" 
                                                and v.specid == uca]
            oframesselect = [o for o in oframes if o.specid == uca] 
            dframesselect = [d for d in dframes if d.specid == uca] 
            cframesselect = [c for c in cframes if c.specid == uca] 
            lframesselect = [l for l in lframes if l.specid == uca] 
            fframesselect = [f for f in fframes if f.specid == uca] 
            tframesselect = [f for f in tframes if f.specid == uca] 
    
    
            commandfile_fn = op.join(redux_dir, "commands_{:s}".format(uca))
            print("Creating {:s}".format(commandfile_fn))
            commandfile.append(open(commandfile_fn, 'w'))
            cmd = []
            for amp in SPEC:
                # Make Error Frame
                cmd = CC.mkerrorframe(vframesselect, amp, cmd, run=args.run_insitu) 
    
                # Subtract Overscan
                cmd = CC.subtractoverscan(biassec, args.suboverscan_options,
                                          vframesselect, amp, cmd, 
                                          run=args.run_insitu)  
                
                # Trim amplifier
                cmd = CC.extractfits(trimsec, vframesselect, amp, cmd,
                                     run=args.run_insitu)
                                          
                # Make Master Bias
                if args.make_masterbias:
                    cmd = CC.meanbiasfits(amp, uca, redux_dir, 'masterbias', 
                                          args.masterbias_options, zframesselect, 
                                          cmd, run=args.run_insitu) 
                                
                # Subtract Master Bias
                if args.use_other_masterbias:
                    masterbiasname = op.join(args.use_other_masterbias, 
                                            'masterbias' 
                                            + '_' 
                                            + uca 
                                            + '_' 
                                            + amp 
                                            + '.fits') 
                else:
                    masterbiasname = op.join(redux_dir, 
                                            'masterbias' 
                                            + '_' 
                                            + uca 
                                            + '_' 
                                            + amp 
                                            + '.fits')
                if args.make_mastertrace or args.make_masterarc or args.make_masterdark:                            
                    cmd = CC.subtractbias(oframesselect, masterbiasname, amp, cmd, 
                                      run=args.run_insitu)
        
            if args.make_mastertrace or args.make_masterarc or args.make_masterdark:                            
                # CCDCombine (multiplies by gain and joins amplifiers)
                cmd = CC.ccdcombine(oframesselect, cmd, run=args.run_insitu)
                '''Change from amps to side convention'''
                for o in oframesselect:
                    o.actionbase["L"] = o.actionbase["LL"]
                    o.actionbase["R"] = o.actionbase["RL"]
                # Clean the "cmp","flt","drk" files that aren't joined
                if args.run_insitu:
                    cmd = clean_folder(oframesselect, cmd, amp="LL")
                    
                # Add photon noise
                for side in SPECBIG:    
                    cmd = CC.addphotonnoise(oframesselect, side, cmd, 
                                            run=args.run_insitu)
                
                if args.drk_file_loc is not None:
                    # Create Master Dark Frames
                    for side in SPECBIG:
                        cmd = CC.meandarkfits(side, uca, redux_dir, 'masterdark',
                                              args.masterdark_options, dframesselect,
                                              cmd, run=args.run_insitu)
                    # Clean "drk" folder
                    if args.run_insitu:        
                        cmd = clean_folder(dframesselect, cmd, side=True)
                                              
                # Normalizing Calibration Frames by Exposure Time
                for cal in cframesselect:
                    for side in SPECBIG:
                        opt     = '-c {:0.1f}'.format(cal.exptime)
                        filename = build_name(cal, side=side)
                        cmd = CC.dividefits(filename, opt, cmd, run=args.run_insitu)
                        rmcmd = CC.build_rmcmd([cal], 'd', side=side)
                        cmd.append(rmcmd)
                        filename = build_name(cal, side=side)
                        cmd = CC.headfits(filename, args.headfits_options, cmd, 
                                          run=args.run_insitu) 
                                               
                # Making Master Arc Frames
                if args.make_masterarc:
                    for side in SPECBIG:
                        for lamp in LAMP_DICT.values():
                            lframesselectl = [l for l in lframes if l.specid == uca 
                                                                and l.object == lamp]
                            if len(lframesselect)>1:
                                cmd = CC.meanlampfits(side, uca, lamp, redux_dir, 
                                                      'masterarc', args.masterarc_options, 
                                                      lframesselectl, cmd, 
                                                      run=args.run_insitu) 
                            else:
                                filename = [build_name(f,side=side) 
                                            for f in lframesselectl]
                                mastername = build_mastername('masterarc', side, uca, 
                                                              redux_dir, lamp)
                                efilename = [build_name(f,side=side,error=True) 
                                             for f in lframesselectl]
                                emastername = build_mastername('e.masterarc', side, uca, 
                                                              redux_dir, lamp)
                                if args.run_insitu:
                                    shutil.copy(filename[0], mastername)
                                    shutil.copy(efilename[0], emastername)
                                cmd.append("cp {:s} {:s}".format(filename[0], mastername))
                                cmd.append("cp {:s} {:s}".format(efilename[0], 
                                                                 emastername))
                        if len(LAMP_DICT.values()) > 1:
                            opt = "--file {:s}".format(build_mastername('masterarc', 
                                                                        side, uca, 
                                                                        redux_dir, 
                                                                        LAMP_DICT[0]))
                            filename = build_mastername('masterarc', side, uca, 
                                                        redux_dir, LAMP_DICT[1])
                            # Adding lamps together (only works for two right now)
                            cmd = CC.addfits(filename, opt, cmd, run=args.run_insitu)
                            smastername = build_mastername('smasterarc', side, uca, 
                                                           redux_dir, LAMP_DICT[1])
                            mastername = build_mastername('masterarc', side, uca, 
                                                          redux_dir)
                            esmastername = build_mastername('e.smasterarc', side, uca, 
                                                            redux_dir, LAMP_DICT[1])
                            emastername = build_mastername('e.masterarc', side, uca, 
                                                           redux_dir)
                            if args.run_insitu:
                                shutil.copy(smastername, mastername)
                                shutil.copy(esmastername, emastername)
                                os.remove(smastername)
                                os.remove(esmastername)
                            cmd.append("cp {:s} {:s}".format(smastername, mastername))
                            cmd.append("cp {:s} {:s}".format(esmastername, emastername))
                            cmd.append("rm {:s} {:s}".format(smastername, esmastername))                                        
                        else:
                            mastername = build_mastername('masterarc', side, uca, 
                                                          redux_dir)
                            emastername = build_mastername('e.masterarc', side, uca, 
                                                           redux_dir)
                            smastername = build_mastername('smasterarc', side, uca, 
                                                           redux_dir, LAMP_DICT[0])
                            esmastername = build_mastername('e.smasterarc', side, uca, 
                                                            redux_dir, LAMP_DICT[0])
                            if args.run_insitu:
                                shutil.copy(smastername, mastername)
                                shutil.copy(esmastername, emastername)
                            cmd.append("cp {:s} {:s}".format(smastername, mastername))
                            cmd.append("cp {:s} {:s}".format(esmastername, emastername))
                            cmd.append("rm {:s} {:s}".format(smastername, esmastername)) 
                    # Clean "cmp" folder
                    if args.run_insitu:
                        cmd = clean_folder(lframesselect, cmd, side=True)
                
                
                # Making Master Trace Frames
                if args.make_mastertrace:
                    for side in SPECBIG:
                        if len(fframesselect) > 1:
                            cmd = CC.meantracefits(side, uca, redux_dir, 'mastertrace', 
                                                   args.mastertrace_options, fframesselect, 
                                                   cmd, run=args.run_insitu)
                        else:
                            
                            filename = [build_name(f, side) for f in fframesselect]
                            mastername = build_mastername('mastertrace', side, uca, 
                                                          redux_dir)
                            efilename = [build_name(f, side, error=True) 
                                         for f in fframesselect]
                            emastername = build_mastername('e.mastertrace', side, uca, 
                                                          redux_dir)
                            if args.run_insitu:
                                shutil.copy(filename[0], mastername)
                                shutil.copy(efilename[0], emastername)
                    # Clean "flt" folder
                    if args.run_insitu:
                        cmd = clean_folder(fframesselect, cmd, side=True)
            
                    # Making Master Twighlight Trace Frames
                    if args.twi_file_loc is not None:
                        for tw in tframesselect:
                            for side in SPECBIG:
                                filename = build_name(tw, side=side)
                                p = pyfits.open(filename)
                                N, D = p[0].data.shape
                                norm = np.median(p[0].data[:,(D/2-100):(D/2+100)])
                                opt     = '-c {:0.1f}'.format(norm)
                                cmd = CC.dividefits(filename, opt, cmd, 
                                                    run=args.run_insitu)
                                rmcmd = CC.build_rmcmd([tw], 'd', side=side)
                                cmd.append(rmcmd)
                                filename = build_name(tw, side=side)
                                cmd = CC.headfits(filename, args.headfits_options, cmd, 
                                                  run=args.run_insitu) 
                        for side in SPECBIG:
                            if len(tframesselect) > 1:
                                cmd = CC.meantracefits(side, uca, redux_dir, 
                                                       'mastertrace_twi', 
                                                       args.mastertrace_options, 
                                                       tframesselect, 
                                                       cmd, run=args.run_insitu)
                            else:
                                
                                filename = [build_name(f, side) for f in tframesselect]
                                mastername = build_mastername('mastertrace_twi', side, uca, 
                                                              redux_dir)
                                efilename = [build_name(f, side, error=True) 
                                             for f in tframesselect]
                                emastername = build_mastername('e.mastertrace_twi', side, 
                                                               uca, redux_dir)
                                if args.run_insitu:
                                    shutil.copy(filename[0], mastername)
                                    shutil.copy(efilename[0], emastername)
                        # Clean "twi" folder
                        if args.run_insitu:
                            cmd = clean_folder(tframesselect, cmd, side=True)
                            
            cmd = flatten(cmd,[])                             
            commandinfo.writecommand(commandfile[-1],cmd)
            commandfile[-1].close()
            
    # Run Deformer
    if args.run_deformer:
        for uca in ucam:
            cmd = [] #DUMMY FOR WORKING PURPOSES
            for side in SPECBIG:
                lines_dir = op.join(configdir, 'lines_files')
                lines_file = 'lines' + '_' + side + '_' + uca +  '.par'
                shutil.copy(op.join(lines_dir, lines_file), 
                            op.join(redux_dir, lines_file))
                cmd.append('cp {:s} {:s}'.format(op.join(lines_dir, lines_file), 
                                             op.join(redux_dir, lines_file)))
                if args.twi_file_loc is not None:
                    if args.twi_deformer:
                        mastertrace = build_mastername('mastertrace_twi', side, 
                                                       uca, redux_dir)
                        if args.use_dist_file:
                            dist_file = args.use_dist_file
                        else:
                            masterarc = build_mastername('masterarc', side, uca, redux_dir)          
                            cmd = CC.deformer(mastertrace, redux_dir, masterarc,  
                                              args.deformer_options, cmd,
                                              op.join(redux_dir, lines_file),
                                              run=args.run_insitu)
                            dist_file = mastertrace[:-5] + '.dist'
                            os.system(('python distplot.py -d %s ' 
                                    '-o distplot_cam%s.pdf %s %s' 
                                    %(dist_file, uca, masterarc, mastertrace)))
                        deformer_opts = ("-d %s --template-spectrum %s" 
                                          %(dist_file,
                                            args.temp_spec))
                        masterarc = build_mastername('mastertrace_twi', side, 
                                                     uca, redux_dir)
                        cmd = CC.deformer(mastertrace, redux_dir, masterarc, 
                                          deformer_opts, cmd, 
                                          run=args.run_insitu)
                        os.system(('python distplot.py -d %s ' 
                                    '-o distplot_cam%s_adj.pdf %s %s' 
                                    %(mastertrace[:-5]+'_adj.dist', uca,
                                      masterarc, mastertrace)))
                else:
                    mastertrace = build_mastername('mastertrace', side, uca, redux_dir)
                    masterarc = build_mastername('masterarc', side, uca, redux_dir)
                    cmd = CC.deformer(mastertrace, redux_dir, masterarc,  
                                      args.deformer_options, cmd, 
                                      op.join(redux_dir, lines_file), 
                                      run=args.run_insitu) 
                    os.system(('python distplot.py -d %s ' 
                                '-o distplot_cam%s.pdf %s %s' 
                                %(dist_file, uca, masterarc, mastertrace)))
                 

        
if __name__ == '__main__':
    main()  
