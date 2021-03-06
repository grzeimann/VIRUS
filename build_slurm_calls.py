# -*- coding: utf-8 -*-
"""
@author: gregz
"""
from __future__ import print_function
import argparse as ap
import os
import os.path as op
import textwrap
import numpy as np

# Default set of spectrographs for reduction
SPECID = ["004","008","012","013","016","017","020","024","025","027","032",
          "037","038","041","047","051"]

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
    description = textwrap.dedent('''Build Slurm Files.''')
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.RawTextHelpFormatter)
                        
    parser.add_argument("--specid", nargs='?', type=str, 
                        help='''List of SPECID's for processing. 
                        Ex: "020,008".''', default = None)

    parser.add_argument("--num_per_slurm", nargs='?', type=int, 
                        help='''Number of calls per slurm''', default=16)

    parser.add_argument("--argstr", nargs='?', type=str, 
                        help='''Argument String.
                        "-r -rc --cal_dir /work/00115/gebhardt/maverick/cal/20160905 --rootdir /work/03946/hetdex/maverick/ -sd 20160903 -so 8 -t -f -s -m -ud"''', 
                        default = None)

    args = parser.parse_args(args=argv)

    if args.specid:
        args.specid = args.specid.replace(" ", "").split(',')
    else:
        args.specid = SPECID  
              
    return args 
    
def build_slurm_file(num, fn):
    with open(fn+'.slurm','w') as f:
        s = ('''#!/bin/bash
#
# Simple SLURM script for submitting multiple serial
# jobs (e.g. parametric studies) using a script wrapper
# to launch the jobs.
#
# To use, build the launcher executable and your
# serial application(s) and place them in your WORKDIR
# directory.  Then, edit the CONTROL_FILE to specify 
# each executable per process.
#-------------------------------------------------------
#-------------------------------------------------------
# 
#------------------Scheduler Options--------------------
#SBATCH -J HETDEX              # Job name
#SBATCH -n {:d}                  # Total number of tasks
#SBATCH -p gpu                 # Queue name
#SBATCH -o HETDEX.o%j          # Name of stdout output file (%j expands to jobid)
#SBATCH -t 01:20:00            # Run time (hh:mm:ss)
#SBATCH -A Hobby-Eberly-Telesco
#------------------------------------------------------
#
# Usage:
#	#$ -pe <parallel environment> <number of slots> 
#	#$ -l h_rt=hours:minutes:seconds to specify run time limit
# 	#$ -N <job name>
# 	#$ -q <queue name>
# 	#$ -o <job output file>
#	   NOTE: The env variable $JOB_ID contains the job id. 
#
#------------------------------------------------------

#------------------General Options---------------------
module load launcher
export EXECUTABLE=$TACC_LAUNCHER_DIR/init_launcher
export WORKDIR=.
export CONTROL_FILE={:s}

# Variable descriptions:
#
#  TACC_LAUNCHER_PPN = number of simultaneous processes per host
#                      - if this variable is not set, value is
#                        determined by the process density/wayness
#                        specified in 'Scheduler Options'
#  EXECUTABLE        = full path to the job launcher executable
#  WORKDIR           = location of working directory
#  CONTROL_FILE      = text input file which specifies
#                      executable for each process
#                      (should be located in WORKDIR)
#------------------------------------------------------

#--------- Intel Xeon Phi Options (EXPERIMENTAL) -------------
export TACC_LAUNCHER_NPHI=0
export TACC_LAUNCHER_PHI_PPN=8
export PHI_WORKDIR=.
export PHI_CONTROL_FILE=phiparamlist

# Variable descriptions:
#  TACC_LAUNCHER_NPHI    = number of Intel Xeon Phi cards to use per node
#                          (use 0 to disable use of Xeon Phi cards)
#  TACC_LAUNCHER_PHI_PPN = number of simultaneous processes per Xeon Phi card
#  PHI_WORKDIR           = location of working directory for Intel Xeon Phi jobs
#  PHI_CONTROL_FILE      = text input file which specifies executable
#                          for each process to be run on Intel Xeon Phi
#                          (should be located in PHI_WORKDIR)
#------------------------------------------------------

#------------ Task Scheduling Options -----------------
export TACC_LAUNCHER_SCHED=interleaved

# Variable descriptions:
#  TACC_LAUNCHER_SCHED = scheduling method for lines in CONTROL_FILE
#                        options (k=process, n=num. lines, p=num. procs):
#                          - interleaved (default): 
#                              process k executes every k+nth line
#                          - block:
#                              process k executes lines [ k(n/p)+1 , (k+1)(n/p) ]
#                          - dynamic:
#                              process k executes first available unclaimed line
#--------------------------------------------------------

#----------------
# Error Checking
#----------------

if [ ! -d $WORKDIR ]; then
        echo " "
	echo "Error: unable to change to working directory."
	echo "       $WORKDIR"
	echo " "
	echo "Job not submitted."
	exit
fi

if [ ! -x $EXECUTABLE ]; then
	echo " "
	echo "Error: unable to find launcher executable $EXECUTABLE."
	echo " "
	echo "Job not submitted."
	exit
fi

if [ ! -e $WORKDIR/$CONTROL_FILE ]; then
	echo " "
	echo "Error: unable to find input control file $CONTROL_FILE."
	echo " "
	echo "Job not submitted."
	exit
fi

#----------------
# Job Submission
#----------------

cd $WORKDIR/
echo " WORKING DIR:   $WORKDIR/"

$TACC_LAUNCHER_DIR/paramrun SLURM $EXECUTABLE $WORKDIR $CONTROL_FILE $PHI_WORKDIR $PHI_CONTROL_FILE

echo " "
echo " Parameteric Job Complete"
echo " "
'''.format(num, fn))
        f.write(s)
        f.flush()

def build_sci_files(fn, calls):
    with open(fn, 'w') as f:
        f.write('\n'.join(calls) + "\n")
        f.flush()
   
def main():
    args = parse_args()
    init_str = 'python science_script.py %s' %args.argstr
    dv = np.linspace(0.7,1.3,7)
    b1v = np.linspace(-0.3,0.3,7)
    calls = []
    for i in dv:
        for j in b1v:
            dirn = 'd_%04.2f_b1_%04.2f' %(i, j)
            if not op.exists(dirn):
                os.mkdir(dirn)
            for spec in args.specid:
                ast = ' --output %s/c%s' %(dirn, spec)
                bst = ' --specid %s' %spec
                cst = ' --dark_mult_val %0.2f' %i
                dst = ' --bias1_mult_val %0.2f' %j
                nstr = init_str + ast + bst + cst + dst
                calls.append(nstr)
    print(calls)
    i = 0
    while (i*args.num_per_slurm) < len(calls):
        fn = 'rsci_d'+str(i)
        ca = calls[(i*args.num_per_slurm):((i+1)*args.num_per_slurm)]
        num = int(len(ca))
        build_sci_files(fn, ca)
        build_slurm_file(num, fn)
        i+=1
    
if __name__ == '__main__':
    main()  