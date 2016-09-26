# -*- coding: utf-8 -*-
"""
Created on Tue May  3 10:15:15 2016

Cure commands that can be run through python scripts.
The format of the commands requires the class VirusFrame.

@author: gregz
"""

import os
from os import environ
import os.path as op
import sys


''' Set your curebin here.  If you have an environment variable called 
'CUREBIN' already set then you can put 'CUREBIN = None' and the program will 
use your environment path'''
CUREBIN = None

if not CUREBIN:
    CUREBIN = environ.get('CUREBIN')
    
if not CUREBIN:
    print("Please set CUREBIN as  environment variable or in the script")
    sys.exit(1)
    
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
    if amp:
        filename = op.join(frame.origloc, frame.actionbase[amp] 
                                          + frame.basename 
                                          + '_' 
                                          + frame.ifuslot 
                                          + amp 
                                          + '_' 
                                          + frame.type 
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
    
def build_rmcmd(frames, action, side=None, amp=None):
    if frames:
        if frames[0].clean:
            if amp:
                x = ['rm ']
                x.append([f.addbase(action, amp=amp) for f in frames])
                rmcmd = ''.join(flatten(x, []))
            if side:
                x = ['rm ']
                x.append([f.addbase(action, side=side) for f in frames]) 
                rmcmd = ''.join(flatten(x, []))
    return rmcmd


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
        return os.system(op.join(CUREBIN, command) 
                         + ' 1>>output.log  2>> error.log')
    else:
        return os.system(op.join(CUREBIN, command))
        
        
def mkerrorframe(frames, amp, cmd, run=True):
    
    filenames = [build_name(f, amp=amp) for f in frames]

    command = 'mkerrorframe' 
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
        
    if run:        
        run_cure_command(command, 0)
    
    cmd.append(command)
    
    return cmd
    
        
def subtractoverscan(biassec, opts, frames, amp, cmd, run=True):
    
    filenames = [build_name(f, amp=amp) for f in frames]

    command = 'subtractfits %s -o %s' %(opts, biassec)

    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
    
    if run:    
        run_cure_command(command, 0)
       
    cmd.append(command)
    
    rmcmd = build_rmcmd(frames, 's', amp=amp)
    
    cmd.append(rmcmd)
    
    return cmd    
    
def subtractbias(frames, masterbiasname, amp, cmd, run=True):
        
    filenames = [build_name(f, amp=amp) for f in frames]

    command = 'subtractfits -f %s' % (masterbiasname) 
        
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
    
    if run:
        run_cure_command(command, 0)

    cmd.append(command)
        
    rmcmd = build_rmcmd(frames, 's', amp=amp)
    
    cmd.append(rmcmd)
    
    return cmd    
    
def subtractdark(frames, masterdarkname, side, cmd, run=True):
        
    filenames = [build_name(f, side=side) for f in frames]

    command = 'subtractfits -f %s' % (masterdarkname) 
        
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
    
    if run:
        run_cure_command(command, 0)

    cmd.append(command)
        
    rmcmd = build_rmcmd(frames, 's', side=side)
    
    cmd.append(rmcmd)
    
    return cmd  

def meanbiasfits(amp, specid, dest_dir, basename, opts, frames, cmd, run=True):
    
    filenames = [build_name(f, amp=amp) for f in frames]
                 
    mastername = build_mastername(basename, amp, specid, dest_dir)
    
    command = 'meanfits %s -o %s' % (opts, mastername)  
    
    for i in xrange(len(filenames)):       
        command = command + ' ' + filenames[i]
    
    if run:    
        run_cure_command(command, 0)
    
    cmd.append(command)

    rmcmd = build_rmcmd(frames, '', amp=amp)
    
    cmd.append(rmcmd)
    
    return cmd    
    
def meandarkfits(side, specid, dest_dir, basename, opts, frames, cmd, 
                 run=True):
    
    filenames = [build_name(f, side=side) for f in frames]

    mastername = build_mastername(basename, side, specid, dest_dir)
    
    command = 'meanfits %s -o %s' % (opts, mastername)  
    
    for i in xrange(len(filenames)):   
        command = command + ' ' + filenames[i]
     
    if run: 
        run_cure_command(command, 0)
        
    cmd.append(command)    
    
    rmcmd = build_rmcmd(frames, '', side=side)
    
    cmd.append(rmcmd)
    
    return cmd
    

def extractfits(trimsec, frames, amp, cmd, run=True):

    filenames = [build_name(f, amp=amp) for f in frames]
    
    command = 'extractfits -r %s' % (trimsec)  
    
    for i in xrange(len(filenames)):
        command = command + ' ' + filenames[i]
    
    if run:    
        run_cure_command (command, 0)
        
    cmd.append(command)
    
    rmcmd = build_rmcmd(frames, 'e', amp=amp)
    
    cmd.append(rmcmd)
    
    return cmd
    
    
def ccdcombine(frames, cmd, run=True):
    
    command = []    
    for i in xrange(len(frames)):
        
        filename = op.join(frames[i].origloc, frames[i].actionbase["LL"] 
                                              + frames[i].basename 
                                              + '_' 
                                              + frames[i].ifuslot)
        
        command.append('ccdcombine ' + filename + ' ' + frames[i].type)   
        
        if run:
            run_cure_command(command[-1], 0)
    
    cmd.append(command)

    return cmd
    

def addphotonnoise(frames, side, cmd, run=True):

    command = 'addphotonnoise'
    
    filenames = [build_name(f, side=side) for f in frames]

    for i in xrange(len(filenames)): 
        command = command + ' ' + filenames[i]
    
    if run:    
        run_cure_command(command, 0)
    
    cmd.append(command)

    rmcmd = build_rmcmd(frames, 'p', side=side)
    
    cmd.append(rmcmd)
    
    return cmd
 
    
def dividepixelflat(frames, opts, side, cmd):

    filenames = [build_name(f, side=side) for f in frames]

    for i in filenames:
        
        command = 'dividefits %s %s' % (opts, i)

        run_cure_command(command, 0)
    
    cmd.append(command)
    
    rmcmd = build_rmcmd(frames, 'd', side=side)
    
    cmd.append(rmcmd)
    
    return cmd 
    
    
def dividefits(filename, opts, cmd, run=True):

    command = 'dividefits %s %s' % (opts,filename)
    
    if run:
        run_cure_command(command, 0)
    
    cmd.append(command)
    
    return cmd
    
    
def addfits(filename, opts, cmd, run=True):

    command = 'addfits %s %s' % (opts,filename)
    
    if run:
        run_cure_command(command, 0)

    cmd.append(command)
        
    return cmd
    
    
def headfits(filename, opts, cmd, run=True):

    command = 'headfits %s %s' % (opts,filename)
    
    if run:    
        run_cure_command(command, 0)

    cmd.append(command)
    
    return cmd
    
    
def meanlampfits(side, specid, lamp, dest_dir, basename , opts, frames, cmd, 
                 run=True):
    
    filenames = [build_name(f, side=side) for f in frames]

    mastername = build_mastername(basename, side, specid, dest_dir, lamp)
    
    command = 'meanfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):  
        command = command + ' ' + filenames[i]
    
    if run:    
        run_cure_command(command, 0)

    cmd.append(command)
    
    return cmd
    
    
def meantracefits(side, specid, dest_dir, basename, opts, frames, cmd, 
                  run=True):
    
    filenames = [build_name(f, side=side) for f in frames]

    mastername = build_mastername(basename, side, specid, dest_dir)
    
    command = 'meanfits %s -o %s' % (opts,mastername)  
    
    for i in xrange(len(filenames)):  
        command = command + ' ' + filenames[i]
    
    if run:
        run_cure_command(command, 0)

    cmd.append(command)
    
    return cmd
    

def deformer(mastertrace, cal_dir, opts, cmd, masterarc=None, linesfile=None, 
             dist=None, fmod=None, rerun=False, run=True, side=None):
  
    if rerun:
        filenames = [build_name(f, side=side) for f in mastertrace]
        for f in filenames:
            command = 'deformer %s -o \"%s\" -d %s -f %s %s' %(opts, cal_dir,
                                                           dist, fmod,
                                                           f)
            if run:
                run_cure_command( command, 0 )
            
            cmd.append(command)
            
    else:
        command = 'deformer %s -o \"%s\" -l %s -a %s %s' %(opts, 
                                                       cal_dir, 
                                                       linesfile,
                                                       masterarc, 
                                                       mastertrace)  

        if run:
            run_cure_command( command, 0 )

        cmd.append(command)

    return cmd
    
    
def subtractsky(frames, side, distmodel, fibermodel, opts, cmd, run=True, 
                skymaster=""):
    
    filenames = [build_name(f, side=side) for f in frames]

    command = []    

    for f in filenames:
        command.append('subtractsky %s %s -d %s -f %s -D %s -F %s %s' %(opts,
                                                                   skymaster,
                                                                   distmodel,
                                                                   fibermodel,
                                                                   distmodel,
                                                                   fibermodel,
                                                                   f))  
        if run:                                                               
            run_cure_command(command[-1], 0)
        
    cmd.append(command)
    
    for f in frames:
        f.actionbase[side] = 'S' + f.actionbase[side]
           
    return cmd
    
    
def fiberextract(frames, side, distmodel, fibermodel, opts, cmd, run=True, 
                 new_deformer=False, cal_folder=None):
    
    filenames = [build_name(f, side=side) for f in frames]
    
    command = []        
    
    for i in xrange(len(filenames)):
        if new_deformer:
            base = op.basename(filenames[i])[:-5]
            distmodel = op.join(cal_folder, base + '_adj.dist')
            fibermodel = op.join(cal_folder, base + '.fmod')
            command.append('fiberextract %s -d %s -f %s %s' %(opts, 
                                                     distmodel, 
                                                     fibermodel, 
                                                     filenames[i]))            
        else:
            command.append('fiberextract %s -d %s -f %s %s' %(opts, 
                                                     distmodel, 
                                                     fibermodel, 
                                                     filenames[i]))
        if run:                                                               
            run_cure_command(command[-1], 0)

    cmd.append(command)
              
    return cmd    
    
def detect(ifufile, ditherfile, outfile, opts, cmd, run=True):
    
    command = 'detect %s -i %s -o %s %s' %(opts,
                                           ifufile,
                                           outfile,
                                           ditherfile)
        
    if run:    
        run_cure_command(command, 0)

    cmd.append(command)
    
    return cmd   
    
def mkcube(ifufile, ditherfile, opts, cmd, run=True):
    
    command = 'mkcube %s -i %s %s' %(opts,
                                     ifufile,
                                     ditherfile)
        
    if run:    
        run_cure_command(command, 0)

    cmd.append(command)
    
    return cmd