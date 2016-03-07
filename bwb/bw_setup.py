#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2015-2016 James Clark <james.clark@ligo.org>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import numpy as np
import time
import sys, json
import os, shutil
import subprocess

from optparse import OptionParser
import ConfigParser


def parser():
    """
    Parser for input (command line and ini file)
    """

    # --- cmd line
    parser = OptionParser()
    parser.add_option("-t", "--user-tag", default="TEST", type=str)
    parser.add_option("-o", "--workdir", type=str, default=None)
    parser.add_option("--trigger-time", type=float, default=None)
    parser.add_option("--server", type=str, default=None)
    parser.add_option("--copy-frames", default=False, action="store_true")

    (opts,args) = parser.parse_args()

    if opts.workdir is None:
        print >> sys.stderr, "ERROR: must specify --output-dir"
        sys.exit()


    if len(args)==0:
        print >> sys.stderr, "ERROR: require config file"
        sys.exit()


    # --- Read config file
    configparser = ConfigParser.ConfigParser()
    configparser.read(args[0])


    return opts, args, configparser



# ----------------
# Set Parameters
# ----------------

# Read inputs
opts, args, configparser = parser()

workdir = opts.workdir 

# Get trigger time(s)

gps = opts.trigger_time
if gps is None:
    print >> sys.stderr, "ERROR: must provide --trigger-time for the time being"
    sys.exit()

#
# --- Params from config file
#

lag = configparser.getfloat('analysis', 'lag')

ifoList = configparser.get('analysis', 'ifoList')
channelList = configparser.get('analysis', 'channelList')
frtypeList = configparser.get('analysis', 'frtypeList')

# parse channels etc
ifoList=ifoList.split(',')
channelList=channelList.split(',')
frtypeList=frtypeList.split(',')

h1_channel=np.array(channelList)[np.array(ifoList)=='H1'][0]
h1_type=np.array(frtypeList)[np.array(ifoList)=='H1'][0]
l1_channel=np.array(channelList)[np.array(ifoList)=='L1'][0]
l1_type=np.array(frtypeList)[np.array(ifoList)=='L1'][0]

# executable
bwb = configparser.get('paths', 'bw_executable')

# Time-Freq config
flow = configparser.getfloat('analysis', 'flow')
seglen = configparser.getfloat('analysis', 'seglen')
psdlen = configparser.getfloat('analysis', 'psdlen')
srate = configparser.getint('analysis', 'sample-rate')



#############################################

# ----------------------------------------
# Setup analysis directory for deployment
# ----------------------------------------

name = workdir

datafind = os.path.join(workdir, 'datafind')
if not os.path.exists(datafind): os.makedirs(datafind)

shutil.copy(bwb, workdir)



#############################################

# ------------------
# Call LIGO Data find
# -------------------


# Set min/max gps times for LIGO data find:
min_gps = gps - 25.0 # XXX: why 25???
max_gps = gps + 25.0

start = min_gps 
end   = max_gps

for ifo, frtype in zip(ifoList,frtypeList):
    cachefilefmt = os.path.join(datafind, '{0}.cache')
    if opts.server is not None:
        ldfcmd = "gw_data_find --observatory {o} --type {frtype} -s {start} -e\
{end} --lal-cache --server={server} | grep file > {cachefile}".format(
                o=ifo[0], frtype=frtype, cachefile = cachefilefmt.format(ifo),
                start=start, end=end, server=opts.server)
    else:
        ldfcmd = "gw_data_find --observatory {o} --type {frtype} -s {start} -e {end} --lal-cache | grep file > {cachefile}".format(
                o=ifo[0], frtype=frtype, cachefile = cachefilefmt.format(ifo),
                start=start, end=end)
    print "Calling LIGO data find ..."
    print ldfcmd
    try:
        subprocess.call(ldfcmd, shell=True)
    except:
        print >> sys.stderr, "ERROR: datafind failed, continuing for testing"


#############################################
# Get frame files from cache

if opts.copy_frames:
    
    print >> sys.stdout, "Copying frame files for input"

    # XXX Having found these files, we now want to copy them to the working
    # directory and make fresh, local cache files

    # 1) read cache file
    # 2) identify unique frames
    # 3) copy unique frames to datafind directory
    # 4) add to transfer files

    frames_to_copy=[]
    for ifo in ifoList:
        cache_entries = np.loadtxt('{workdir}/datafind/{ifo}.cache'.format(
            workdir=workdir, ifo=ifo), dtype=str)

        frame_paths=[]
        frame_files=[]
        for c,cache_entry in enumerate(cache_entries):
            frame_paths.append(cache_entry[-1].split('localhost')[-1])
            frame_files.append(frame_paths[c].split('/')[-1])

        #unique_frames, unique_idx = np.unique(frame_files, return_index=True)

        unique_frames = list(set(frame_files))
        unique_idx = []
        for f,frame_file in enumerate(frame_files):
            current_idx=unique_frames.index(frame_file)
            if f==0:
                unique_idx.append(current_idx)
            else:
                if unique_idx[f-1] !=current_idx:
                    unique_idx.append(current_idx)

        unique_idx=np.array(unique_idx)
        frame_paths=np.array(frame_paths)
        if len(unique_frames)>1:
            for frame_path in frame_paths[unique_idx]:
                frames_to_copy.append(frame_path)
        else:
            frames_to_copy.append(frame_paths[unique_idx][0])

    print frames_to_copy

    for frame in frames_to_copy:
        print >> sys.stdout, "Copying %s"%frame
        shutil.copy(frame, os.path.join(workdir, 'datafind'))


    # Now we need to make a new, local cache file
    shutil.copy('{ifo}.cache'.format(ifo=ifo), '{ifo}.cache.bk'.format(ifo=ifo))


sys.exit()

#############################################

# ----------------------------------
# Write Job submission files and DAG
# ----------------------------------

# -----------------
# BWB cmdline 
# -----------------
bwbcmdline = """--ifo H1 --H1-flow $(macroflow) --H1-channel $(macroh1channel)   \
--ifo L1 --L1-flow $(macroflow) --L1-channel $(macrol1channel)  \
--H1-cache ./datafind/H1.cache \
--L1-cache ./datafind/L1.cache \
--trigtime $(macrotrigtime) --srate $(macrosrate) --seglen $(macroseglen)} \
--bayesLine --PSDstart $(macropsdstart) --PSDlength $(macropsdlen) \
--outputDir $(macrooutdir)  \
--L1-timeslide $(macrol1timeslide) 
"""

# ----------
# Template for bwb submit file
# ----------
submit_str = """
executable=BayesWaveBurst
universe=standard
arguments={bwbcmdline}
output=condorOut.out
error=condorError.err
log=condorLog.log
notification=never
should_transfer_files=YES
when_to_transfer_output = ON_EXIT
stream_error=False
stream_output=False
WantRemoteIO=False
accounting_group=ligo.prod.o1.burst.paramest.bayeswave
transfer_input_files=BayesWaveBurst,datafind,{outdir}
transfer_output_files={outdir}
queue 1
"""

dagfile = open(os.path.join(workdir, 'dagfile.dag'), 'w')

outdir  = 'job_' + str(int(gps)) + '_' + str(lag)
trigdir = os.path.join(workdir, 'job_' + str(int(gps)) + '_' + str(lag))

if os.path.exists(trigdir): print str(int(gps)) + '_' + str(lag)+' exists'

if not os.path.exists(trigdir): os.makedirs(trigdir)
   

# -- Create BWB submit file

submitname = os.path.join(workdir, 'submitBWB.sub')
submitfile = open(submitname, 'w')
submitfile.write(submit_str.format(bwbcmdline=bwbcmdline, outdir=outdir))
submitfile.close()


# ---- write the dag file

dagfile.write("JOB {jobname} submitBWB.sub\n".format(jobname=outdir))

# -----------------
# BWB arguments 
# -----------------
bwbargsfmt = """macroflow=\"{flow}\" macroh1channel=\"{h1_channel}\" \
macrol1channel=\"{l1_channel}\" macrotrigtime=\"{gps}\" macrosrate=\"{srate}\" \
macroseglen=\"{seglen}\" macropsdstart=\"{gps}\" macropsdlen=\"{psdlen}\" \
macrooutdir=\"{outdir}\" macrol1timeslide=\"{lag}\" 
"""

bwbvars = bwbargsfmt.format(flow=flow, h1_channel=h1_channel,
        l1_channel=l1_channel, gps=gps, srate=srate, seglen=seglen,
        psdlen=psdlen, outdir=outdir, lag=lag )

dagfile.write("VARS {jobname} {bwbvars}".format(jobname=outdir,
    bwbvars=bwbvars))
dagfile.write("RETRY {jobname} 1 \n\n".format(jobname=outdir))

dagfile.close()
