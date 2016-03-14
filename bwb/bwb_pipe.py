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
import sys
import os, shutil
import subprocess
import uuid

#import bwb_pipe_utils as pipe_utils
from glue import pipeline

from optparse import OptionParser
import ConfigParser

import bwb_pipe_utils as pipe_utils


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
    parser.add_option("--skip-datafind", default=False, action="store_true")

    (opts,args) = parser.parse_args()

    if opts.workdir is None:
        print >> sys.stderr, "ERROR: must specify --output-dir"
        sys.exit()


    if len(args)==0:
        print >> sys.stderr, "ERROR: require config file"
        sys.exit()


    # --- Read config file
    cp = ConfigParser.ConfigParser()
    cp.read(args[0])


    return opts, args, cp



# ----------------
# Set Parameters
# ----------------

# Read inputs
opts, args, cp = parser()

workdir = opts.workdir 

# Get trigger time(s)

gps = opts.trigger_time
if gps is None:
    print >> sys.stderr, "ERROR: must provide --trigger-time for the time being"
    sys.exit()

#
# --- Params from config file
#

ifoList = cp.get('datafind', 'ifoList')
channelList = cp.get('datafind', 'channelList')
frtypeList = cp.get('datafind', 'frtypeList')

# parse channels etc
ifoList=ifoList.split(',')
channelList=channelList.split(',')
frtypeList=frtypeList.split(',')



#############################################

# ----------------------------------------
# Setup analysis directory for deployment
# ----------------------------------------

datafind = os.path.join(workdir, 'datafind')
if not os.path.exists(datafind): os.makedirs(datafind)

shutil.copy(cp.get('paths','bw_executable'), workdir)



#############################################

# ------------------
# Call LIGO Data find
# -------------------
cacheFiles = {}
if not opts.skip_datafind:

    # We'll do this here so that we can, if desired, ship frames with jobs


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

        cacheFiles[ifo]=os.path.join('datafind', '{0}.cache'.format(ifo))

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

            if cache_entries.ndim==1: cache_entries=[cache_entries]

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


        for frame in frames_to_copy:
            print >> sys.stdout, "Copying %s"%frame
            shutil.copy(frame, os.path.join(workdir, 'datafind'))


        #
        # Now we need to make a new, local cache file
        # - do this by manipulating the path string in the cache file to be relative 
        for ifo in ifoList:
            cache_file = os.path.join(workdir, 'datafind/{ifo}.cache'.format(ifo=ifo))
            shutil.copy(cache_file, cache_file.replace('cache','cache.bk'))

            cache_entries = np.loadtxt(cache_file, dtype=str)
            if cache_entries.ndim>1:
                cache_entries = cache_entries[unique_idx]

            if cache_entries.ndim==1:
                cache_entries = [cache_entries]

            new_cache = open(cache_file, 'w')
            for entry in cache_entries:

                local_path=os.path.join('datafind',entry[4].split('/')[-1])

                new_cache.writelines('{ifo} {type} {gps} {length} {path}\n'.format(
                    ifo=ifo, type=entry[1], gps=entry[2], length=entry[3],
                    path=local_path))

            new_cache.close()

else:
    # user-specified cache files
    for i,ifo in enumerate(ifoList):
        cacheFiles[ifo] = \
                cp.get('datafind','cacheFiles').split(',')[i]


#############################################


# ----------------------------------
# Directory Structure
# ----------------------------------

logdir = os.path.join(workdir, 'logs')

if not os.path.exists(logdir): os.makedirs(logdir)

# -----------------------------------------------------------------------
# DAG Writing
# -----------------------------------------------------------------------

#
# Initialise DAG and Jobs
#

# ---- Create a dag to which we can add jobs.
dag = pipeline.CondorDAG(log=opts.user_tag+'.log')

# ---- Set the name of the file that will contain the DAG.
dag.set_dag_file( os.path.join(workdir,'BayesWave_{0}'.format(opts.user_tag)) )

# ---- Make instance of BayesWaveBurstJob.
bwb_job = pipe_utils.BayesWaveBurstJob(cp, workdir, cacheFiles)

#
# Build Nodes
#

# XXX: ultimately want 1 output dir per job!

outputDir  = 'BayesWaveBurst_' + str(int(gps)) + '_' + str(uuid.uuid1())

if not os.path.exists(os.path.join(workdir,outputDir)): os.makedirs(os.path.join(workdir,outputDir))

bwb_node = pipe_utils.BayesWaveBurstNode(bwb_job)

# add options
bwb_node.set_trigtime(opts.trigger_time)
bwb_node.set_PSDstart(opts.trigger_time)
bwb_node.set_outputDir(outputDir)


# Add Nodes to DAG
dag.add_node(bwb_node)

#
# Finalise DAG
#
# ---- Write out the submit files needed by condor.
dag.write_sub_files()

# ---- Write out the DAG itself.
dag.write_dag()
dag.write_script()

sys.exit()
# -----------------
# BWB cmdline 
# -----------------
bwbcmdline = """--ifo H1 --H1-flow $(macroflow) --H1-channel $(macroh1channel)   \
--ifo L1 --L1-flow $(macroflow) --L1-channel $(macrol1channel)  \
--H1-cache ./datafind/H1.cache \
--L1-cache ./datafind/L1.cache \
--trigtime $(macrotrigtime) --srate $(macrosrate) --seglen $(macroseglen) \
--bayesLine --PSDstart $(macropsdstart) --PSDlength $(macropsdlen) \
--outputDir $(macrooutputDir)  \
--L1-timeslide $(macrol1timeslide) 
"""

# ----------
# Template for bwb submit file
# ----------
submit_str = """
executable=BayesWaveBurst
universe=standard
arguments={bwbcmdline}
output={outputDir}/{outputDir}.out
error={outputDir}/{outputDir}.err
log={outputDir}/{outputDir}.log
notification=never
should_transfer_files=YES
when_to_transfer_output = ON_EXIT
stream_error=False
stream_output=False
WantRemoteIO=False
accounting_group=ligo.prod.o1.burst.paramest.bayeswave
transfer_input_files=BayesWaveBurst,datafind,{outputDir},logs
transfer_output_files={outputDir},logs
queue 1
"""

dagfile = open(os.path.join(workdir, 'dagfile.dag'), 'w')

   

# -- Create BWB submit file

submitname = os.path.join(workdir, 'submitBWB.sub')
submitfile = open(submitname, 'w')
submitfile.write(submit_str.format(bwbcmdline=bwbcmdline, outputDir=outputDir))
submitfile.close()


# ---- write the dag file

dagfile.write("JOB {jobname} submitBWB.sub\n".format(jobname=outputDir))

# -----------------
# BWB arguments 
# -----------------
bwbargsfmt = """macroflow=\"{flow}\" macroh1channel=\"{h1_channel}\" \
macrol1channel=\"{l1_channel}\" macrotrigtime=\"{gps}\" macrosrate=\"{srate}\" \
macroseglen=\"{seglen}\" macropsdstart=\"{gps}\" macropsdlen=\"{psdlen}\" \
macrooutputDir=\"{outputDir}\" macrol1timeslide=\"{lag}\" 
"""

bwbvars = bwbargsfmt.format(flow=flow, h1_channel=h1_channel,
        l1_channel=l1_channel, gps=gps, srate=srate, seglen=seglen,
        psdlen=psdlen, outputDir=outputDir, lag=lag )

dagfile.write("VARS {jobname} {bwbvars}".format(jobname=outputDir,
    bwbvars=bwbvars))
dagfile.write("RETRY {jobname} 1 \n\n".format(jobname=outputDir))

dagfile.close()

# -----------------
# BWB shell script 
# -----------------
fullcmdline = """./BayesWaveBurst \
--ifo H1 --H1-flow {flow} --H1-channel {h1_channel}   \
--ifo L1 --L1-flow {flow} --L1-channel {l1_channel}  \
--H1-cache ./datafind/H1.cache \
--L1-cache ./datafind/L1.cache \
--trigtime {gps} --srate {srate} --seglen {seglen} \
--bayesLine --PSDstart {gps} --PSDlength {psdlen} \
--outputDir {outputDir}  \
--L1-timeslide {lag}
""".format(flow=flow, h1_channel=h1_channel, l1_channel=l1_channel, gps=gps,
        srate=srate, seglen=seglen,
        psdlen=psdlen, outputDir=outputDir, lag=lag )

shellname = os.path.join(workdir, 'runBWB.sh')
shellfile = open(shellname, 'w')
shellfile.write(fullcmdline)
shellfile.close()
os.chmod(shellname,0755)





 
