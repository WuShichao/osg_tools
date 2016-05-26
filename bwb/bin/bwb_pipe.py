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
import fileinput

from glue import pipeline

from glue.ligolw import lsctables
from glue.ligolw import table
from glue.ligolw import utils
from pylal.SimInspiralUtils import ExtractSimInspiralTableLIGOLWContentHandler
lsctables.use_in(ExtractSimInspiralTableLIGOLWContentHandler)

from optparse import OptionParser
import ConfigParser

import bwb_pipe_utils as pipe_utils

def localize_xml(xmlfile, old_path, new_path):
    """
    Modify absolute paths in xml files to relative paths
    """

    f = open(xmlfile,'r')
    filedata = f.read()
    f.close()

    newdata = filedata.replace(old_path,new_path)

    shutil.move(xmlfile, xmlfile+'.bak')

    f = open(xmlfile,'w')
    f.write(newdata)
    f.close()

    return 0

def parser():
    """
    Parser for input (command line and ini file)
    """

    # --- cmd line
    parser = OptionParser()
    parser.add_option("-t", "--user-tag", default="TEST", type=str)
    parser.add_option("-o", "--workdir", type=str, default=None)
    parser.add_option("--trigger-time", type=float, default=None)
    parser.add_option("--trigger-list", type=str, default=None)
    parser.add_option("--server", type=str, default=None)
    parser.add_option("--copy-frames", default=False, action="store_true")
    parser.add_option("--skip-datafind", default=False, action="store_true")

# XXX: putting this in the config.ini for now
#   parser.add_option("--inj", default=None)
#   parser.add_option("--nrhdf5", default=None)
#   parser.add_option("--events", default="all")
 
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

def hyphen_range(s):
    """
    yield each integer from a complex range string like "1-9,12, 15-20,23"

    Stolen from:
    http://code.activestate.com/recipes/577279-generate-list-of-numbers-from-hyphenated-and-comma/
    """

    for x in s.split(','):
        elem = x.split('-')
        if len(elem) == 1: # a number
            yield int(elem[0])
        elif len(elem) == 2: # a range inclusive
            start, end = map(int, elem)
            for i in xrange(start, end+1):
                yield i
        else: # more than one hyphen
            raise ValueError('format error in %s' % x)

# TODO
# 1) set retry=1 or 2


# ----------------
# Set Parameters
# ----------------

# --- Parse options, arguments and ini file
opts, args, cp = parser()

workdir = opts.workdir 
if not os.path.exists(workdir): os.makedirs(workdir)

# --- Make local copies of necessary input files
shutil.copy(args[0], os.path.join(workdir, 'config.ini'))

# Injection file (e.g., sim-inspiral table)
try:
    injfile=cp.get('injections', 'injfile')
except:
    injfile=None

if injfile is not None:
    shutil.copy(injfile, workdir)
    injfile=os.path.basename(injfile)


# NR HDF5 data
try:
    nrdata=cp.get('injections', 'nrdata')
    nr_full_path=cp.get('injections', 'nrdata')
except:
    nrdata=None
if nrdata is not None:
    shutil.copy(nrdata, workdir)
    nrdata=os.path.basename(nrdata)

    # Make sure normal permissions on hdf5
    os.chmod(os.path.join(workdir, nrdata), 0644)

    # Modify xml IN WORKDIR to point to local hdf5
    localize_xml(os.path.join(workdir, injfile), nr_full_path, nrdata)

#
# --- Params from config file
#

ifoList=cp.get('datafind', 'ifoList')
channelList=cp.get('datafind', 'channelList')
frtypeList=cp.get('datafind', 'frtypeList')

# parse channels etc
ifoList=ifoList.split(',')
channelList=channelList.split(',')
frtypeList=frtypeList.split(',')


#############################################
#
# Get trigger time(s)
#
if opts.trigger_time is not None:
    trigtimes = [opts.trigger_time]

if opts.trigger_list is not None:
    trigtimes = np.loadtxt(opts.trigger_list)

if injfile is not None:

    #
    # Read inspinj file
    #
    xmldoc=utils.load_filename(os.path.join(workdir,injfile), contenthandler=
            ExtractSimInspiralTableLIGOLWContentHandler, verbose=True)
    table=table.get_table(xmldoc, lsctables.SimInspiralTable.tableName)

    # Get gps list from sim_inspiral; for some reason we need both the trigtime
    # and the event number
    trigtimes=np.array([sim.geocent_end_time+1e-9*sim.geocent_end_time_ns \
            for sim in table])
    
    # reduce to specified values
    events=cp.get('injections', 'events')

    if events!='all':

        injevents=list(hyphen_range(events))

        trigtimes=trigtimes[injevents]


#############################################

# ----------------------------------------
# Setup analysis directory for deployment
# ----------------------------------------

topdir=os.getcwd()
os.chdir(workdir)

datafind = 'datafind'
if not os.path.exists(datafind): os.makedirs(datafind)

shutil.copy(cp.get('paths','bwb_executable'), '.')
shutil.copy(cp.get('paths','bwp_executable'), '.')


# -------------------
# Call LIGO Data find
# -------------------
cacheFiles = {}

if "LALSimAdLIGO" in channelList: opts.skip_datafind=True

if not opts.skip_datafind:


    # Set min/max gps times for LIGO data find:
    min_gps = min(trigtimes) - 25.0 # XXX: why 25???
    max_gps = max(trigtimes) + 25.0

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
        
        #print >> sys.stdout, "Copying frame files for input IS DUMB, EXITING"

        # XXX THIS NEEDS TO BE OVERHAULED SO WE SHIP FRAMES WITH JOBS, NOT COPY
        # FIRST!  DO this by adding the unique frames to a list whose elements
        # are then associated to individual jobs
        #sys.exit()

        # XXX Having found these files, we now want to copy them to the working
        # directory and make fresh, local cache files

        # 1) read cache file
        # 2) identify unique frames
        # 3) copy unique frames to datafind directory
        # 4) add to transfer files

        frames_to_copy=[]
        for ifo in ifoList:
            cache_entries = np.loadtxt('datafind/{ifo}.cache'.format(ifo=ifo),
                    dtype=str)

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
            shutil.copy(frame, 'datafind')


        #
        # Now we need to make a new, local cache file
        # - do this by manipulating the path string in the cache file to be relative 
        for ifo in ifoList:
            cache_file = 'datafind/{ifo}.cache'.format(ifo=ifo)
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


    for i,ifo in enumerate(ifoList):


        if "LALSimAdLIGO" not in channelList:
            # skipping datafind but using user-specified cache files

            cacheFile = \
                    os.path.abspath(cp.get('datafind','cacheFiles').split(',')[i])
            try:
                shutil.copy(cacheFile, 'datafind')
            except IOError:
                print >> sys.stderr, "Warning: cachefiles not found, remember to copy by hand"
            cacheFiles[ifo] = os.path.join('datafind',os.path.basename(cacheFile))

        else:

            cacheFiles[ifo] = cp.get('datafind','cacheFiles').split(',')[i]


#############################################



# -----------------------------------------------------------------------
# DAG Writing
# -----------------------------------------------------------------------

#
# Initialise DAG and Jobs
#

# ---- Create a dag to which we can add jobs.
dag = pipeline.CondorDAG(log=opts.user_tag+'.log')

# ---- Set the name of the file that will contain the DAG.
dag.set_dag_file( 'bayeswave_{0}'.format(opts.user_tag) )

# ---- Make instance of bayeswaveJob.
bwb_job = pipe_utils.bayeswaveJob(cp, cacheFiles, injfile=injfile,
        nrdata=nrdata)
bwp_job = pipe_utils.bayeswave_postJob(cp, cacheFiles, injfile=injfile,
        nrdata=nrdata)

#
# Build Nodes
#
if "LALSimAdLIGO" in channelList:
    try:
        dataseed=cp.getint('datafind', 'dataseed')
    except ConfigParser.NoOptionError:
        print >> sys.stderr, "datafind section requires dataseed for sim data"
        sys.exit()

for g,gps in enumerate(trigtimes):

    outputDir  = 'bayeswave_' + str(int(gps)) + '_' + str(uuid.uuid4())

    if not os.path.exists(outputDir): os.makedirs(outputDir)

    bwb_node = pipe_utils.bayeswaveNode(bwb_job)
    bwp_node = pipe_utils.bayeswave_postNode(bwp_job)

    # add options for bayeswave node
    bwb_node.set_trigtime(gps)
    bwb_node.set_PSDstart(gps)
    bwb_node.set_retry(1)
    bwb_node.set_outputDir(outputDir)

    if "LALSimAdLIGO" in channelList:
        bwb_node.set_dataseed(dataseed)
        bwp_node.set_dataseed(dataseed)
        dataseed+=1

    # add options for bayeswave_post node
    bwp_node.set_trigtime(gps)
    bwp_node.set_PSDstart(gps)
    bwp_node.set_retry(1)
    bwp_node.set_outputDir(outputDir)

    if injfile is not None:

        # STILL TO SUPPORT:
        # 1) xml, hdf5 file transfer
        # 2) xml needs to point to hdf5 correctly after transfer

        bwb_node.set_injevent(injevents[g])
        bwp_node.set_injevent(injevents[g])

    bwp_node.add_parent(bwb_node)

    # Add Nodes to DAG
    dag.add_node(bwb_node)
    dag.add_node(bwp_node)

#
# Finalise DAG
#
# ---- Write out the submit files needed by condor.
dag.write_sub_files()

# ---- Write out the DAG itself.
dag.write_dag()
dag.write_script()

# move back
os.chdir(topdir)



























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
executable=bayeswave
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
transfer_input_files=bayeswave,datafind,{outputDir},logs
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
fullcmdline = """./bayeswave \
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





 
