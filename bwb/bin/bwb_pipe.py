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
import ast

from glue import pipeline

from glue.ligolw import lsctables
from glue.ligolw import table
from glue.ligolw import utils
from pylal.SimInspiralUtils import ExtractSimInspiralTableLIGOLWContentHandler
lsctables.use_in(ExtractSimInspiralTableLIGOLWContentHandler)

from lalapps import inspiralutils
from glue import segmentsUtils, segments

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

def job_times(trigtime, psdlen, padding):
    """
    Compute the gps times corresponding to a given trigger time

    start = min(trigtime - (psdlen + padding), trigtime-0.5*seglen)
    stop  = max(start+psdlen, trigtime+0.5*Sseglen)
    """

    start = min(trigtime - (psdlen + padding), trigtime-0.5*seglen)
    stop = max(start+psdlen, trigtime+0.5*Sseglen)

    return segments.segment(start,stop)

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
        print >> sys.stderr, "ERROR: must specify --workdir"
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
    # Copy injfile locally
    shutil.copy(injfile, workdir)
    injfile=os.path.basename(injfile)

    # Also locate and copy rom chirp time data
    rom_data_path = os.path.join(os.environ.get('LAL_DATA_PATH'),
            'SEOBNRv2ChirpTimeSS.dat')
    shutil.copy(rom_data_path, workdir)


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


#
# --- Determine min/max times for data coverage
#
padding = cp.getint('input','padding')
psdlen = cp.getint('input','PSDlength')

if cp.has_option('input','gps-start-time'):
    gps_start_time = cp.getint('input','gps-start-time')
else:
    gps_start_time = job_times(min(trigtimes), psdlen, padding)[0]
    cp.set('input','gps-start-time',gps_start_time)

if cp.has_option('input','gps-end-time'):
    gps_end_time = cp.getint('input','gps-end-time')
else:
    gps_end_time = job_times(max(trigtimes), psdlen, padding)[1]
    cp.set('input','gps-end-time',gps_end_time)


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


################################################

# ----------------------------------------------
# Data Acquisition: gw_data_find & segdb queries
# ----------------------------------------------


#
# --- datafind params from config file
#

ifoList=ast.literal_eval(cp.get('datafind', 'ifoList'))
channelList=ast.literal_eval(cp.get('datafind', 'channelList'))
frtypeList=ast.literal_eval(cp.get('datafind', 'frtypeList'))

cacheFiles = {}
segments = {}

if not opts.skip_datafind:

    for ifo in ifoList:

        if frtypeList[ifo] == "LALSimAdLIGO": 
            cacheFiles[ifo]= "LALSimAdLIGO"
            continue
        else:

            #
            # --- Run DataFind query to produce cache files for frames
            #
            cachefilefmt = os.path.join(datafind, '{0}.cache')

            if opts.server is not None:
                ldfcmd = "gw_data_find --observatory {o} --type {frtype} -s {gps_start_time} -e\
        {gps_end_time} --lal-cache --server={server} | grep file > {cachefile}".format(
                        o=ifo[0], frtype=frtype, cachefile = cachefilefmt.format(ifo),
                        gps_start_time=gps_start_time, gps_end_time=gps_end_time, server=opts.server)
            else:
                ldfcmd = "gw_data_find --observatory {o} --type {frtype} -s {gps_start_time} -e {gps_end_time} --lal-cache | grep file > {cachefile}".format(
                        o=ifo[0], frtype=frtype, cachefile = cachefilefmt.format(ifo),
                        gps_start_time=gps_start_time, gps_end_time=gps_end_time)
            print >> sys.stdout, "Calling LIGO data find ..."
            print >> sys.stdout, ldfcmd

            subprocess.call(ldfcmd, shell=True)

            cacheFiles[ifo]=os.path.join('datafind', '{0}.cache'.format(ifo))

            #
            # --- Run segdb query
            #
            # 0) move into datafind directory
            # 1) Find segments \in [gps-start-time, gps-end-time]
            # 2) Loop over trigtimes and identify analyzeable jobs
            # 3) move out of datafind directory

            if config.has_option('datafind','veto-categories'):
              veto_categories=ast.literal_eval(config.get('datafind','veto-categories'))
            else: veto_categories=[]

            (segFileName,dqVetoes)=inspiralutils.findSegmentsToAnalyze(config, ifo,
                    veto_categories, generate_segments=True,
                    use_available_data=False, data_quality_vetoes=False)

            segfile=open(segFileName)
            segments[ifo]=segmentsUtils.fromsegwizard(segfile)


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


    for ifo in ifoList:


        if "LALSimAdLIGO" not in channelList:
            # skipping datafind but using user-specified cache files

            cacheFile = \
                    os.path.abspath(ast.literal_eval(cp.get('datafind','cacheFiles')[ifo]))
            try:
                shutil.copy(cacheFile, 'datafind')
            except IOError:
                print >> sys.stderr, "Warning: cachefiles not found, remember to copy by hand"
            cacheFiles[ifo] = os.path.join('datafind',os.path.basename(cacheFile))

        else:

            cacheFiles[ifo] = \
                    ast.literal_eval(cp.get('datafind','cacheFiles')[ifo])


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

