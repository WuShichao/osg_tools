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
import socket
import subprocess
import uuid
import fileinput
import ast

from glue import pipeline

import lalinspiral, lalburst
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

def job_times(trigtime, seglen, psdlen, padding):
    """
    Compute the gps times corresponding to a given trigger time

    psdstart = trigtime - (psdlen + padding)
    start = floor(min(psdstart, trigtime-0.5*seglen))
    stop  = ceil(max(start+psdlen, trigtime+0.5*Sseglen))

    returns segment(start,stop), psdstart

    so that start can be used easily as a psd start
    """

    psdstart=trigtime - (psdlen + padding)
    start = np.floor(min(psdstart, trigtime-0.5*seglen))
    stop = np.ceil(max(start+psdlen, trigtime+0.5*seglen))

    return segments.segment(start,stop), psdstart

def parser():
    """
    Parser for input (command line and ini file)
    """

    # --- cmd line
    parser = OptionParser()
    parser.add_option("-t", "--user-tag", default="", type=str)
    parser.add_option("-o", "--workdir", type=str, default=None)
    parser.add_option("--trigger-time", type=float, default=None)
    parser.add_option("--trigger-list", type=str, default=None)
    parser.add_option("--server", type=str, default=None)
    parser.add_option("--copy-frames", default=False, action="store_true")
    parser.add_option("--skip-datafind", default=False, action="store_true")

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

def read_injection_table(filename):
    try:
        sim_inspiral_table = lalinspiral.SimInspiralTableFromLIGOLw(filename,0,0)
        if sim_inspiral_table == -1:
            raise AttributeError("%s has no sim-inspiral table")
    
        print >> sys.stdout, "Reading trigger times from sim_inspiral table"
        trigger_times=[]
        while True:
            trigger_times.append(float(sim_inspiral_table.geocent_end_time))

            if sim_inspiral_table.next is None: break
            sim_inspiral_table = sim_inspiral_table.next
    
        return trigger_times

    except:
        pass
    
    try:
        sim_burst_table = lalburst.SimBurstTableFromLIGOLw(filename,None,None)
    
        print >> sys.stdout, "Reading trigger times from sim_burst table"
        trigger_times=[]
        while True:
            trigger_times.append(float(sim_burst_table.time_geocent_gps))
    
            if sim_burst_table.next is None: break
            sim_burst_table = sim_burst_table.next
    
        return trigger_times
    except:
        pass
    
    try:
        sim_ringdown_table = lalinspiral.SimRingdownTableFromLIGOLw(filename,0,0)
        if sim_ringdown_table == -1:
            raise AttributeError("%s has no sim-ringdown table")
    
        print >> sys.stdout, "Reading trigger times from sim_ringdown table"
        trigger_times=[]
        while True:
            trigger_times.append(float(sim_ringdown_table.geocent_start_time))
    
            if sim_ringdown_table.next is None: break
            sim_ringdown_table = sim_ringdown_table.next
    
        return trigger_times
    
    except:
        print >> sys.stderr, "Error: cannot read injection file %s"%injfile
        sys.exit()

# ----------------
# Set Parameters
# ----------------

# --- Parse options, arguments and ini file
opts, args, cp = parser()
cp.set('condor','copy-frames',str(opts.copy_frames))

workdir = opts.workdir 
if not os.path.exists(workdir): 
    print >> sys.stdout, "making work-directory: %s"%workdir
    os.makedirs(workdir)
else:
    print >> sys.stderr, "WARNING: work-directory %s exists"%workdir


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
    try:
        rom_data_path = os.path.join(os.environ['LAL_DATA_PATH'],
                'SEOBNRv2ChirpTimeSS.dat')
        shutil.copy(rom_data_path, workdir)
    except KeyError:
        print >> sys.stderr, "Warning: LAL_DATA_PATH not set"
        print >> sys.stderr, \
                "CBC injections require SEOBNRv2ChirpTimeSS.dat in \
LAL_DATA_PATH"


# NR HDF5 data
try:
    nrdata=cp.get('injections', 'nrhdf5')
    nr_full_path=cp.get('injections', 'nrhdf5')
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
#   Use cases: single trigger, ascii trigger list, sim_inspiral table
#   Not yet supported: sim_burst, time-slides
#
if opts.trigger_time is not None:
    trigger_times = [opts.trigger_time]

if opts.trigger_list is not None:
    trigger_times = np.loadtxt(opts.trigger_list)

if injfile is not None:

    #
    # Read injection file
    #
    filename=os.path.join(workdir,injfile)
    trigger_times = read_injection_table(filename)
    
    # reduce to specified values
    events=cp.get('injections', 'events')

    if events!='all':
        injevents=list(hyphen_range(events))
        trigger_times=np.array(trigger_times)[injevents]


#
# --- Determine min/max times for data coverage
#
seglen = cp.getfloat('input','seglen')
psdlen = cp.getfloat('input','PSDlength')
padding = cp.getfloat('input','padding')

if cp.has_option('input','gps-start-time'):
    gps_start_time = cp.getint('input','gps-start-time')
else:
    seg, _ = job_times(min(trigger_times), seglen, psdlen, padding)
    gps_start_time = seg[0]
    cp.set('input','gps-start-time',str(int(gps_start_time)))

if cp.has_option('input','gps-end-time'):
    gps_end_time = cp.getint('input','gps-end-time')
else:
    seg,_ = job_times(max(trigger_times), seglen, psdlen, padding)
    gps_end_time = seg[1]
    cp.set('input','gps-end-time',str(int(gps_end_time)))

#############################################

# ----------------------------------------
# Setup analysis directory for deployment
# ----------------------------------------

topdir=os.getcwd()
os.chdir(workdir)

datafind_dir = 'datafind'
if not os.path.exists(datafind_dir): os.makedirs(datafind_dir)
if cp.has_option('injections', 'mdc-cache'):
    shutil.copy(cp.get('injections', 'mdc-cache'),
            os.path.join('datafind','MDC.cache'))


segment_dir = 'segments'
if not os.path.exists(segment_dir): os.makedirs(segment_dir)

shutil.copy(cp.get('bw_paths','bwb_executable'), '.')
shutil.copy(cp.get('bw_paths','bwp_executable'), '.')


################################################

# ----------------------------------------------
# Data Acquisition: gw_data_find & segdb queries
# ----------------------------------------------


#
# --- datafind params from config file
#

ifo_list=ast.literal_eval(cp.get('input','ifo-list'))
channel_list=ast.literal_eval(cp.get('datafind', 'channel-list'))
frtype_list=ast.literal_eval(cp.get('datafind', 'frtype-list'))

cache_files = {}
segmentList = {}
framePaths={}
frameSegs={}

if not opts.skip_datafind:

    for ifo in ifo_list:

        if frtype_list[ifo] == "LALSimAdLIGO": 
            cache_files[ifo] = "LALSimAdLIGO"
            segmentList[ifo] = \
                    segments.segmentlist([segments.segment(gps_start_time,
                        gps_end_time)])
            continue
        else:

            #
            # --- Run DataFind query to produce cache files for frames
            #
            cachefilefmt = os.path.join(datafind_dir, '{0}.cache')

            if opts.server is not None:
                ldfcmd = "gw_data_find --observatory {o} --type {frtype} \
    -s {gps_start_time} -e {gps_end_time} --lal-cache\
    --server={server} -u {url_type} > {cachefile}".format(
                        o=ifo[0], frtype=frtype_list[ifo],
                        cachefile=cachefilefmt.format(ifo),
                        gps_start_time=gps_start_time,
                        gps_end_time=gps_end_time, server=opts.server,
                        url_type=cp.get('datafind','url-type'))
            else:
                ldfcmd = "gw_data_find --observatory {o} --type {frtype} -s \
{gps_start_time} -e {gps_end_time} --lal-cache -u {url_type}>\
{cachefile}".format( o=ifo[0], frtype=frtype_list[ifo],
    cachefile=cachefilefmt.format(ifo), gps_start_time=gps_start_time,
    gps_end_time=gps_end_time, url_type=cp.get('datafind','url-type'))
            print >> sys.stdout, "Calling LIGO data find ..."
            print >> sys.stdout, ldfcmd

            subprocess.call(ldfcmd, shell=True)

            cache_files[ifo]=os.path.join('datafind', '{0}.cache'.format(ifo))

            #
            # Record frame segments so we can identify frames for OSG transfers
            #
            frameSegs[ifo] = segmentsUtils.fromlalcache(open(cache_files[ifo]))


            #
            # --- Run segdb query
            #
            # 1) Find segments \in [gps-start-time, gps-end-time]
            # 2) Loop over trigger_times and identify analyzeable jobs

            if cp.has_option('datafind','veto-categories'):
              veto_categories=ast.literal_eval(cp.get('datafind','veto-categories'))
            else: veto_categories=[]

            curdir=os.getcwd()
            os.chdir(segment_dir)

            (segFileName,dqVetoes)=inspiralutils.findSegmentsToAnalyze(cp, ifo,
                    veto_categories, generate_segments=True,
                    use_available_data=False, data_quality_vetoes=False)

            segfile=open(segFileName)
            segmentList[ifo]=segmentsUtils.fromsegwizard(segfile)
            segmentList[ifo].coalesce()
            segfile.close()

            if segmentList[ifo] == []:
                print >> sys.stderr, "No matching segments for %s"%ifo
                sys.exit()

            os.chdir(curdir)


        # --------------------------------------------------------------------
        # Set up cache files to point to local copies of frames in the working
        # directory

        if opts.copy_frames:

            #
            # Now we need to make a new, local cache file
            # - do this by manipulating the path string in the cache file to be relative 
            cache_file = 'datafind/{ifo}.cache'.format(ifo=ifo)
            shutil.copy(cache_file, cache_file.replace('cache','cache.bk'))

            cache_entries = np.loadtxt(cache_file, dtype=str)
            if cache_entries.ndim==1: cache_entries = [cache_entries]
            
            framePaths[ifo]=[]
            new_cache = open(cache_file, 'w')
            for c,cache_entry in enumerate(cache_entries):
                frame = cache_entry[-1].split('localhost')[-1]
                framePaths[ifo].append(frame)

                #local_path=os.path.join('datafind',cache_entry[4].split('/')[-1])
                local_path=cache_entry[4].split('/')[-1]

                new_cache.writelines('{ifo} {type} {gps} {length} {path}\n'.format(
                    ifo=ifo, type=cache_entry[1], gps=cache_entry[2],
                    length=cache_entry[3], path=local_path))

            new_cache.close()

else:

    print "SKIPPING DATAFIND & SEGDB"

    for ifo in ifo_list:

        if  frtype_list[ifo] == "LALSimAdLIGO":
            cache_files[ifo] == "LALSimAdLIGO"
        else:
            cache_files[ifo] = os.path.join('datafind','%s.cache'%ifo)

        segmentList[ifo] = segments.segment(gps_start_time, gps_end_time)


#############################################



# -----------------------------------------------------------------------
# DAG Writing
# -----------------------------------------------------------------------

#
# Initialise DAG and Jobs
#

# ---- Create a dag to which we can add jobs.
dag = pipeline.CondorDAG(log=opts.workdir+'.log')

# ---- Set the name of the file that will contain the DAG.
dag.set_dag_file( 'bayeswave_{0}'.format(opts.workdir) )

# ---- Make instance of bayeswaveJob.
bwb_job = pipe_utils.bayeswaveJob(cp, cache_files, injfile=injfile,
        nrdata=nrdata)
bwp_job = pipe_utils.bayeswave_postJob(cp, cache_files, injfile=injfile,
        nrdata=nrdata)

#
# Build Nodes
#
if "LALSimAdLIGO" in cache_files.values():
    try:
        dataseed=cp.getint('input', 'dataseed')
    except ConfigParser.NoOptionError:
        print >> sys.stderr, "[input] section requires dataseed for sim data"
        print >> sys.stderr, "...removing %s"%workdir
        os.chdir(topdir)
        shutil.rmtree(workdir)
        sys.exit()

unanalyzeable_jobs = []

# XXX: Testing times (one of these has no data)
#trigger_times = [1126252133, 1126259365]

transferFrames={}
for t, trigger_time in enumerate(trigger_times):

    print >> sys.stdout, "---------------------------------------"

    # -------------------------------------------
    # Check job times fall within available data
    job_segment, psd_start = job_times(trigger_time, seglen, psdlen, padding)

    for ifo in ifo_list:

        job_in_segments = [seg.__contains__(job_segment) \
                for seg in segmentList[ifo]]

        if not any(job_in_segments):

            bad_job={}
            bad_job['ifo']=ifo
            bad_job['trigger_time']=trigger_time
            bad_job['seglen']=seglen
            bad_job['psdlen']=psdlen
            bad_job['padding']=padding
            bad_job['job_segment']=job_segment
            bad_job['data_segments']=segmentList[ifo]

            unanalyzeable_jobs.append(bad_job)
            
            print >> sys.stderr, "Warning: No matching %s segments for job %d of %d"%(
                    ifo, t+1, len(trigger_times))
            print >> sys.stderr, bad_job
            break

    else:

        print >> sys.stdout, "Adding node for GPS %d (%d of %d)"%(trigger_time, t+1,
                len(trigger_times))


        if "LALSimAdLIGO" not in cache_files.values():
            #
            # Identify frames associated with this job
            if opts.copy_frames:
                for ifo in ifo_list:
                    frame_idx = [seg.intersects(job_segment) for seg in frameSegs[ifo]]
                    transferFrames[ifo] = [frame for f,frame in
                            enumerate(framePaths[ifo]) if frame_idx[f]] 

        outputDir  = 'bayeswave_' + str(int(trigger_time)) + '_' + str(uuid.uuid4())

        if not os.path.exists(outputDir): os.makedirs(outputDir)

        bwb_node = pipe_utils.bayeswaveNode(bwb_job)
        bwp_node = pipe_utils.bayeswave_postNode(bwp_job)


        # add options for bayeswave node
        bwb_node.set_trigtime(trigger_time)
        bwb_node.set_PSDstart(psd_start)
        bwb_node.set_retry(1)
        bwb_node.set_outputDir(outputDir)
        if transferFrames: bwb_node.add_frame_transfer(transferFrames)

        if "LALSimAdLIGO" in cache_files.values():
            bwb_node.set_dataseed(dataseed)
            bwp_node.set_dataseed(dataseed)
            #gpsNow = int(os.popen('lalapps_tconvert now').readline())
            #dataseed+=np.random.randint(1,gpsNow)
            dataseed+=1

        # add options for bayeswave_post node
        bwp_node.set_trigtime(trigger_time)
        bwp_node.set_PSDstart(psd_start)
        bwp_node.set_retry(1)
        bwp_node.set_outputDir(outputDir)

        if injfile is not None:

            bwb_node.set_injevent(injevents[t])
            bwp_node.set_injevent(injevents[t])

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

# print some summary info:
if len(trigger_times)-len(unanalyzeable_jobs)>0:
    print """
    Total number of requested trigger times: {ntrigs_desired}
    Number of triggers successfully added to DAG: {ntrigs_added}
    Number of triggers failing data criteria: {ntrigs_failed}

    To submit:
        cd {workdir}
        condor_submit_dag {dagfile}
    """.format(ntrigs_desired=len(trigger_times),
            ntrigs_added=len(trigger_times)-len(unanalyzeable_jobs),
            ntrigs_failed=len(unanalyzeable_jobs),
            workdir=workdir, dagfile=dag.get_dag_file())
else:
    print ""
    print "No analyzeable jobs in requested time"





