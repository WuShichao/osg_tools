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

from lalapps import inspiralutils
from glue import segmentsUtils, segments

from optparse import OptionParser
import ConfigParser

import bwb_pipe_utils as pipe_utils

#############################################
#
# Local function defs

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

def dump_job_info(job_dir, trigger):
    """
    Writes a text file with job info to outputDir:

    GPS time, lag or GraceID, frequency, and cWB’s rho
    """
    f=open(os.path.join(job_dir, 'job_info.txt'), 'w')

    f.write('# rho gps lag freq veto1 veto2 graceID\n')
    f.write('{rho} {gps_time} {time_lag} {trig_frequency} {veto1} \
{veto2} {graceID}\n'.format(
        gps_time=trigger.trigger_time,
        time_lag=trigger.time_lag,
        trig_frequency=trigger.trigger_frequency,
        rho=trigger.rho,
        veto1=trigger.veto1,
        veto2=trigger.veto2,
        graceID=trigger.graceID))
    f.close()



def parser():
    """
    Parser for input (command line and ini file)
    """

    # --- cmd line
    parser = OptionParser()
    parser.add_option("-t", "--user-tag", default="", type=str)
    parser.add_option("-r", "--workdir", type=str, default=None)
    parser.add_option("--trigger-time", type=float, default=None)
    parser.add_option("--trigger-list", type=str, default=None)
    parser.add_option("--cwb-trigger-list", type=str, default=None)
    parser.add_option("--server", type=str, default=None)
    parser.add_option("--copy-frames", default=False, action="store_true")
    parser.add_option("--skip-datafind", default=False, action="store_true")
    parser.add_option("--sim-data", , default=False, action="store_true")
    parser.add_option("-I", "--injfile", default=None)
    parser.add_option("-G", "--graceID", default=None)
    parser.add_option("--graceID-list", default=None)
    parser.add_option("--bw-inject", default=False, action="store_true")
    parser.add_option("--condor-submit", default=False, action="store_true")
    parser.add_option("--submit-to-gracedb", default=False, action="store_true")
    parser.add_option("--html-root", default=None)
    parser.add_option("--skip-megapy", default=False, action="store_true")
    parser.add_option("--tidy-up", default=False, action="store_true")


    (opts,args) = parser.parse_args()

    if opts.workdir is None:
        print >> sys.stderr, "ERROR: must specify --workdir"
        sys.exit()


    if len(args)==0:
        print >> sys.stderr, "ERROR: require config file"
        sys.exit()
    if not os.path.isfile(args[0]):
        print >> sys.stderr, "ERROR: config file %s does not exist"%args[0]
        sys.exit()


    # --- Read config file
    cp = ConfigParser.ConfigParser()
    cp.optionxform = str
    cp.read(args[0])


    return opts, args, cp

# END --- Local function defs
#############################################

#############################################
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

# Injection file (e.g., sim-inspiral table).  Try commandline first, if none,
# try config file
injfile=opts.injfile
if injfile is None:
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
# Get Trigger Info
#

# XXX: Careful, there's nothing here to handle the non-exclusivity of these
# options other than common sense
if opts.trigger_time is not None and not\
    cp.has_option('bayeswave_options','BW-inject'):
    #
    # Read trigger from commandline
    #
    trigger_list = pipe_utils.triggerList(cp, [opts.trigger_time])

if opts.trigger_list is not None:
    #
    # Read triggers from ascii list 
    #
    trigger_list = pipe_utils.triggerList(cp, trigger_file=opts.trigger_list)

if opts.cwb_trigger_list is not None:
    #
    # Read triggers from ascii list 
    #
    trigger_list = pipe_utils.triggerList(cp, cwb_trigger_file=opts.cwb_trigger_list)

if injfile is not None:
    #
    # Read injection file
    #
    filename=os.path.join(workdir,injfile)
    trigger_list = pipe_utils.triggerList(cp, injection_file=filename)

if cp.has_option('bayeswave_options','BW-inject'):
    # Check the option is valid:
    if cp.get('bayeswave_options','BW-inject') not in ['signal','glitch']:
        print >> sys.stderr, "Error: BW-inject must be in ", ['signal','glitch']
        sys.exit()
    #
    # Perform internal injections drawn from the signal or glitch model
    #
    if opts.trigger_time is None:
        opts.trigger_time=1126259462.392 
    print >> sys.stdout, "Setting trigger time to %f"%opts.trigger_time
    trigger_list = pipe_utils.triggerList(cp, gps_times=opts.trigger_time,
            internal_injections=True)
    
# GraceDB support
if opts.graceID is not None:

    graceIDs = [opts.graceID]
    trigger_list = pipe_utils.triggerList(cp, graceIDs=graceIDs)


if opts.graceID_list is not None:

    graceIDs = np.loadtxt(opts.graceID_list)
    trigger_list = pipe_utils.triggerList(cp, graceIDs=graceIDs)

if opts.submit_to_gracedb:
    if opts.html_root is None:
        html_root = cp.get('bayeswave_paths', 'html-root')
    else:
        html_root = opts.html_root
    if html_root is None:
        print >> sys.stder, "demanding submit to gdb but no html-root"
        sys.exit()



    if not os.path.exists(html_root):
        os.makedirs(html_root)
    else:
        print >> sys.stderr, "Warning: html-root %s exists"%html_root


# Extract trigger times for readability
trigger_times = [trig.trigger_time for trig in trigger_list.triggers]
lag_times = [trig.time_lag for trig in trigger_list.triggers]

#
# --- Determine min/max times for data coverage
#
seglen = cp.getfloat('input','seglen')
psdlen = cp.getfloat('input','PSDlength')
padding = cp.getfloat('input','padding')

if cp.has_option('input','gps-start-time'):
    gps_start_time = cp.getint('input','gps-start-time')
else:
    trigtime = min(trigger_times) - (max(np.absolute(lag_times))+25.0)
    seg, _ = job_times(trigtime, seglen, psdlen, padding)
    gps_start_time = seg[0]

if cp.has_option('input','gps-end-time'):
    gps_end_time = cp.getint('input','gps-end-time')
else:
    trigtime = max(trigger_times) + (max(np.absolute(lag_times))+25.0)
    seg,_ = job_times(trigtime, seglen, psdlen, padding)
    gps_end_time = seg[1]

# Timelag adjustment
#gps_start_time = min(trigger_times) - (max(np.absolute(lag_times))+25.0)
#gps_end_time   = max(trigger_times) + (max(np.absolute(lag_times))+25.0)

# Update config parser
cp.set('input','gps-start-time',str(int(gps_start_time)))
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

# Decide whether simulating data
if not cp.has_option('input','sim-data'):
    cp.set('input', 'sim-data', opts.sim_data)
elif cp.has_option('input','sim-data') and opts.sim_data:
    # Override the config file with the command line
    cp.set('input', 'sim-data', opts.sim_data)

cache_files = {}
segmentList = {}
framePaths={}
frameSegs={}

#
# --- Handle special cases for segdb
#

if (opts.cwb_trigger_list is not None) \
        or (opts.trigger_list is not None) \
        or (opts.graceID is not None) \
        or (opts.graceID_list is not None):

    # Assume triggers lie in analyzeable segments
    cp.set('datafind','ignore-science-segments', True)

for ifo in ifo_list:

    if cp.get('input','sim-data'):
        # Get the type of simulated data from the frame type list
        # E.g., to simulate from LALSimAdLIGO put this in the config.ini:
        #   frtype-list={'H1':'LALSimAdLIGO','L1':'LALSimAdLIGO'}
        # To simulate from arbitraray *A*SD:
        #   frtype-list={'H1':'interp:/home/tyson/O2/review/bayesline/IFO0_asd.dat','L1':'interp:/home/tyson/O2/review/bayesline/IFO0_asd.dat'}

        cache_files[ifo] = frtype_list[ifo]

        segmentList[ifo] = \
                segments.segmentlist([segments.segment(gps_start_time,
                    gps_end_time)])

    else:

        #
        # --- Run DataFind query to produce cache files for frames
        #
        cachefilefmt = os.path.join(datafind_dir, '{0}.cache')

        if opts.skip_datafind:
            continue
        else:

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

        # Record frame segments so we can identify frames for OSG transfers
        if opts.skip_datafind:
            # XXX: if no datafind, assume frames include jobs.  But should
            # change to take cache file locations
            frameSegs[ifo] = \
                    segments.segmentlist([segments.segment(gps_start_time,
                        gps_end_time)])
        else:
            frameSegs[ifo] = segmentsUtils.fromlalcache(open(cache_files[ifo]))

        if cp.has_option('datafind','ignore-science-segments'):
            segmentList[ifo] = \
                    segments.segmentlist([segments.segment(gps_start_time,
                        gps_end_time)])
        else:

            #
            # --- Run segdb query
            #

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
            print "Setting up frame copying"

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

#########################################################################


#########################################################################
# -----------------------------------------------------------------------
# DAG Writing
# -----------------------------------------------------------------------

#
# Initialise DAG and Jobs
#

# ---- Create a dag to which we can add jobs.
dag = pipeline.CondorDAG(log=opts.workdir+'.log')

# ---- Set the name of the file that will contain the DAG.
dag.set_dag_file( 'bayeswave_{0}'.format(os.path.basename(opts.workdir)) )

# ---- Create DAG jobs
#   bayeswave: main bayeswave analysis
#   bayeswave_post: bayeswave_post
#   megasky: skymap job
#   megaplot: remaining plots & webpage generation
#   submitToGraceDB: upload skymap & PE to graceDB (optional)
bayeswave_job = pipe_utils.bayeswaveJob(cp, cache_files, injfile=injfile,
        nrdata=nrdata)
bayeswave_post_job = pipe_utils.bayeswave_postJob(cp, cache_files, injfile=injfile,
        nrdata=nrdata)
megasky_job = pipe_utils.megaskyJob(cp)
megaplot_job = pipe_utils.megaplotJob(cp)

if opts.submit_to_gracedb:
    submitToGraceDB_job = pipe_utils.submitToGraceDB(cp)

if opts.tidy_up:
    archiver_job = pipe_utils.archiverJob(cp)

#
# Build Nodes
#
#if "LALSimAdLIGO" in cache_files.values():
# XXX: post jobs currently require a data seed for the dummy LALSimAdLIGO data
try:
    dataseed=cp.getint('input', 'dataseed')
except ConfigParser.NoOptionError:
    print >> sys.stderr, "[input] section requires dataseed for sim data"
    print >> sys.stderr, " (you need this in bayeswave_post, even if real data"
    print >> sys.stderr, "...removing %s"%workdir
    os.chdir(topdir)
    shutil.rmtree(workdir)
    sys.exit()

unanalyzeable_jobs = []

# XXX: Testing times (one of these has no data)
#trigger_times = [1126252133, 1126259365]

transferFrames={}
totaltrigs=0
for t,trigger in enumerate(trigger_list.triggers):

    print >> sys.stdout, "---------------------------------------"

    # -------------------------------------------
    # Check job times fall within available data
    job_segment, psd_start = job_times(trigger.trigger_time, seglen, psdlen, padding)

    for ifo in ifo_list:

        job_in_segments = [seg.__contains__(job_segment) \
                for seg in segmentList[ifo]]

        if not any(job_in_segments):

            bad_job={}
            bad_job['ifo']=ifo
            bad_job['trigger_time']=trigger.trigger_time
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

        print >> sys.stdout, "Adding node for GPS %d, L1-timeslide %f (%d of %d)"%(
                trigger.trigger_time, trigger.time_lag, totaltrigs+1,
                len(trigger_times))


        if "LALSimAdLIGO" not in cache_files.values():
            #
            # Identify frames associated with this job
            if opts.copy_frames:
                for ifo in ifo_list:
                    frame_idx = [seg.intersects(job_segment) for seg in frameSegs[ifo]]
                    transferFrames[ifo] = [frame for f,frame in
                            enumerate(framePaths[ifo]) if frame_idx[f]] 

        outputDir  = 'bayeswave_' + str(int(trigger.trigger_time)) + '_' + \
                str(float(trigger.time_lag)) + '_' + str(uuid.uuid4())

        if not os.path.exists(outputDir): os.makedirs(outputDir)

        # XXX Dump job info file
        #GPS time, lag or GraceID, frequency, and cWB’s rho
        dump_job_info(outputDir, trigger) 

        # Create DAG nodes
        #   bayeswave: main bayeswave analysis
        #   bayeswave_post: bayeswave_post
        #   megasky: skymap job
        #   megaplot: distribution plots & webpage generation
        #   submitToGraceDB: upload skymap & PE to graceDB (optional)
        bayeswave_node = pipe_utils.bayeswaveNode(bayeswave_job)
        bayeswave_post_node = pipe_utils.bayeswave_postNode(bayeswave_post_job)
        megasky_node = pipe_utils.megaskyNode(megasky_job, outputDir)
        megaplot_node = pipe_utils.megaplotNode(megaplot_job, outputDir)

        if opts.submit_to_gracedb:
            htmlDir=os.path.join(html_root, outputDir)
            if not os.path.exists(htmlDir):
                os.makedirs(htmlDir)
            gracedb_node = pipe_utils.submitToGraceDBNode(submitToGraceDB_job,
                    outputDir, htmlDir)

        if opts.tidy_up:
            archiver_node = pipe_utils.archiverNode(archiver_job, outputDir)

        #
        # --- Add options for bayeswave node
        #
        bayeswave_node.set_trigtime(trigger.trigger_time)
        bayeswave_node.set_segment_start(trigger.trigger_time -
                trigger.seglen/2.)
        bayeswave_node.set_srate(trigger.srate)
        bayeswave_node.set_seglen(trigger.seglen)
        bayeswave_node.set_window(trigger.window)
        bayeswave_node.set_PSDstart(psd_start)
        bayeswave_node.set_outputDir(outputDir)
        if transferFrames: bayeswave_node.add_frame_transfer(transferFrames)

        if "LALSimAdLIGO" in cache_files.values():
            bayeswave_node.set_dataseed(dataseed)
        bayeswave_post_node.set_dataseed(dataseed)
        dataseed+=1


        if cp.has_option('bayeswave_options','BW-inject'):
            bayeswave_node.set_BW_event(trigger.BW_event)

        #
        # --- Add options for bayeswave_post node
        #
        bayeswave_post_node.set_trigtime(trigger.trigger_time)
        bayeswave_post_node.set_segment_start(trigger.trigger_time -
                trigger.seglen/2.)
        bayeswave_post_node.set_srate(trigger.srate)
        bayeswave_post_node.set_seglen(trigger.seglen)
        bayeswave_post_node.set_window(trigger.window)
        bayeswave_post_node.set_PSDstart(psd_start)
        bayeswave_post_node.set_outputDir(outputDir)

        if injfile is not None:
            bayeswave_node.set_injevent(trigger.injevent)
            bayeswave_post_node.set_injevent(trigger.injevent)

        bayeswave_node.set_L1_timeslide(trigger.time_lag)
        bayeswave_post_node.set_L1_timeslide(trigger.time_lag)

        if cp.has_option('bayeswave_options','BW-inject'):
            bayeswave_post_node.set_BW_event(trigger.BW_event)

        #
        # --- Add options for mega-scripts
        #
        megasky_node.set_outputDir(outputDir)
        megaplot_node.set_outputDir(outputDir)


        #
        # --- Add parent/child relationships
        #
        bayeswave_post_node.add_parent(bayeswave_node)
        if not opts.skip_megapy:
            megasky_node.add_parent(bayeswave_post_node)
            megaplot_node.add_parent(bayeswave_post_node) 
        if opts.submit_to_gracedb:
            gracedb_node.add_parent(megaplot_node) 
            gracedb_node.add_parent(megasky_node) 
        if opts.tidy_up:
            if opts.skip_megapy:
                archiver_node.add_parent(bayeswave_post_node) 
                archiver_node.add_parent(bayeswave_post_node) 
            else:
                archiver_node.add_parent(megaplot_node)
                archiver_node.add_parent(megasky_node) 

            archiver_node.set_post_script(cp.get("bayeswave_paths","cleaner"))
            archiver_node.add_post_script_arg(outputDir)

        # Add Nodes to DAG
        dag.add_node(bayeswave_node)
        dag.add_node(bayeswave_post_node)
        if not opts.skip_megapy:
            dag.add_node(megasky_node)
            dag.add_node(megaplot_node)
        if opts.submit_to_gracedb:
            dag.add_node(gracedb_node)
        if opts.tidy_up:
            dag.add_node(archiver_node)


        # --- Add

        totaltrigs+=1


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


if opts.condor_submit:
    # Auto-submit dag by cd-ing into the work-directory and submitting
    # chdir is useful with the OSG-friendly relative paths

    print "Submitting DAG..."
     
    os.chdir(workdir)
    x = subprocess.Popen(['condor_submit_dag',dag.get_dag_file()])
    x.wait()
    if x.returncode==0:
        print 'Submitted DAG file: ',dag.get_dag_file()
    else:
        print 'Unable to submit DAG file'
    os.chdir(topdir)




