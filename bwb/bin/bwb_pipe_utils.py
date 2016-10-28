#!/usr/bin/env python     
# -*- coding: utf-8 -*-     
# Copyright (C) 2016-2017 James Clark <james.clark@ligo.org>     
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

# DAG Class definitions for bayeswave

import sys,os,subprocess
# Test BAYESWAVE_PREFIX
try:
    BAYESWAVE_PREFIX = os.environ['BAYESWAVE_PREFIX']
except KeyError:
    print >> sys.stdout, "BAYESWAVE_PREFIX is unset"
    sys.exit()

from glue import pipeline
from glue.ligolw import ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables
#import lalinspiral, lalburst

import ConfigParser
import itertools
import socket
import ast
import numpy as np
import random

# define a content handler
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass
lsctables.use_in(LIGOLWContentHandler)

#
# Convenience Defs
#
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

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

    xmldoc = ligolw_utils.load_filename(filename, contenthandler =
            LIGOLWContentHandler, verbose = True)

    sim_inspiral_table = lsctables.SimInspiralTable.get_table(xmldoc)

    return ( sim_inspiral_table.get_column('geocent_end_time') + \
            1e-9*sim_inspiral_table.get_column('geocent_end_time') )



#    try:
#       sim_inspiral_table = lalinspiral.SimInspiralTableFromLIGOLw(filename,0,0)
#       if sim_inspiral_table == -1:
#           raise AttributeError("%s has no sim-inspiral table")
#
#       trigger_times=[]
#       while True:
#           trigger_times.append(float(sim_inspiral_table.geocent_end_time))
#
#           if sim_inspiral_table.next is None: break
#           sim_inspiral_table = sim_inspiral_table.next
#
#       return trigger_times
#
#    except:
#        pass
    
#   try:
#       sim_burst_table = lalburst.SimBurstTableFromLIGOLw(filename,None,None)
#   
#       print >> sys.stdout, "Reading trigger times from sim_burst table"
#       trigger_times=[]
#       while True:
#           trigger_times.append(float(sim_burst_table.time_geocent_gps))
#   
#           if sim_burst_table.next is None: break
#           sim_burst_table = sim_burst_table.next
#   
#       return trigger_times
#   except:
#       pass
#   
#   try:
#       sim_ringdown_table = lalinspiral.SimRingdownTableFromLIGOLw(filename,0,0)
#       if sim_ringdown_table == -1:
#           raise AttributeError("%s has no sim-ringdown table")
#   
#       print >> sys.stdout, "Reading trigger times from sim_ringdown table"
#       trigger_times=[]
#       while True:
#           trigger_times.append(float(sim_ringdown_table.geocent_start_time))
#   
#           if sim_ringdown_table.next is None: break
#           sim_ringdown_table = sim_ringdown_table.next
#   
#       return trigger_times
#   
#   except:
#       print >> sys.stderr, "Error: cannot read injection file %s"%injfile
#       sys.exit()

class eventTrigger:
    """
    Stores event characteristics and determines run configuration for this event
    """
    def __init__(self, cp, trigger_time=None, time_lag=0.0,
            trigger_frequency=None, rho=None, graceID=None, injevent=None,
            frequency_threshold=200., default_srate=1024., min_srate=1024.,
            max_srate=4096., default_seglen=4., max_seglen=4., min_seglen=2.,
            default_window=1.0, min_window=0.5, max_window=1.0, veto1=None,
            veto2=None, BW_event=None):


        #
        # Get run configuration
        #
        try:
            self.frequency_threshold = cp.getfloat('input', 'frequency-threshold')
        except:
            self.frequency_threshold = frequency_threshold

        try:
            self.default_srate = cp.getfloat('input', 'srate')
        except:
            self.default_srate = default_srate

        try:
            self.default_seglen = cp.getfloat('input', 'seglen')
        except:
            self.default_seglen = default_seglen

        try:
            self.default_window = cp.getfloat('input', 'window')
        except:
            self.default_window = default_window

        try:
            self.min_srate = cp.getfloat('input', 'min-srate')
        except:
            self.min_srate = min_srate

        try:
            self.max_srate = cp.getfloat('input', 'max-srate')
        except:
            self.max_srate = max_srate

        try:
            self.max_seglen = cp.getfloat('input', 'max-seglen')
        except:
            self.max_seglen = max_seglen

        try:
            self.min_seglen = cp.getfloat('input', 'min-seglen')
        except:
            self.min_seglen = min_seglen

        try:
            self.max_window = cp.getfloat('input', 'max-window')
        except:
            self.max_window = max_window

        try:
            self.min_window = cp.getfloat('input', 'min-window')
        except:
            self.min_window = min_window

        #
        # Add trigger properties
        #
        self.trigger_time = trigger_time
        self.time_lag = time_lag
        self.trigger_frequency = trigger_frequency

        # Variable sample rate / window length [fixed TF volume]

        # XXX new variables: seglen, segment-start, window

        if trigger_frequency is not None:
            # Adjust sample rate for this trigger
            # - min srate => max_seglen
            # - max srate => min_seglen
            if trigger_frequency < self.frequency_threshold:
               self.srate = self.min_srate
               self.seglen = self.max_seglen
               self.window = self.max_window
            else:
               self.srate = self.max_srate
               self.seglen = self.min_seglen
               self.window = self.min_window
        else:
            self.srate = self.default_srate
            self.seglen = self.default_seglen
            self.window = self.default_window

        self.rho = rho
        self.injevent = injevent

        self.veto1=veto1
        self.veto2=veto2

        self.BW_event=BW_event

        #
        # GraceDB Support
        #

        # If graceID is given, override other trigger values
        self.graceID = graceID
        if graceID is not None:
            self.query_graceDB(graceID)

    #
    # Update trigger properties
    #
    def set_injevent(self, injevent):
        self.injevent = injevent


    def query_graceDB(self,graceid):

        from ligo.gracedb.rest import GraceDb 

        # Instantiate graceDB event
        gracedb = GraceDb()
        event = gracedb.event(graceid)
        event_info = event.json()

        # Get loudness (for informational, not analysis, purposes)
        try:
            self.rho = event_info['extra_attributes']['MultiBurst']['snr']
        except KeyError:
            print >> sys.stderr, \
                    "graceDB UID %s has no MultiBurst snr attribute"%(graceid)

        # Set time
        self.trigger_time = event_info['gpstime']

        # Set frequency
        try:
            self.trigger_frequency = \
                    event_info['extra_attributes']['MultiBurst']['central_freq']

            if self.trigger_frequency < self.frequency_threshold:
               self.srate = self.min_srate
               self.seglen = self.max_seglen
               self.window = self.max_window
            else:
               self.srate = self.max_srate
               self.seglen = self.min_seglen
               self.window = self.min_window

        except KeyError:
            print >> sys.stderr, \
                    "graceDB UID %s has no MultiBurst central_freq attribute"%(graceid)
            print >> sys.stderr, "...using default sample rate"
            self.srate = self.default_srate




class triggerList:
    """
    Object to store trigger properties and associated configuration

    Allowed formats:
        trigger_gps 
        trigger_gps | time_lag
        trigger_gps | time_lag | trigger_frequency
        trigger_gps | time_lag | trigger_frequency | rho
    """

    def __init__(self, cp, gps_times=None, trigger_file=None,
            injection_file=None, followup_injections=None,
            cwb_trigger_file=None, rho_threshold=-1.0,
            internal_injections=False, graceIDs=None):

        #
        # Assign trigger data
        #
        if gps_times is not None and not internal_injections:
            # Create trigger list from gps times
            self.triggers=list()
            for gps_time in gps_times:
                self.triggers.append(eventTrigger(cp, trigger_time=gps_time))

        elif trigger_file is not None:
            # Create trigger list from ascii file
            self.triggers = self.parse_trigger_list(cp, trigger_file)

        elif injection_file is not None:
            # Create trigger list from sim* LIGOLW-XML table
            self.triggers = self.parse_injection_file(cp, injection_file,
                    followup_injections=followup_injections)

        elif cwb_trigger_file is not None:
            # Create trigger list from cwb triggers
            self.triggers = self.parse_cwb_trigger_list(cp, cwb_trigger_file)

        elif internal_injections:
            # Set up || runs to sample from the prior
            self.triggers = self.build_internal_injections(cp, gps_times)

        elif graceIDs is not None:
            # Create trigger list from graceDB queries
            self.triggers = self.parse_graceDB_triggers(cp, graceIDs)

        else:
            # Fail
            print >> sys.stdout, "don't know what to do."
            sys.exit()

    def parse_graceDB_triggers(self, cp, graceIDs):

        triggers=[]
        for graceid in graceIDs:
            triggers.append(eventTrigger(cp, graceID=graceid))

        return triggers

    def build_internal_injections(self, cp, gps_time):

        BW_Nsamples = cp.getint('bayeswave_options', 'BW-Nsamples')

        # Determine chain length
        injtype=cp.get('bayeswave_options', 'BW-inject')
        try:
            injname=cp.get('bayeswave_options', 'BW-injName')
            injname+='_'
        except:
            injname=''


        try:
            BW_chainLength=cp.getint('bayeswave_options','BW-chainLength')
        except ConfigParser.NoOptionError:

            print >> sys.stdout, "Reading chainlength from files in %s"%(
                    cp.get('bayeswave_options','BW-path'))

            # O1 names:
            if injtype=='glitch':
                filename=injname+'glitch_glitchchain_ifo0.dat.0'
            else: filename=injname+'signal_wavechain.dat.0' 
            filename=os.path.join(cp.get('bayeswave_options','BW-path'), filename)

            try:
                # O1 names
                o1=os.path.exists(filename)
                if not o1: 
                    raise ValueError(
                        "o1 style chain-names not found,trying o2-style")
            except:
                # O2 names:
                if injtype=='glitch':
                    filename=injname+'glitch_params_ifo0.dat.0'
                else: filename=injname+'signal_params.dat.0' 
                filename=os.path.join(cp.get('bayeswave_options','BW-path'), filename)

            BW_chainLength = file_len(filename)

        try:
            BW_seed = cp.getint('bayeswave_options', 'BW-seed')
        except:
            BW_seed = None

        random.seed(BW_seed)
        BW_events = random.sample(xrange(0,BW_chainLength), BW_Nsamples)

        triggers=[]
        for BW_event in BW_events:
            triggers.append(eventTrigger(cp, trigger_time=gps_time,
                BW_event=BW_event))

        return triggers



    def parse_injection_file(self, cp, injection_file, followup_injections=None,
            injwindow=2.0):


        xmldoc = ligolw_utils.load_filename(injection_file, contenthandler =
                LIGOLWContentHandler, verbose = True)
        sim_inspiral_table = lsctables.SimInspiralTable.get_table(xmldoc)

        injection_times = sim_inspiral_table.get_column('geocent_end_time') + \
                1e-9*sim_inspiral_table.get_column('geocent_end_time_ns')

        print "..read %d injections"%len(injection_times)

        triggers=[]
        if followup_injections is None:

            print 'downsampling to requested injections using events= in config'

            # reduce to specified values
            events=cp.get('injections', 'events')

            if events!='all':
                injevents=list(hyphen_range(events))
            else:
                injevents=range(len(injection_times))

            for i in injevents:
                triggers.append(eventTrigger(cp, trigger_time=injection_times[i],
                    injevent=i))

        else:

            # Parse the detected injections

            print "downsampling to events listed in %s"%followup_injections
            trigger_list_from_file = triggerList(cp,
                    trigger_file=followup_injections)

            # Find corresponding injection events
            for trigger in trigger_list_from_file.triggers:

                injevent = np.concatenate(np.argwhere(
                    abs(trigger.trigger_time - injection_times) < injwindow))[0]

                trigger.set_injevent(injevent)

                triggers.append(trigger)


        return triggers


    def parse_cwb_trigger_list(self, cp, cwb_trigger_file, rho_threshold=-1.0,
            keep_frac=1.0):

        # Get rho threshold
        try:
            rho_threshold = cp.getfloat('input', 'rho-threshold')
        except:
            rho_threshold = rho_threshold

        print >> sys.stdout, "Discarding rho<=%f"%rho_threshold
        
        names = ['veto1', 'veto2', 'rho', 'cc1', 'cc2', 'cc3', 'amp', 'tshift',
                'tsupershift', 'like', 'penalty', 'disbalance', 'f',
                'bandwidth', 'duration', 'pixels', 'resolution', 'runnumber',
                'Lgps', 'Hgps', 'sSNRL', 'sSNRH', 'hrssL', 'hrssH', 'phi',
                'theta', 'psi']

        data = np.recfromtxt(cwb_trigger_file,names=names)

        Hgps = data['Hgps']
        Lgps = data['Lgps']
        rhoList = data['rho']
        freqList = data['f']

        plusveto = data['veto1']
        minusveto = data['veto2']

        lagList = []

        for h,l in zip(Hgps,Lgps):
           lagList.append(round(h-l))

        gpsList = Hgps

        triggers=[]
        for gps, lag, freq, rho, veto1, veto2 in zip(gpsList, lagList, freqList,
                rhoList, plusveto, minusveto):
            # Apply rho threshold
            if rho < rho_threshold: continue

            triggers.append(eventTrigger(cp,
                trigger_time=gps,
                time_lag=lag,
                trigger_frequency=freq,
                rho=rho,
                veto1=veto1,
                veto2=veto2))

        # Finally, downsample to a smaller fraction of triggers
        try:
            keep_frac = cp.getfloat('input', 'keep-frac')
        except:
            keep_frac = keep_frac

        nall=len(triggers)
        nkeep=int(np.ceil(keep_frac*nall))
        keepidx=random.sample(range(0,len(triggers)), nkeep)
        triggers_out = [ triggers[i] for i in sorted(keepidx) ]

        print >> sys.stdout, "Read %d triggers, following up %d"%(
                nall, len(triggers_out))

        return triggers_out


    def parse_trigger_list(self, cp, trigger_file, rho_threshold=-1.0,
            keep_frac=1.0):

        trigger_data = np.loadtxt(trigger_file)
        try:
            nrows, ncols = trigger_data.shape
        except ValueError:
            nrows = len(trigger_data)
            ncols = 1

        triggers = list()

        if ncols==1:
            # Just have trigger time
            for i in xrange(nrows):
                triggers.append(eventTrigger(cp,
                    trigger_time=trigger_data[i]))

        elif ncols==2:
            # Trigger time, lag
            for i in xrange(nrows):
                triggers.append(eventTrigger(cp,
                    trigger_time=trigger_data[i,0],
                    time_lag=trigger_data[i,1]))

        elif ncols==3:
            # Trigger time, lag, frequency
            for i in xrange(nrows):
                triggers.append(eventTrigger(cp,
                    trigger_time=trigger_data[i,0],
                    time_lag=trigger_data[i,1],
                    trigger_frequency=trigger_data[i,2]))

        elif ncols==4:
            # Trigger time, lag, frequency, rho
            try:
                rho_threshold = cp.getfloat('input', 'rho-threshold')
            except:
                rho_threshold = rho_threshold

            print >> sys.stdout, "Discarding rho<=%f"%rho_threshold

            for i in xrange(nrows):
                # Apply rho threshold
                if trigger_data[i,3] < rho_threshold: continue
                triggers.append(eventTrigger(cp,
                    trigger_time=trigger_data[i,0],
                    time_lag=trigger_data[i,1],
                    trigger_frequency=trigger_data[i,2],
                    rho=trigger_data[i,3]))


        # Finally, downsample to a smaller fraction of triggers
        try:
            keep_frac = cp.getfloat('input', 'keep-frac')
        except:
            keep_frac = keep_frac

        nkeep=int(np.ceil(keep_frac*len(triggers)))
        keepidx=random.sample(range(0,len(triggers)), nkeep)
        triggers_out = [ triggers[i] for i in sorted(keepidx) ]

        print >> sys.stdout, "Read %d triggers, following up %d"%(
                nrows, len(triggers_out))

        return triggers_out



    # -- END trigger_list class



#
# Condor Definitions
#

class bayeswaveJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, cacheFiles, injfile=None, nrdata=None, dax=False):

        universe=cp.get('condor','universe')

        bayeswave=os.path.join(os.environ['BAYESWAVE_PREFIX','src/bayeswave')

        pipeline.CondorDAGJob.__init__(self,universe,bayeswave)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)

        if cp.has_option('condor', 'accounting_group'):
            self.add_condor_cmd('accounting_group', cp.get('condor', 'accounting_group'))   

        self.set_stdout_file('$(macrooutputDir)/bayeswave_$(cluster)-$(process)-$(node).out')
        self.set_stderr_file('$(macrooutputDir)/bayeswave_$(cluster)-$(process)-$(node).err')
        self.set_log_file('$(macrooutputDir)/bayeswave_$(cluster)-$(process)-$(node).log')

        # --- Allow desired sites
        if cp.has_option('condor','desired-sites'):
            self.add_condor_cmd('+DESIRED_Sites',cp.get('condor','desired-sites'))

        if cp.has_option('condor', 'bayeswave-request-memory'):
            self.add_condor_cmd('request_memory',
                    cp.get('condor', 'bayeswave-request-memory'))   
        #
        # Identify osg vs ldg site
        #
        if cp.getboolean('condor', 'osg-jobs'):
            print >> sys.stdout, "Configuring file transfers for OSG deployment"


            # --- Perform file transfers
            self.add_condor_cmd('should_transfer_files', 'YES')
            self.add_condor_cmd('when_to_transfer_output', 'ON_EXIT_OR_EVICT')
            self.add_condor_cmd('transfer_output_files', '$(macrooutputDir)')

            # --- Files to include in transfer
            transferstring='datafind,$(macrooutputDir)'

            if cp.has_option('condor','transfer-files'):
                # allow specification of additional files to transfer
                transferstring+=',%s'%cp.get('condor','transfer-files')

            if cp.getboolean('condor','copy-frames'): transferstring+=',$(macroframes)'

            self.add_condor_cmd('transfer_input_files', transferstring)

        self.add_condor_cmd('getenv', 'True')


        # ----------------------------------------------------------------------------------
        # --- Required options  ------------------------------------------------------------
        # ----------------------------------------------------------------------------------
        ifo_list = ast.literal_eval(cp.get('input', 'ifo-list'))
        if not cp.getboolean('datafind','sim-data'):
            channel_list = ast.literal_eval(cp.get('datafind', 'channel-list'))

        # XXX: hack to repeat option for --ifo H1 --ifo L1 etc
        ifo_list_opt = ifo_list[0]
        for ifo in ifo_list[1:]: ifo_list_opt += ' --ifo {0}'.format(ifo)
        self.add_opt('ifo', ifo_list_opt)

#        self.add_opt('srate', cp.get('input', 'srate'))
#        self.add_opt('seglen', cp.get('input', 'seglen'))
        self.add_opt('PSDlength', cp.get('input', 'PSDlength'))
 
        flow = ast.literal_eval(cp.get('input','flow'))
        for ifo in ifo_list:
            self.add_opt('{ifo}-flow'.format(ifo=ifo), str(flow[ifo]))
            self.add_opt('{ifo}-cache'.format(ifo=ifo), cacheFiles[ifo])

            if not cp.getboolean('datafind','sim-data'):
                # only specify channels for real data
                self.add_opt('{ifo}-channel'.format(ifo=ifo), channel_list[ifo])


        # dataseed
        if cp.has_option('bayeswave_options', 'dataseed'):
            self.add_opt('dataseed', cp.get('bayeswave_options', 'dataseed'))

        # ----------------------------------------------------------------------------------
        # --- Run parameters   -------------------------------------------------------------
        # ----------------------------------------------------------------------------------

        # segment-start
        if cp.has_option('bayeswave_options', 'segment-start'):
            self.add_opt('segment-start', cp.get('bayeswave_options',
                'segment-start'))

        # Niter
        if cp.has_option('bayeswave_options', 'Niter'):
            self.add_opt('Niter', cp.get('bayeswave_options', 'Niter'))

        # Nchain
        if cp.has_option('bayeswave_options', 'Nchain'):
            self.add_opt('Nchain', cp.get('bayeswave_options', 'Nchain'))
            
        # Ncycle
        if cp.has_option('bayeswave_options', 'Ncycle'):
            self.add_opt('Ncycle', cp.get('bayeswave_options', 'Ncycle'))

        # Nburnin
        if cp.has_option('bayeswave_options', 'Nburnin'):
            self.add_opt('Nburnin', cp.get('bayeswave_options', 'Nburnin'))

        # maxLogL
        if cp.has_option('bayeswave_options', 'maxLogL'):
            self.add_opt('maxLogL', cp.get('bayeswave_options', 'maxLogL'))

        # chainseed
        if cp.has_option('bayeswave_options', 'chainseed'):
            self.add_opt('chainseed', cp.get('bayeswave_options', 'chainseed'))

        # runName
        if cp.has_option('bayeswave_options', 'runName'):
            self.add_opt('runName', cp.get('bayeswave_options', 'runName'))

        # 0noise
        if cp.has_option('bayeswave_options', '0noise'):
            self.add_opt('0noise', cp.get('bayeswave_options', '0noise'))

        # zeroLogL
        if cp.has_option('bayeswave_options', 'zeroLogL'):
            self.add_opt('zeroLogL', cp.get('bayeswave_options', 'zeroLogL'))

        # restart
        if cp.has_option('bayeswave_options', 'restart'):
            self.add_opt('restart', cp.get('bayeswave_options', 'restart'))

        # gnuplot
        if cp.has_option('bayeswave_options', 'gnuplot'):
            self.add_opt('gnuplot', cp.get('bayeswave_options', 'gnuplot'))

        # verbose
        if cp.has_option('bayeswave_options', 'verbose'):
            self.add_opt('verbose', cp.get('bayeswave_options', 'verbose'))

        # window
        #if cp.has_option('bayeswave_options', 'window'):
        #    self.add_opt('window', cp.get('bayeswave_options', 'window'))

        # self-checkpointing
        if cp.has_option('condor', 'checkpoint'):
            self.add_opt('checkpoint', cp.get('condor', 'checkpoint'))

        # version
        if cp.has_option('bayeswave_options', 'version'):
            self.add_opt('version', cp.get('bayeswave_options', 'version'))

        # help
        if cp.has_option('bayeswave_options', 'help'):
            self.add_opt('version', cp.get('bayeswave_options', 'help'))

        # ----------------------------------------------------------------------------------
        # --- Run parameters   -------------------------------------------------------------
        # ----------------------------------------------------------------------------------

        # fullOnly
        if cp.has_option('bayeswave_options', 'fullOnly'):
            self.add_opt('fullOnly', cp.get('bayeswave_options', 'fullOnly'))

        # noClean
        if cp.has_option('bayeswave_options', 'noClean'):
            self.add_opt('noClean', cp.get('bayeswave_options', 'noClean'))

        # noSignal
        if cp.has_option('bayeswave_options', 'noSignal'):
            self.add_opt('noSignal', cp.get('bayeswave_options', 'noSignal'))

        # cleanOnly
        if cp.has_option('bayeswave_options', 'cleanOnly'):
            self.add_opt('cleanOnly', cp.get('bayeswave_options', 'cleanOnly'))

        # noiseOnly
        if cp.has_option('bayeswave_options', 'noiseOnly'):
            self.add_opt('noiseOnly', cp.get('bayeswave_options', 'noiseOnly'))

        # signalOnly
        if cp.has_option('bayeswave_options', 'signalOnly'):
            self.add_opt('signalOnly', cp.get('bayeswave_options', 'signalOnly'))

        # glitchOnly
        if cp.has_option('bayeswave_options', 'glitchOnly'):
            self.add_opt('glitchOnly', cp.get('bayeswave_options', 'glitchOnly'))

        # noPSDfit
        if cp.has_option('bayeswave_options', 'noPSDfit'):
            self.add_opt('noPSDfit', cp.get('bayeswave_options', 'noPSDfit'))

        # bayesLine
        if cp.has_option('bayeswave_options', 'bayesLine'):
            self.add_opt('bayesLine', cp.get('bayeswave_options', 'bayesLine'))

        # stochastic
        if cp.has_option('bayeswave_options', 'stochastic'):
            self.add_opt('stochastic', cp.get('bayeswave_options', 'stochastic'))

        # ----------------------------------------------------------------------------------
        # --- Priors & Proposasl  ----------------------------------------------------------
        # ----------------------------------------------------------------------------------

        # Dimensions
        if cp.has_option('bayeswave_options', 'Dmin'):
            self.add_opt('Dmin', cp.get('bayeswave_options', 'Dmin'))
        if cp.has_option('bayeswave_options', 'Dmax'):
            self.add_opt('Dmax', cp.get('bayeswave_options', 'Dmax'))

        # fixD
        if cp.has_option('bayeswave_options', 'fixD'):
            self.add_opt('fixD', cp.get('bayeswave_options', 'fixD'))

        # Quality factor
        if cp.has_option('bayeswave_options', 'Qmin'):
            self.add_opt('Qmin', cp.get('bayeswave_options', 'Qmin'))
        if cp.has_option('bayeswave_options', 'Qmax'):
           self.add_opt('Qmax', cp.get('bayeswave_options', 'Qmax'))

        # waveletPrior
        if cp.has_option('bayeswave_options', 'waveletPrior'):
             self.add_opt('waveletPrior', cp.get('bayeswave_options',
                 'waveletPrior'))

        # clusterPrior
        if cp.has_option('bayeswave_options', 'clusterPrior'):
             self.add_opt('clusterPrior', cp.get('bayeswave_options',
                 'clusterPrior'))

        # clusterPath
        if cp.has_option('bayeswave_options', 'clusterPath'):
             self.add_opt('clusterPath', cp.get('bayeswave_options',
                 'clusterPath'))

        # clusterAlpha
        if cp.has_option('bayeswave_options', 'clusterAlpha'):
             self.add_opt('clusterAlpha', cp.get('bayeswave_options',
                 'clusterAlpha'))
        
        # clusterBeta
        if cp.has_option('bayeswave_options', 'clusterBeta'):
             self.add_opt('clusterBeta', cp.get('bayeswave_options',
                 'clusterBeta'))

        # clusterGamma
        if cp.has_option('bayeswave_options', 'clusterGamma'):
             self.add_opt('clusterGamma', cp.get('bayeswave_options',
                 'clusterGamma'))

        # backgroundPrior
        if cp.has_option('bayeswave_options', 'backgroundPrior'):
            if cp.get('bayeswave_options', 'backgroundPrior') == None:
                print >> sys.stderr, "must specifiy name of 2-column bkg frequency distribution file"
                sys.exit()
            self.add_opt('backgroundPrior', cp.get('bayeswave_options','backgroundPrior'))


        # noOrientationProposal
        if cp.has_option('bayeswave_options', 'noOrientationProposal'):
             self.add_opt('noOrientationProposal', cp.get('bayeswave_options',
                 'noOrientationProposal'))


        # uniformAmplitudePrior 
        if cp.has_option('bayeswave_options', 'uniformAmplitudePrior'):
             self.add_opt('uniformAmplitudePrior', cp.get('bayeswave_options',
                 'uniformAmplitudePrior'))

        # noSignalAmplitudePrior
        if cp.has_option('bayeswave_options', 'noSignalAmplitudePrior'):
             self.add_opt('noSignalAmplitudePrior', cp.get('bayeswave_options',
                 'noSignalAmplitudePrior'))

        # noAmplitudeProposal
        if cp.has_option('bayeswave_options', 'noAmplitudeProposal'):
             self.add_opt('noAmplitudeProposal', cp.get('bayeswave_options',
                 'noAmplitudeProposal'))

        # varyExtrinsicAmplitude
        if cp.has_option('bayeswave_options', 'varyExtrinsicAmplitude'):
             self.add_opt('varyExtrinsicAmplitude', cp.get('bayeswave_options',
                 'varyExtrinsicAmplitude'))

        # noClusterProposal
        if cp.has_option('bayeswave_options', 'noClusterProposal'):
             self.add_opt('noClusterProposal', cp.get('bayeswave_options',
                 'noClusterProposal'))

        # clusterWeight
        if cp.has_option('bayeswave_options', 'clusterWeight'):
             self.add_opt('clusterWeight', cp.get('bayeswave_options', 'clusterWeight'))

        # ampPriorPeak
        if cp.has_option('bayeswave_options', 'ampPriorPeak'):
             self.add_opt('ampPriorPeak', cp.get('bayeswave_options', 'ampPriorPeak'))

        # signalPriorPeak
        if cp.has_option('bayeswave_options', 'signalPriorPeak'):
             self.add_opt('signalPriorPeak', cp.get('bayeswave_options', 'signalPriorPeak'))

        # dimensionDecayRate
        if cp.has_option('bayeswave_options', 'dimensionDecayRate'):
             self.add_opt('dimensionDecayRate', cp.get('bayeswave_options', 'dimensionDecayRate'))

        # fixIntrinsicParams
        if cp.has_option('bayeswave_options', 'fixIntrinsicParams'):
             self.add_opt('fixIntrinsicParams', cp.get('bayeswave_options', 'fixIntrinsicParams'))

        # fixExtrinsicParams
        if cp.has_option('bayeswave_options', 'fixExtrinsicParams'):
             self.add_opt('fixExtrinsicParams', cp.get('bayeswave_options', 'fixExtrinsicParams'))

        # ----------------------------------------------------------------------------------
        # --- Parallel Tempering parameters  -----------------------------------------------
        # ----------------------------------------------------------------------------------

        # tempMin
        if cp.has_option('bayeswave_options', 'tempMin'):
             self.add_opt('tempMin', cp.get('bayeswave_options', 'tempMin'))

        # noAdaptTemperature
        if cp.has_option('bayeswave_options', 'noAdaptTemperature'):
             self.add_opt('noAdaptTemperature', cp.get('bayeswave_options', 'noAdaptTemperature'))
             k
        # tempSpacing
        if cp.has_option('bayeswave_options', 'tempSpacing'):
             self.add_opt('tempSpacing', cp.get('bayeswave_options', 'tempSpacing'))

        # noSplineEvidence
        if cp.has_option('bayeswave_options', 'noSplineEvidence'):
             self.add_opt('noSplineEvidence', cp.get('bayeswave_options', 'noSplineEvidence'))

        # ----------------------------------------------------------------------------------
        # --- LALInference Injection Options  ----------------------------------------------
        # ----------------------------------------------------------------------------------

        # Injection file
        if injfile is not None:
            injfile=os.path.join('..',injfile)
            self.add_opt('inj', injfile)

        # NR file
        if nrdata is not None:
            nrdata=os.path.join('..',nrdata)
            self.add_opt('inj-numreldata', nrdata)

        # ----------------------------------------------------------------------------------
        # --- Burst MDC injection ----------------------------------------------------------
        # ----------------------------------------------------------------------------------

        # mdc-cache
        if cp.has_option('injections', 'mdc-cache'):
            mdc_cache_list=str(['../datafind/MDC.cache' for ifo in
                ifo_list]).replace("'",'')
            mdc_cache_list=mdc_cache_list.replace(' ','')
            self.add_opt('MDC-cache', mdc_cache_list)

        # mdc-channels
        if cp.has_option('injections', 'mdc-channels'):
            mdc_channel_list=ast.literal_eval(cp.get('injections','mdc-channels'))
            mdc_channel_str=str(mdc_channel_list.values()).replace("'",'')
            mdc_channel_str=mdc_channel_str.replace(' ','')
            self.add_opt('MDC-channel', mdc_channel_str)

        # mdc-prefactor
        if cp.has_option('injections', 'mdc-prefactor'):
            self.add_opt('MDC-prefactor', cp.get('injections', 'mdc-prefactor'))

        # ----------------------------------------------------------------------------------
        # --- BayesWave internal injection options -----------------------------------------
        # ----------------------------------------------------------------------------------
        # BW-inject
        if cp.has_option('bayeswave_options', 'BW-inject'):
             self.add_opt('BW-inject', cp.get('bayeswave_options', 'BW-inject'))

        # BW-injName
        if cp.has_option('bayeswave_options', 'BW-injName'):
             self.add_opt('BW-injName', cp.get('bayeswave_options', 'BW-injName'))

        # BW-path
        if cp.has_option('bayeswave_options', 'BW-path'):
             self.add_opt('BW-path', cp.get('bayeswave_options', 'BW-path'))


        # XXX: where is this?
        # NC
        if cp.has_option('bayeswave_options','NC'):
            self.add_opt('NC', cp.get('bayeswave_options', 'NC'))
        # NCmin
        if cp.has_option('bayeswave_options','NCmin'):
            self.add_opt('NCmin', cp.get('bayeswave_options', 'NCmin'))
        # NCmax
        if cp.has_option('bayeswave_options','NCmax'):
            self.add_opt('NCmax', cp.get('bayeswave_options', 'NCmax'))


        self.set_sub_file('bayeswave.sub')

class bayeswaveNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

    def __init__(self, bayeswave_job):

        pipeline.CondorDAGNode.__init__(self,bayeswave_job)
        pipeline.AnalysisNode.__init__(self)

    def set_trigtime(self, trigtime):
        self.add_var_opt('trigtime', '%.9f'%trigtime)
        self.__trigtime = trigtime

    def set_segment_start(self, segment_start):
        self.add_var_opt('segment-start', '%.9f'%segment_start)
        self.__segment_start = segment_start

    def set_srate(self, srate):
        self.add_var_opt('srate', srate)
        self.__srate = srate

    def set_seglen(self, seglen):
        self.add_var_opt('seglen', seglen)
        self.__seglen = seglen

    def set_window(self, window):
        self.add_var_opt('window', window)
        self.__window = window

    def set_rolloff(self, rolloff):
        self.add_var_opt('padding', rolloff)
        self.__rolloff = rolloff

    def set_PSDstart(self, PSDstart):
        self.add_var_opt('PSDstart', '%.9f'%PSDstart)
        self.__PSDstart = PSDstart

    def set_outputDir(self, outputDir):
        self.add_var_opt('outputDir', outputDir)
        self.__outputDir = outputDir

    def set_injevent(self, event):
        self.add_var_opt('event', event)
        self.__event = event

    def set_dataseed(self, dataseed):
        self.add_var_opt('dataseed', dataseed)
        self.__dataseed = dataseed

    def add_frame_transfer(self, framedict):
        """
        Add a list of frames to transfer
        """
        self.__frames=""
        for ifo in framedict.keys():
            for frame in framedict[ifo]:
                self.__frames += frame + ','
        self.__frames.strip(',')
        self.add_var_opt('frames', self.__frames)
  
    def set_L1_timeslide(self, L1_timeslide):
        self.add_var_opt('L1-timeslide', L1_timeslide)
        self.__L1_timeslide = L1_timeslide

    def set_BW_event(self, BW_event):
        self.add_var_opt('BW-event', BW_event)
        self.__BW_event = BW_event

#
# Post-processing
#

class bayeswave_postJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, cacheFiles, injfile=None, nrdata=None, dax=False):

        universe=cp.get('condor','universe')

        bayeswave_post=os.path.join(os.environ['BAYESWAVE_PREFIX','src/bayeswave_post')

        pipeline.CondorDAGJob.__init__(self,universe,bayeswave_post)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)

        if cp.has_option('condor', 'accounting_group'):
            self.add_condor_cmd('accounting_group', cp.get('condor', 'accounting_group'))   

        if cp.has_option('condor', 'bayeswave_post-request-memory'):
            self.add_condor_cmd('request_memory',
                    cp.get('condor', 'bayeswave_post-request-memory'))   

        self.set_stdout_file('$(macrooutputDir)/bayeswave_post_$(cluster)-$(process)-$(node).out')
        self.set_stderr_file('$(macrooutputDir)/bayeswave_post_$(cluster)-$(process)-$(node).err')
        self.set_log_file('$(macrooutputDir)/bayeswave_post_$(cluster)-$(process)-$(node).log')

        # Request 4GB of RAM for pp jobs
        #self.add_condor_cmd('request_memory', '4000')

        # --- Allow desired sites
        if cp.has_option('condor','desired-sites'):
            self.add_condor_cmd('+DESIRED_Sites',cp.get('condor','desired-sites'))

        #
        # Identify osg vs ldg site
        #
        if cp.getboolean('condor', 'osg-jobs'):
            print >> sys.stdout, "Configuring file transfers for OSG deployment"


            # --- Perform file transfers
            self.add_condor_cmd('should_transfer_files', 'YES')
            self.add_condor_cmd('when_to_transfer_output', 'ON_EXIT_OR_EVICT')
            self.add_condor_cmd('transfer_output_files', '$(macrooutputDir)')

            # --- Files to include in transfer
            # FIXME: PostProc doesn't currently need frame transfer
            transferstring='datafind,$(macrooutputDir)'

            if cp.has_option('condor','transfer-files'):
                # allow specification of additional files to transfer
                transferstring+=',%s'%cp.get('condor','transfer-files')

            self.add_condor_cmd('transfer_input_files', transferstring)

        self.add_condor_cmd('getenv', 'True')

        # --- Required options
        ifo_list = ast.literal_eval(cp.get('input', 'ifo-list'))
        if not cp.get('datafind','sim-data'):
            channel_list = ast.literal_eval(cp.get('datafind', 'channel-list'))

        # XXX: hack to repeat option
        ifo_list_opt = ifo_list[0]
        for ifo in ifo_list[1:]:
            ifo_list_opt += ' --ifo {0}'.format(ifo)
        self.add_opt('ifo', ifo_list_opt)

        #self.add_opt('srate', cp.get('input', 'srate'))
        #self.add_opt('seglen', cp.get('input', 'seglen'))
        self.add_opt('PSDlength', cp.get('input', 'PSDlength'))
 
        flow = ast.literal_eval(cp.get('input','flow'))

        for i,ifo in enumerate(ifo_list):
            self.add_opt('{ifo}-flow'.format(ifo=ifo), str(flow[ifo]))


        # --- Optional options
        # bayesLine
        if cp.has_option('bayeswave_post_options', 'bayesLine'):
            self.add_opt('bayesLine', cp.get('bayeswave_post_options', 'bayesLine'))

        # 0noise
        if cp.has_option('bayeswave_post_options', '0noise'):
            self.add_opt('0noise', cp.get('bayeswave_post_options', '0noise'))

        #
        # Injection file
        #
        if injfile is not None:
            # XXX: note that bayeswave works within the outputDir, so point to
            # injection
            injfile=os.path.join('..',injfile)
            self.add_opt('inj', injfile)

        if nrdata is not None:
            nrdata=os.path.join('..',nrdata)
            self.add_opt('inj-numreldata', nrdata)

        #
        # MDC Setup
        #
        if cp.has_option('injections', 'mdc-cache'):
            mdc_cache_list=str(['../datafind/MDC.cache' for ifo in
                ifo_list]).replace("'",'')
            mdc_cache_list=mdc_cache_list.replace(' ','')
            self.add_opt('MDC-cache', mdc_cache_list)

        if cp.has_option('injections', 'mdc-channels'):
            #mdc_channel_list=ast.literal_eval(cp.get('injections','mdc-channels'))
            mdc_channel_list=ast.literal_eval(cp.get('injections','mdc-channels'))
            mdc_channel_str=str(mdc_channel_list.values()).replace("'",'')
            mdc_channel_str=mdc_channel_str.replace(' ','')
            self.add_opt('MDC-channel', mdc_channel_str)

        if cp.has_option('injections', 'mdc-prefactor'):
            self.add_opt('MDC-prefactor', cp.get('injections', 'mdc-prefactor'))

        # BW-inject
        if cp.has_option('bayeswave_options', 'BW-inject'):
             self.add_opt('BW-inject', cp.get('bayeswave_options', 'BW-inject'))

        # BW-injName
        if cp.has_option('bayeswave_options', 'BW-injName'):
             self.add_opt('BW-injName', cp.get('bayeswave_options', 'BW-injName'))

        # BW-path
        if cp.has_option('bayeswave_options', 'BW-path'):
             self.add_opt('BW-path', cp.get('bayeswave_options', 'BW-path'))

        self.set_sub_file('bayeswave_post.sub')


class bayeswave_postNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

    def __init__(self, bayeswave_post_job):

        pipeline.CondorDAGNode.__init__(self, bayeswave_post_job)
        pipeline.AnalysisNode.__init__(self)

    def set_trigtime(self, trigtime):
        self.add_var_opt('trigtime', '%.9f'%trigtime)
        self.__trigtime = trigtime

    def set_segment_start(self, segment_start):
        self.add_var_opt('segment-start', '%.9f'%segment_start)
        self.__segment_start = segment_start

    def set_window(self, window):
        self.add_var_opt('window', window)
        self.__window = window

    def set_rolloff(self, rolloff):
        self.add_var_opt('padding', rolloff)
        self.__rolloff = rolloff

    def set_srate(self, srate):
        self.add_var_opt('srate', srate)
        self.__srate = srate

    def set_seglen(self, seglen):
        self.add_var_opt('seglen', seglen)
        self.__seglen = seglen

    def set_PSDstart(self, PSDstart):
        self.add_var_opt('PSDstart', '%.9f'%PSDstart)
        self.__PSDstart = PSDstart

    def set_outputDir(self, ifo_list, outputDir):
        self.add_var_opt('outputDir', outputDir)
        self.__outputDir = outputDir

        # bayeswave_post now uses PSD estimates straight from bayeswave and
        # no channel name needed.  These estimates lie in the outputDir so add
        # the variable option here
        for i,ifo in enumerate(ifo_list):
            self.add_var_opt('{ifo}-cache'.format(ifo=ifo),
                    "interp:{outputDir}/IFO{i}_asd.dat".format(i=i,
                        outputDir=outputDir))

    def set_injevent(self, event):
        self.add_var_opt('event', event)
        self.__event = event

    def set_dataseed(self, dataseed):
        self.add_var_opt('dataseed', dataseed)
        self.__dataseed = dataseed

    def set_L1_timeslide(self, L1_timeslide):
        self.add_var_opt('L1-timeslide', L1_timeslide)
        self.__L1_timeslide = L1_timeslide

    def set_BW_event(self, BW_event):
        self.add_var_opt('BW-event', BW_event)
        self.__BW_event = BW_event


#
# skymap job
#

class megaskyJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, dax=False):

        # XXX consider local universe?
        universe='vanilla'

        # Point this to the src dir
        megasky = os.path.join(os.environ['BAYESWAVE_PREFIX','postprocess/skymap/megasky.py')

        pipeline.CondorDAGJob.__init__(self,universe,megasky)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)

        if cp.has_option('condor', 'accounting_group'):
            self.add_condor_cmd('accounting_group', cp.get('condor', 'accounting_group'))   

        self.add_condor_cmd('getenv', 'True')

        #
        # Identify osg vs ldg site
        #
        if cp.getboolean('condor', 'osg-jobs'):
            print >> sys.stdout, "Configuring file transfers for OSG deployment"

            # --- Perform file transfers
            self.add_condor_cmd('should_transfer_files', 'YES')
            self.add_condor_cmd('when_to_transfer_output', 'ON_EXIT_OR_EVICT')
            self.add_condor_cmd('transfer_input_files', '$(macroargument0)')
            self.add_condor_cmd('transfer_output_files', '$(macroargument0)')

        self.set_stdout_file('$(macroargument0)/megasky_$(cluster)-$(process)-$(node).out')
        self.set_stderr_file('$(macroargument0)/megasky_$(cluster)-$(process)-$(node).err')
        self.set_log_file('$(macroargument0)/megasky_$(cluster)-$(process)-$(node).log')
        self.set_sub_file('megasky.sub')

        hostname = socket.gethostname()

        # --- Allow desired sites
        if cp.has_option('condor','desired-sites'):
            self.add_condor_cmd('+DESIRED_Sites',cp.get('condor','desired-sites'))


class megaskyNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

    def __init__(self, megasky_job, rundir):

        pipeline.CondorDAGNode.__init__(self, megasky_job)
        pipeline.AnalysisNode.__init__(self)

    # Set work dir
    def set_outputDir(self, outputDir):
        self.add_var_arg(outputDir)
        self.__outputDir = outputDir

#
# megaplot job
#

class megaplotJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, dax=False):

        universe='vanilla'

        # Point this to the src dir
        megaplot = os.path.join(os.environ['BAYESWAVE_PREFIX','postprocess/megaplot.py')

        pipeline.CondorDAGJob.__init__(self,universe, megaplot)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)

        # --- Allow desired sites
        if cp.has_option('condor','desired-sites'):
            self.add_condor_cmd('+DESIRED_Sites',cp.get('condor','desired-sites'))

        if cp.has_option('condor', 'accounting_group'):
            self.add_condor_cmd('accounting_group', cp.get('condor', 'accounting_group'))   

        self.add_condor_cmd('getenv', 'True')

        #
        # Identify osg vs ldg site
        #
        if cp.getboolean('condor', 'osg-jobs'):
            print >> sys.stdout, "Configuring file transfers for OSG deployment"

            css=os.path.join(os.environ['BAYESWAVE_PREFIX'], 'postprocess', 'BWBweb.css')
            ajax=os.path.join(os.environ['BAYESWAVE_PREFIX'], 'postprocess', 'secure_ajax.js')
            nav=os.path.join(os.environ['BAYESWAVE_PREFIX'], 'postprocess', 'navigate.js')

            # --- Perform file transfers
            self.add_condor_cmd('should_transfer_files', 'YES')
            self.add_condor_cmd('when_to_transfer_output', 'ON_EXIT_OR_EVICT')
            self.add_condor_cmd('transfer_input_files',
                '$(macroargument0),{css},{ajax},{nav}'.format(css=css,ajax=ajax,nav=nav))
            self.add_condor_cmd('transfer_output_files', '$(macroargument0)')

        self.set_stdout_file('$(macroargument0)/megaplot_$(cluster)-$(process)-$(node).out')
        self.set_stderr_file('$(macroargument0)/megaplot_$(cluster)-$(process)-$(node).err')
        self.set_log_file('$(macroargument0)/megaplot_$(cluster)-$(process)-$(node).log')
        self.set_sub_file('megaplot.sub')


class megaplotNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

    def __init__(self, megaplot_job, rundir):

        pipeline.CondorDAGNode.__init__(self, megaplot_job)
        pipeline.AnalysisNode.__init__(self)

    # Set work dir
    def set_outputDir(self, outputDir):
        self.add_var_arg(outputDir)
        self.__outputDir = outputDir

#
# submitGraceDB
#

class submitToGraceDB(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, dax=False):

        universe='vanilla'

        # Point this to the src dir
        gdb_submitter = cp.get('bayeswave_paths','gdb-submitter')
        pipeline.CondorDAGJob.__init__(self,universe, gdb_submitter)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)

        hostname = socket.gethostname()

        # --- Allow desired sites
        if cp.has_option('condor','desired-sites'):
            self.add_condor_cmd('+DESIRED_Sites',cp.get('condor','desired-sites'))

        if cp.has_option('condor', 'accounting_group'):
            self.add_condor_cmd('accounting_group', cp.get('condor', 'accounting_group'))   

        self.add_condor_cmd('getenv', 'True')

        self.set_stdout_file('gdb_submitter_$(cluster)-$(process)-$(node).out')
        self.set_stderr_file('gdb_submitter_$(cluster)-$(process)-$(node).err')
        self.set_log_file('gdb_submitter_$(cluster)-$(process)-$(node).log')
        self.set_sub_file('gdb_submitter.sub')


class submitToGraceDBNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

    def __init__(self, gdb_submitter_job, rundir, htmlDir):

        pipeline.CondorDAGNode.__init__(self, gdb_submitter_job)
        pipeline.AnalysisNode.__init__(self)
        # Set job initialdir, so python codes know where to expect input files
        self.add_var_condor_cmd('initialdir', rundir)   
        self.__rundir = rundir

        # Set html directory
        self.add_var_opt('htmlDir', htmlDir)
        self.__htmlDir = htmlDir


#
# Housekeeping
#

class archiverJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, dax=False):

        universe='vanilla'

        # Point this to the src dir
        archiver = cp.get('bayeswave_paths','archiver')
        pipeline.CondorDAGJob.__init__(self,universe, archiver)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)

        # --- Allow desired sites
        if cp.has_option('condor','desired-sites'):
            self.add_condor_cmd('+DESIRED_Sites',cp.get('condor','desired-sites'))

        if cp.has_option('condor', 'accounting_group'):
            self.add_condor_cmd('accounting_group', cp.get('condor', 'accounting_group'))   
        #
        # Identify osg vs ldg site
        #
        # FIXME: currently only associates PACE (GaTech) as an OSG site
        hostname = socket.gethostname()
        if 'pace.gatech.edu' in hostname:
            print >> sys.stdout, "Looks like you're on PACE; configuring file transfers"

            # --- Perform file transfers
            self.add_condor_cmd('should_transfer_files', 'YES')
            self.add_condor_cmd('when_to_transfer_output', 'ON_EXIT_OR_EVICT')
            self.add_condor_cmd('transfer_input_files', '$(macroargument0)')
            self.add_condor_cmd('transfer_output_files',
                    '$(macroargument0).tar.bz2')

        self.add_condor_cmd('getenv', 'True')

        self.set_stdout_file('$(macroargument0)/archiver_$(cluster)-$(process)-$(node).out')
        self.set_stderr_file('$(macroargument0)/archiver_$(cluster)-$(process)-$(node).err')
        self.set_log_file('archiver_$(cluster)-$(process)-$(node).log')
        self.set_sub_file('archiver.sub')


class archiverNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

    def __init__(self, archiver_job, rundir):

        pipeline.CondorDAGNode.__init__(self, archiver_job)
        pipeline.AnalysisNode.__init__(self)

        # Set run directory
        self.add_var_arg(rundir)
        self.__rundir = rundir


