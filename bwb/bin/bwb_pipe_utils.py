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

from glue import pipeline
import itertools
import socket
import sys,os
import ast

#
# Main analysis
#

class bayeswaveJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, cacheFiles, injfile=None, nrdata=None, dax=False):

        universe=cp.get('condor','universe')

        bayeswave_exec = cp.get('bayeswave_paths','bayeswave_executable')

        pipeline.CondorDAGJob.__init__(self,universe,bayeswave_exec)
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)

        if cp.has_option('condor', 'accounting_group'):
            self.add_condor_cmd('accounting_group', cp.get('condor', 'accounting_group'))   

        self.set_stdout_file('$(macrooutputDir)/bayeswave_$(cluster)-$(process)-$(node).out')
        self.set_stderr_file('$(macrooutputDir)/bayeswave_$(cluster)-$(process)-$(node).err')
        self.set_log_file('$(macrooutputDir)/bayeswave_$(cluster)-$(process)-$(node).log')

        #
        # Identify osg vs ldg site
        #
        # FIXME: currently only associates PACE (GaTech) as an OSG site
        hostname = socket.gethostname()
        #hostname = 'pace.gatech.edu'
        if 'pace.gatech.edu' in hostname:
            print >> sys.stdout, "Looks like you're on PACE; configuring file transfers"

            # --- Perform file transfers
            self.add_condor_cmd('should_transfer_files', 'YES')
            self.add_condor_cmd('when_to_transfer_output', 'ON_EXIT_OR_EVICT')
            self.add_condor_cmd('transfer_output_files', '$(macrooutputDir)')

            # --- Files to include in transfer
            transferstring='datafind,$(macrooutputDir)'

            if cp.getboolean('condor','copy-frames'): transferstring+=',$(macroframes)'

            if injfile is not None:
                transferstring+=','+'SEOBNRv2ChirpTimeSS.dat,'+injfile
            if nrdata is not None: transferstring+=','+nrdata
            self.add_condor_cmd('transfer_input_files', transferstring)

            # --- Point to ROM data (which should have been copied
            if injfile is not None:
                self.add_condor_cmd("environment", "LAL_DATA_PATH=./")


        self.add_condor_cmd('getenv', 'True')

        # --- Required options
        ifo_list = ast.literal_eval(cp.get('input', 'ifo-list'))
        channel_list = ast.literal_eval(cp.get('datafind', 'channel-list'))

        # XXX: hack to repeat option for --ifo H1 --ifo L1 etc
        ifo_list_opt = ifo_list[0]
        for ifo in ifo_list[1:]: ifo_list_opt += ' --ifo {0}'.format(ifo)
        self.add_opt('ifo', ifo_list_opt)

        self.add_opt('srate', cp.get('input', 'srate'))
        self.add_opt('seglen', cp.get('input', 'seglen'))
        self.add_opt('PSDlength', cp.get('input', 'PSDlength'))
 
        flow = ast.literal_eval(cp.get('input','flow'))
        for ifo in ifo_list:
            self.add_opt('{ifo}-flow'.format(ifo=ifo), str(flow[ifo]))
            self.add_opt('{ifo}-cache'.format(ifo=ifo), cacheFiles[ifo])
            self.add_opt('{ifo}-channel'.format(ifo=ifo), channel_list[ifo])

        # --- Optional options

        # self-checkpointing
        if cp.has_option('condor', 'checkpoint'):
            self.add_opt('checkpoint', cp.get('condor', 'checkpoint'))

        #
        # Priors
        #
        # Quality factor
        if cp.has_option('bayeswave_args', 'Qmin'):
            self.add_opt('Qmin', cp.get('bayeswave_args', 'Qmin'))
        if cp.has_option('bayeswave_args', 'Qmax'):
            self.add_opt('Qmax', cp.get('bayeswave_args', 'Qmax'))


        # dataseed
        if cp.has_option('bayeswave_args', 'dataseed'):
            self.add_opt('dataseed', cp.get('bayeswave_args', 'dataseed'))

        # Niter
        if cp.has_option('bayeswave_args', 'Niter'):
            self.add_opt('Niter', cp.get('bayeswave_args', 'Niter'))

        # Nchain
        if cp.has_option('bayeswave_args', 'Nchain'):
            self.add_opt('Nchain', cp.get('bayeswave_args', 'Nchain'))
            
        # Ncycle
        if cp.has_option('bayeswave_args', 'Ncycle'):
            self.add_opt('Ncycle', cp.get('bayeswave_args', 'Ncycle'))

        # Nburnin
        if cp.has_option('bayeswave_args', 'Nburnin'):
            self.add_opt('Nburnin', cp.get('bayeswave_args', 'Nburnin'))

        # chainseed
        if cp.has_option('bayeswave_args', 'chainseed'):
            self.add_opt('chainseed', cp.get('bayeswave_args', 'chainseed'))

        # runName
        if cp.has_option('bayeswave_args', 'runName'):
            self.add_opt('runName', cp.get('bayeswave_args', 'runName'))

        # 0noise
        if cp.has_option('bayeswave_args', '0noise'):
            self.add_opt('0noise', cp.get('bayeswave_args', '0noise'))

        # zeroLogL
        if cp.has_option('bayeswave_args', 'zeroLogL'):
            self.add_opt('zeroLogL', cp.get('bayeswave_args', 'zeroLogL'))

        # bayesLine
        if cp.has_option('bayeswave_args', 'BayesLine'):
            self.add_opt('bayesLine', cp.get('bayeswave_args', 'BayesLine'))

        # noClean
        if cp.has_option('bayeswave_args', 'noClean'):
            self.add_opt('noClean', cp.get('bayeswave_args', 'noClean'))

        # fixD
        if cp.has_option('bayeswave_args', 'fixD'):
            self.add_opt('fixD', cp.get('bayeswave_args', 'fixD'))

        # signalOnly
        if cp.has_option('bayeswave_args', 'signalOnly'):
            self.add_opt('signalOnly', cp.get('bayeswave_args', 'signalOnly'))

        # fullOnly
        if cp.has_option('bayeswave_args', 'fullOnly'):
            self.add_opt('fullOnly', cp.get('bayeswave_args', 'fullOnly'))

        # noSignal
        if cp.has_option('bayeswave_args', 'noSignal'):
            self.add_opt('noSignal', cp.get('bayeswave_args', 'noSignal'))

        # cleanOnly
        if cp.has_option('bayeswave_args', 'cleanOnly'):
            self.add_opt('cleanOnly', cp.get('bayeswave_args', 'cleanOnly'))

        # noiseOnly
        if cp.has_option('bayeswave_args', 'noiseOnly'):
            self.add_opt('noiseOnly', cp.get('bayeswave_args', 'noiseOnly'))

        # glitchOnly
        if cp.has_option('bayeswave_args', 'glitchOnly'):
            self.add_opt('glitchOnly', cp.get('bayeswave_args', 'glitchOnly'))

        # noPSDfit
        if cp.has_option('bayeswave_args', 'noPSDfit'):
            self.add_opt('noPSDfit', cp.get('bayeswave_args', 'noPSDfit'))

        # stochastic
        if cp.has_option('bayeswave_args', 'stochastic'):
            self.add_opt('stochastic', cp.get('bayeswave_args', 'stochastic'))

        # verbose
        if cp.has_option('bayeswave_args', 'verbose'):
            self.add_opt('verbose', cp.get('bayeswave_args', 'verbose'))

        # gnuplot
        if cp.has_option('bayeswave_args', 'gnuplot'):
            self.add_opt('gnuplot', cp.get('bayeswave_args', 'gnuplot'))

        # NC
        if cp.has_option('bayeswave_args','NC'):
            self.add_opt('NC', cp.get('bayeswave_args', 'NC'))
        # NCmin
        if cp.has_option('bayeswave_args','NCmin'):
            self.add_opt('NCmin', cp.get('bayeswave_args', 'NCmin'))
        # NCmax
        if cp.has_option('bayeswave_args','NCmax'):
            self.add_opt('NCmax', cp.get('bayeswave_args', 'NCmax'))

        #
        # Proposals
        #
        if cp.has_option('bayeswave_args', 'orientationProposal'):
            self.add_opt('orientationProposal', cp.get('bayeswave_args', 'orientationProposal'))

        #
        # Injection file
        #
        if injfile is not None:
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

        self.set_sub_file('bayeswave.sub')

class bayeswaveNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

    def __init__(self, bayeswave_job):

        pipeline.CondorDAGNode.__init__(self,bayeswave_job)
        pipeline.AnalysisNode.__init__(self)

    def set_trigtime(self, trigtime):
        self.add_var_opt('trigtime', trigtime)
        self.__trigtime = trigtime

    def set_PSDstart(self, PSDstart):
        self.add_var_opt('PSDstart', PSDstart)
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

#
# Post-processing
#

class bayeswave_postJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, cacheFiles, injfile=None, nrdata=None, dax=False):

        universe=cp.get('condor','universe')

        bayeswave_post_exec = cp.get('bayeswave_paths','bayeswave_post_executable')

        pipeline.CondorDAGJob.__init__(self,universe,'bayeswave_post')
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)

        if cp.has_option('condor', 'accounting_group'):
            self.add_condor_cmd('accounting_group', cp.get('condor', 'accounting_group'))   

        self.set_stdout_file('$(macrooutputDir)/bayeswave_post_$(cluster)-$(process)-$(node).out')
        self.set_stderr_file('$(macrooutputDir)/bayeswave_post_$(cluster)-$(process)-$(node).err')
        self.set_log_file('$(macrooutputDir)/bayeswave_post_$(cluster)-$(process)-$(node).log')

        # Request 4GB of RAM for pp jobs
        self.add_condor_cmd('request_memory', '4000')

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
            self.add_condor_cmd('transfer_output_files', '$(macrooutputDir)')

            # --- Files to include in transfer
            # FIXME: PostProc doesn't currently need frame transfer
            transferstring='datafind,$(macrooutputDir)'
            if injfile is not None:
                transferstring+=','+'SEOBNRv2ChirpTimeSS.dat,'+injfile
            if nrdata is not None: transferstring+=','+nrdata
            self.add_condor_cmd('transfer_input_files', transferstring)

            # --- Point to ROM data (which should have been copied
            if injfile is not None:
                self.add_condor_cmd("environment", 'LAL_DATA_PATH="./"')

        self.add_condor_cmd('getenv', 'True')

        # --- Required options
        ifo_list = ast.literal_eval(cp.get('input', 'ifo-list'))
        channel_list = ast.literal_eval(cp.get('datafind', 'channel-list'))

        # XXX: hack to repeat option
        ifo_list_opt = ifo_list[0]
        for ifo in ifo_list[1:]:
            ifo_list_opt += ' --ifo {0}'.format(ifo)
        self.add_opt('ifo', ifo_list_opt)

        self.add_opt('srate', cp.get('input', 'srate'))
        self.add_opt('seglen', cp.get('input', 'seglen'))
        self.add_opt('PSDlength', cp.get('input', 'PSDlength'))
 
        flow = ast.literal_eval(cp.get('input','flow'))
        for ifo in ifo_list:
            self.add_opt('{ifo}-flow'.format(ifo=ifo), str(flow[ifo]))

            # XXX: Postproc currently expects LALSimAdLIGO
            self.add_opt('{ifo}-cache'.format(ifo=ifo), "LALSimAdLIGO")
            self.add_opt('{ifo}-channel'.format(ifo=ifo), "LALSimAdLIGO")


        # --- Optional options
        # bayesLine
        if cp.has_option('bayeswave_args', 'BayesLine'):
            self.add_opt('bayesLine', cp.get('bayeswave_args', 'BayesLine'))

        # 0noise
        if cp.has_option('bayeswave_post_args', '0noise'):
            self.add_opt('0noise', cp.get('bayeswave_post_args', '0noise'))

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


        self.set_sub_file('bayeswave_post.sub')


class bayeswave_postNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

    def __init__(self, bayeswave_post_job):

        pipeline.CondorDAGNode.__init__(self, bayeswave_post_job)
        pipeline.AnalysisNode.__init__(self)

    def set_trigtime(self, trigtime):
        self.add_var_opt('trigtime', trigtime)
        self.__trigtime = trigtime

    def set_PSDstart(self, PSDstart):
        self.add_var_opt('PSDstart', PSDstart)
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

    def set_L1_timeslide(self, L1_timeslide):
        self.add_var_opt('L1-timeslide', L1_timeslide)
        self.__L1_timeslide = L1_timeslide


#
# skymap job
#

class megaskyJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, dax=False):

        # XXX consider local universe?
        universe='vanilla'

        # Point this to the src dir
        megasky_exec = cp.get('bayeswave_paths','megasky')
        pipeline.CondorDAGJob.__init__(self,universe,'megasky')
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)

        if cp.has_option('condor', 'accounting_group'):
            self.add_condor_cmd('accounting_group', cp.get('condor', 'accounting_group'))   

        self.set_stdout_file('megasky_$(cluster)-$(process)-$(node).out')
        self.set_stderr_file('megasky_$(cluster)-$(process)-$(node).err')
        self.set_log_file('megasky_$(cluster)-$(process)-$(node).log')

        self.set_sub_file('megasky.sub')


class megaskyNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

    def __init__(self, megasky_job, rundir):

        pipeline.CondorDAGNode.__init__(self, megasky_job)
        pipeline.AnalysisNode.__init__(self)

    def set_initialdir(self, rundir):
        #self.add_var_opt('trigtime', trigtime)
        self.__rundir = rundir

#
# megaplot job
#

class megaplotJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, dax=False):

        universe='vanilla'

        # Point this to the src dir
        pipeline.CondorDAGJob.__init__(self,universe,'megaplot')
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)

        if cp.has_option('condor', 'accounting_group'):
            self.add_condor_cmd('accounting_group', cp.get('condor', 'accounting_group'))   

        self.set_stdout_file('megaplot_$(cluster)-$(process)-$(node).out')
        self.set_stderr_file('megaplot_$(cluster)-$(process)-$(node).err')
        self.set_log_file('megaplot_$(cluster)-$(process)-$(node).log')

        self.set_sub_file('megaplot.sub')


class megaplotNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

    def __init__(self, megaplot_job, rundir):

        pipeline.CondorDAGNode.__init__(self, megaplot_job)
        pipeline.AnalysisNode.__init__(self)

#   def set_trigtime(self, trigtime):
#       self.add_var_opt('trigtime', trigtime)
#       self.__trigtime = trigtime

