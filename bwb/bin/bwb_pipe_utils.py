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

#
# Main analysis
#

class bayeswaveJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, cacheFiles, injfile=None, nrdata=None, dax=False):

        universe=cp.get('condor','universe')

        pipeline.CondorDAGJob.__init__(self,universe,'bayeswave')
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
        if 'pace.gatech.edu' in hostname:
            print >> sys.stdout, "Looks like you're on PACE; configuring file transfers"

            # --- Perform file transfers
            self.add_condor_cmd('should_transfer_files', 'YES')
            self.add_condor_cmd('when_to_transfer_output', 'ON_EXIT_OR_EVICT')
            self.add_condor_cmd('transfer_output_files', '$(macrooutputDir)')

            # --- Files to include in transfer
            transferstring='datafind,$(macrooutputDir)'
            if injfile is not None: transferstring+=','+injfile
            if nrdata is not None: transferstring+=','+nrdata
            self.add_condor_cmd('transfer_input_files', transferstring)

        self.add_condor_cmd('getenv', 'True')

        # --- Required options
        ifoList = cp.get('datafind', 'ifoList').split(',')
        channelList = cp.get('datafind', 'channelList').split(',')

        # XXX: hack to repeat option
        ifo_list_opt = ifoList[0]
        for ifo in ifoList[1:]: ifo_list_opt += ' --ifo {0}'.format(ifo)
        self.add_opt('ifo', ifo_list_opt)

        self.add_opt('srate', cp.get('bwb_args', 'srate'))
        self.add_opt('seglen', cp.get('bwb_args', 'seglen'))
        self.add_opt('PSDlength', cp.get('bwb_args', 'PSDlength'))
 
        flow = cp.get('bwb_args','flow')
        for i,ifo in enumerate(ifoList):
            self.add_opt('{ifo}-flow'.format(ifo=ifo), flow)
            self.add_opt('{ifo}-cache'.format(ifo=ifo), cacheFiles[ifo])
            self.add_opt('{ifo}-channel'.format(ifo=ifo), channelList[i])

        # --- Optional options

        # self-checkpointing
        if cp.has_option('condor', 'checkpoint'):
            self.add_opt('checkpoint', cp.get('condor', 'checkpoint'))

        #
        # Priors
        #
        # Quality factor
        if cp.has_option('bwb_args', 'Qmin'):
            self.add_opt('Qmin', cp.get('bwb_args', 'Qmin'))
        if cp.has_option('bwb_args', 'Qmax'):
            self.add_opt('Qmax', cp.get('bwb_args', 'Qmax'))


        # dataseed
        if cp.has_option('bwb_args', 'dataseed'):
            self.add_opt('dataseed', cp.get('bwb_args', 'dataseed'))

        # Niter
        if cp.has_option('bwb_args', 'Niter'):
            self.add_opt('Niter', cp.get('bwb_args', 'Niter'))

        # Nchain
        if cp.has_option('bwb_args', 'Nchain'):
            self.add_opt('Nchain', cp.get('bwb_args', 'Nchain'))
            
        # Ncycle
        if cp.has_option('bwb_args', 'Ncycle'):
            self.add_opt('Ncycle', cp.get('bwb_args', 'Ncycle'))

        # Nburnin
        if cp.has_option('bwb_args', 'Nburnin'):
            self.add_opt('Nburnin', cp.get('bwb_args', 'Nburnin'))

        # chainseed
        if cp.has_option('bwb_args', 'chainseed'):
            self.add_opt('chainseed', cp.get('bwb_args', 'chainseed'))

        # runName
        if cp.has_option('bwb_args', 'runName'):
            self.add_opt('runName', cp.get('bwb_args', 'runName'))

        # 0noise
        if cp.has_option('bwb_args', '0noise'):
            self.add_opt('0noise', cp.get('bwb_args', '0noise'))

        # zeroLogL
        if cp.has_option('bwb_args', 'zeroLogL'):
            self.add_opt('zeroLogL', cp.get('bwb_args', 'zeroLogL'))

        # bayesLine
        if cp.has_option('bwb_args', 'BayesLine'):
            self.add_opt('bayesLine', cp.get('bwb_args', 'BayesLine'))

        # noClean
        if cp.has_option('bwb_args', 'noClean'):
            self.add_opt('noClean', cp.get('bwb_args', 'noClean'))

        # fixD
        if cp.has_option('bwb_args', 'fixD'):
            self.add_opt('fixD', cp.get('bwb_args', 'fixD'))

        # signalOnly
        if cp.has_option('bwb_args', 'signalOnly'):
            self.add_opt('signalOnly', cp.get('bwb_args', 'signalOnly'))

        # fullOnly
        if cp.has_option('bwb_args', 'fullOnly'):
            self.add_opt('fullOnly', cp.get('bwb_args', 'fullOnly'))

        # noSignal
        if cp.has_option('bwb_args', 'noSignal'):
            self.add_opt('noSignal', cp.get('bwb_args', 'noSignal'))

        # cleanOnly
        if cp.has_option('bwb_args', 'cleanOnly'):
            self.add_opt('cleanOnly', cp.get('bwb_args', 'cleanOnly'))

        # noiseOnly
        if cp.has_option('bwb_args', 'noiseOnly'):
            self.add_opt('noiseOnly', cp.get('bwb_args', 'noiseOnly'))

        # glitchOnly
        if cp.has_option('bwb_args', 'glitchOnly'):
            self.add_opt('glitchOnly', cp.get('bwb_args', 'glitchOnly'))

        # noPSDfit
        if cp.has_option('bwb_args', 'noPSDfit'):
            self.add_opt('noPSDfit', cp.get('bwb_args', 'noPSDfit'))

        # stochastic
        if cp.has_option('bwb_args', 'stochastic'):
            self.add_opt('stochastic', cp.get('bwb_args', 'stochastic'))

        # verbose
        if cp.has_option('bwb_args', 'verbose'):
            self.add_opt('verbose', cp.get('bwb_args', 'verbose'))

        # gnuplot
        if cp.has_option('bwb_args', 'gnuplot'):
            self.add_opt('gnuplot', cp.get('bwb_args', 'gnuplot'))

        #
        # Proposals
        #
        if cp.has_option('bwb_args', 'orientationProposal'):
            self.add_opt('orientationProposal', cp.get('bwb_args', 'orientationProposal'))

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
            self.add_opt('MDC-cache', cp.get('injections', 'mdc-cache'))

        if cp.has_option('injections', 'mdc-channel'):
            self.add_opt('MDC-channel', cp.get('injections', 'mdc-channel'))

        if cp.has_option('injections', 'mdc-prefactor'):
            self.add_opt('MDC-prefactor', cp.get('injections', 'mdc-prefactor'))

        self.set_sub_file('bayeswave.sub')

class bayeswaveNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

    def __init__(self, bwb_job):

        pipeline.CondorDAGNode.__init__(self,bwb_job)
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

  

#
# Post-processing
#

class bayeswave_postJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, cacheFiles, injfile=None, nrdata=None, dax=False):

        universe=cp.get('condor','universe')

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
            transferstring='datafind,$(macrooutputDir)'
            if injfile is not None: transferstring+=','+injfile
            if nrdata is not None: transferstring+=','+nrdata
            self.add_condor_cmd('transfer_input_files', transferstring)

        self.add_condor_cmd('getenv', 'True')

        # --- Required options
        ifoList = cp.get('datafind', 'ifoList').split(',')
        channelList = cp.get('datafind', 'channelList').split(',')

        # XXX: hack to repeat option
        ifo_list_opt = ifoList[0]
        for ifo in ifoList[1:]:
            ifo_list_opt += ' --ifo {0}'.format(ifo)
        self.add_opt('ifo', ifo_list_opt)

        self.add_opt('srate', cp.get('bwb_args', 'srate'))
        self.add_opt('seglen', cp.get('bwb_args', 'seglen'))
        self.add_opt('PSDlength', cp.get('bwb_args', 'PSDlength'))
 
        flow = cp.get('bwb_args','flow')
        for i,ifo in enumerate(ifoList):
            self.add_opt('{ifo}-flow'.format(ifo=ifo), flow)
            self.add_opt('{ifo}-cache'.format(ifo=ifo), cacheFiles[ifo])
            self.add_opt('{ifo}-channel'.format(ifo=ifo), channelList[i])


        # --- Optional options
        # bayesLine
        if cp.has_option('bwb_args', 'BayesLine'):
            self.add_opt('bayesLine', cp.get('bwb_args', 'BayesLine'))

        # 0noise
        if cp.has_option('bwp_args', '0noise'):
            self.add_opt('0noise', cp.get('bwp_args', '0noise'))

        #
        # Injection file
        #
        if injfile is not None:
            self.add_opt('inj', injfile)

        if nrdata is not None:
            self.add_opt('inj-numreldata', nrdata)

        #
        # MDC Setup
        #
        if cp.has_option('inj', 'mdc-cache'):
            self.add_opt('MDC-cache', cp.get('inj', 'mdc-cache'))

        if cp.has_option('inj', 'mdc-channel'):
            self.add_opt('MDC-channel', cp.get('inj', 'mdc-channel'))


        self.set_sub_file('bayeswave_post.sub')


class bayeswave_postNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

    def __init__(self, bwp_job):

        pipeline.CondorDAGNode.__init__(self, bwp_job)
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


