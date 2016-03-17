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

# DAG Class definitions for BayesWaveBurst

from glue import pipeline
import itertools

#
# Main analysis
#

class BayesWaveBurstJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, cacheFiles, dax=False):


        universe=cp.get('condor','universe')

        pipeline.CondorDAGJob.__init__(self,universe,'BayesWaveBurst')
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)

        if cp.has_option('bwb_args', 'accounting_group'):
            self.add_condor_cmd('accounting_group', cp.get('condor', 'accounting_group'))   

        self.set_stdout_file('logs/BayesWaveBurst_$(cluster)-$(process)-$(node).out')
        self.set_stderr_file('logs/BayesWaveBurst_$(cluster)-$(process)-$(node).err')
        self.set_log_file('logs/BayesWaveBurst_$(cluster)-$(process)-$(node).log')

        self.add_condor_cmd('should_transfer_files', 'YES')
        self.add_condor_cmd('when_to_transfer_output', 'ON_EXIT')
        self.add_condor_cmd('transfer_input_files',
                'BayesWaveBurst,datafind,$(macrooutputDir),logs')
        self.add_condor_cmd('transfer_output_files', '$(macrooutputDir),logs')

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

        # noClean
        if cp.has_option('bwb_args', 'noClean'):
            self.add_opt('noClean', cp.get('bwb_args', 'noClean'))

        # fixD
        if cp.has_option('bwb_args', 'fixD'):
            self.add_opt('fixD', cp.get('bwb_args', 'fixD'))

        self.set_sub_file('BayesWaveBurst.sub')


class BayesWaveBurstNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

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
  

#
# Post-processing
#

class BayesWavePostJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, cp, cacheFiles, dax=False):


        universe=cp.get('condor','universe')

        pipeline.CondorDAGJob.__init__(self,universe,'BayesWavePost')
        pipeline.AnalysisJob.__init__(self,cp,dax=dax)

        if cp.has_option('bwb_args', 'accounting_group'):
            self.add_condor_cmd('accounting_group', cp.get('condor', 'accounting_group'))   

        self.set_stdout_file('logs/BayesWavePost_$(cluster)-$(process)-$(node).out')
        self.set_stderr_file('logs/BayesWavePost_$(cluster)-$(process)-$(node).err')
        self.set_log_file('logs/BayesWavePost_$(cluster)-$(process)-$(node).log')

        self.add_condor_cmd('should_transfer_files', 'YES')
        self.add_condor_cmd('when_to_transfer_output', 'ON_EXIT')
        self.add_condor_cmd('transfer_input_files',
                'BayesWavePost,datafind,$(macrooutputDir),logs')
        self.add_condor_cmd('transfer_output_files', '$(macrooutputDir),logs')

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

        self.set_sub_file('BayesWavePost.sub')


class BayesWavePostNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

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


