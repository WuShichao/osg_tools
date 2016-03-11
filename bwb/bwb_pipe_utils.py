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

class BayesWaveBurstJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):

    def __init__(self, configparser, workdir, cacheFiles, dax=False):


        universe=configparser.get('general','universe')

        pipeline.CondorDAGJob.__init__(self,universe,'BayesWaveBurst')
        pipeline.AnalysisJob.__init__(self,configparser,dax=dax)

        self.set_stdout_file('logs/BayesWaveBurst_$(cluster)-$(process)-$(node).out')
        self.set_stderr_file('logs/BayesWaveBurst_$(cluster)-$(process)-$(node).err')
        self.set_log_file('logs/BayesWaveBurst_$(cluster)-$(process)-$(node).log')

        self.add_condor_cmd('should_transfer_files', 'YES')
        self.add_condor_cmd('when_to_transfer_output', 'ON_EXIT')
        self.add_condor_cmd('transfer_input_files',
                'BayesWaveBurst,datafind,$(macrooutdir),logs')

        # --- Common options
        ifoList = configparser.get('datafind', 'ifoList').split(',')
        channelList = configparser.get('datafind', 'channelList').split(',')

        self.add_ini_opts(configparser, 'bwb_args')
        # XXX: hack to repeat option on purpose...
        ifo_list_opt = ifoList[0]
        for ifo in ifoList[1:]:
            ifo_list_opt += ' --ifo {0}'.format(ifo)
        self.add_opt('ifo', ifo_list_opt)
 
        flow = configparser.get('general','flow')
        for i,ifo in enumerate(ifoList):
            self.add_opt('{ifo}-flow'.format(ifo=ifo), flow)
            self.add_opt('{ifo}-cache'.format(ifo=ifo), cacheFiles[ifo])
            self.add_opt('{ifo}-channel'.format(ifo=ifo), channelList[i])

        self.set_sub_file('{workdir}/BayesWaveBurst.sub'.format(workdir=workdir))


class BayesWaveBurstNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):

    new_id = itertools.count().next

    def __init__(self, bwb_job):

        pipeline.CondorDAGNode.__init__(self,bwb_job)
        pipeline.AnalysisNode.__init__(self)

    def set_trigtime(self, trigtime):
        self.add_var_opt('trigtime', trigtime)
        self.__trigtime = trigtime

    def set_psdtime(self, psdtime):
        self.add_var_opt('psdtime', psdtime)
        self.__psdtime = psdtime

    def set_outputDir(self, outputDir):
        self.add_var_opt('outputDir', outputDir)
        self.__outputDir = outputDir
  


#class BayesWavePostJob():
