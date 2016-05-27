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



import sys
from optparse import OptionParser
import ConfigParser
from lalapps import inspiralutils
from glue import segmentsUtils

parser=OptionParser()
opts,args = parser.parse_args()
config=ConfigParser.ConfigParser()
config.read(args[0])

ifo='H1'
veto_categories=[1,2]


(segFileName,dqVetoes)=inspiralutils.findSegmentsToAnalyze(config, ifo,
        veto_categories, generate_segments=True,
        use_available_data=True, data_quality_vetoes=True)

segfile=open(segFileName)
segs=segmentsUtils.fromsegwizard(segfile)



