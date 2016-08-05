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

import sys
import numpy as np
from glue import lal

# Get frame lists
cit_list = np.loadtxt(sys.argv[1], dtype=str)
gatech_list = np.loadtxt(sys.argv[2], dtype=str)

# Convert to cache, then segment lists
cit_cache=lal.Cache.from_urls(cit_list)
gatech_cache=lal.Cache.from_urls(gatech_list)

cit_segmentlist=cit_cache.to_segmentlistdict()
gatech_segmentlist=gatech_cache.to_segmentlistdict()


for segment in cit_segmentlist:
    if segment not in gatech_segmentlist:
        print "missing segment ", segment

print "*** Data Set Summary ***"
print "CIT: Start: {start} {end} [{duration}]".format(
        start
