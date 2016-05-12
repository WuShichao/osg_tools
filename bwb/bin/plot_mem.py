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
from matplotlib import pyplot as pl


data = np.loadtxt(sys.argv[1], dtype=str)
time_stamps = data[:,0]
memsize = data[:,1]
memsize = np.array([float(val) for val in memsize])

f, ax = pl.subplots()
ax.bar(range(len(memsize)), memsize)
ax.set_xticks(range(len(memsize)))
ax.set_xticklabels(data[:,0], rotation='vertical')
ax.set_xlabel('Time Stamp (local @ CIT)')
ax.set_ylabel('Image size (kb)')
ax.set_title(sys.argv[1])
f.tight_layout()

pl.show()
