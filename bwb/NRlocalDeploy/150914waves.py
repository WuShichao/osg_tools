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
"""
bhex_pca.py
"""

import sys, os
import os.path
import subprocess
import cPickle as pickle

import timeit
import numpy as np

import nrburst_pca_utils as nrbu_pca
import nrburst_utils as nrbu

#
# --- catalog Definition
#
bounds = dict()
bounds['Mchirpmin30Hz'] = [-np.inf, 27.0]

#
# --- Generate initial catalog
#
print >> sys.stdout,  '~~~~~~~~~~~~~~~~~~~~~'
print >> sys.stdout,  'Selecting Simulations'
print >> sys.stdout,  ''
then = timeit.time.time()

# Select simulations
simulations = nrbu.simulation_details(param_bounds=bounds,
        catdir='/data/lvc_nr/GaTech')


f=open('allowed_waves.txt','w')
for sim in simulations.simulations:
    wavefile=sim['wavefile'].split('/')[-1].replace('.h5','')
    f.writelines("%s\n"%wavefile)
f.close()


