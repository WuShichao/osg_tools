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
This script selects GaTech waveforms with sufficient cycles to model chirp
masses down to 28 Msun with a 30 Hz cutoff.

The output is a text file with the names of the permitted waveforms and the
total mass to which they should be scaled such that the chirp mass of that
waveform matches the best estimate from CBC PE of GW150914
"""

import sys, os
import os.path
import subprocess
import cPickle as pickle

import timeit
import numpy as np

from pycbc import pnutils
import nrburst_pca_utils as nrbu_pca
import nrburst_utils as nrbu

# XXX Useful info for GW150914
#detector frame chirp mass: 30.3 +2.1±0.4 (−1.9±0.4)
pe_mchirp=30.3


#
# --- catalog Definition
#
bounds = dict()
bounds['Mchirpmin30Hz'] = [-np.inf, 28.0]

#
# --- Generate initial catalog
#
print >> sys.stdout,  '~~~~~~~~~~~~~~~~~~~~~'
print >> sys.stdout,  'Selecting Simulations'
print >> sys.stdout,  ''
then = timeit.time.time()

# Select simulations
simulations = nrbu.simulation_details(param_bounds=bounds,
        catdir='/home/jclark308/lvc_nr/GaTech')


f=open('GW150914_masses.txt','w')
for sim in simulations.simulations:

    wavefile=sim['wavefile'].split('/')[-1].replace('.h5','')

    m1,m2 = pnutils.mchirp_eta_to_mass1_mass2(pe_mchirp, sim['eta'])
    print m1,m2

    f.writelines("%s %f\n"%(wavefile,m1+m2))

f.close()


