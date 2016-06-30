#!/usr/bin/python

import sys, json
from ligo.gracedb.rest import GraceDb
import os
import subprocess
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math



# -----------
# Main
# ----------

np.set_printoptions(precision=2)

# Get GraceID
dirname = os.getcwd()
graceid = dirname.split('_')[-1]
print graceid

IFO = [0,1]

freq = []
bandwidth = []
duration = []

dur_low = [0,0]
dur_high = [0,0]
dur_c = [0,0]

# Get frequency, bandwidth, duration info (duration on logscale)
for ifo in IFO:
   data = np.genfromtxt('tables/mode_{ifo}.txt'.format(ifo=ifo))
   freq.append(data[12])
   bandwidth.append(data[14])
   dur_low[ifo] = np.around(10**data[10,1],3)
   dur_high[ifo] = np.around(10**data[10,2],3)
   dur_c[ifo] = np.around(10**data[10,0],3)

freq = np.around(freq,2)
bandwidth = np.around(bandwidth,2)

print freq
print bandwidth
print dur_low, dur_high, dur_c

# Get Bayes factor +error info
data = np.genfromtxt('evidence_stacked.dat')

BSG = data[2]-data[1]
BSN = data[2]-data[0]

err_SG = math.sqrt(data[5]+data[4])
err_SN = math.sqrt(data[5]+data[3])

err_SG = round(err_SG,2)
err_SN = round(err_SN,2)

print BSG, err_SG
print BSN, err_SN


# Format information to be sent to gracedb
gdbtable = '<table>\
<tr><th colspan=2>BWB parameter estimation</th></tr> \
<tr><td>frequency (Hz)</td><td align=right>{freq}</td></tr> \
<tr><td>&nbsp;&nbsp;&nbsp;90% CI lower bound</td><td align=right>{freq}</td></tr> \
<tr><td>&nbsp;&nbsp;&nbsp;90% CI lower bound</td><td align=right>{freq}</td></tr> \
<tr><td>bandwidth (Hz)</td><td align=right>{bw}</td></tr> \
<tr><td>Duration (s)</td><td align=right>{dur}</td></tr> \
<tr><td>lnBSG</td><td align=right>{BSG}</td></tr> \
<tr><td>lnBSN</td><td align=right>{BSN}</td></tr> \
</table>'.format(freq=freq,bw=bandwidth,dur=duration,BSG=BSG,BSN=BSN)


gdbtable = '<table> \
<tr><th colspan=2>BWB parameter estimation</th></tr> \
<tr><th>Median</th><th>90% CI lower bound</th>\
<tr><td>frequency (Hz)</td><td align=right>{freq}</td></tr> \
<tr><td>&nbsp;&nbsp;&nbsp;90% CI lower bound</td><td align=right>{freq}</td></tr> \
<tr><td>&nbsp;&nbsp;&nbsp;90% CI lower bound</td><td align=right>{freq}</td></tr> \
<tr><td>bandwidth (Hz)</td><td align=right>{bw}</td></tr> \
<tr><td>Duration (s)</td><td align=right>{dur}</td></tr> \
<tr><td>lnBSG</td><td align=right>{BSG}</td></tr> \
<tr><td>lnBSN</td><td align=right>{BSN}</td></tr> \
</table>'.format(freq=freq,bw=bandwidth,dur=duration,BSG=BSG,BSN=BSN)

paramtable = '<table> \
<tr><th colspan=4>BWB parameter estimation</th></tr> \
<tr><th colspan=2>Param</th><th>Median</th><th>&nbsp;&nbsp;&nbsp;90%CI lower bound</th><th>&nbsp;&nbsp;&nbsp;90%CI upper bound</th></tr> \
<tr><td rowspan=2>frequency (Hz)</td><td>H1</td><td align=right>{freqH}</td><td align=right>{freqHlow}</td><td align=right>{freqHhigh}</td></tr> \
<tr><td>L1</td><td align=right>{freqL}</td><td align=right>{freqLlow}</td><td align=right>{freqLhigh}</td></tr> \
<tr><td rowspan=2>bandwidth (Hz)</td><td>H1</td><td align=right>{bwH}</td><td align=right>{bwHlow}</td><td align=right>{bwHhigh}</td></tr> \
<tr><td>L1</td><td align=right>{bwL}</td><td align=right>{bwLlow}</td><td align=right>{bwLhigh}</td></tr> \
<tr><td rowspan=2>duration (s)</td><td>H1</td><td align=right>{durH}</td><td align=right>{durHlow}</td><td align=right>{durHhigh}</td></tr> \
<tr><td>L1</td><td align=right>{durL}</td><td align=right>{durLlow}</td><td align=right>{durLhigh}</td></tr></table> \
'.format(freqH=freq[0][0],freqHlow=freq[0][1],freqHhigh=freq[0][2],freqL=freq[1][0],freqLlow=freq[1][1],freqLhigh=freq[1][2],bwH=bandwidth[0][0],bwHlow=bandwidth[0][1],bwHhigh=bandwidth[0][2],bwL=bandwidth[1][0],bwLlow=bandwidth[1][1],bwLhigh=bandwidth[1][2],durH=dur_c[0],durL=dur_c[1],durHlow=dur_low[0],durLlow=dur_low[1],durHhigh=dur_high[0],durLhigh=dur_high[1])

BFtable = '<table> \
<tr><th colspan=2>BWB Bayes Factors</th></tr> \
<tr><td>lnBSG</td><td align=right>{BSG}+/-{errBSG}</td></tr> \
<tr><td>lnBSN</td><td align=right>{BSN}+/-{errBSN}</td></tr> \
</table>'.format(BSG=BSG,BSN=BSN,errBSG=err_SG,errBSN=err_SN)

# Sky map
skyname = glob.glob('skymap*.fits')[0]

os.system('cp {sky} BW_skymap.fits'.format(sky=skyname)) # (change name so it's clear on GraceDB which skymap is ours)

#skytag = ["sky_loc","lvem"]
skytag = ["sky_loc"]

# Actually send info to gracedb
gracedb = GraceDb()
#gracedb.writeLog(graceid, "BayesWave Skymap image", filename='plots/skymap.png', tagname='sky_loc')
gracedb.writeLog(graceid, "BayesWave skymap FITS", filename='BW_skymap.fits', tagname=skytag)
gracedb.writeLog(graceid, "<a href='https://ldas-jobs.ligo.caltech.edu/~meg.millhouse/O1/zero_lag/job_{0}'>BWB Follow-up results</a>".format(graceid), tagname='pe')
gracedb.writeLog(graceid,paramtable,tagname='pe')
gracedb.writeLog(graceid,BFtable,tagname='pe')

os.chdir('..')
os.system('cp -r '+dirname+' /home/meg.millhouse/public_html/O1/zero_lag/job_'+graceid)
