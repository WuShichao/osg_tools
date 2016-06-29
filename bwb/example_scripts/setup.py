#!/usr/bin/python

import numpy as np
import time
import sys, json
from ligo.gracedb.rest import GraceDb
import os
import subprocess

# ----------------
# Set Parameters
# ----------------

########### change me ########################
topdir = '/home/jclark/O1/BG/C01/unconstrained_e' # CHANGE ME TO YOUR OUTPUT DIRECTORY
triglist ='/home/meg.millhouse/O1/misc_stuff/CO1_large_data_set/C01_unconstrained_e.dat' # CHANGE ME TO YOUR TRIGLIST
#############################################

cachedir = os.path.join(topdir, 'cachedir')
bwb = '/home/meg.millhouse/bayeswave/branches/O2_dev/src/BayesWaveBurst' # change to your executable if you want--- but can keep it to point to Meg's
svn = '/home/meg.millhouse/bayeswave/trunk/'
ifoList = ['H1', 'L1']
frtypeList = ['H1_HOFT_C01', 'L1_HOFT_C01']
name = topdir

# ------------------------------------------
# Read list of GPS times and trigger numbers
# ------------------------------------------


names = ['veto1','veto2','rho','cc1','cc2','cc3','amp','tshift','tsupershift','like','penalty','disbalance','f','bandwidth','duration','pixels','resolution','runnumber','Lgps','Hgps','sSNRL','sSNRH','hrssL','hrssH','phi','theta','psi']

data = np.recfromtxt(triglist,names=names)

Hgps = data['Hgps']
Lgps = data['Lgps']
rhoList = data['rho']
freqList = data['f']

plusveto = data['veto1']
minusveto = data['veto2']

lagList = []

for h,l in zip(Hgps,Lgps):
   lagList.append(round(h-l))

gpsList = Hgps

# -----------------
# BWB arguments 
# -----------------
bwbargsfmt = """--ifo H1 --H1-flow 32 --H1-channel H1:DCS-CALIB_STRAIN_C01   \
       --ifo L1 --L1-flow 32 --L1-channel L1:DCS-CALIB_STRAIN_C01   \
       --H1-cache {cachedir}/H1.cache \
       --L1-cache {cachedir}/L1.cache \
       --trigtime {gps} --srate {srate} --seglen 4 \
       --bayesLine --PSDstart {gps} --PSDlength 4 \
       --outputDir {trigdir}  \
       --L1-timeslide {lag} \

"""
#       --Niter 2000000 --NCmin 35 --NCmax 35  \
#       --MDC-channel [H1:GW-H,L1:GW-H]   \
#       --MDC-cache [{cachefile},{cachefile}] \
#        --MDC-prefactor {scale} \

# ----------
# Template for bwb submit file
# ----------
submit_str = """
executable={bwb}
universe=standard
arguments={bwbargs}
output={top}/condorOut.txt
error={top}/condorError.txt
log={top}/condorLog.txt
notification=never
should_transfer_files=YES
when_to_transfer_output = ON_EXIT
stream_error=False
stream_output=False
WantRemoteIO=False
accounting_group = ligo.prod.o1.burst.paramest.bayeswave
queue 1
"""

# ----------------------------
# Template for BWP submit file
# ----------------------------
post_submitstr = """
executable={top}/pp.sh
universe=vanilla
arguments=""
output={top}/ppOut.txt
error={top}/ppError.txt
log={top}/ppLog.txt
notification=never
should_transfer_files=YES
when_to_transfer_output = ON_EXIT
RequestMemory = 2000
stream_error=False
stream_output=False
WantRemoteIO=False
accounting_group = ligo.prod.o1.burst.paramest.bayeswave
queue 1
"""



# ------------------
# Call LIGO Data find
# -------------------


# Set min/max gps times for LIGO data find:
min_gps = min(gpsList) - (max(np.absolute(lagList))+25.0)
max_gps = max(gpsList) + (max(np.absolute(lagList))+25.0)

start = min_gps # For start and end GPS times- pad high/low GPS times by largest magnitude lag
end   = max_gps

print min_gps, max_gps

if not os.path.exists(cachedir): os.makedirs(cachedir)

for ifo, frtype in zip(ifoList,frtypeList):
    cachefilefmt = os.path.join(cachedir, '{0}.cache')
    ldfcmd = "gw_data_find --observatory {o} --type {frtype} -s {start} -e {end} --lal-cache | grep file > {cachefile}".format(o=ifo[0], frtype=frtype, cachefile = cachefilefmt.format(ifo), start=start, end=end)
    print "Calling LIGO data find ..."
    print ldfcmd
    subprocess.call(ldfcmd, shell=True)



# ---------------
# Open dag file
# ---------------

dagfile = open( os.path.join(name, 'dagfile.dag'), 'w')

counter = 0

for gps, lag, freq, rho, veto1, veto2 in zip(gpsList, lagList, freqList, rhoList, plusveto, minusveto):
   if rho < 8.0: continue

   counter += 1
   number = float(gps)
   trigdir = os.path.join(name, 'job_' + str(int(number)) + '_' + str(lag))
   if os.path.exists(trigdir): print str(int(number)) + '_' + str(lag)+' exists'
   if not os.path.exists(trigdir): os.makedirs(trigdir)
   if freq > 200.0:
      srate = 2048
   else:
      srate = 1024
   
   # --- make cache file

# -- Create info file
   infoname = os.path.join(trigdir, 'job_info.txt')
   infofile = open(infoname,'w')
   infofile.write("{rho} {gps} {lag} {freq} {veto1} {veto2}\n".format(rho=rho,gps=gps,lag=lag,freq=freq,veto1=veto1,veto2=veto2))
   #infofile.write(str(rho)+" "+str(gps)+" "+str(lag)+" "+str(freq)+'\n')
   infofile.close()

# -- Create BWB submit file
   bwbargs = bwbargsfmt.format(gps=gps, cachedir=cachedir, trigdir=trigdir, lag=lag, srate=srate)
   submitname = os.path.join(trigdir, 'submitBWB.txt')
   submitfile = open(submitname, 'w')
   submitfile.write(submit_str.format(top=trigdir, bwb=bwb, bwbargs=bwbargs))
   submitfile.close()

# -- Create BWP submit file
   post_submitname = os.path.join(trigdir, 'submitBWP.txt')
   submitfile = open(post_submitname, 'w')
   submitfile.write(post_submitstr.format(top=trigdir))
   submitfile.close()

# -- Create BWP shell script
   ppname = os.path.join(trigdir, 'pp.sh')
   doitfile = open(ppname, 'w')
   doitfile.write("""#! /bin/sh
    cd {directory}
    source /home/meg.millhouse/.bashrc
    python {svn}/postprocess/s6/run_bwp.py
    python {svn}/postprocess/megaplot.py
    """.format(directory=trigdir, svn=svn))
   doitfile.close()


# ---- write the dag gile
   graceid = 'job_' + str(number) + '_' + str(lag)

   dagfile.write("JOB {graceid} {submitname}\n".format(graceid=graceid, submitname=submitname))
   dagfile.write("RETRY {graceid} 1 \n\n".format(graceid=graceid))

dagfile.close()
