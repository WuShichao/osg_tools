#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal.xlal.date import XLALGreenwichSiderealTime
import find_trig as ft
import healpy as hp
from subprocess import call
# from scipy.optimize import curve_fit
import acor
from pylal import bayespputils
import sky_area.sky_area_clustering as sky
import getopt
from readbwb import BwbParams

import lalinference.fits as lf
import lalinference.plot as lp
import lalinference.cmap

# -------------------------------------------
# Module to make skymaps, skyview webpage, etc.
# Works with single output directory from 
# BayesWaveBurst
# -------------------------------------------
def make_skyview(directory='.', mdc=None, NSIDE=128, ra=None, dec=None, results=None, npost=5000):

    # -- Close any open figures
    plt.close('all')

    # ----------------
    # Read S6 MDC log file
    # -----------------
    if mdc is None:
        print "No MDC log provided."

    # -- Save PWD for future reference
    topdir = os.getcwd()

    # -- Enter requested directory
    os.chdir(directory)
    print "Entered", os.getcwd()
    jobname = glob.glob('*bayeswave.run')[0]
    
    # ---------------------------------------
    # Extract information from bayeswave.run
    # Note: We may only need the trigger time from this
    # ---------------------------------------
    bayeswave = open(jobname, 'r')
    ifoNames = []
    for line in bayeswave:
        spl = line.split()
        if len(spl) < 1: continue
        
        # Determine names of IFOs and MDC scale factor
        if spl[0] == 'Command' and spl[1] == 'line:':
            for index, splElement in enumerate(spl):
                if splElement == '--ifo':
                    ifoNames.append(spl[index+1]) 
                    
        del spl[-1]
        # Determine number of IFOs
        if ' '.join(spl) == 'number of detectors':
            spl = line.split()
            ifoNum = spl[3]
            continue         
        # Determine gps trigger time
        if ' '.join(spl) == 'trigger time (s)':
            spl = line.split()
            gps = spl[3]
            break
    bayeswave.close()

    # --------------------------------
    # Create skymap summary statistics
    # --------------------------------
    # -- Get run name
    params = BwbParams()


    jobname = params.jobname

    # -- Input skymap data
    # num, L, ralist, sin_dec, psi, e, dA, dphi, dt = np.loadtxt(filename, unpack=True)
    print "Extracting RA/DEC samples"
    filename = './chains/' + jobname + 'signal_params.dat.0'
    data = np.loadtxt(filename, unpack=True,usecols=(0,1,2))
    ralist = data[1]
    sin_dec = data[2]
    print "Total samples are {0}".format(ralist.size)

    # -- Remove burn in samples
    burnin = ralist.size/4
    ralist = ralist[burnin:]
    sin_dec = sin_dec[burnin:]
    print "After removing burn-in samples are {0}".format(ralist.size)

    declist = np.arcsin(sin_dec)
    thetalist = np.pi/2.0 - declist

    radec = np.column_stack((ralist, declist))

    if radec.shape[0] > npost:
        radec = np.random.permutation(radec)[:npost,:]

    kde = sky.ClusteredSkyKDEPosterior(radec)

    # -- Get the sidereal time
    print "Finding sidereal time"
    print "Using GPS {0}".format(gps)
    trigtime = float(gps)
    lalgps = LIGOTimeGPS(trigtime)
    sidtime = XLALGreenwichSiderealTime(lalgps, 0)
    sidtime = sidtime % (np.pi*2)

    # -- Get the injection location
    print "GOT MDC?"
    print mdc
    if mdc is None:
        injtheta = 0
        injphi   = 0
        injdec   = 0
        injra    = 0
    else:
        injtheta, injphi = mdc.get_theta_phi(trigtime)
        injra = injphi + sidtime
        injdec = np.pi/2 - injtheta
        print "GOT INJECTION PARAMTERS"
        print injtheta, injra

    # -- Special handling for injections drawn from prior
    if params.ra is not None:
        injdec = params.dec 
        injra =  params.ra
        injtheta = np.pi/2 - injdec
        mdc = True
    
    # -- Make plots directory, if needed
    plotsDir = './plots'
    if not os.path.exists(plotsDir):
        os.makedirs(plotsDir)

    # Add markers (e.g., for injections or external triggers).
    if (ra is not None and dec is not None):
        # Convert the right ascension to either a reference angle from -pi to pi
        # or a wrapped angle from 0 to 2 pi, depending on the version of Matplotlib.
        #for rarad, decrad in [np.deg2rad([ra, dec])]:
        (rarad, decrad) = np.deg2rad([ra, dec])
        if geo:
            dlon = -sidtime # + or -?
            #rarad = lalinference.plot.reference_angle(rarad + dlon)
            rarad = rarad + dlon
            while rarad > np.pi:
                rarad = rarad - 2*np.pi
            while rarad < -np.pi:
                rarad = rarad + 2*np.pi
        else:
            #rarad = lalinference.plot.wrapped_angle(rarad)
            while rarad > 2*np.pi:
                rarad = rarad - 2*np.pi
            while rarad < 0:
                rarad = rarad + 2*np.pi
        # Shift the declination to match hp.projscatter conventions
        decrad = np.pi*0.5 - decrad
   
    # -- Plot the skymap and injection location
    skymap = kde.as_healpix(NSIDE, nest=True)
   
    #fig = plt.figure(frameon=False)
    fig = plt.figure(figsize=(8,6), frameon=False)
    ax = plt.subplot(111, projection='astro mollweide')
    ax.cla()
    ax.grid()
    lp.healpix_heatmap(skymap, nest=True, vmin=0.0, vmax=np.max(skymap), cmap=plt.get_cmap('cylon'))
	
    if (ra is not None and dec is not None):
        plt.plot(rarad, dec, 'kx', ms=30, mew=1)
    if mdc is not None:
        plt.plot(injra, injdec, 'kx', ms=30, mew=1)
    plt.savefig(plotsDir+'/skymap.png')
    plt.close()

    lf.write_sky_map('skymap_{0}.fits'.format(gps), skymap, nest=True, gps_time=trigtime)

    # -- Calculate the 50 and 90% credible intervals
    sq_deg = (180.0/np.pi)**2
    print "Calculating sky area ..."
    (area50, area90) = kde.sky_area( [0.5, 0.9] )*sq_deg
    print "Got area50 of {0}".format(area50)

    # -- Get the found contour and searched areas
    if mdc is not None:
        print "Calculating p value of injection"
        injcontour = kde.p_values( np.array([[injra, injdec]]) )[0]
        print "Got contour of {0}".format(injcontour)
    
        print "Calculating searched area..."
        searcharea = (kde.searched_area( np.array([[injra, injdec]]) )*sq_deg)[0]
        print "Got searched area of {0}".format(searcharea)
    else:
        injcontour = 0
        searcharea = 0

    pixarea = hp.nside2pixarea(NSIDE, degrees=True)
    print("here")
    if results is not None:
        for name in ['injcontour', 'area50', 'area90', 'searcharea']:
            if name not in results: results[name] = []
            results[name].append(eval(name)) 


    # -- Make an output file with key statistics
    outfile = open('skystats.txt', 'w')
    outfile.write("# area50  area90  searcharea  injcontour \n")
    outfile.write("{0} {1} {2} {3}".format(area50, area90, searcharea, injcontour))
    outfile.close()


# -- Write main script for command line running
if __name__ == "__main__":

    # Allow navigation into specified working directory
    topdir=os.getcwd()
    try:
        workdir=sys.argv[1]
    except IndexError:
        # No work directory specified, workdir=./
        workdir=os.getcwd()
    os.chdir(workdir)

    opts, args = getopt.getopt(sys.argv[1:], "", ['directory=', 'mdc=', 'NSIDE=', 'ra=', 'dec=', 'geo'])

    # -- Set default argument values
    directory='.'
    mdc=None
    NSIDE=128
    ra=None
    dec=None
    geo=False
    results=None
    
    print "Got these arguments"
    print opts

    print "With these values:"
    print args

    for opt, arg in opts:
        print opt
        if opt=='--mdc':
            # -- mdc argument should be the name of the MDC log
            mdc = ft.Mdc(arg)
            print "Reading MDC log {0}".format(arg)
            try:
                mdc = ft.Mdc(arg)
            except: 
                print "WARNING!  Failed to read MDC log file"
                #mdc = None
        if opt == '--directory':
            directory = arg
        if opt == '--NSIDE':
            NSIDE = int(arg)
        if opt == '--sim':
            sim = 'True' == arg
        if opt == '--ra':
            ra = float(arg)
        if opt == '--dec':
            dec = float(arg)
        if opt == '--geo':
            geo = 'True'
	
    make_skyview(directory, mdc, NSIDE, ra, dec, results)

    # Move back to original dir
    os.chdir(topdir)

