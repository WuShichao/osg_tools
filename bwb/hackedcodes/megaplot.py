#!/usr/bin/env python
"""
This script generates a webpage to display the results of a BayesWave run.

This script was tested on ldas-pcdev1.ligo.caltech.edu only.
"""

__author__ = "Jonah Kanner, Francesco Pannarale"
__email__ = "jkanner@caltech.edu, francesco.pannarale@ligo.org"
__version__ = "1.3"
__date__ = "07.10.2014"

######################################################################################################################
#
# Import modules
#
######################################################################################################################
import argparse
import glob
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import os
import pwd
import subprocess
import sys
import re,subprocess

print 'Path to megaplot: '
print sys.argv[0]
print '\n'

# Allow navigation into specified working directory
topdir=os.getcwd()
try:
    workdir=sys.argv[1]
except IndexError:
    # No work directory specified, workdir=./
    workdir=os.getcwd()
os.chdir(workdir)

# See: http://pyinstaller.readthedocs.io/en/latest/runtime-information.html
if getattr(sys, 'frozen', False):
        # we are running in a bundle
        print 'Running as compiled binary'
        postprocesspath = sys._MEIPASS

        svn_info = open(os.path.join(postprocesspath, 'svn_info.txt'),'r')
        print "---- SVN INFO ---- "
        for info_line in svn_info.readlines():
            print info_line
        svn_info.close()
        print "---- END SVN INFO ---- "

else:
        # we are running in a normal Python environment
        print 'Running as interpreted script'
        postprocesspath = os.path.dirname(os.path.abspath(__file__))

        print 'svn revision: '
        try:
            os.system("svn info "+postprocesspath+" | grep \"Revision\" | awk '{print $2}'")
        except:
            print 'Couldn\'t find svn version from directory: ', postprocesspath
        print '\n'

#from cycler import cycler
from scipy import integrate

######################################################################################################################
#
# Dictionaries and other global variables
#
######################################################################################################################
# Item labels must coincide with moment names given in BayesWavePost
# 1st entry = x-axis label
# 2nd entry = file name tag
moments_dict = {}
moments_dict['dur_rec']               = ['log10(Duration)', 'duration']
moments_dict['t_energy_rec']          = ['log10(Signal Energy)', 'energy']
moments_dict['f0_rec']                = ['Central Freq (Hz)', 'freq']
moments_dict['band_rec']              = ['Bandwidth (Hz)', 'band']
moments_dict['t0_rec']                = ['Central Time (s)', 'tzero']
moments_dict['overlap']               = ['Overlap - recovered signal and injection', 'overlap']
moments_dict['network_overlap']       = ['Network overlap - recovered signal and injection', 'network_overlap']
moments_dict['overlap_plus']          = ['Plus Polarization Overlap - recovered signal and injection', 'overlap_plus']
moments_dict['network_overlap_plus']  = ['Plus Polarization Network overlap - recovered signal and injection', 'network_overlap_plus']
moments_dict['overlap_cross']         = ['Cross Polarization Overlap - recovered signal and injection', 'overlap_cross']
moments_dict['network_overlap_cross'] = ['Cross Polarization Network overlap - recovered signal and injection', 'network_overlap_cross']
moments_dict['h_max']                 = ['|h(t)| maximum value', 'hmax']
moments_dict['t_at_h_max']            = ['Time which |h(t)| is maximum (s)', 't_at_hmax']
moments_dict['snr']                   = ['Recovered SNR','snr']

# 1st entry = html file name
# 2nd entry = subpage header
# 3rd entry = plot file name tag
html_dict = {}
html_dict['dur_rec']           = ['duration', 'Histogram of duration', moments_dict['dur_rec'][1]]
html_dict['t_energy_rec']      = ['energy', 'Histogram of signal energy, (aka ||<i>h(t)</i>||<sup>2</sup>)', moments_dict['t_energy_rec'][1]]
html_dict['f0_rec']            = ['f0', 'Histogram of central frequency', moments_dict['f0_rec'][1]]
html_dict['band_rec']          = ['band', 'Histogram of bandwidth', moments_dict['band_rec'][1]]
html_dict['t0_rec']            = ['t0', 'Histogram of central time', moments_dict['t0_rec'][1]]
html_dict['overlap']           = ['overlap', 'Overlap Histogram between recovered signal and injection', moments_dict['overlap'][1]]
html_dict['overlap_plus']      = ['overlap', 'Overlap Histogram between recovered signal and injection for the plus polarization', moments_dict['overlap_plus'][1]]
html_dict['overlap_cross']     = ['overlap', 'Overlap Histogram between recovered signal and injection for the cross polarization', moments_dict['overlap_cross'][1]]
html_dict['h_max']             = ['hmax', 'Histogram of the maximum of ||<i>h(t)</i>||', moments_dict['h_max'][1]]
html_dict['t_at_h_max']        = ['t_at_hmax', 'Histogram of the time at the maximum of ||<i>h(t)</i>||', moments_dict['t_at_h_max'][1]]
#html_dict['waveform']          = ['waveform', 'Median Model Waveforms and 1-sigma Uncertainties', 'waveform']
html_dict['signal_waveform']   = ['', 'Median Signal Model Waveforms and 1-sigma Uncertainties', 'signal_waveform']
html_dict['glitch_waveform']   = ['', 'Median Glitch Model Waveforms and 1-sigma Uncertainties', 'glitch_waveform']
html_dict['bayesogram']        = ['bayesogram', 'Bayesogram (time): 2-D Histogram of Posterior', 't_bgram']
html_dict['signal_bayesogram'] = ['', 'Bayesogram for the signal model waveforms', 'signal_t_bgram']
html_dict['glitch_bayesogram'] = ['', 'Bayesogram for the glitch model waveforms', 'glitch_t_bgram']
html_dict['spec']              = ['spec', 'Spectrogram of median reconstructed model waveform', 'self_spec']
html_dict['signal_spec']       = ['', 'Spectrogram of median reconstructed signal model waveform', 'signal_self_spec']
html_dict['glitch_spec']       = ['', 'Spectrogram of median reconstructed glitch model waveform', 'glitch_self_spec']
html_dict['diagnostics']       = ['diagnostics', 'Diagnostic plots']
html_dict['skymap']            = ['skymap', 'Skymap']
html_dict['inj_spec']          = ['injections', 'Injection spectrograms', 'inj_self_spec']
html_dict['snr']               = ['snr','SNR','snr']
html_dict['time_domain']       = ['timedomain','Time Domain Information','waveform']
html_dict['signal_time_domain']       = ['timedomain','Recovered signal time domain waveform','signal_waveform']
html_dict['glitch_time_domain']       = ['timedomain','Recovered glitch time domain waveform','glitch_waveform']
html_dict['frequency_domain']  = ['freqdomain','Frequency Domain Information','powerspec']
html_dict['signal_frequency_domain']  = ['freqdomain','Median PSD and recovered signal in frequency domain','signal_frequence_domain']
html_dict['glitch_frequency_domain']  = ['freqdomain','Median PSD and recovered glitch in frequency domain','glitch_frequence_domain']
html_dict['time_frequency']    = ['timefreq','Time-Frequency Information','freq_domain']

modelList = ('signal', 'glitch')
upperCaseModel_dict = {'glitch': 'Glitch', 'signal': 'Signal', 'noise': 'Noise', 'clean': 'Clean'}

postDir   = 'post/'
plotsDir  = 'plots/'
tablesDir = 'tables/'
htmlDir   = 'html/'

#Adopt common color scheme for different models
ncolor = 'darkgrey'
gcolor = 'darkgoldenrod'
scolor = 'darkorchid'

H1color = 'darkgoldenrod'
L1color = 'darkkhaki'
V1color = 'dkarsage'
ifoColors = [H1color,L1color,V1color]

injcolor = 'teal'

######################################################################################################################
#
# Read in run data and info
#
######################################################################################################################

# ----------------------------------
# Extract information about the run
# ----------------------------------
def readbwb():
    # -- Initialize variables with info about the job
    info = ""
    ifoNames = []
    ifoList = []
    jobName = ''
    restrictModel = ''
    injFlag = False
    lalsimFlag = False
    mdc = False
    snrList = []
    # -- Open *bayeswave.run 
    bayeswaverunfile = glob.glob('*bayeswave.run')
    if not len(bayeswaverunfile):
        sys.exit("\n *bayeswave.run not found! \n")
    bayeswaverunfile = bayeswaverunfile[0]
    bayeswave = open(bayeswaverunfile, 'r')
    # -- Parse command line arguments
    lines = bayeswave.readlines()
    cmdline = lines[7].split(' ')

    for index, arg in enumerate(cmdline):
        # Determine the IFOs involved
        if arg=='--ifo':
            ifoNames.append(cmdline[index+1])
        # Determine the trigger time
        elif arg=='--trigtime':
            gps = float(cmdline[index+1])
            info = info + "The trigger time is GPS: {0}\n".format(gps)
        # Determine the job name
        elif arg=='--runName':
            jobName = cmdline[index+1]
            print "The job name is: {0}".format(jobName)
            jobName = jobName+'_'
        # Determine if the job was glitchOnly, noiseOnly, or signalOnly
        elif arg=='--glitchOnly':
            print '\nThis run was executed with the --glitchOnly flag\n'
            restrictModel = 'glitch'
        elif arg=='--noiseOnly':
            print '\nThis run was executed with the --noiseOnly flag\n'
            restrictModel = 'noise'
        elif arg=='--signalOnly':
            print '\nThis run was executed with the --signalOnly flag\n'
            restrictModel = 'signal'
        elif arg=='--inj':
            injFlag = True
            mdc = True
        elif arg=='--MDC-cache':
            mdc = True
        elif arg=='--lalinspiralinjection':
            lalsimFlag = True

    if lalsimFlag: injFlag = False
    # -- Determine number of IFOs
    ifoNum = len(ifoNames)
    # -- Convert the IFO name list into a list of numbers useful for opening datafiles
    if ifoNum == 1:
        ifoList = ['0']
    elif ifoNum == 2:
        ifoList = ['0', '1']
    elif ifoNum == 3: 
        ifoList = ['0', '1', '2']
    info = info + "Detectors(s) found: {0}\n".format(len(ifoList))
    info = info + "The detectors are: {0}\n".format(', '.join(ifoNames))
    # -- If MDC, read SNR info
    if mdc:
        for ifo in ifoList:
            snrFile = jobName+'post/injection_whitened_moments.dat.'+ifo
            if not os.path.exists(snrFile):
                sys.exit("\n {0} not found! \n".format(snrFile))
            snrdata = open(snrFile, 'r')
            snrdata.readline()
            snr = float(snrdata.readline().split(' ')[0])
            snrList.append(snr)
            snrdata.close()
            #TODO: thse are different from the values in condorOut*.txt
            #Check numbers in injection_colored_moments.dat.*
            info = info + 'Injected SNR in detector {0} = {1}\n'.format(ifoNames[int(ifo)],snrList[-1])
    bayeswave.close()
    # -- Report to user
    print "{0}".format(info)
    return(jobName, restrictModel, mdc, injFlag, bayeswaverunfile, ifoList, ifoNames, gps, snrList, info)

# --------------------------------------------------
# Define helper function to extract median waveform
# --------------------------------------------------
def get_median(strain_master, ifo):    
    median_waveform = np.median(strain_master, axis=0)
    down_vec = np.zeros(median_waveform.size)
    up_vec = np.zeros(median_waveform.size)
    for sample in range(median_waveform.size):
        vector = np.transpose(strain_master)[sample]
        sort_vec = np.sort(vector)
        up_quart = sort_vec[ int(0.95*sort_vec.size) ]
        down_quart = sort_vec[ int(0.05*sort_vec.size) ]
        down_vec[sample] = down_quart
        up_vec[sample] = up_quart            
    return (median_waveform, up_vec, down_vec)

# ---------------------------------------------
# Extract median waveform from big output file
# ---------------------------------------------
def get_waveform_bigfile(filename, ifo):
    strain_master = np.loadtxt(filename)
    if len(strain_master.shape) == 1:
        median_waveform = strain_master
        down_vec = strain_master 
        up_vec = strain_master
    else:
        median_waveform, up_vec, down_vec = get_median(strain_master, ifo)
    return (median_waveform, up_vec, down_vec)

# ---------------------------------------------
# Get Median waveform with CIs
# ---------------------------------------------
def get_waveform(filename):
    names = ['samples','median_waveform','50low','50high','90low','90high']
    data = np.recfromtxt(filename,names=names)
    return (data['samples'],data['median_waveform'],data['50low'],data['50high'],data['90low'],data['90high'])



# ------------------------------------------------------------
# Calculate the mode and confidence interval from a histogram
# ------------------------------------------------------------
def get_mode(n, bins):
    # -- Take average of bins to get rid of size mismatch
    if (bins.size) == (n.size + 1):
        new_bins = [0.5*(bins[i] + bins[i+1]) for i in range(0,n.size)]
        bins = new_bins
    # -- Set the index of the mode to the tallest peak
    modeindx = n.argmax()
    samps_90 = 0.90*n.sum()
    
    for i in np.arange(n.size):
        lowindx = np.max([0, modeindx-i])
        highindx = np.min([n.size-1, modeindx+i])
        n_inrange = n[lowindx:highindx+1]
        if n_inrange.sum() >= samps_90:
            low = bins[lowindx]
            high = bins[highindx]
            mode = bins[modeindx]
            break
    return(mode, low, high)
    

######################################################################################################################
#
# Plotting functions
#
######################################################################################################################
#....................................................
# New function to get SNR^2(t) for x-axis
#....................................................
def snrfunctime(median_waveform):
    powerList = []
    wave_integral = median_waveform
    time_integrate = time
    #integral at each
    dt = 0.0009765
    h_t_2 = []
    for line, line_1, time_i in zip(wave_integral, wave_integral[1:], time_integrate):
        t = time_i
        w0 = line
        h_t_2.append(line**2)
        w1  = line_1
        snr_t = (((w0**2)+(w1**2))/2)*dt
        powerList.append([snr_t,t])
    
    #plot snr^2(t)
#    power, time_int = zip(*powerList)
#    plt.plot(time_int, power)
#    plt.title("SNR^2(t) for axes determination")
#    plt.xlabel('Time (s)')
#    plt.ylabel('{0} SNR^2'.format(model))
#    plt.title('IFO {0}'.format(ifo))
#    plt.axis([-1.0, 1.0, -0.1, np.max(power)+15])
#    plt.savefig(plotsDir+'{1}_SNR_TIME_{0}.png'.format(ifo, model))

    return(powerList)

# -------------------------------------------------
# Read in data and median waveform to get plot axis
# Now using SNR^2 for time axis determination
# -------------------------------------------------
def get_axes(jobName, postDir, ifoList, worc, model, time):
    axisList = []
    for ifo in ifoList:
        # -- Read Signal model
        try:
            filename = str(jobName)+postDir+'{1}_recovered_{2}_waveform.dat.{0}'.format(ifo, model, worc)
            median_waveform, up_vec, down_vec = get_waveform_bigfile(filename, ifo)
        except:
            filename = str(jobName)+postDir+'{1}_median_time_domain_waveform.dat.{0}'.format(ifo, model)
            timesamp, median_waveform, up_vec, down_vec = get_waveform(filename, ifo)

#        filename = str(jobName)+postDir+'{1}_median_time_domain_waveform.dat.{0}'.format(ifo, model)
#        timesamp, median_waveform, up_vec, down_vec, up_vec2, down_vec2 = get_waveform(filename)

        # -- Get axis info
        snr_t = snrfunctime(median_waveform)
        #y
        wave = up_vec
        wave_max = wave.max()
        #x
        wave = median_waveform
        power, and_time = zip(*snr_t)
        power = list(power)
        sig_times = []
        #Fixed ~90$ rep. from 5% to 95%
        for row in snr_t:
            if (row[0] >= 0.001*np.max(power) and row[0] <= 0.999*np.max(power)):
                sig_times.append(row[1])
        axisList.append([np.min(sig_times), np.max(sig_times)+0.1*(np.max(sig_times)-np.min(sig_times)), -wave_max*1.1, wave_max*1.1])
    
    
    # -- Select Axis with biggest y-value
    ymax = 0
    for axcand in axisList:
        if axcand[3] > ymax:
            ymax = axcand[3]
            axwinner = axcand
    return(axwinner)


def get_axes_fdom(jobName, postDir, ifoList, worc, model, time):
    
    ymin = np.zeros(len(ifoList))
    ymax = np.zeros(len(ifoList))
    
    for ifo in ifoList:
        # -- Read Signal model
    
        names = ['f','dpower','rpower','psd']
    
        filename = str(jobName)+postDir+'gaussian_noise_model_ifo{0}.dat'.format(ifo)
        
        data = np.recfromtxt(filename,names=names)
        
        # figure out what fmin is-- psd padded with 1's until fmin
        imin = 0
        while data['psd'][imin] == 1:
            imin += 1
        
        ymin[int(ifo)] = min(data['rpower'][imin:])
        ymax[int(ifo)] = max(data['rpower'][imin:])

        xmin = data['f'][imin]
        xmax = data['f'][-1]
    

    axwinner = [xmin, xmax, min(ymin), max(ymax)]
    return(axwinner)


# --------------
# Plot Evidence
# --------------
def plot_evidence(jobName, plotsDir):

    sig_noise=0
    sig_gl=0
    sig_si=0
    err_sig_noise=0
    err_sig_gl=0
    err_sig_si=0

    # -- Read evidence data
    try:
        infile = open(str(jobName)+'evidence.dat', 'r')
        maxodds = 20
        for line in infile:
            spl = line.split()
            if spl[0] == 'noise':
                sig_noise = float(spl[1])
                err_sig_noise = float(spl[2])
                #evidence.dat files have alread store variance (not std dev)
                #err_sig_noise *= err_sig_noise
            if spl[0] == 'glitch':
                sig_gl = float(spl[1])
                err_sig_gl = float(spl[2])
                #err_sig_gl *= err_sig_gl
            if spl[0] == 'signal':
                sig_si = float(spl[1])
                err_sig_si = float(spl[2])
                #err_sig_si *= err_sig_si
        infile.close()
    except:
        sig_noise=0
        sig_gl=0
        sig_si=0
        err_sig_noise=0
        err_sig_gl=0
        err_sig_si=0

    sig_noise = sig_si - sig_noise
    sig_gl    = sig_si - sig_gl

    err_sig_noise += err_sig_si
    err_sig_noise = math.sqrt(err_sig_noise)
    err_sig_gl += err_sig_si
    err_sig_gl = math.sqrt(err_sig_gl)
    # -- Report to the user 
    print
    print '   log( E_signal / E_noise ) =', sig_noise
    print '   log( E_signal / E_glitch ) =', sig_gl
    # -- Plot the data point 
    plt.figure()
    plt.errorbar(sig_gl, sig_noise, 2*err_sig_gl, 2*err_sig_noise, color='black')
    # -- Store maxima and establish axes 
    maxodds = 1.1*np.array( [np.abs(sig_gl)+2*err_sig_gl, np.abs(sig_noise)+2*err_sig_noise, 20] ).max()
    xmaxodds = 1.1*np.maximum(np.abs(sig_gl)+2*err_sig_gl, 20)
    ymaxodds = 1.1*np.maximum(np.abs(sig_noise)+2*err_sig_noise, 20)
    plt.axis([-xmaxodds, xmaxodds, -ymaxodds, ymaxodds])
    # -- Color in the plot 
    plt.fill_between([0,maxodds], [0, 0], [maxodds, maxodds], facecolor=scolor, interpolate=True, alpha=0.3)
    plt.fill_between([-maxodds,0,maxodds], [-maxodds,0,0], [-maxodds, -maxodds, -maxodds], facecolor=ncolor, interpolate=True, alpha=0.3)
    plt.fill_between([-maxodds,0], [-maxodds, 0], [maxodds, maxodds], facecolor=gcolor, interpolate=True, alpha=0.3)
    plt.grid()
    # -- Labels on the plot
    plt.text(0.9*xmaxodds, 0.9*ymaxodds, 'Signal', horizontalalignment='right', verticalalignment='top')
    plt.text(-0.9*xmaxodds, 0.9*ymaxodds, 'Glitch', horizontalalignment='left', verticalalignment='top')
    plt.text(0.9*xmaxodds, -0.9*ymaxodds, 'Noise', horizontalalignment='right', verticalalignment='bottom')
    # -- Final touches 
    plt.xlabel(r'LN($B_{signal/glitch}$)')
    plt.ylabel(r'LN($B_{signal/noise}$)')
    plt.savefig(plotsDir+'odds.png')
    plt.close()
    
    
    # -- Return values to be used later
    return(sig_gl, sig_noise, err_sig_gl, err_sig_noise)

# -----------------------------------------------
# Plot the median waveform and injected waveform
# -----------------------------------------------
def plot_waveform(jobName, postDir, ifo, plotsDir, worc, mdc, model, axwinner, time, low_50, high_50, low_90, high_90, median_waveform):
    plt.figure()

    try:
        filename = str(jobName)+postDir+'{1}_data.dat.{0}'.format(ifo, worc)
        inj_median_waveform, inj_up_vec, inj_down_vec = get_waveform_bigfile(filename, ifo)
        plt.plot(time, inj_median_waveform, color = '0.75', linewidth=2, alpha=0.25)
    except:
        print "I couldn't find the file {0}".format(filename)

    if model == 'glitch':
       colour = gcolor
    elif model == 'signal':
       colour = scolor

    plt.fill_between(time, low_50, high_50, facecolor=colour, edgecolor=colour, alpha=0.5)
    plt.fill_between(time, low_90, high_90, facecolor=colour, edgecolor=colour, alpha=0.3)

    if mdc:
        # -- Read in the injected waveform
        try:
            filename = str(jobName)+postDir+'injected_{1}_waveform.dat.{0}'.format(ifo, worc)
            inj_median_waveform, inj_up_vec, inj_down_vec = get_waveform_bigfile(filename, ifo)    
        except:
            filename = str(jobName)+postDir+'injected_waveform.dat.{0}'.format(ifo)
            inj_median_waveform, inj_up_vec, inj_down_vec = get_waveform_bigfile(filename, ifo)    
        inj_time = time
        plt.plot(inj_time, inj_median_waveform, injcolor, linewidth=1)

    plt.plot(time, median_waveform, color=colour, linewidth=1, alpha=1)

    #axisList.append([sig_times.min(), sig_times.max(), -wave_max*1.1, wave_max*1.1])
    plt.xlabel('Time (s)')
    plt.ylabel('Whitened strain'.format(model,colour))
    plt.title('Reconstructed {0} model in {1}'.format(model,ifoNames[int(ifo)]))
    
    # -- Save the full versions of the plot
    if worc == 'whitened':
        plt.savefig(plotsDir+'{1}_waveform_{0}_full.png'.format(ifo, model))
    else:
        plt.savefig(plotsDir+'c{1}_waveform_{0}_full.png'.format(ifo, model))

    plt.axis(axwinner)
    # -- Save the plot
    if worc == 'whitened':
        plt.savefig(plotsDir+'{1}_waveform_{0}.png'.format(ifo, model))
    else:
        plt.savefig(plotsDir+'c{1}_waveform_{0}.png'.format(ifo, model))

    plt.close()


def plot_power_spectrum(jobName, postDir, ifo, plotsDir, worc, mdc, model, axwinner, freq, low_50, high_50, low_90, high_90, median_waveform,type):
    plt.figure()
    
    # To do: data and injected waveforms if possible
    
    if model == 'glitch':
        colour = gcolor
    elif model == 'signal':
        colour = scolor


    if mdc:
        filename = st(jobName)+postDir+'injected_whitened_spectrum.dat.{0}'.format(ifo)

        injected_spectrum = np.genfromtxt(filename)

        plt.semilogy(injected_spectrum[:,0],injected_spectrum[:,1],color = '0.75', linewidth=2, alpha=0.25)



    plt.fill_between(freq, low_50, high_50, facecolor=colour, edgecolor=colour, alpha=0.5)
    plt.fill_between(freq, low_90, high_90, facecolor=colour, edgecolor=colour, alpha=0.3)
    
    if mdc:
        # -- Read in the injected waveform
        try:
            filename = str(jobName)+postDir+'injected_{1}_waveform.dat.{0}'.format(ifo, worc)
            inj_median_waveform, inj_up_vec, inj_down_vec = get_waveform_bigfile(filename, ifo)
        except:
            filename = str(jobName)+postDir+'injected_waveform.dat.{0}'.format(ifo)
            inj_median_waveform, inj_up_vec, inj_down_vec = get_waveform_bigfile(filename, ifo)
        inj_time = time
        plt.plot(inj_time, inj_median_waveform, injcolor, linewidth=1)

    if type == 'psd':
        filename = str(jobName)+postDir+'gaussian_noise_model_ifo{0}.dat'.format(ifo)
        data = np.recfromtxt(filename, names = ['f', 'dpower', 'rpower', 'Sn'])
        plt.semilogy(data['f'],data['rpower'], ncolor, alpha=0.6)
        plt.ylim(min(data['rpower']),max(data['rpower']))


    plt.semilogy(freq, median_waveform, color=colour, linewidth=1, alpha=1)
    
    #axisList.append([sig_times.min(), sig_times.max(), -wave_max*1.1, wave_max*1.1])
    plt.xlabel('Frequency (Hz)')
    if type == 'psd':
        plt.ylabel('PSD and data (grey)'.format(model))
    if type == 'powerspec':
         plt.ylabel('Whitened power spectrum of median waveform'.format(model))
    plt.title('{0}'.format(ifoNames[int(ifo)]))



    # -- Save the full versions of the plot
    if worc == 'whitened':
        plt.savefig(plotsDir+'{1}_{2}_{0}_full.png'.format(ifo, model, type))
    else:
        plt.savefig(plotsDir+'c{1}_{2}_{0}_full.png'.format(ifo, model, type))


    plt.xlim(32,512)
    
#plt.axis(axwinner)
    # -- Save the plot
    if worc == 'whitened':
        plt.savefig(plotsDir+'{1}_{2}_{0}.png'.format(ifo, model, type))
    else:
        plt.savefig(plotsDir+'c{1}_{2}_{0}.png'.format(ifo, model, type))
    
    plt.close()

def plot_full_spectro(jobName, postDir, ifo, plotsDir, worc, mdc, model, axwinner, psd_info, powerspec_info, median_waveform):
    
    plt.figure()

    if model == 'glitch':
        colour = gcolor
    elif model == 'signal':
        colour = scolor


    if mdc:
        filename = str(jobName)+postDir+'injected_whitened_spectrum.dat.{0}'.format(ifo)

        injected_spectrum = np.genfromtxt(filename)

        plt.semilogy(injected_spectrum[:,0],injected_spectrum[:,1],injcolor, linewidth=1)


    names = ['f','dpower','rpower','psd']
    filename = str(jobName)+postDir+'gaussian_noise_model_ifo{0}.dat'.format(ifo)
    data = np.recfromtxt(filename,names=names)

    # plot data
    plt.semilogx(data['f'],data['rpower'],color=ncolor,alpha=0.5)

    # plot psd
    plt.fill_between(psd_info[0],psd_info[4],psd_info[5],color='grey',alpha=0.8)
    plt.semilogy(psd_info[0],psd_info[1],color='k',ls='-')

    # plot powerspec
    plt.fill_between(powerspec_info[0],powerspec_info[2],powerspec_info[3],color=colour,alpha=0.5)
    plt.fill_between(powerspec_info[0],powerspec_info[4],powerspec_info[5],color=colour,alpha=0.3)
    plt.plot(powerspec_info[0],powerspec_info[1],color=colour)

    plt.xscale('log')
    plt.yscale('log')

    plt.ylim(axwinner[2],axwinner[3])
    plt.xlim(axwinner[0],axwinner[1])

    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Strain")
    plt.title("Power Spectra for {0} model in {1}".format(model,ifoNames[int(ifo)]))


    plt.savefig(plotsDir+'{0}_frequence_domain_{1}.png'.format(model,ifo))

    plt.close()



def plot_tf_tracks(jobName, postDir, ifo, plotsDir, worc, mdc, model, axwinner, f_axwinner,time, low_50, high_50, low_90, high_90, median_waveform):
    plt.figure()
    
    # To do: data and injected waveforms if possible
    
    if model == 'glitch':
        colour = gcolor
    elif model == 'signal':
        colour = scolor

    t = time - 2.0


    plt.fill_between(t, low_50, high_50, facecolor=colour, edgecolor=colour, alpha=0.5)
    plt.fill_between(t, low_90, high_90, facecolor=colour, edgecolor=colour, alpha=0.3)
    
    if mdc:
        # -- Read in the injected waveform
        try:
            filename = str(jobName)+postDir+'injected_{1}_waveform.dat.{0}'.format(ifo, worc)
            inj_median_waveform, inj_up_vec, inj_down_vec = get_waveform_bigfile(filename, ifo)
        except:
            filename = str(jobName)+postDir+'injected_waveform.dat.{0}'.format(ifo)
            inj_median_waveform, inj_up_vec, inj_down_vec = get_waveform_bigfile(filename, ifo)
        inj_time = time
        plt.plot(inj_time, inj_median_waveform, injcolor, linewidth=1)
    
    
    plt.plot(t, median_waveform, color=colour, linewidth=1, alpha=1)
    
    #axisList.append([sig_times.min(), sig_times.max(), -wave_max*1.1, wave_max*1.1])
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.title('Time vs. Frequency for {0} model in {1}'.format(model,ifoNames[int(ifo)]))
    
    
    
    # -- Save the full versions of the plot
    if worc == 'whitened':
        plt.savefig(plotsDir+'{1}_tf_{0}_full.png'.format(ifo, model))
    else:
        plt.savefig(plotsDir+'c{1}_tf_{0}_full.png'.format(ifo, model))
    
    
    plt.xlim(axwinner[0],axwinner[1])
    plt.ylim(f_axwinner[0],f_axwinner[1])
    
    #plt.axis(axwinner)
    # -- Save the plot
    if worc == 'whitened':
        plt.savefig(plotsDir+'{1}_tf_{0}.png'.format(ifo, model))
    else:
        plt.savefig(plotsDir+'c{1}_tf_{0}.png'.format(ifo, model))
    
    plt.close()






def plot_cross_waveform(jobName, postDir, ifo, plotsDir, worc, mdc, model, axwinner, time, down_vec, up_vec, median_waveform):
    plt.figure()

    plt.fill_between(time, down_vec, up_vec, facecolor=scolor, edgecolor=scolor, alpha=0.4)
    if mdc:
        # -- Read in the injected waveform
        try:
            filename = str(jobName)+postDir+'injected_cross_{1}_waveform.dat.{0}'.format(ifo, worc)
            inj_median_waveform, inj_up_vec, inj_down_vec = get_waveform_bigfile(filename, ifo)    
        except:
            filename = str(jobName)+postDir+'injected_cross_waveform.dat.{0}'.format(ifo)
            inj_median_waveform, inj_up_vec, inj_down_vec = get_waveform_bigfile(filename, ifo)    
        inj_time = time
        plt.plot(inj_time, inj_median_waveform, injcolor, linewidth=1)

    plt.plot(time, median_waveform, scolor, linewidth=1, alpha=0.75)

    #axisList.append([sig_times.min(), sig_times.max(), -wave_max*1.1, wave_max*1.1])
    plt.xlabel('Time (s)')
    plt.ylabel('Best Fit {0} Waveform (blue) and injected waveform (red)'.format(model))
    plt.title('{0}'.format(ifoNames[int(ifo)]))
    
    # -- Save the full versions of the plot
    if worc == 'whitened':
        plt.savefig(plotsDir+'{1}_cross_waveform_{0}_full.png'.format(ifo, model))
    else:
        plt.savefig(plotsDir+'c{1}_cross_waveform_{0}_full.png'.format(ifo, model))

    plt.axis(axwinner)
    # -- Save the plot
    if worc == 'whitened':
        plt.savefig(plotsDir+'{1}_cross_waveform_{0}.png'.format(ifo, model))
    else:
        plt.savefig(plotsDir+'c{1}_cross_waveform_{0}.png'.format(ifo, model))

    plt.close()

def plot_plus_waveform(jobName, postDir, ifo, plotsDir, worc, mdc, model, axwinner, time, down_vec, up_vec, median_waveform):
    plt.figure()

    plt.fill_between(time, down_vec, up_vec, facecolor=scolor, edgecolor=scolor, alpha=0.4)
    if mdc:
        # -- Read in the injected waveform
        try:
            filename = str(jobName)+postDir+'injected_plus_{1}_waveform.dat.{0}'.format(ifo, worc)
            inj_median_waveform, inj_up_vec, inj_down_vec = get_waveform_bigfile(filename, ifo)    
        except:
            filename = str(jobName)+postDir+'injected_plus_waveform.dat.{0}'.format(ifo)
            inj_median_waveform, inj_up_vec, inj_down_vec = get_waveform_bigfile(filename, ifo)    
        inj_time = time
        plt.plot(inj_time, inj_median_waveform, injcolor, linewidth=1)

    plt.plot(time, median_waveform, scolor, linewidth=1, alpha=0.75)

    #axisList.append([sig_times.min(), sig_times.max(), -wave_max*1.1, wave_max*1.1])
    plt.xlabel('Time (s)')
    plt.ylabel('Best Fit {0} Waveform (blue) and injected waveform (red)'.format(model))
    plt.title('{0}'.format(ifoNames[int(ifo)]))
    
    # -- Save the full versions of the plot
    if worc == 'whitened':
        plt.savefig(plotsDir+'{1}_plus_waveform_{0}_full.png'.format(ifo, model))
    else:
        plt.savefig(plotsDir+'c{1}_plus_waveform_{0}_full.png'.format(ifo, model))

    plt.axis(axwinner)
    # -- Save the plot
    if worc == 'whitened':
        plt.savefig(plotsDir+'{1}_plus_waveform_{0}.png'.format(ifo, model))
    else:
        plt.savefig(plotsDir+'c{1}_plus_waveform_{0}.png'.format(ifo, model))

    plt.close()

# ----------------------
# Create the Bayesogram
# ----------------------
def bayesogram(filename, plotsDir, worc, time, twodrange, median_waveform, axwinner, ifo, model):
    strain_master = np.loadtxt(filename)
    multiply = strain_master.size / time.size
    X = list(time)*multiply
    Y = strain_master.flatten()
    # twodrange = [timerange, [Y.min(), Y.max()]]
    timerange = twodrange[0]
    Nx = np.abs( (timerange[1]-timerange[0])/(time[1] - time[0]) )
    H, xedges, yedges = np.histogram2d(X,Y,bins=[Nx, 100],range=twodrange)

    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
  
    # Plot 2D histogram using pcolor
    plt.figure()
    plt.pcolormesh(xedges,yedges,H,cmap='OrRd')
    # cbar = plt.colorbar(orientation='horizontal')
    # cbar.ax.set_ylabel('Counts')
    plt.plot(time, median_waveform, 'w', linewidth=1, alpha=0.5)
    plt.xlabel('Time (s)')
    plt.ylabel('h(t)')
    plt.title('{0} - {1} model'.format(ifoNames[int(ifo)],model))
    plt.axis(axwinner)
    if worc == 'whitened':
        plt.savefig(plotsDir+'{1}_t_bgram_{0}.png'.format(ifo, model))
    else:
        plt.savefig(plotsDir+'c{1}_t_bgram_{0}.png'.format(ifo, model))

    plt.close()

# --------------------------------------
# Plot Q scans of waveforms
# --------------------------------------
def Q_scan(model, Q, t, f, ifo, axwinner, f_axwinner, climwinner=[1,1]):
    data = np.genfromtxt(postDir+'/{0}_spectrogram_{1}.dat.{2}'.format(model, Q, ifo))
    data = data[...,::-1] # gets printed out time reversed, need to fix it

#plt.clim(np.min(data)/np.sum(data),np.max(data)/np.sum(data))

    fig = plt.figure()
    plt.imshow(data,aspect='auto',origin='lower',extent=[t[0],t[-1],f[0],f[-1]],cmap='OrRd')
    
    if model == 'data':
        plt.clim(np.min(data),np.max(data))
    
    if model == 'residual':
        plt.clim(climwinner)
    
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.title('Spectrogram of median {0} waveform in {1}, Q={2}'.format(model,ifoNames[int(ifo)], Q))
    
    plt.xlim(axwinner[0],axwinner[1])
    plt.ylim(f_axwinner[0],f_axwinner[1])

    plt.savefig(plotsDir+'{0}_spectrogram_Q{1}_{2}.png'.format(model,Q,ifo))

    plt.close()

    return [np.min(data),np.max(data)]



# --------------------------------------
# Make a spectrogram of median waveform
# --------------------------------------
def waveform_spectrogram(ifo, plotsDir, worc, time, median_waveform, model):
    ts = time[1] - time[0]
    fs = int(np.round(1.0/ts))
    NFFT = fs / 16
    window = np.blackman(NFFT)
    holder = np.zeros(NFFT)
    spec_power, freqs, bins = mlab.specgram(median_waveform, NFFT=NFFT, Fs=fs,
                                            window=window,
                                            noverlap=int(0.8*NFFT))

    # -- Try eliminating quiet pixels
    sort_power = np.sort(spec_power.flatten())[::-1]
    sum_power = np.cumsum(sort_power)
    indx = np.nonzero( sum_power < 0.95*sum_power.max() )
    try:
        point90 = sort_power[indx[0].max()]
        spec_power[ spec_power < point90 ] = spec_power.min()
        vmin = np.log10(spec_power[ spec_power > spec_power.min()].min())
        vmax = np.log10(spec_power.max())
    except: # h(t) = 0
        point90 = sum_power.max() 
        spec_power[ spec_power >= point90 ] = 0.1
        vmin = -1.0 
        vmax = 0.0

    # -- Try do it yourself spectrogram
    fig = plt.figure()
    ax = plt.subplot(111)
    binwidth = (bins[1] - bins[0]) / 2.0
    # Adjust time axis: trigger ends up at segment length-2.0s, so we must shift back to place the trigger at t=0s
    tshift = np.abs(time[0])+binwidth
    cax = ax.pcolormesh(bins-tshift, freqs, np.log10(spec_power), vmin=vmin, vmax=vmax, cmap='OrRd')
    plt.title('{0} - {1} model'.format(ifoNames[int(ifo)],model))
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.axis([time[0], time[-1], 0, fs/2])
    cb = plt.colorbar(cax, orientation='horizontal')
    cb.set_label('Logrithmic Signal Power')
    if worc == 'whitened':
        plt.savefig(plotsDir+'{1}_self_spec_{0}_full.png'.format(ifo, model))
    else:
        plt.savefig(plotsDir+'c{1}_self_spec_{0}_full.png'.format(ifo, model))

    plt.axis([axwinner[0], axwinner[1], 0, fs/2])
    if worc == 'whitened':
        plt.savefig(plotsDir+'{1}_self_spec_{0}.png'.format(ifo, model))
    else:
        plt.savefig(plotsDir+'c{1}_self_spec_{0}.png'.format(ifo, model))

    plt.close()

# -----------
# Histograms
# -----------
def make_histogram(moments, moment, ifo, plotsDir, worc, outfile, xtype, ytype):
  colour    = scolor
  injcolour = injcolor

  # -- Drop '_plus' and '_cross' as these do not appear in the moments' name list
  mom = moment.replace('_plus','')
  mom = mom.replace('_cross','')
  # -- y scale is linear by default, otherwise log10
  ylog = False
  if ytype == 'log10':
    ylog = True
  # -- Start figure
  plt.figure()
  # -- Data from recovered signal
  if xtype == 'log10':
    rn,rbins,patches = plt.hist(np.log10(moments[mom]), bins=50, label='Recovered', alpha=0.5, log=ylog, linewidth=0, color=colour)
  if xtype == 'lin':
    if moment == 't0_rec':
      rn,rbins,patches = plt.hist(moments[mom] - 2, bins=50, label='Recovered', alpha=0.5, log=ylog, linewidth=0, color=colour)
    elif 'overlap' in moment:
      rbins = np.arange(-1.02,1.02,0.01)
      rn,rbins,patches = plt.hist(moments[mom], rbins, alpha=0.5, log=ylog, linewidth=0, color=colour)
      plt.xlim(-1.1, 1.1)
    else:
      rn,rbins,patches = plt.hist(moments[mom], bins=50, label='Recovered', alpha=0.5, linewidth=0, log=ylog, color=colour)
  mode, low, high = get_mode(rn, rbins)
  # -- Overlay data relative to the injected signal
  if 'overlap' not in moment:
    if xtype == 'log10':
      inmode = np.median(np.log10(injmoments[mom])) #TODO: fix header in BayesWavePost
    if xtype == 'lin':
      if moment == 't0_rec':
        inmode = np.median(injmoments[mom] - 2) # Correcting for the 2s added by BWB to the data time window
      else:
        inmode = np.median(injmoments[mom])
    # -- Plot the histogram
    if mdc:
      plt.axvline(x=inmode, label='Injected', color=injcolour)
      plt.legend()
    plt.ylim(1, 2*rn.max())
  if not mom == 'network_overlap':
    plt.title('Whitened {0} data'.format(ifoNames[int(ifo)]))
  plt.xlabel(moments_dict[moment][0])
  plt.grid()
  if worc == 'whitened':
    if not mom == 'network_overlap':
      plt.savefig(plotsDir+moments_dict[moment][1]+'_{1}_{0}.png'.format(ifo,ytype))
    else:
      plt.savefig(plotsDir+moments_dict[moment][1]+'_{0}.png'.format(ytype))
  else:
    if not mom == 'network_overlap':
      plt.savefig(plotsDir+'c'+moments_dict[moment][1]+'_{1}_{0}.png'.format(ifo,ytype))
    else:
      plt.savefig(plotsDir+'c'+moments_dict[moment][1]+'_{0}.png'.format(ytype))
  plt.close()
  # -- Write modes to file and report to the user (this is done only once)
  if ylog==False:
    outfile.write("{0} {1} {2} \n".format(mode, low, high))
    if 'overlap' not in moment:
      outfile.write("{0} {0} {0}\n".format(inmode))

# -------------------------------------
# Return recovered and injected modes
# -------------------------------------
def mode_values(moment, ifo, mdc, momentsPlus, momentsCross, ifoList, tablesDir, worc):
    # -- Initialize variables to be returned
    mode, low, high, inmode = ['-', '-', '-', '-']
    # -- Open and read modes table
    if worc == 'whitened':
        tablesrc = './'+tablesDir+'mode_{0}.txt'.format(ifo)
    elif worc == 'colored':
        tablesrc = './'+tablesDir+'cmode_{0}.txt'.format(ifo)
    table_data = np.loadtxt(tablesrc)
    # -- Single detector overlap
    if moment == 'overlap':
        mode, low, high = table_data[0]
    # -- Single detector h+ overlap 
    elif moment == 'overlap_plus':
        if len(momentsPlus) > 0:
            mode, low, high = table_data[1]
    # -- Single detector hx overlap 
    elif moment == 'overlap_cross':
        if len(momentsCross) > 0:
            mode, low, high = table_data[2]
    # -- Detector network overlap
    elif moment == 'network_overlap':
        if len(ifoList)>1:
            mode, low, high = table_data[3]
    # -- Detector network h+ overlap 
    elif moment == 'network_overlap_plus':
        if len(ifoList)>1 and len(momentsPlus) > 0:
            mode, low, high = table_data[4]
    # -- Detector network hx overlap 
    elif moment == 'network_overlap_cross':
        if len(ifoList)>1 and len(momentsCross) > 0:
            mode, low, high = table_data[5]
    # -- Signal energy
    elif moment == 't_energy_rec':
        mode, low, high = table_data[6]
        if mdc:
            inmode = table_data[7][0]
    # -- Central time
    elif moment == 't0_rec':
        mode, low, high = table_data[8]
        if mdc:
            inmode = table_data[9][0]
    # -- Duration
    elif moment == 'dur_rec':
        mode, low, high = table_data[10]
        if mdc:
            inmode = table_data[11][0]
    # -- Central frequency
    elif moment == 'f0_rec':
        mode, low, high = table_data[12]
        if mdc:
            inmode = table_data[13][0]
    # -- Bandwidth
    elif moment == 'band_rec':
        mode, low, high = table_data[14]
        if mdc:
            inmode = table_data[15][0]
    # -- Maximum value of ||h(t)|| 
    elif moment == 'h_max':
        mode, low, high = table_data[16]
        if mdc:
            inmode = table_data[17][0]
    # -- Time at maximum value of ||h(t)|| 
    elif moment == 't_at_h_max':
        mode, low, high = table_data[18]
        if mdc:
            inmode = table_data[19][0]
    elif moment == 'snr':
        mode, low, high = table_data[20]

    # -- Return recovered mode value, 90% credible interval limites, and injected mode
    return(mode, low, high, inmode)

# -----------------
# Plot Temp Chains
# -----------------
def plot_likelihood_1(modelList, plotsDir):
    plt.figure()
    plt.title('Likelihood')
    plt.xlim([1e-3, 1])
    for mod in modelList:
        if mod == 'glitch':
            colour = gcolor
        elif mod == 'signal':
            colour = scolor
        else:
            colour = ncolor
        try:
            names = ['temp','likehood','error','acl']
            data = np.recfromtxt(str(jobName)+"{0}_evidence.dat".format(mod), names=names)
        except:
            try:
               names = ['temp','likehood','error']
               data = np.recfromtxt(str(jobName)+"{0}_evidence.dat".format(mod), names=names)
            except:
                continue
        #error = np.zeros(likehood.shape)
        plt.semilogx(data['temp'], data['likehood'], label=mod, linewidth=2, color=colour)
        plt.errorbar(data['temp'], data['likehood'], 2*data['error'], color=colour)

    plt.xlabel( '1/Temp' )
    plt.ylabel('log(L)')
    plt.grid()
    plt.legend(loc=2)
    plt.savefig(plotsDir+'likelihood.png')
    plt.close()

def plot_likelihood_2(modelList, plotsDir):
    plt.figure()
    plt.title('Likelihood')
    for mod in modelList:
        if mod == 'glitch':
            colour = gcolor
        elif mod == 'signal':
            colour = scolor
        else:
            colour = ncolor
        try:
            names = ['temp','likehood','error','acl']
            data = np.recfromtxt(str(jobName)+"{0}_evidence.dat".format(mod), names=names)
        except:
            try:
               names = ['temp','likehood','error']
               data = np.recfromtxt(str(jobName)+"{0}_evidence.dat".format(mod), names=names)
            except:
                continue
        plt.semilogx(data['temp'], data['likehood']*data['temp'], label=mod, linewidth=2, color=colour)
        plt.errorbar(data['temp'], data['likehood']*data['temp'], 2*data['error'], color=colour)
    plt.grid()
    plt.xlabel('1/Temp')
    plt.ylabel('log(L) X 1/Temp')
    plt.legend(loc=2)
    plt.savefig(plotsDir+'TL.png')
    plt.close()

# -----------
# Plot PSD
# -----------
def plot_psd(jobName, postDir, ifoList, ifoNames, plotsDir):
  colorList = [H1color,L1color,V1color]
  for ifo in ifoList:
    filename = str(jobName)+postDir+'gaussian_noise_model_ifo{0}.dat'.format(ifo)
    data = np.recfromtxt(filename, names=['f', 'dpower', 'rpower', 'Sn'])
    f = data['f']
    plt.close('all')
    plt.figure()
    plt.semilogy(f, data['dpower'], ncolor, alpha=0.6)
    plt.semilogy(f, data['rpower'], colorList[int(ifo)], alpha=0.6)
    plt.semilogy(f, data['Sn'], color='black', alpha=1, linewidth=2.5)
    plt.title('PSD for {0}'.format(ifoNames[int(ifo)]))
    plt.xlim(f.min(), f.max())
    plt.xlabel('Freq (Hz)')
    plt.ylabel('PSD')
    plt.grid()
    plt.savefig(os.path.join(plotsDir, 'psd{0}.png'.format(ifo)))

# -------------------------
# Plot dimensions of model
# -------------------------
def plot_model_dims(modelList, ifoList, ifoNames, plotsDir):
    
  lineStyles = ['-', '--', ':']
  lineColors = [H1color, L1color, V1color]
  glitchChains = []
  signalChains = []
  # -- Read in data
  for mod in modelList:

    intChainFile = 'chains/'+str(jobName)+"{0}_model.dat.0".format(mod)
    N = len(open(intChainFile).readlines())
    chains = np.loadtxt(intChainFile,skiprows=N/2)
    chains = np.transpose(chains)
    if mod == 'glitch':
      glitchChains = chains
    elif mod == 'signal':
      signalChains = chains

  # -- Variation of the model dimension over MCMC iterations. 2 subplots for signal and glitch models

  if len(modelList) == 2:
        fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
        ax1.set_title('Model Dimension')
        if not signalChains == []:
            ax1.plot(signalChains[2], linewidth=2, color=scolor, label='Signal')
        if not glitchChains == []:
            for ifo in ifoList:
                ifoint=int(ifo)
                ax2.plot(glitchChains[3+ifoint], lineStyles[ifoint], color=lineColors[ifoint], linewidth=2, label=ifoNames[ifoint])

        # -- Legend placement outside the plot
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width*0.8, box.height])
        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # -- Make subplots close to each other and hide x ticks for all but bottom plot
        fig.subplots_adjust(hspace=0.1, right=0.8)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        plt.xlabel( 'MC step/10' )
        ax1.set_ylabel('Signal Model')
        ax1.grid()
        ax2.set_ylabel('Glitch Model')
        ax2.grid()
        plt.savefig(plotsDir+'model_dimensions.png')
        plt.close()

  elif len(modelList) == 1:
        fig = plt.figure()
        ax1 = plt.subplot(111)
        ax1.set_title('Model Dimension')
        if not signalChains == []:
            ax1.plot(signalChains[2], linewidth=2, color=scolor, label='Signal')
        if not glitchChains == []:
            for ifo in ifoList:
                ifoint=int(ifo)
                ax1.plot(glitchChains[3+ifoint], lineStyles[ifoint], color=lineColors[ifoint], linewidth=2, label=ifoNames[ifoint])
        
        # -- Legend placement outside the plot
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # -- Make subplots close to each other and hide x ticks for all but bottom plot
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        plt.xlabel( 'MC step/10' )
        ax1.set_ylabel('{0} model'.format(modelList[0]))
        ax1.grid()
        plt.savefig(plotsDir+'model_dimensions.png')
        plt.close()



  if len(modelList) == 2:
      # -- Histogram of the model dimension. 2 subplots for signal and glitch models
      if not signalChains == [] and not glitchChains == []:
        fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
      else:
        fig, (ax1, ax2) = plt.subplots(2, sharex=True)
      ax1.set_title('Model Dimension Histogram')
      if not signalChains == []:
        n,bins,patches = ax1.hist(signalChains[2], bins=range(int(min(signalChains[2])), int(max(signalChains[2]))+1, 1), histtype='bar', color=scolor, log=False)
      if not glitchChains == []:
        data = np.dstack(glitchChains[3:3+len(ifoList)])[0]
        #ax2.set_prop_cycle(cycler('color',['darkgoldenrod','darkkhaki','dkarsage']))
        n,bins,patches = ax2.hist(data, bins=range(int(data.min()), int(data.max())+1, 1), label=ifoNames, histtype='bar', log=False, color=[lineColors[int(i)] for i in ifoList])
      # -- Legend placement outside the plot
      box = ax1.get_position()
      ax1.set_position([box.x0, box.y0, box.width*0.8, box.height])
      box = ax2.get_position()
      ax2.set_position([box.x0, box.y0, box.width*0.8, box.height])
      ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
      # -- Make subplots close to each other and hide x ticks for all but bottom plot
      fig.subplots_adjust(hspace=0.1, right=0.8)
      plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
      plt.xlabel('Model Dimension')
      ax1.set_ylabel('Signal Model')
      ax1.grid()
      ax2.set_ylabel('Glitch Model')
      ax2.grid()
      plt.savefig(plotsDir+'model_dimensions_histo.png')

  elif len(modelList) == 1:
        # -- Histogram of the model dimension. 2 subplots for signal and glitch models
        fig = plt.figure()
        ax1 = plt.subplot(111)
        ax1.set_title('Model Dimension Histogram')
        if not signalChains == []:
            n,bins,patches = ax1.hist(signalChains[2], bins=range(int(min(signalChains[2])), int(max(signalChains[2]))+1, 1), histtype='bar', color=scolor, log=False)
        if not glitchChains == []:
            data = np.dstack(glitchChains[3:3+len(ifoList)])[0]
            #ax2.set_prop_cycle(cycler('color',['darkgoldenrod','darkkhaki','dkarsage']))
            n,bins,patches = ax1.hist(data, bins=range(int(data.min()), int(data.max())+1, 1), label=ifoNames, histtype='bar', log=False, color=[lineColors[int(i)] for i in ifoList])
        # -- Legend placement outside the plot
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # -- Make subplots close to each other and hide x ticks for all but bottom plot
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        plt.xlabel('Model Dimension')
        ax1.set_ylabel('{0} model'.format(modelList[0]))
        ax1.grid()
        plt.savefig(plotsDir+'model_dimensions_histo.png')



######################################################################################################################
#
# BWB webpage production
#
######################################################################################################################
def make_mode_table(page, subpage, imdc, momentsPlus, momentsCross, ifoList, tablesDir, ifoNames, worc):
    if page in ['dur_rec', 't_energy_rec', 'f0_rec', 'band_rec', 't0_rec', 'h_max', 't_at_h_max', 'snr']: # All histogram pages, excluding overlaps
        # -- Table header
        subpage.write('  <br/>\n')
        subpage.write('  <table id="tab1">\n')
        subpage.write('    <tr>\n')
        subpage.write('      <th>IFO</th>\n')
        subpage.write('      <th>Recovered value</th>\n')
        subpage.write('      <th>90% C.I. lower limit</th>\n')
        subpage.write('      <th>90% C.I. upper limit</th>\n')
        subpage.write('      <th>Injected value</th>')
        subpage.write('    </tr>\n')
        # -- Loop over IFOs 
        for ifo in ifoList:
            mode, low, high, inmode = mode_values(page, ifo, mdc, momentsPlus, momentsCross, ifoList, tablesDir, worc)
            subpage.write('    <tr>\n')
            subpage.write('      <td>{0}</td>\n'.format(ifoNames[int(ifo)]))
            for value in [mode, low, high, inmode]:
                subpage.write('      <td>{0}</td>\n'.format(value))
            subpage.write('    </tr>\n')
        subpage.write('  </table>\n')
    elif page == 'overlap': # Overlap histogram page
        # -- Table header
        subpage.write('  <br/>\n')
        subpage.write('  <table id="tab1">\n')
        subpage.write('    <tr>\n')
        subpage.write('      <th>Overlap</th>\n')
        subpage.write('      <th>Recovered value</th>\n')
        subpage.write('      <th>90% C.I. lower limit</th>\n')
        subpage.write('      <th>90% C.I. upper limit</th>\n')
        subpage.write('    </tr>\n')
        # -- Loop over whole signal, + and x polarizations
        for polarization in ['', '_plus', '_cross']:
            for ifo in ifoList: # Single detector
                plotsrc = './'+plotsDir+html_dict[page][2]+polarization+'_{1}_{0}.png'.format(ifo,'lin')
                if os.path.exists(plotsrc): # Indirect way of making sure this overlap was calculated 
                    mode, low, high, inmode = mode_values('overlap'+polarization, ifo, mdc, momentsPlus, momentsCross, ifoList, tablesDir, worc)
                    subpage.write('    <tr>\n')
                    if polarization == '':
                        subpage.write('      <td>{0} - Signal </td>\n'.format(ifoNames[int(ifo)]))
                    elif polarization == '_plus':
                        subpage.write('      <td>{0} - + Polarization </td>\n'.format(ifoNames[int(ifo)]))
                    else:
                        subpage.write('      <td>{0} - x Polarization </td>\n'.format(ifoNames[int(ifo)]))
                    for value in [mode, low, high]:
                        subpage.write('      <td>{0}</td>\n'.format(value))
                    subpage.write('    </tr>\n')
            if len(ifoList)>1: # Network overlap
                mode, low, high, inmode = mode_values('network_overlap'+polarization, '0', mdc, momentsPlus, momentsCross, ifoList, tablesDir, worc)
                subpage.write('    <tr>\n')
                if polarization == '':
                    subpage.write('      <td>Network - Signal </td>\n')
                elif polarization == '_plus':
                    subpage.write('      <td>Network - + Polarization </td>\n')
                else:
                    subpage.write('      <td>Network - x Polarization </td>\n')
                for value in [mode, low, high]:
                    subpage.write('      <td>{0}</td>\n'.format(value))
                subpage.write('    </tr>\n')
        subpage.write('  </table>\n')

def make_subpage(htmlDir, plotsDir, page, model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel):
    #if not os.path.exists(plotsDir):
    #if 'waveform' in page: 
    #    subpage = open(htmlDir+'waveform.html', 'w')
    #else:
    #    subpage = open(htmlDir+html_dict[page][0]+'.html', 'w')
    subpage = open(htmlDir+html_dict[page][0]+'.html', 'w')
    # -- Start of html subpage
    subpage.write('<html>\n')
    subpage.write('<head>\n')
    subpage.write('</head>\n')
    subpage.write('<body>\n')
    # -- Whitened data plots
    #if page == 'spec':
    #    subpage.write('    <h3>{0}</h3>\n'.format(html_dict[page][1].format(model)))
    #else:
    subpage.write('    <h2>{0}</h2>\n'.format(html_dict[page][1]))
    if page == 'waveform':
    #if 'waveform' in page:
        subpage.write('    Reconstructed waveforms in whitened data.  Red is injected waveform, blue is recovered waveform<br/>\n')
    if page == 'diagnostics':
        for plot in ['likelihood', 'TL', 'model_dimensions_histo', 'model_dimensions']:
            plotsrc = './'+plotsDir+plot+'.png'
            subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
    elif page == 'skymap':
        # Display the skymap
        plotsrc = './'+plotsDir+'skymap.png'
        subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=650></a><br/>\n')
        # Display all the megasky results: remember to run megasky first, then megaplot
        if os.path.exists(str(jobName)+'megasky_results.dat'):
            megasky_results = np.loadtxt(str(jobName)+'megasky_results.dat')
            pixarea, area50, area90, searcharea, injcontour = megasky_results 
            subpage.write('    The area per pixel is {0:.2f} sq. degrees.<br/>\n'.format(pixarea))
            subpage.write('    The area of the 50% region is {0:.2f} sq. degrees.<br/>\n'.format(area50))
            subpage.write('    The area of the 90% region is {0:.2f} sq. degrees.<br/>\n'.format(area90))
            subpage.write('    The searched area is {0:.2f} sq. degrees.<br/>\n'.format(searcharea))
            subpage.write('    The injection was found at the {0:.2f} percent area contour.<br/>\n'.format(injcontour))
        # Display the other megaksy plots, if present
        plotsrc = './'+plotsDir+'ra_acf.png'
        if os.path.exists(plotsrc):
            subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
        plotsrc = './'+plotsDir+'ra_chain.png'
        if os.path.exists(plotsrc):
            subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
    elif page in ['waveform']:
        subpage.write('  <br/>\n')
        for mod in modelList:
            subpage.write('    <center><h3>{0}</h3></center>\n'.format(html_dict[mod+'_'+page][1]))
            for ifo in ifoList:
                plotsrc = './'+plotsDir+html_dict[mod+'_'+page][2]+'_{0}.png'.format(ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            # -- Toggle for full plots
            if page != 'bayesogram':
                subpage.write('<br/>See full axis: <a id="displayFull'+page+mod+'" href="javascript:toggle(\'divFull'+page+mod+'\',\'displayFull'+page+mod+'\');">show</a>\n')
                subpage.write('  <div id="divFull'+page+mod+'" style="display: none">\n')
                for ifo in ifoList:
                    plotsrc = './'+plotsDir+html_dict[mod+'_'+page][2]+'_{0}_full.png'.format(ifo)
                    subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
                subpage.write("</div>\n")
    elif page == 'inj_spec':
        subpage.write('  <br/>\n')
        for ifo in ifoList:
            plotsrc = './'+plotsDir+html_dict[page][2]+'_{0}.png'.format(ifo)
            subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
        # -- Toggle for full plots
        subpage.write('<br/>See full axis: <a id="displayFull'+page+'" href="javascript:toggle(\'divFull'+page+'\',\'displayFull'+page+'\');">show</a>\n')
        subpage.write('  <div id="divFull'+page+'" style="display: none">\n')
        for ifo in ifoList:
            plotsrc = './'+plotsDir+html_dict[page][2]+'_{0}_full.png'.format(ifo)
            subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
        subpage.write("</div>\n")
    elif page in ['dur_rec', 't_energy_rec', 'f0_rec', 'band_rec', 't0_rec', 'h_max', 't_at_h_max','snr']: # All histogram pages, excluding overlaps
        subpage.write('  <br/>\n')
        for ytype in ['lin', 'log10']:
            if ytype == 'log10':
                subpage.write('  Histograms with logarithmic y-axis: <a id="displayLog10'+html_dict[page][0]+'" href="javascript:toggle(\'divLog10'+html_dict[page][0]+'\',\'displayLog10'+html_dict[page][0]+'\');">show</a>\n')
                subpage.write('  <div id="divLog10'+html_dict[page][0]+'" style="display: none">\n')
            for ifo in ifoList:
                plotsrc = './'+plotsDir+html_dict[page][2]+'_{1}_{0}.png'.format(ifo,ytype)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            subpage.write('    </br>\n')
            if ytype == 'log10':
                subpage.write('  </div>\n')
        # Make table with recovered and injected values of the mode
        make_mode_table(page, subpage, mdc, momentsPlus, momentsCross, ifoList, tablesDir, ifoNames, 'whitened')
        # -- Time domain page
    elif page == 'time_domain':
        subpage.write('  <br/>\n')
        for mod in modelList:
            subpage.write('    <center><h3>{0}</h3></center>\n'.format(html_dict[mod+'_'+page][1]))
            for ifo in ifoList:
                plotsrc = './'+plotsDir+html_dict[mod+'_'+page][2]+'_{0}.png'.format(ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            subpage.write('<br/>See full axis: <a id="displayFull'+page+mod+'" href="javascript:toggle(\'divFull'+page+mod+'\',\'displayFull'+page+mod+'\');">show</a>\n')
            subpage.write('  <div id="divFull'+page+mod+'" style="display: none">\n')
            for ifo in ifoList:
                plotsrc = './'+plotsDir+html_dict[mod+'_'+page][2]+'_{0}_full.png'.format(ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            subpage.write("</div>\n")
    elif page == 'frequency_domain':
        subpage.write('  <br/>\n')
        for mod in modelList:
            subpage.write('    <center><h3>{0}</h3></center>\n'.format(html_dict[mod+'_'+page][1]))
            for ifo in ifoList:
                plotsrc = './'+plotsDir+html_dict[mod+'_'+page][2]+'_{0}.png'.format(ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')


    elif page == 'time_frequency':
        subpage.write('  <br/>\n')
        for mod in modelList:
            subpage.write('    <center><h3>Median {0} model Spectrograms</h3></center>\n'.format(mod))
            for ifo in ifoList:
                plotsrc = './'+plotsDir+'{0}_spectrogram_Q8_{1}.png'.format(mod,ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            subpage.write('<br/>Other Q resolutions: <a id="displayFull'+page+mod+'" href="javascript:toggle(\'divFull'+page+mod+'\',\'displayFull'+page+mod+'\');">show</a>\n')
            subpage.write('  <div id="divFull'+page+mod+'" style="display: none">\n')
            for ifo in ifoList:
                plotsrc = './'+plotsDir+'{0}_spectrogram_Q4_{1}.png'.format(mod,ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            for ifo in ifoList:
                plotsrc = './'+plotsDir+'{0}_spectrogram_Q16_{1}.png'.format(mod,ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            subpage.write("</div>\n")

            

            subpage.write('    <center><h3>Median {0} time vs. frequency tracks</h3></center>\n'.format(mod))
            for ifo in ifoList:
                plotsrc = './'+plotsDir+'{0}_tf_{1}.png'.format(mod,ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')


        if mdc:
            subpage.write('    <center><h3>Injected waveform spectrograms</h3></center>\n'.format(mod))
            for ifo in ifoList:
                plotsrc = './'+plotsDir+'injected_spectrogram_Q8_{0}.png'.format(ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            subpage.write('<br/>Other Q resolutions: <a id="displayFull'+page+'data'+'" href="javascript:toggle(\'divFull'+page+'injected'+'\',\'displayFull'+page+'injected'+'\');">show</a>\n')
            subpage.write('  <div id="divFull'+page+'injected'+'" style="display: none">\n')
            for ifo in ifoList:
                plotsrc = './'+plotsDir+'injected_spectrogram_Q4_{0}.png'.format(ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            for ifo in ifoList:
                plotsrc = './'+plotsDir+'injected_spectrogram_Q16_{0}.png'.format(ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            subpage.write("</div>\n")



        subpage.write('    <center><h3>Data spectrograms</h3></center>\n'.format(mod))
        for ifo in ifoList:
            plotsrc = './'+plotsDir+'data_spectrogram_Q8_{0}.png'.format(ifo)
            subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
        subpage.write('<br/>Other Q resolutions: <a id="displayFull'+page+'data'+'" href="javascript:toggle(\'divFull'+page+'data'+'\',\'displayFull'+page+'data'+'\');">show</a>\n')
        subpage.write('  <div id="divFull'+page+'data'+'" style="display: none">\n')
        for ifo in ifoList:
            plotsrc = './'+plotsDir+'data_spectrogram_Q4_{0}.png'.format(ifo)
            subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
        for ifo in ifoList:
            plotsrc = './'+plotsDir+'data_spectrogram_Q16_{0}.png'.format(ifo)
            subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
        subpage.write("</div>\n")
            
            
        if not restrictModel == 'glitch':
            subpage.write('    <center><h3>Residual (data-signal) spectrograms</h3></center>\n'.format(mod))
            for ifo in ifoList:
                plotsrc = './'+plotsDir+'residual_spectrogram_Q8_{0}.png'.format(ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            subpage.write('<br/>Other Q resolutions: <a id="displayFull'+page+'data'+'" href="javascript:toggle(\'divFull'+page+'residuals'+'\',\'displayFull'+page+'residuals'+'\');">show</a>\n')
            subpage.write('  <div id="divFull'+page+'residuals'+'" style="display: none">\n')
            for ifo in ifoList:
                plotsrc = './'+plotsDir+'residual_spectrogram_Q4_{0}.png'.format(ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            for ifo in ifoList:
                plotsrc = './'+plotsDir+'residual_spectrogram_Q16_{0}.png'.format(ifo)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            subpage.write("</div>\n")



    else: # Overlap histograms
        subpage.write('  <br/>\n')
        for ytype in ['lin', 'log10']:
            if ytype == 'log10':
                subpage.write('  Histograms with logarithmic y-axis: <a id="displayLog10'+html_dict[page][0]+'" href="javascript:toggle(\'divLog10'+html_dict[page][0]+'\',\'displayLog10'+html_dict[page][0]+'\');">show</a>\n')
                subpage.write('  <div id="divLog10'+html_dict[page][0]+'" style="display: none">\n')
                subpage.write('    <h2>{0}</h2>\n'.format(html_dict[page][1]))
            for ifo in ifoList:
                plotsrc = './'+plotsDir+html_dict[page][2]+'_{1}_{0}.png'.format(ifo,ytype)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            if page == 'overlap' and len(ifoList)>1: # Network overlap
                plotsrc = './'+plotsDir+'network_overlap_{0}.png'.format(ytype)
                subpage.write('    <h3>Nework overlap between recovered signal and injection</h3>\n')
                subpage.write('    <a href="'+plotsrc+'"><img src="./'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            for polarization in ['plus', 'cross']: # Repeat for + and x polarization results
                for ifo in ifoList:
                    plotsrc = './'+plotsDir+html_dict[page][2]+'_'+polarization+'_{1}_{0}.png'.format(ifo,ytype)
                    if os.path.exists(plotsrc):
                        if ifo == '0':
                            subpage.write('    <h3>{0}</h3>\n'.format(html_dict[page+'_'+polarization][1]))
                        subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
                if page == 'overlap' and len(ifoList)>1:
                    plotsrc = './'+plotsDir+'network_overlap_'+polarization+'_{0}.png'.format(ytype)
                    if os.path.exists(plotsrc):
                        subpage.write('    <h3>Nework overlap between recovered signal and injection for the '+polarization+' polarization</h3>\n')
                        subpage.write('    <a href="'+plotsrc+'"><img src="./'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            subpage.write('  <br/>\n')
            if ytype == 'log10':
                subpage.write('  </div>\n')
            else:
                # Make table with overlaps 
                make_mode_table('overlap', subpage, mdc, momentsPlus, momentsCross, ifoList, tablesDir, ifoNames, 'whitened')
    # -- Coloured data plots
    subpage.write('  <br/>\n')
    if page in ['dur_rec', 't_energy_rec', 'f0_rec', 'band_rec', 't0_rec', 'h_max', 't_at_h_max']: # All histogram pages, excluding overlaps
        subpage.write('  Results for coloured data: <a id="displayC'+html_dict[page][0]+'" href="javascript:toggle(\'divC'+html_dict[page][0]+'\',\'displayC'+html_dict[page][0]+'\');">show</a>\n')
        subpage.write('  <div id="divC'+html_dict[page][0]+'" style="display: none">\n')
        for ytype in ['lin', 'log10']:
            if ytype == 'log10':
                subpage.write('  <br/>\n')
                subpage.write('  Histograms for coloured data with logarithmic y-axis: <a id="displayLog10C'+html_dict[page][0]+'" href="javascript:toggle(\'divLog10C'+html_dict[page][0]+'\',\'displayLog10C'+html_dict[page][0]+'\');">show</a>\n')
                subpage.write('  <div id="divLog10C'+html_dict[page][0]+'" style="display: none">\n')
            for ifo in ifoList:
                plotsrc = './'+plotsDir+'c'+html_dict[page][2]+'_{1}_{0}.png'.format(ifo,ytype)
                subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
            if ytype == 'log10':
                subpage.write('  </div>\n')
            else:
                # Make table with recovered and injected values of the mode
                make_mode_table(page, subpage, mdc, momentsPlus, momentsCross, ifoList, tablesDir, ifoNames, 'colored')
    elif page == 'inj_spec': # Spectrocram injection
        subpage.write('  Results for coloured data: <a id="displayC'+html_dict[page][0]+'" href="javascript:toggle(\'divC'+html_dict[page][0]+'\',\'displayC'+html_dict[page][0]+'\');">show</a>\n')
        subpage.write('  <div id="divC'+html_dict[page][0]+'" style="display: none">\n')
        for ifo in ifoList:
            plotsrc = './'+plotsDir+'c'+html_dict[page][2]+'_{0}.png'.format(ifo)
            subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
        # -- Toggle for full plots
        subpage.write('<br/>See full axis: <a id="displayCFull'+page+'" href="javascript:toggle(\'divCFull'+page+'\',\'displayCFull'+page+'\');">show</a>\n')
        subpage.write('  <div id="divCFull'+page+'" style="display: none">\n')
        for ifo in ifoList:
            plotsrc = './'+plotsDir+'c'+html_dict[page][2]+'_{0}_full.png'.format(ifo)
            subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
        subpage.write("</div>")
    elif page in ['waveform', 'bayesogram', 'spec', 'overlap']: # All other subpages
        subpage.write('  Results for coloured data: <a id="displayC'+html_dict[page][0]+'" href="javascript:toggle(\'divC'+html_dict[page][0]+'\',\'displayC'+html_dict[page][0]+'\');">show</a>\n')
        subpage.write('  <div id="divC'+html_dict[page][0]+'" style="display: none">\n')
        if page in ['waveform', 'bayesogram', 'spec']:
            for mod in modelList:
                subpage.write('    <h3><center>{0}</h3></center>\n'.format(html_dict[mod+'_'+page][1]))
                for ifo in ifoList:
                    plotsrc = './'+plotsDir+'c'+html_dict[mod+'_'+page][2]+'_{0}.png'.format(ifo)
                    subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
        else: # Overlap plots
            for ytype in ['lin', 'log10']:
                if ytype == 'log10':
                    subpage.write('  <br/>\n')
                    subpage.write('  Histograms for coloured data with logarithmic y-axis: <a id="displayLog10C'+html_dict[page][0]+'" href="javascript:toggle(\'divLog10C'+html_dict[page][0]+'\',\'displayLog10C'+html_dict[page][0]+'\');">show</a>\n')
                    subpage.write('  <div id="divLog10C'+html_dict[page][0]+'" style="display: none">\n')
                    subpage.write('    <h2>{0}</h2>\n'.format(html_dict[page][1]))
                for ifo in ifoList:
                    plotsrc = './'+plotsDir+'c'+html_dict[page][2]+'_{1}_{0}.png'.format(ifo,ytype)
                    subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
                if page == 'overlap' and len(ifoList)>1: # Network overlap
                    plotsrc = './'+plotsDir+'cnetwork_overlap_{0}.png'.format(ytype)
                    if os.path.exists(plotsrc):
                        subpage.write('    <h3>Nework overlap between recovered signal and injection for the '+polarization+' polarization</h3>\n')
                        subpage.write('    <a href="'+plotsrc+'"><img src="./'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
                for polarization in ['plus', 'cross']: # Repeat for + and x polarization results
                    for ifo in ifoList:
                        plotsrc = './'+plotsDir+'c'+html_dict[page][2]+'_'+polarization+'_{1}_{0}.png'.format(ifo,ytype)
                        if os.path.exists(plotsrc):
                            if ifo == '0':
                                subpage.write('    <h3>{0}</h3>\n'.format(html_dict[page+'_'+polarization][1]))
                            subpage.write('    <a href="'+plotsrc+'"><img src="'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
                    if page == 'overlap' and len(ifoList)>1:
                        plotsrc = './'+plotsDir+'cnetwork_overlap_'+polarization+'_{0}.png'.format(ytype)
                        if os.path.exists(plotsrc):
                            subpage.write('    <h3>Nework overlap between recovered signal and injection for the '+polarization+' polarization</h3>\n')
                            subpage.write('    <a href="'+plotsrc+'"><img src="./'+plotsrc+'" alt="'+plotsrc+'" width=500></a>\n')
                if ytype == 'log10':
                    subpage.write('  </div>\n')
                else:
                    # Make table with overlaps 
                    make_mode_table('overlap', subpage, mdc, momentsPlus, momentsCross, ifoList, tablesDir, ifoNames, 'colored')
    subpage.write('  </div>\n')
    # -- End of html nubpage
    subpage.write('</body>\n')
    subpage.write('</html>\n')
    subpage.close()

def make_homepage(htmlDir, plotsDir, tablesDir, model, ifoNames, snrList, bayeswaverunfile, sig_gl, sig_noise):
    summary = open(htmlDir+'summary.html', 'w')
    # -- Start of html page with basic info about the run
    summary.write('<html>\n')
    summary.write('<head>\n')
    summary.write('</head>\n')
    summary.write('<body>\n')
    summary.write('  <p>\n')
    summary.write('\n')
    # -- Odds plot
    summary.write('    <a href="./'+plotsDir+'odds.png"><img src="./'+plotsDir+'odds.png" style="float: right;" width=600 alt="odds.png"></a>\n')
    # -- Other information
    #summary.write('    <h3>{0} Model</h3>\n'.format(upperCaseModel_dict[model]))
    summary.write('    <h3>Detector Names: {0}</h3>\n'.format(', '.join(ifoNames)))
    if mdc:
      summary.write('    <h3>Matched Filter SNRs of Injections</h3>\n')
      for ifo, snr in zip(ifoNames, snrList):
        summary.write('    {0} injected with SNR {1:.1f}<br/>\n'.format(ifo, snr))
    summary.write('    <h3>Log Info</h3>\n')
    summary.write('    <a href="./'+bayeswaverunfile+'">See full log file</a>\n')
    summary.write('    <h3>Evidence for Signal</h3>\n')
    summary.write('    log(Evidence_signal / Evidence_glitch) = {0:.1f} &plusmn; {1:.1f} <br/>\n'.format(sig_gl,err_sig_gl))
    summary.write('    log(Evidence_signal / Evidence_noise) = {0:.1f} &plusmn; {1:.1f} <br/>\n'.format(sig_noise,err_sig_noise))
    summary.write('  </p>\n')
    # -- End of html page with basic info about the run
    summary.write('</body>\n')
    summary.write('</html>\n')
    summary.close()

def make_index(htmlDir, plotsDir, tablesDir, model, gps, ifoList, ifoNames, snrList, bayeswaverunfile, sig_gl, sig_noise):
    index = open('index.html', 'w')
    # -- Start of the BWB webpage
    index.write('<!DOCTYPE html/css>\n')
    index.write('<html>\n')
    index.write('<head>\n')
    index.write('<link rel="stylesheet" type="text/css" href="./html/BWBweb.css">\n')
    index.write('<!--<script src="http://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>--!>\n')
    index.write('<script src="./html/secure_ajax.js"></script>\n')
    index.write('<script src="./html/navigate.js"></script>\n')
    index.write('</head>\n')
    # -- Title
    index.write('<body>\n')
    index.write('<div class="container wrapper">\n')
    index.write('  <div id="top">\n')
    index.write('    <h1>BayesWave Output Page</h1>\n')
    index.write('       <p>Results for trigger at GPS {0}</p>\n'.format(gps))
    index.write('  </div>\n')
    # -- Navigation bar
    index.write('  <div class="wrapper">\n')
    index.write('    <div id="menubar">\n')
    index.write('      <ul id="menulist">\n')
    index.write('        <li class="menuitem" id="summary">Model selection\n')
    index.write('        <li class="menuitem" id="timedomain">Time domain waveforms\n')
    index.write('        <li class="menuitem" id="freqdomain">Frequency domain information\n')
    index.write('        <li class="menuitem" id="timefreq">Time vs Frequency information\n')
    index.write('        <li class="menuitem" id="snr">SNR\n')
    index.write('        <li class="menuitem" id="energy">Signal energy\n')
    index.write('        <li class="menuitem" id="t0">Central time\n')
    index.write('        <li class="menuitem" id="duration">Duration\n')
    index.write('        <li class="menuitem" id="f0">Central frequency\n')
    index.write('        <li class="menuitem" id="band">Bandwidth\n')
    index.write('        <li class="menuitem" id="hmax">Max Amplitude\n')
    #index.write('        <li class="menuitem" id="t_at_hmax">Time at Max Amplitude\n')
    index.write('        <li class="menuitem" id="skymap">Skymap\n')
    if mdc:
        index.write('        <li class="menuitem" id="overlap">Overlap\n')


    index.write('        <li class="menuitem" id="diagnostics">Diagnostics\n')
    index.write('      </ul>\n')
    index.write('    </div>\n')
    # -- Main part
    index.write('    <div id="main">\n')
    index.write('      <p>\n')
    index.write('\n')
    index.write('        <a href="./'+plotsDir+'odds.png"><img src="./'+plotsDir+'odds.png" style="float: right;" width=600 alt="odds.png"></a>\n')
    #index.write('        <h3>{0} Model</h3>\n'.format(upperCaseModel_dict[model]))
    index.write('        <h3>Detector Names: {0}</h3>\n'.format(', '.join(ifoNames)))
    if mdc:
        index.write('        <h3>Matched Filter SNRs of Injections</h3>\n')
        for ifo, snr in zip(ifoNames, snrList):
            index.write('        {0} injected with SNR {1:.3f}<br/>\n'.format(ifo, snr))
    index.write('        <h3>Log Info</h3>\n')
    index.write('    <a href="./'+bayeswaverunfile+'">See full log file</a>\n')
    index.write('        <h3>Evidence for Signal</h3>\n')
    index.write('    log(Evidence_signal / Evidence_glitch) = {0:.2f} &plusmn; {1:.2f} <br/>\n'.format(sig_gl,err_sig_gl))
    index.write('    log(Evidence_signal / Evidence_noise) = {0:.2f} &plusmn; {1:.2f} <br/>\n'.format(sig_noise,err_sig_noise))
    index.write('      </p>\n')
    index.write('    </div>\n')
    # -- End of the BWB webpage
    index.write('</body>\n')
    index.write('</html>\n')
    index.close()


def make_navpage(htmlDir):
    
    
    navpage = open(htmlDir+'navigate.js','w')
    
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#summary").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/summary.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('  $("#timedomain").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/timedomain.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('  $("#freqdomain").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/freqdomain.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('  $("#timefreq").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/timefreq.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#overlap").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/overlap.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#energy").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/energy.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#t0").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/t0.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#duration").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/duration.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#f0").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/f0.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#band").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/band.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#hmax").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/hmax.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#t_at_hmax").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/t_at_hmax.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#spec").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/spec.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#diagnostics").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/diagnostics.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#skymap").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/skymap.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#injections").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/injections.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('$(document).ready(function(){')
    navpage.write('\n')
    navpage.write('    $("#snr").click(function(){')
    navpage.write('\n')
    navpage.write('    $("#main").load("./html/snr.html");')
    navpage.write('\n')
    navpage.write('  });')
    navpage.write('\n')
    navpage.write('});')
    navpage.write('\n')
    navpage.write('function toggle(showHideDiv, switchTextDiv) {')
    navpage.write('\n')
    navpage.write('	var ele = document.getElementById(showHideDiv);')
    navpage.write('\n')
    navpage.write('	var text = document.getElementById(switchTextDiv);')
    navpage.write('\n')
    navpage.write('	if(ele.style.display == "block") {')
    navpage.write('\n')
    navpage.write('    		ele.style.display = "none";')
    navpage.write('\n')
    navpage.write('		text.innerHTML = "show";')
    navpage.write('\n')
    navpage.write('  	}')
    navpage.write('\n')
    navpage.write('	else {')
    navpage.write('\n')
    navpage.write('		ele.style.display = "block";')
    navpage.write('\n')
    navpage.write('		text.innerHTML = "hide";')
    navpage.write('\n')
    navpage.write('	}')
    navpage.write('\n')
    navpage.write('} ')
    navpage.write('\n')
    

def make_css(htmlDir):
    csspage = open(htmlDir+'BWBweb.css')

    csspage.write('body {\n')
    csspage.write('    font: 90% Lucida Sans;\n')
    csspage.write('    margin: 20px;\n')
    csspage.write('    line-height: 26px;\n')
    csspage.write('}\n')
    csspage.write('\n')
    csspage.write('.container {\n')
    csspage.write('    min-width: 1220px;\n')
    csspage.write('}\n')
    csspage.write('\n')
    csspage.write('.wrapper {\n')
    csspage.write('    position: relative;\n')
    csspage.write('    overflow: auto;\n')
    csspage.write('}\n')
    csspage.write('\n')
    csspage.write('/*#top, #sidebar, #bottom, .menuitem { */\n')
    csspage.write('#top, #bottom, .menuitem {\n')
    csspage.write('    border-radius: 4px;\n')
    csspage.write('    margin: 4px;\n')
    csspage.write('}\n')
    csspage.write('\n')
    csspage.write('#top {\n')
    csspage.write('    background-color: #71B09F;\n')
    csspage.write('    color: #ffffff;\n')
    csspage.write('    border: 1px solid #374040;\n')
    csspage.write('    padding: 15px;\n')
    csspage.write('}\n')
    csspage.write('\n')
    csspage.write('#menubar {\n')
    csspage.write('    width: 200px;\n')
    csspage.write('    float: left\n')
    csspage.write('}\n')
    csspage.write('\n')
    csspage.write('#main {\n')
    csspage.write('    padding: 10px;\n')
    csspage.write('    margin: 0 210px;\n')
    csspage.write('    /*padding-right: 0px;\n')
    csspage.write('    margin-right: 0px;*/\n')
    csspage.write('    min-width: 1010px;\n')
    csspage.write('}\n')
    csspage.write('\n')
    csspage.write('/* #sidebar {\n')
    csspage.write('    background-color: #32a4e7;\n')
    csspage.write('    color: #ffffff;\n')
    csspage.write('    padding: 10px;\n')
    csspage.write('    width: 180px;\n')
    csspage.write('    bottom: 0;\n')
    csspage.write('    top: 0;\n')
    csspage.write('    right: 0;\n')
    csspage.write('    position: absolute;\n')
    csspage.write('} */\n')
    csspage.write('\n')
    csspage.write('#bottom {\n')
    csspage.write('    border: 1px solid #d4d4d4;\n')
    csspage.write('    background-color: #f1f1f1;\n')
    csspage.write('    text-align: center;\n')
    csspage.write('    padding: 10px;\n')
    csspage.write('    font-size: 70%;\n')
    csspage.write('    line-height: 14px;\n')
    csspage.write('}\n')
    csspage.write('\n')
    csspage.write('#top h1, #top p, #menulist {\n')
    csspage.write('    margin: 0;\n')
    csspage.write('    padding: 0;\n')
    csspage.write('}\n')
    csspage.write('\n')
    csspage.write('.menuitem {\n')
    csspage.write('    background-color: #B1C7C7;\n')
    csspage.write('    border: 1px solid #333333;\n')
    csspage.write('    list-style-type: none;\n')
    csspage.write('    padding: 2px;\n')
    csspage.write('    cursor: pointer;\n')
    csspage.write('}\n')
    csspage.write('\n')
    csspage.write('.menuitem:hover {\n')
    csspage.write('  background-color: #7A8787;\n')
    csspage.write('}\n')
    csspage.write('\n')
    csspage.write('a {\n')
    csspage.write('    color: #000000;\n')
    csspage.write('    text-decoration: underline;\n')
    csspage.write('}\n')
    csspage.write('\n')
    csspage.write('a:hover {\n')
    csspage.write('    text-decoration: none;\n')
    csspage.write('}\n')
    csspage.write('\n')
    csspage.write('table {\n')
    csspage.write('    width:100%;\n')
    csspage.write('}\n')
    csspage.write('table, th, td {\n')
    csspage.write('    border: 1px solid #333333;\n')
    csspage.write('    border-collapse: collapse;\n')
    csspage.write('}\n')
    csspage.write('th, td {\n')
    csspage.write('    padding: 5px;\n')
    csspage.write('    text-align: left;\n')
    csspage.write('}\n')
    csspage.write('table#tab1 tr:nth-child(even) {\n')
    csspage.write('    background-color: #d4d4d4;\n')
    csspage.write('}\n')
    csspage.write('table#tab1 tr:nth-child(odd) {\n')
    csspage.write('   background-color:#fff;\n')
    csspage.write('}\n')
    csspage.write('table#tab1 th {\n')
    csspage.write('    background-color: #71B09F;\n')
    csspage.write('    color: white;\n')
    csspage.write('}\n')



def make_webpage(htmlDir, model, mdc, gps, ifoList, ifoNames, modelList, snrList, sig_gl, sig_noise, momentsPlus, momentsCross, postprocesspath):
    # TODO: move in separate function
    # -- Find out the path to the BayesWave executable

    os.system('cp '+postprocesspath+'/BWBweb.css '+htmlDir+'.')
    os.system('cp '+postprocesspath+'/secure_ajax.js '+htmlDir+'.')
    os.system('cp '+postprocesspath+'/navigate.js '+htmlDir+'.')

    # -- Write the index
    make_index(htmlDir, plotsDir, tablesDir, model, gps, ifoList, ifoNames, snrList, bayeswaverunfile, sig_gl, sig_noise)
    # -- Write summary page (works as a homepage) #TODO: index.html and summary.html have code in common
    make_homepage(htmlDir, plotsDir, tablesDir, model, ifoNames, snrList, bayeswaverunfile, sig_gl, sig_noise)
    # -- Write all subpages
    make_subpage(htmlDir, plotsDir, 'dur_rec', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    make_subpage(htmlDir, plotsDir, 't_energy_rec', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    make_subpage(htmlDir, plotsDir, 'f0_rec', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    make_subpage(htmlDir, plotsDir, 'band_rec', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    make_subpage(htmlDir, plotsDir, 't0_rec', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    make_subpage(htmlDir, plotsDir, 'h_max', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    make_subpage(htmlDir, plotsDir, 't_at_h_max', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
#make_subpage(htmlDir, plotsDir, 'waveform', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross)
    make_subpage(htmlDir, plotsDir, 'time_domain', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    make_subpage(htmlDir, plotsDir, 'frequency_domain', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    make_subpage(htmlDir, plotsDir, 'time_frequency', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    make_subpage(htmlDir, plotsDir, 'bayesogram', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    make_subpage(htmlDir, plotsDir, 'spec', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    make_subpage(htmlDir, plotsDir, 'diagnostics', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    make_subpage(htmlDir, plotsDir, 'skymap', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    make_subpage(htmlDir, plotsDir, 'snr', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)
    if mdc:
        make_subpage(htmlDir, plotsDir, 'overlap', model, modelList, ifoList, ifoNames, momentsPlus, momentsCross,restrictModel)



######################################################################################################################
#
# Main
#
######################################################################################################################
# -- Parse command line arguments
parser = argparse.ArgumentParser(description='Produce html page for a bayeswave run.')


# -- Get basic info on the run
jobName, restrictModel, mdc, injFlag, bayeswaverunfile, ifoList, ifoNames, gps, snrList, info = readbwb()
if not restrictModel == '':
    model = restrictModel  #TODO: think about the --noiseOnly case
    modelList = [restrictModel]
    restrictModel = modelList
else:
    model = 'signal'
    restrictModel = ('signal', 'glitch', 'noise')
print "The mdc status is: {0}\n".format(mdc)
if injFlag:
    print "The injection was performed via an xml table\n"

if len(restrictModel) == 1:
    restrictModel = restrictModel[0]


# -- Create directory that will contain all plots
if not os.path.exists(plotsDir):
    os.makedirs(plotsDir)

# -- Create directory that will contain all the tables 
if not os.path.exists(tablesDir):
    os.makedirs(tablesDir)
    
# -- Read in time vector
time = np.loadtxt(str(jobName)+postDir+'timesamp.dat')
freq = np.loadtxt(str(jobName)+postDir+'freqsamp.dat')

for mod in modelList: 
    print "Analyzing the {0} model\n".format(mod)
    
    # -- Loop over interferometers to create plots
    for ifo in ifoList:
        # --------------------
        # Make waveform plots
        # --------------------
        for worc in ['whitened']:

            # -- Determine axes for waveform plots
            axwinner = get_axes(jobName, postDir, ifoList, worc, model, time)

            # -- Read in the reconstructed waveform
            try:
                filename = str(jobName)+postDir+'{1}_recovered_{2}_waveform.dat.{0}'.format(ifo, mod, worc)
                #median_waveform, up_vec, down_vec = get_waveform_bigfile(filename, ifo)
                timesamp,median_waveform, high_50, low_50, high_90, low_90 = get_waveform(filename)
            except:
                filename = str(jobName)+postDir+'{1}_median_time_domain_waveform.dat.{0}'.format(ifo, mod)
                #median_waveform, up_vec, down_vec = get_waveform_bigfile(filename, ifo)
                timesamp,median_waveform, high_50, low_50, high_90, low_90 = get_waveform(filename)
    
            # -- Plot the median waveform and injected waveform (if present)
            plot_waveform(jobName, postDir, ifo, plotsDir, worc, mdc, mod, axwinner, time, low_50, high_50, low_90, high_90, median_waveform)
            
            
            filename = str(jobName)+postDir+'{1}_median_frequency_domain_waveform.dat.{0}'.format(ifo,mod)
            powerspec_info = get_waveform(filename)
            
            filename = str(jobName)+postDir+'{1}_median_PSD.dat.{0}'.format(ifo,mod)
            psd_info = get_waveform(filename)
            
            f_axwinner = get_axes_fdom(jobName, postDir, ifoList, worc, model, time)
            
            plot_full_spectro(jobName, postDir, ifo, plotsDir, worc, mdc, mod, f_axwinner, psd_info, powerspec_info, median_waveform)
            

            
            # ---- tf track
            filename = str(jobName)+postDir+'{1}_median_tf.dat.{0}'.format(ifo,mod)
            freqsamp, median_powerspec, high_50, low_50, high_90, low_90 = get_waveform(filename)
            
            plot_tf_tracks(jobName, postDir, ifo, plotsDir, worc, mdc, mod, axwinner, f_axwinner, freqsamp, low_50, high_50, low_90, high_90, median_powerspec);
            
            
              #if injFlag:
              # plot_plus_waveform(jobName, postDir, ifo, plotsDir, worc, injFlag, mod, axwinner, time, down_vec, up_vec, median_waveform)
              # plot_cross_waveform(jobName, postDir, ifo, plotsDir, worc, injFlag, mod, axwinner, time, down_vec, up_vec, median_waveform)
    
        #---------
        # Qscan
        #---------
        climwinner4 = [1,1]
        climwinner8 = [1,1]
        climwinner16 = [1,1]
        
        Q_scan(mod, 4, time, freq, ifo, axwinner, f_axwinner)
        Q_scan(mod, 8, time, freq, ifo, axwinner, f_axwinner)
        Q_scan(mod, 16, time, freq, ifo, axwinner, f_axwinner)
        

        climwinner4 = Q_scan('data', 4, time, freq, ifo, axwinner, f_axwinner)
        climwinner8 = Q_scan('data', 8, time, freq, ifo, axwinner, f_axwinner)
        climwinner16 = Q_scan('data', 16, time, freq, ifo, axwinner, f_axwinner)

        if not restrictModel == 'glitch':
            Q_scan('residual', 4, time, freq, ifo, axwinner, f_axwinner, climwinner4)
            Q_scan('residual', 8, time, freq, ifo, axwinner, f_axwinner, climwinner8)
            Q_scan('residual', 16, time, freq, ifo, axwinner, f_axwinner, climwinner16)
        
        if mdc:
            Q_scan('injected', 4, time, freq, ifo, axwinner, f_axwinner)
            Q_scan('injected', 8, time, freq, ifo, axwinner, f_axwinner)
            Q_scan('injected', 16, time, freq, ifo, axwinner, f_axwinner)

        
        # ---------------------------------
        # Make histograms of PEC parameters
        # ---------------------------------
        if mod == model:
            for worc in ['whitened', 'colored']:
                # -- Read moments file and open mode output file
                moments = np.recfromtxt(str(jobName)+postDir+'{1}_{2}_moments.dat.{0}'.format(ifo, model, worc), names=True)
                momentsPlus = [] # Moments relative to the + polarization only
                momentsCross = [] # Moments relative to the x polarization only
                if mdc: 
                    injmoments = np.recfromtxt(str(jobName)+postDir+'injection_{1}_moments.dat.{0}'.format(ifo, worc), names=True)
                      #if injFlag:
                      # momentsPlus = np.recfromtxt(str(jobName)+postDir+'plus{1}_{2}_moments.dat.{0}'.format(ifo, upperCaseModel_dict[model], worc), names=True)
                      # momentsCross = np.recfromtxt(str(jobName)+postDir+'cross{1}_{2}_moments.dat.{0}'.format(ifo, upperCaseModel_dict[model], worc), names=True)
                else: 
                    injmoments = moments
                if worc == 'whitened':
                    outfile = open(tablesDir+'mode_{0}.txt'.format(ifo), 'w')
                else:
                    outfile = open(tablesDir+'cmode_{0}.txt'.format(ifo), 'w')

                # -- Write output file header 
                outfile.write("# Rows are: 1  Overlap of the signal\n")
                outfile.write("#           2  Overlap of the plus polarized signal          (0's if this info is not available)\n")
                outfile.write("#           3  Overlap of the cross polarized signal         (0's if this info is not available)\n")
                outfile.write("#           4  Network Overlap of the signal                 (0's if this info is not available)\n")
                outfile.write("#           5  Network Overlap of the plus polarized signal  (0's if this info is not available)\n")
                outfile.write("#           6  Network Overlap of the cross polarized signal (0's if this info is not available)\n")
                outfile.write("#           7  Signal Energy               (mode, 90% lower bound, 90% upper bound)\n")
                outfile.write("#           8  Injected Signal Energy      (repeated three times, 0's if this info is not available)\n")
                outfile.write("#           9  Central Time                (mode, 90% lower bound, 90% upper bound)\n")
                outfile.write("#           10 Injected Central Time       (repeated three times, 0's if this info is not available)\n")
                outfile.write("#           11 Duration                    (mode, 90% lower bound, 90% upper bound)\n")
                outfile.write("#           12 Injected Duration           (repeated three times, 0's if this info is not available)\n")
                outfile.write("#           13 Central Frequency           (mode, 90% lower bound, 90% upper bound)\n")
                outfile.write("#           14 Injected Central Freq       (repeated three times, 0's if this info is not available)\n")
                outfile.write("#           15 Bandwidth                   (mode, 90% lower bound, 90% upper bound)\n")
                outfile.write("#           16 Injected Bandwidth          (repeated three times, 0's if this info is not available)\n")
                outfile.write("#           17 Max ||h(t)||                (mode, 90% lower bound, 90% upper bound)\n")
                outfile.write("#           18 Injected Max ||h(t)||       (repeated three times, 0's if this info is not available)\n")
                outfile.write("#           19 t at Max ||h(t)||           (mode, 90% lower bound, 90% upper bound)\n")
                outfile.write("#           20 Injected t at Max ||h(t)||  (repeated three times, 0's if this info is not available)\n")
                outfile.write("#           21 SNR                         (mode, 90% lower bound, 90% upper bound)\n")

                # -- Histogram of the overlaps: use both linear and log10 y-axis
                if mdc:
                    for ytype in ['lin', 'log10']:
                        # Single detector
                        make_histogram(moments, 'overlap', ifo, plotsDir, worc, outfile, 'lin', ytype)
                        if len(momentsPlus) > 0:
                            make_histogram(momentsPlus, 'overlap_plus', ifo, plotsDir, worc, outfile, 'lin', ytype)
                        elif ytype == 'lin':
                            outfile.write("0 0 0\n")
                        if len(momentsCross) > 0:
                            make_histogram(momentsCross, 'overlap_cross', ifo, plotsDir, worc, outfile, 'lin', ytype)
                        elif ytype == 'lin':
                            outfile.write("0 0 0\n")
                        # Detector network
                        if len(ifoList)>1:
                            make_histogram(moments, 'network_overlap', ifo, plotsDir, worc, outfile, 'lin', ytype)
                            if len(momentsPlus) > 0:
                                make_histogram(momentsPlus, 'network_overlap_plus', ifo, plotsDir, worc, outfile, 'lin', ytype)
                            elif ytype == 'lin':
                                outfile.write("0 0 0\n")
                            if len(momentsCross) > 0:
                                make_histogram(momentsCross,'network_overlap_cross', ifo, plotsDir, worc, outfile, 'lin', ytype)
                            elif ytype == 'lin':
                                outfile.write("0 0 0\n")
                        elif ytype == 'lin':
                            i = 0
                            while i < 3:
                                outfile.write("0 0 0\n")
                                i = i+1

                else:
                   # Since this is not an MDC, overlaps will be ignored and filled with 0's
                    i = 0
                    while i < 6:
                        outfile.write("0 0 0\n")
                        i = i+1
                # -- Make histograms with linear and log10 y-axis
                for ytype in ['lin', 'log10']:
                        # Histogram of signal energy
                    make_histogram(moments, 't_energy_rec', ifo, plotsDir, worc, outfile, 'log10', ytype)
                        # Histogram central time
                    make_histogram(moments, 't0_rec', ifo, plotsDir, worc, outfile, 'lin', ytype)
                        # Histogram of duration
                    make_histogram(moments, 'dur_rec', ifo, plotsDir, worc, outfile, 'log10', ytype)
                        # Histogram of central frequency
                    make_histogram(moments, 'f0_rec', ifo, plotsDir, worc, outfile, 'lin', ytype)
                        # Histogram of bandwidth
                    make_histogram(moments, 'band_rec', ifo, plotsDir, worc, outfile, 'lin', ytype)
                    # Histogram of the maximum ||h(t)||
                    make_histogram(moments, 'h_max', ifo, plotsDir, worc, outfile, 'lin', ytype)
                        # Histogram of t at maximum ||h(t)||
                    make_histogram(moments, 't_at_h_max', ifo, plotsDir, worc, outfile, 'lin', ytype)
                        # Histogram of snr
                    make_histogram(moments,'snr', ifo, plotsDir, worc, outfile, 'lin', ytype)
                                                
            # -- End the loop over white/colored moments
                outfile.close()

                        


# -- Plot evidence
sig_gl, sig_noise, err_sig_gl, err_sig_noise = plot_evidence(jobName, plotsDir)

# -- Plot Temp Chains
# WARNING: This feature appears to be broken because of a new file format!

try:
    plot_likelihood_1(modelList, plotsDir)
    plot_likelihood_2(modelList, plotsDir)
except:
    print "Failed to plot temp vs. likelihood.  This may be due to a change in file format for files like signal_evidence.dat"

# -- Plot the evolution of the dimensions of each model; do it only if signal and/or glitch model are present


if 'signal' in modelList or 'glitch' in modelList:
    plot_model_dims(modelList, ifoList, ifoNames, plotsDir)


# -- Plot PSD
plot_psd(jobName, postDir, ifoList, ifoNames, plotsDir)

# -- Create directory that will contain all webpage code
if not os.path.exists(htmlDir):
    os.makedirs(htmlDir)

# -- Make webpage
make_webpage(htmlDir, model, mdc, gps, ifoList, ifoNames, modelList, snrList, sig_gl, sig_noise, momentsPlus, momentsCross, postprocesspath)

# Move back to original dir
os.chdir(topdir)
