#!/bin/bash

##################################################################################
#                                                                                #
# Example BayesWave run for BurstMDC  SineGaussian injection in SIMULATED data.  #   
#                                                                                # 
# RUN TIME ~10 hours                                                             #  
#  - Processor: 2 GHz Intel Core i7                                              #
#  - Memory: 4 GB 1600 Mhz                                                       #
#                                                                                #
##################################################################################


##################################################################################
#                                                                                #  
# NOTE: The burstMDC frame file is large so the transfer will take some time.    #
#                                                                                #
##################################################################################

source /usr/local/etc/master_rc
source /Users/tyson/bayeswave/trunk/utils/setup_paths.sh

# Get S6 BurstMDC injection frames with the data we need off of CIT
#gsiscp ldas-grid.ligo.caltech.edu:/mnt/qfs3/tyson/public_html/BayesWave/examples/burst_mdc_injection/HL-SineGaussian-968046516-4095.gwf .

# Make cache file for the injection frames
find HL-SineGaussian-968046516-4095.gwf | lalapps_path2cache > MDC.cache

# Main BayesWave analysis. Note that this injection is in simulated data (LALSimAdLIGO)
bayeswave \
	--ifo H1 --H1-flow 16 \
	--ifo L1 --L1-flow 16 \
	--H1-cache LALSimAdLIGO --H1-channel H1:LALSimAdLIGO \
	--L1-cache LALSimAdLIGO --L1-channel L1:LALSimAdLIGO \
	--srate 1024 --seglen 4 \
	--trigtime 968046586.000000 --PSDstart 968046586.000000 --PSDlength 8 \
	--MDC-channel [H1:Science,L1:Science] \
       	--MDC-cache [MDC.cache,MDC.cache] \
       	--MDC-prefactor 1 \
	--dataseed 1234 --bayesLine --Niter 2000000 --gnuplot

# BayesWave's post processing
bayeswave_post \
        --ifo H1 --H1-flow 16 \
        --ifo L1 --L1-flow 16 \
        --H1-cache LALSimAdLIGO --H1-channel H1:LALSimAdLIGO \
        --L1-cache LALSimAdLIGO --L1-channel L1:LALSimAdLIGO \
        --srate 1024 --seglen 4 \
        --trigtime 968046586.000000 --PSDstart 968046586.000000 --PSDlength 8 --bayesLine \
	--MDC-channel [H1:Science,L1:Science] \
        --MDC-cache [MDC.cache,MDC.cache] \
        --MDC-prefactor 1 \
	--dataseed 1234 \
	--0noise 

# Make skymap
python /Users/tyson/bayeswave/trunk/postprocess/skymap/megasky.py 

# Make the output page
python /Users/tyson/bayeswave/trunk/postprocess/megaplot.py 
