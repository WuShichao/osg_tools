#!/bin/bash

###############################################################
#                                                             #
# Example BayesWave run sampling from prior                   #   
#                                                             # 
# RUN TIME ~20 minutes                                        #  
#  - Processor: 2 GHz Intel Core i7                           #
#  - Memory: 4 GB 1600 Mhz                                    #
#                                                             #
###############################################################

# Setup paths for LAL install
#source /usr/local/etc/master_rc
#source /Users/tyson/bayeswave/trunk/utils/setup_paths.sh

# Main BayesWave analysis
bayeswave \
	--ifo H1 --H1-flow 16 \
	--ifo L1 --L1-flow 16 \
	--H1-cache LALSimAdLIGO --H1-channel H1:LALSimAdLIGO \
	--L1-cache LALSimAdLIGO --L1-channel L1:LALSimAdLIGO \
	--srate 1024 --seglen 4 \
	--trigtime 900000000.00 --PSDstart 900000000 --PSDlength 1024 \
	--dataseed 1234 --NC 5 --zeroLogL --orientationProposal 

# Check bayeswave exited normally before heading to post processing
rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

# BayesWave post processing
bayeswave_post \
	--ifo H1 --H1-flow 16 \
	--ifo L1 --L1-flow 16 \
       	--H1-cache LALSimAdLIGO --H1-channel H1:LALSimAdLIGO \
       	--L1-cache LALSimAdLIGO --L1-channel L1:LALSimAdLIGO \
       	--srate 1024 --seglen 4 \
       	--trigtime 900000000.00 --PSDstart 900000000 --PSDlength 1024 \
       	--dataseed 1234 \
	--0noise --noClean

# Generate skymap
python /Users/tyson/bayeswave/trunk/postprocess/skymap/megasky.py

# Generate output webpage
python /Users/tyson/bayeswave/trunk/postprocess/megaplot.py 
