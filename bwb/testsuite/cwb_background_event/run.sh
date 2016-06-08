#!/bin/bash

###############################################################
#                                                             #
# Example BayesWave run for cWB background event  in S6 data. #   
#                                                             # 
# RUN TIME ~6 hours                                           #  
#  - Processor: 2 GHz Intel Core i7                           #
#  - Memory: 4 GB 1600 Mhz                                    #
#                                                             #
###############################################################

source ~/lscsoft/master/etc/lscsoftrc 

# Get S6 frames that contain the bkgd event off of CIT
#gsiscp ldas-grid.ligo.caltech.edu:/data/node196/frames/S6/LDAShoftC02/LHO/H-H1_LDAS_C02_L2-9650/H-H1_LDAS_C02_L2-965060608-128.gwf .
#gsiscp ldas-grid.ligo.caltech.edu:/data/node197/frames/S6/LDAShoftC02/LLO/L-L1_LDAS_C02_L2-9650/L-L1_LDAS_C02_L2-965060608-128.gwf .

# Make LIGO cache file for frames
find H-H1_LDAS_C02_L2-965060608-128.gwf | lalapps_path2cache > H1.cache
find L-L1_LDAS_C02_L2-965060608-128.gwf | lalapps_path2cache > L1.cache

# The main BayesWave analysis
bayeswave \
        --ifo H1 --H1-flow 16 \
        --ifo L1 --L1-flow 16 \
        --H1-cache H1.cache --H1-channel H1:LDAS-STRAIN \
        --L1-cache L1.cache --L1-channel L1:LDAS-STRAIN \
        --L1-timeslide -26 \
        --srate 512 --seglen 4 \
        --trigtime 965060683.3113 --PSDstart 965060684.3113 --PSDlength 32 \
        --bayesLine --Niter 1000000 --gnuplot 

# BayesWave's post processing.  Notice that the IFO-cache and IFO-channel arguments are different
bayeswave_post \
        --ifo H1 --H1-flow 16 \
        --ifo L1 --L1-flow 16 \
        --H1-cache LALSimAdLIGO --H1-channel LALSimAdLIGO \
        --L1-cache LALSimAdLIGO --L1-channel LALSimAdLIGO \
        --srate 512 --seglen 4 \
        --trigtime 965060683.3113 --PSDstart 965060683.3113 --PSDlength 8 \
        --bayesLine \
        --dataseed 1234 \
        --0noise

# Make skymap
python ~/bayeswave/trunk/postprocess/skymap/megasky.py

# Make the output page
python ~/bayeswave/trunk/postprocess/megaplot.py

