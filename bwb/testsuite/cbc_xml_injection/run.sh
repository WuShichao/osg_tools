#!/bin/bash


###############################################################
#                                                             #
# Example BayesWave run for EOB injection in SIMULATED data.  #   
#                                                             # 
# RUN TIME ~16 hours                                          #  
#  - Processor: 2 GHz Intel Core i7                           #
#  - Memory: 4 GB 1600 Mhz                                    #
#                                                             #
###############################################################

source /usr/local/etc/master_rc
source /Users/tyson/bayeswave/trunk/utils/setup_paths.sh

# Copy injection XML to your working directory
#gsiscp ldas-grid.ligo.caltech:/home/tyson/examples/cbc_xml_injection/HL-INJECTIONS-EOB-ALIGO.xml .

# The main BayesWave run
bayeswave \
        --ifo H1 --H1-flow 16 \
        --ifo L1 --L1-flow 16 \
        --H1-cache LALSimAdLIGO --H1-channel LALSimAdLIGO \
        --L1-cache LALSimAdLIGO --L1-channel LALSimAdLIGO \
        --srate 1024 --seglen 4 \
        --trigtime 952402816.00 --PSDstart 952402818.00 --PSDlength 64 \
        --inj HL-INJECTIONS-EOB-ALIGO.xml --event 20 \
        --dataseed 1234 --Niter 2000000 --bayesLine --gnuplot 

rc=$?; if [[ $rc != 0 ]]; then exit $rc; fi

# BayesWave's post processing
bayeswave_post \
        --ifo H1 --H1-flow 16 \
        --ifo L1 --L1-flow 16 \
        --H1-cache LALSimAdLIGO --H1-channel LALSimAdLIGO \
        --L1-cache LALSimAdLIGO --L1-channel LALSimAdLIGO \
        --srate 1024 --seglen 4 \
        --trigtime 952402816.00 --PSDstart 952402816.00 --PSDlength 4 \
        --inj HL-INJECTIONS-EOB-ALIGO.xml --event 20 --noClean \
        --dataseed 1234 \
        --0noise

# Generate sky map
#python /Users/tyson/bayeswave/trunk/postprocess/skymap/megasky.py --mdc=HL-INJECTIONS-EOB-ALIGO.xml
python /Users/tyson/bayeswave/trunk/postprocess/skymap/megasky.py 

# Generate output webpages
python /Users/tyson/bayeswave/trunk/postprocess/megaplot.py




