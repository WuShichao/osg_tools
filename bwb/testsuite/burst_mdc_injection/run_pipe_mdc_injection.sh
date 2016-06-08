#!/bin/bash
#
###########################################################################################
#                                                                                         #
# Example BayesWave pipeline run for BurstMDC  SineGaussian injection in SIMULATED data.  #   
#                                                                                         # 
# RUN TIME ~10 hours                                                                      #  
#  - Processor: 2 GHz Intel Core i7                                                       #
#  - Memory: 4 GB 1600 Mhz                                                                #
#                                                                                         #
###########################################################################################

#source ~jclark/etc/source_master.sh
source ~jclark/etc/bayeswave_ldg-user-env.sh

# Make cache file for the injection frames
find HL-SineGaussian-968046516-4095.gwf | lalapps_path2cache > MDC.cache

# Setup the pipeline.  Can handle multiple trig-times with ascii list
# (--triger-list <list>).
bwb_pipe.py burst_mdc_injection.ini \
    --workdir burst_mdc_injection  \
    --trigger-time 968046586.000000
