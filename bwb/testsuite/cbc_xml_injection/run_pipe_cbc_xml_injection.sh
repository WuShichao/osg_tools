#!/bin/bash
#
########################################################################
#                                                                      #
# Example BayesWave pipeline run for EOB injection in SIMULATED data.  #   
#                                                                      # 
# RUN TIME ~16 hours                                                   #  
#  - Processor: 2 GHz Intel Core i7                                    #
#  - Memory: 4 GB 1600 Mhz                                             #
#                                                                      #
########################################################################

source ~jclark/etc/source_master.sh

#
# Generate sim-inspiral table using lalapps_inspinj
#

seed=`lalapps_tconvert now`
gpsstart=`lalapps_tconvert now`
gpsend=$((${gpsstart} + 100))
outfile="HL-INJECTIONS-SEOB-ALIGO.xml"


lalapps_inspinj \
    --seed ${seed} --f-lower 20 --gps-start-time ${gpsstart} \
    --gps-end-time ${gpsend} --waveform SEOBNRv2threePointFivePN \
    --amp-order 0 \
    --time-step 10 --time-interval 5 --l-distr random \
    --i-distr uniform --disable-spin \
    --m-distr totalMass --min-mass1 50 --max-mass1 100\
    --min-mass2 50 --max-mass2 100\
    --output ${outfile} \
    --snr-distr volume \
    --min-snr 15 --max-snr 30 \
    --ligo-fake-psd LALAdLIGO \
    --ligo-start-freq 16 \
    --ifos H1,L1 --verbose

source ~jclark/etc/bayeswave_ldg-user-env.sh


#
# Setup the pipeline.  
#
bwb_pipe.py cbc_xml_injection.ini \
    --workdir cbc_xml_injection  
