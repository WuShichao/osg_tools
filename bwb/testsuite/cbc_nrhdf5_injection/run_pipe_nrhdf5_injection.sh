#!/bin/bash
#
########################################################################
#                                                                      #
# Example BayesWave pipeline run for NR injection in SIMULATED data.   #   
#                                                                      # 
# RUN TIME XX hours                                                    #  
#  - Processor: 2 GHz Intel Core i7                                    #
#  - Memory: 4 GB 1600 Mhz                                             #
#                                                                      #
########################################################################

source ~jclark/etc/source_master.sh

#
# Generate sim-inspiral table using lalapps_inspinj and NR data
#

catalogfile="GT0887.xml.gz"
seed=`lalapps_tconvert now`
gpsstart=${seed} 
gpsend=$((${gpsstart}+100))

outfile="HL-INJECTIONS-NR_`echo ${catalogfile} | sed 's/.xml.gz/.xml/g'`"

lalapps_inspinj \
    --seed ${seed} --f-lower 30 --gps-start-time ${gpsstart} \
    --gps-end-time ${gpsend} --waveform NR_hdf5threePointFivePN \
    --amp-order 0 \
    --ninja2-mass --nr-file ${catalogfile} \
    --time-step 10 --time-interval 5 --l-distr random \
    --i-distr uniform \
    --m-distr nrwaves --disable-spin \
    --min-mtotal 100 --max-mtotal 500 \
    --taper-injection start --output ${outfile} \
    --dchirp-distr uniform  --min-distance 50000 --max-distance 100000 #\

#
# Setup the pipeline.  
#
source ~jclark/etc/bayeswave_ldg-user-env.sh

bwb_pipe.py cbc_nrhdf5_injection.ini \
    --workdir cbc_nrhdf5_injection  
