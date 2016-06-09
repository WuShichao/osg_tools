#!/bin/bash

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

