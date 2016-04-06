#!/bin/bash


catalogdir="/home/jclark308/lvc_nr_xml/GaTech"
catalogfile="${1}.xml.gz"

seed=`lalapps_tconvert now`
gpsstart=1126621184 # first CO2 frame for O1
gpsend=`python -c "print ${gpsstart} + 1000"`
outfile=`echo ${catalogfile} | sed 's/.xml.gz/_siminsp.xml/g'`


lalapps_inspinj \
    --seed ${seed} --f-lower 30 --gps-start-time ${gpsstart} \
    --gps-end-time ${gpsend} --waveform NR_hdf5threePointFivePN \
    --amp-order 0 \
    --ninja2-mass --nr-file ${catalogdir}/${catalogfile} \
    --time-step 1 --time-interval 0 --l-distr random \
    --dchirp-distr uniform --i-distr uniform \
    --min-distance 50000 --max-distance 50000 \
    --m-distr nrwaves --disable-spin \
    --min-mtotal 70 --max-mtotal 70\
    --taper-injection start --output ${outfile}

