#!/bin/bash


catalogdir="/home/jclark308/lvc_nr_xml/GaTech"
catalogfile="${1}.xml.gz"
totalmass=${2}

: ${catalogfile:?"Need to specify NR waveform (e.g., GT0001)"}
: ${totalmass:?"Need to specify total mass (e.g., 100)"}

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
    --i-distr uniform \
    --m-distr nrwaves --disable-spin \
    --min-mtotal ${totalmass} --max-mtotal ${totalmass}\
    --taper-injection start --output ${outfile} \
    --dchirp-distr uniform  --min-distance 50000 --max-distance 50000 #\
#       --snr-distr volume \
#       --min-snr 15 --max-snr 15 \
#       --ligo-psd aligopsd.txt \
#       --ligo-start-freq 30 \
#       --ifos H1,L1 \
#       --ninja-snr 

