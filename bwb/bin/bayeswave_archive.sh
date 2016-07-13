#!/bin/bash
#
#
# This script archives important output directories to reduce disk usage
# footpring
#
# Things to remove:
#
# 1. chains
# 2. checkpoint
# 3. waveforms
#
# Archive (bz2) the post directory
#
# Large runs should --skip-megapy
#
# This can be run as a DAG node

echo "Performing housekeeping"

outputDir=${1}

pushd ${outputDir}
echo "Removing: chains checkpoint waveforms"
for d in chains checkpoint waveforms snr evidence.gpi
do
    if [ -d "${d}" ]
    then
        rm -r ${d}
    fi
done

popd

echo "Compressing results directory ${outputDir}"
tar -cjf ${outputDir}.tar.bz2 ${outputDir}

echo "Housekeeping complete, time for tea"
popd
