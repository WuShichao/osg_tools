#!/bin/bash

# Things to remove:
#
#
# 1. chains
# 2. checkpoint
# 3. waveforms
#
# Archive (bz2) the post directory
#
# Large runs should --skip-megapy

echo "Performing housekeeping"

outputDir=${1}

pushd ${outputDir}
echo "Removing: chains checkpoint waveforms"
for d in chains checkpoint waveforms
do
    rm -r ${d}
done

