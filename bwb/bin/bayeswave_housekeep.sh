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

echo "Compressing post directory"
tar -cjf post.tar.bz2 post

echo "Removing old post directory"
rm -r post

echo "Housekeeping complete, time for tea"
popd
