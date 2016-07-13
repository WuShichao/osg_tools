#!/bin/bash
#
#
# This POST script REMOVES the work directory, which is now compressed
#
# Large runs should --skip-megapy

echo "Performing housekeeping"

outputDir=${1}

rm -rf ${outputDir}

echo "Tidying complete, time for tea"
