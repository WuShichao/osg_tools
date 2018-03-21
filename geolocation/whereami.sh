#!/bin/bash

echo "Reading input file:"
cat $1

echo "Putting output in:"
outfile=$2

echo "Geolocating.."
curl ipinfo.io > ${outfile}
echo "Geolocation success, exiting"

