#!/bin/bash

echo 'universe = vanilla
executable = /home/jclark/Projects/pmnsdata/compute_peaks.sh
arguments = "$(macrooutputDir)"
getenv = True
log = bayeswave_peaks.log
accounting_group = ligo.dev.o2.burst.paramest.bayeswave
error = $(macrooutputDir)/bayeswave_peaks_$(macrooutputDir).err
output = $(macrooutputDir)/bayeswave_peaks_$(macrooutputDir).out
notification = never
queue 1
' > bayeswave_peaks.sub


jobcounter=0
for outputdir in bayeswave_11*
do

    echo "JOB bayeswave_peaks_${jobcounter} bayeswave_peaks.sub" >> bayeswave_peaks.dag
    echo "RETRY bayeswave_peaks_${jobcounter} 0">> bayeswave_peaks.dag
    echo "VARS bayeswave_peaks_${jobcounter} macrooutputDir=\"${outputdir}\"" >> bayeswave_peaks.dag

    jobcounter=$((${jobcounter}+1))
done

