#!/bin/bash

echo 'universe = vanilla
executable = '`which megasky.py`'
arguments = "$(macrooutputDir)"
getenv = True
log = megasky_$(macrooutputDir).log
accounting_group = ligo.dev.o2.burst.paramest.bayeswave
error = $(macrooutputDir)/megasky_$(macrooutputDir).err
output = $(macrooutputDir)/megasky_$(macrooutputDir).out
should_transfer_files = YES
transfer_input_files $(macrooutputDir)
transfer_output_files $(macrooutputDir)
when_to_transfer_output = ON_EXIT_OR_EVICT
notification = never
queue 1
' > megasky.sub

echo 'universe = vanilla
executable = '`which megaplot.py`'
arguments = "$(macrooutputDir)"
getenv = True
log = megaplot_$(macrooutputDir).log
accounting_group = ligo.dev.o2.burst.paramest.bayeswave
error = $(macrooutputDir)/megaplot_$(macrooutputDir).err
output = $(macrooutputDir)/megaplot_$(macrooutputDir).out
should_transfer_files = YES
transfer_input_files $(macrooutputDir)
transfer_output_files $(macrooutputDir)
when_to_transfer_output = ON_EXIT_OR_EVICT
notification = never
queue 1
' > megaplot.sub

jobcounter=0
for outputdir in bayeswave_11*
do

    echo "JOB megasky_${jobcounter} megasky.sub" >> megasky.dag
    echo "RETRY megasky_${jobcounter} 0">> megasky.dag
    echo "VARS megasky_${jobcounter} macrooutputDir=\"${outputdir}\"" >> megasky.dag

    echo "JOB megaplot_${jobcounter} megaplot.sub" >> megaplot.dag
    echo "RETRY megaplot_${jobcounter} 0">> megaplot.dag
    echo "VARS megaplot_${jobcounter} macrooutputDir=\"${outputdir}\"" >> megaplot.dag

    jobcounter=$((${jobcounter}+1))
done

