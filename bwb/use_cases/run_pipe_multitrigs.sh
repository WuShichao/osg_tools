#!/bin/bash
#
# Example call to the pipeline, using a list of trigger times in an ascii file.
# This would be useful for e.g., following up a list of CWB triggers. 
#
# Trigger times are (arbirarily chosen for example useage) 100 quasi-random
# times around GW150914.  
#
# To run quickly, the config file specifies --signalOnly and --noClean


bwb_pipe.py usecase_multitrigs.ini \
    --workdir usecase_multitrigs \
    --trigger-list usecase_multitrigs.asc
