#!/bin/bash
#
# Example call to the pipeline, using a list of trigger times in an ascii file.
# This would be useful for e.g., following up a list of CWB triggers. 
#
# Trigger times are (arbirarily chosen for example useage) 100 quasi-random
# times around GW150914.  We therefore use the same config file as the example
# GW150914 analysis.


bwb_pipe.py usecase_GW150914-C02.ini \
    --workdir usecase_GW150914-C02 \
    --trigger-list usecase_multitrigs.asc
