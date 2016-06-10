#!/bin/bash
#
# Example call to the pipeline, using trigger time from offline CWB CED and
# timeslides (see config file).  Timeslides are specified in the config file as a python literal.
#
# E.g., L1-timeslides=[0.1, 0.5, 1]
#       L1-timeslides=np.arange(0.1,10, 0.1) 
# etc

bwb_pipe.py cwb_GW150914_timeslides.ini \
    --workdir cwb_GW150914_timeslides \
    --trigger-time 1126259462.392
