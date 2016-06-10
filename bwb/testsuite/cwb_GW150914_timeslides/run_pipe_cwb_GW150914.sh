#!/bin/bash
#
# Example call to the pipeline, using trigger time from offline CWB CED and
# timeslides (see config file)

bwb_pipe.py cwb_GW150914_timeslides.ini \
    --workdir cwb_GW150914_timeslides \
    --trigger-time 1126259462.392
