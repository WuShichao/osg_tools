#!/bin/bash
#
# Example call to the pipeline, using trigger list and sim data
bwb_pipe.py simulated_noise.ini \
    --workdir simulated_noise_single \
    --trigger-time 1126258467.576778 \
    --sim-data
    #--condor-submit
    #--trigger-list trigger_list.asc \
