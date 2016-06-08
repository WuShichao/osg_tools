#!/bin/bash
#
# Example call to the pipeline, using trigger list and sim data
bwb_pipe.py simulated_noise.ini \
    --workdir simulated_noise \
    --trigger-list trigger_list.asc
