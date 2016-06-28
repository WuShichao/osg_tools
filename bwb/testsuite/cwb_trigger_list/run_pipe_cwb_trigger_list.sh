#!/bin/bash
#
# Example call to the pipeline, using a CWB trigger list
bwb_pipe.py config_triggers.ini \
    --workdir external_triggers \
    --cwb-trigger-list cwb_background.asc

# Example call to the pipeline, using one of the other trigger lists
bwb_pipe.py config_triggers.ini \
    --workdir external_triggers \
    --trigger-list trigger_times-and-lags-and-freqs.asc
