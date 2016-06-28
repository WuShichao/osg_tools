#!/bin/bash
#
# Example call to the pipeline, using trigger list and sim data
bwb_pipe.py config_triggers.ini \
    --workdir external_triggers \
    --cwb-trigger-list $1
