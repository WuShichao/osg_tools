#!/bin/bash
#
# Example call to the pipeline, using a CWB trigger list
rm -r HL
bwb_pipe.py config_triggers.ini \
    --workdir HL --sim-data\
    --cwb-trigger-list HL-EVENTS.txt

rm -r HLV
bwb_pipe.py config_triggers-Vslides.ini \
    --workdir HLV --sim-data\
    --cwb-trigger-list HLV-EVENTS.txt

