#!/bin/bash
#
# Example call to the pipeline, using trigger list and sim data
bwb_pipe.py cwb_trigger_list.ini \
    --workdir cwb_trigger_list \
    --trigger-list $1
