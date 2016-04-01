#!/bin/bash

export PATH=${PATH}:${HOME}/Projects/osg_tools/bwb
export PYTHONPATH=${PYTHONPATH}:${HOME}/Projects/osg_tools/bwb

bwb_pipe.py gw150914_bw.ini \
    --user-tag multijobs \
    --workdir multijobs 
    #--trigger-list trigtimes.txt
