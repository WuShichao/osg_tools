#!/bin/bash

export PATH=${PATH}:${HOME}/Projects/osg_tools/bwb
export PYTHONPATH=${PYTHONPATH}:${HOME}/Projects/osg_tools/bwb

bwb_pipe.py signalonlyinjections.ini \
    --user-tag GT0420 \
    --workdir GT0420 
