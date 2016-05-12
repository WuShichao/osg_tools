#!/bin/bash

bwb_pipe.py ${1}.ini \
    --user-tag ${1} \
    --workdir ${1} \
    --trigger-time 1126259462.39 \
    --skip-datafind
