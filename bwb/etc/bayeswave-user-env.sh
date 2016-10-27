#!/bin/bash


export LAL_DATA_PATH=${HOME}/data/ROM_data

# -------------------------------------------
#                 BAYESWAVE                 #
# -------------------------------------------

# BayesWave
export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/src
export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess
export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess/dist
export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess/skymap
export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess/skymap/dist
export PYTHONPATH=${PYTHONPATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess
export PYTHONPATH=${PYTHONPATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess/skymap
export PYTHONPATH=${PYTHONPATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/utils

# BayesWave Pipeline
export PATH=${PATH}:${HOME}/Projects/osg_tools/bwb/bin
export PYTHONPATH=${PYTHONPATH}:${HOME}/Projects/osg_tools/bwb/bin

