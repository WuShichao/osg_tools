#!/bin/bash

LIGO_LIBS="/cvmfs/oasis.opensciencegrid.org/ligo/pipeline/bayeswave"

lscsoft="${LIGO_LIBS}/lscsoft/"
lalsuite="${lscsoft}lalsuite-v6.38"

source ${lscsoft}/lscsoft-user-env.sh
source ${lalsuite}/etc/lalsuiterc
source ${lalsuite}/pylal/etc/pylal-user-env.sh
source ${lalsuite}/glue/etc/glue-user-env.sh

export LAL_DATA_PATH=${HOME}/data/ROM_data

export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/src
export PATH=${PATH}:${HOME}/Projects/osg_tools/bwb/bin
export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess
export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess/dist
export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess/skymap
export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess/skymap/dist

#export PYTHONPATH=${PYTHONPATH}:${LIGO_LIBS}/ligo-gracedb-1.20/lib/python2.6/site-packages
export PYTHONPATH=${PYTHONPATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess
export PYTHONPATH=${PYTHONPATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess/skymap
export PYTHONPATH=${PYTHONPATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/utils
export PYTHONPATH=${PYTHONPATH}:${HOME}/Projects/osg_tools/bwb/bin

export PYTHONPATH=${PYTHONPATH}:/nv/hp11/jclark308/opt/python-cjson-1.1.0/lib64/python2.6/site-packages

#source /nv/hp11/jclark308/.local/etc/dqsegdb-user-env.sh


