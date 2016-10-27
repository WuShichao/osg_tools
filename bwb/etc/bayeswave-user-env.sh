#!/bin/bash

# -------------------------------------------
#                   LAL                     #
# -------------------------------------------

if [ -e /cvmfs/oasis.opensciencegrid.org/ligo/pipeline/bayeswave/bayeswave_deps-user-env.sh ]
then 
    echo "BayesWave dependency setup detected in /cvmfs. Sourcing bayeswave_deps-user-env.sh"
elif [ -z ${LALSUITE_PREFIX} ]
    echo "LALSUITE_PREFIX is set to ${LALSUITE_PREFIX}, using that installation"
    source ${LALSUITE_PREFIX}/etc/lalsuiterc
    source ${LALSUITE_PREFIX}/pylal/etc/pylal-user-env.sh
    source ${LALSUITE_PREFIX}/glue/etc/glue-user-env.sh
else
    echo "LALSUITE_PREFIX is unset, better hope there's a global installation"
fi

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
export BWB_PIPE_PATH=${HOME}/Projects/osg_tools/bwb
export PATH=${PATH}:${BWB_PIPE_PATH}/bin
export PYTHONPATH=${PYTHONPATH}:${BWB_PIPE_PATH}/bin

