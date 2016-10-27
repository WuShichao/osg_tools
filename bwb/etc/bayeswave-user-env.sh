#!/bin/bash

# ----------------------------------------- #
#                   LAL                     #
# ----------------------------------------- #

if [ -d /cvmfs/oasis.opensciencegrid.org/ligo/pipeline ]
then

    echo "CVMFS detected: OSG deployment"
    export CVMFS_PIPEDIR=/cvmfs/oasis.opensciencegrid.org/ligo/pipeline

    # Set up BayesWave dependencies
    source ${CVMFS_PIPEDIR}/bayeswave/bayeswave_osg_deps-user-env.sh

    # Set LAL_DATA_PATH
    if [ -d ${CVMFS_PIPEDIR}/pycbc/pycbc/ROM/lal-data/lalsimulation ]
    then
        echo "Setting LAL_DATA_PATH=${CVMFS_PIPEDIR}/pycbc/pycbc/ROM/lal-data/lalsimulation"
        export LAL_DATA_PATH=${CVMFS_PIPEDIR}/pycbc/pycbc/ROM/lal-data/lalsimulation
    fi

fi


if [ ! -z ${LALSUITE_PREFIX} ]
then
    source ${LALSUITE_PREFIX}/etc/lalsuiterc
    source ${LALSUITE_PREFIX}/pylal/etc/pylal-user-env.sh
    source ${LALSUITE_PREFIX}/glue/etc/glue-user-env.sh
else
    echo "LALSUITE_PREFIX is unset, you better hope there's a global installation"
fi


# ----------------------------------------- #
#                 LAL_DATA_PATH             #
# ----------------------------------------- #

if [ ! -z ${LAL_DATA_PATH} ]
then 
    echo "Using LAL_DATA_PATH=${LAL_DATA_PATH}"
else
    echo "!!! WARNING: LAL_DATA_PATH is unset !!!"
fi

# ----------------------------------------- #
#                 BAYESWAVE                 #
# ----------------------------------------- #

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

