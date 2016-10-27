#!/bin/bash

# ----------------------------------------- #
#                   LAL                     #
# ----------------------------------------- #

if [ -d /cvmfs/oasis.opensciencegrid.org/ligo/pipeline ]

    echo "CVMFS detected, assuming OSG deployment"
    OSG_DEPLOY=1
    CVMFS_PIPEDIR=/cvmfs/oasis.opensciencegrid.org/ligo/pipeline

    # Set up BayesWave dependencies
    if [ -e ${CVMFS_PIPEDIR}/bayeswave/bayeswave_deps-user-env.sh ]
    then 
        echo "BayesWave dependency setup detected in /cvmfs. Sourcing bayeswave_osg_deps-user-env.sh"
        source ${CVMFS_PIPEDIR}/bayeswave/bayeswave_osg_deps-user-env.sh
    else
        echo "!!! WARNING: no bayeswave dependency setup found in ${CVMFS} !!!"
    fi

    # Set LAL_DATA_PATH
    echo "Setting LAL_DATA_PATH=${CVMFS_PIPEDIR}/pycbc/pycbc/ROM/lal-data/lalsimulation"
    export LAL_DATA_PATH=${CVMFS_PIPEDIR}/pycbc/pycbc/ROM/lal-data/lalsimulation

else
    OSG_DEPLOY=0
fi


if [ -z ${LALSUITE_PREFIX} ]
then
    echo "LALSUITE_PREFIX is set to ${LALSUITE_PREFIX}, using that installation"
    source ${LALSUITE_PREFIX}/etc/lalsuiterc
    source ${LALSUITE_PREFIX}/pylal/etc/pylal-user-env.sh
    source ${LALSUITE_PREFIX}/glue/etc/glue-user-env.sh
else
    echo "LALSUITE_PREFIX is unset, you better hope there's a global installation"
fi


# ----------------------------------------- #
#                 LAL_DATA_PATH             #
# ----------------------------------------- #

if [ ${OSG_DEPLOY} ]
then
    export LAL_DATA_PATH=${CVMFS_PIPEDIR}/pycbc/pycbc/ROM/lal-data/lalsimulation
    echo "Setting LAL_DATA_PATH=${LAL_DATA_PATH}"
else
    if [ -z ${LAL_DATA_PATH} ]
    then 
        echo "Using LAL_DATA_PATH=${LAL_DATA_PATH}"
    else
        echo "!!! WARNING: LAL_DATA_PATH is unset !!!"
    fi
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

