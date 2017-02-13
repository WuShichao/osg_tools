#!/bin/bash

# ----------------------------------------- #
#                   LAL                     #
# ----------------------------------------- #


# User should set OSG_DEPLOY=true if desired
if [ ${OSG_DEPLOY} ]
then

    echo "OSG_DEPLOY=${OSG_DEPLOY}"
    export CVMFS_PATH=/cvmfs/oasis.opensciencegrid.org/ligo

    # Set up BayesWave dependencies
    source ${CVMFS_PATH}/pipeline/bayeswave/bayeswave_osg_deps-user-env.sh

    if [ -z ${BAYESWAVE_PIPE_PREFIX} ]
    then
        echo "BAYESWAVE_PIPE_PREFIX is unset"
        return
    else
        echo "BAYESWAVE_PREFIX=${BAYESWAVE_PIPE_PREFIX}"
    fi

    # Set LAL_DATA_PATH
    if [ -z ${LAL_DATA_PATH} ]
    then
        echo "LAL_DATA_PATH is unset, using ${CVMFS_PATH}/sw/pycbc/lalsuite-extra/current/share/lalsimulation"
        export LAL_DATA_PATH=${CVMFS_PATH}/sw/pycbc/lalsuite-extra/current/share/lalsimulation
    else
        echo "LAL_DATA_PATH=${LAL_DATA_PATH}"
    fi

fi


if [ ! -z ${LALSUITE_PREFIX} ]
then
    echo "Using LALSUITE_PREFIX=${LALSUITE_PREFIX}"
    source ${LALSUITE_PREFIX}/etc/lalsuite-user-env.sh
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

if [ -z ${BAYESWAVE_PREFIX} ]
then
    echo "BAYESWAVE_PREFIX is unset"
    echo "set e.g., BAYESWAVE_PREFIX=${HOME}/data/src/lscsoft/bayeswave/trunk"
    return
else
    echo "BAYESWAVE_PREFIX=${BAYESWAVE_PREFIX}"
fi

# BayesWave
export PATH=${BAYESWAVE_PREFIX}/src:${PATH}
export PATH=${BAYESWAVE_PREFIX}/postprocess:${PATH}
export PATH=${BAYESWAVE_PREFIX}/postprocess/dist:${PATH}
export PATH=${BAYESWAVE_PREFIX}/postprocess/skymap:${PATH}
export PATH=${BAYESWAVE_PREFIX}/postprocess/skymap/dist:${PATH}
export PYTHONPATH=${BAYESWAVE_PREFIX}/postprocess:${PYTHONPATH}
export PYTHONPATH=${BAYESWAVE_PREFIX}/postprocess/skymap:${PYTHONPATH}
export PYTHONPATH=${BAYESWAVE_PREFIX}/utils:${PYTHONPATH}

# BayesWave Pipeline
if [ -z ${BAYESWAVE_PIPE_PREFIX} ]
then
    echo "BAYESWAVE_PIPE_PREFIX is unset"
    echo "set e.g., BAYESWAVE_PIPE_PREFIX=${HOME}/Projects/osg_tools/bwb"
    return
else
    echo BAYESWAVE_PIPE_PREFIX=${BAYESWAVE_PIPE_PREFIX}
fi
export PATH=${BAYESWAVE_PIPE_PREFIX}/bin:${PATH}
export PYTHONPATH=${BAYESWAVE_PIPE_PREFIX}/bin:${PYTHONPATH}







