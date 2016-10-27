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

if [ -z ${BAYESWAVE_ROOT} ]
then
    echo "BAYESWAVE_ROOT is unset"
    echo "set e.g., BAYESWAVE_ROOT=${HOME}/data/src/lscsoft/bayeswave/trunk"
    return
else
    echo "BAYESWAVE_ROOT=${BAYESWAVE_ROOT}"
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
export BWB_PIPE_PATH=${HOME}/Projects/osg_tools/bwb
export PATH=${BWB_PIPE_PATH}/bin:${PATH}
export PYTHONPATH=${BWB_PIPE_PATH}/bin:${PYTHONPATH}







