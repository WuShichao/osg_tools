#!/bin/bash

LIGO_LIBS="/cvmfs/oasis.opensciencegrid.org/ligo"

lalsuite="${LIGO_LIBS}/lalsuite_master"
source ${lalsuite}/etc/lalsuiterc

export LAL_DATA_PATH=${HOME}/data/ROM_data


export PATH=${LIGO_LIBS}/swig-3.0.8/bin:${PATH}
export PATH=${LIGO_LIBS}/hdf5-1.8.16/bin:${PATH}
export PATH=${LIGO_LIBS}/gsl-2.0/bin:${PATH}
export PATH=${LIGO_LIBS}/${lalsuite}:${PATH}

export LIBRARY_PATH=${LIGO_LIBS}/FFTW3/fftw3.3.4/lib:${LIBRARY_PATH}
export LIBRARY_PATH=${LIGO_LIBS}/FFTW3F/lib:${LIBRARY_PATH}
export LIBRARY_PATH=${LIGO_LIBS}/gsl-2.0/lib:${LIBRARY_PATH}
export LIBRARY_PATH=${LIGO_LIBS}/hdf5-1.8.16/lib:${LIBRARY_PATH}
export LIBRARY_PATH=${LIGO_LIBS}/libframe-8.20/lib:${LIBRARY_PATH}
export LIBRARY_PATH=${LIGO_LIBS}/metaio-8.3.0/lib:${LIBRARY_PATH}
export LIBRARY_PATH=${lalsuite}/lib:/usr/lib64:/usr/lib:$LIBRARY_PATH

export LD_LIBRARY_PATH=${LIGO_LIBS}/FFTW3/fftw3.3.4/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${LIGO_LIBS}/FFTW3F/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${LIGO_LIBS}/gsl-2.0/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${LIGO_LIBS}/hdf5-1.8.16/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${LIGO_LIBS}/libframe-8.20/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${LIGO_LIBS}/metaio-8.3.0/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${lalsuite}/lib:/usr/lib64:/usr/lib:$LD_LIBRARY_PATH

export LD_RUN_PATH=${LIGO_LIBS}/gsl-2.0/lib:${LD_RUN_PATH}
export LD_RUN_PATH=${LIGO_LIBS}/FFTW3/fftw3.3.4/lib:${LD_RUN_PATH}
export LD_RUN_PATH=${LIGO_LIBS}/FFTW3F/lib:${LD_RUN_PATH}

export C_INCLUDE_PATH=${LIGO_LIBS}/FFTW3/fftw3.3.4/include:${C_INCLUDE_PATH}
export C_INCLUDE_PATH=${LIGO_LIBS}/FFTW3F/include:${C_INCLUDE_PATH}
export C_INCLUDE_PATH=${LIGO_LIBS}/gsl-2.0/include/gsl:${C_INCLUDE_PATH}
export C_INCLUDE_PATH=${LIGO_LIBS}/hdf5-1.8.16/include:${C_INCLUDE_PATH}
export C_INCLUDE_PATH=${LIGO_LIBS}/libframe-8.20/include:${C_INCLUDE_PATH}
export C_INCLUDE_PATH=${LIGO_LIBS}/metaio-8.3.0/include:${C_INCLUDE_PATH}
export C_INCLUDE_PATH=${lalsuite}/include:/usr/include:/usr/include/gsl:$C_INCLUDE_PATH

export PKG_CONFIG_PATH=${LIGO_LIBS}/libframe-8.20/lib/pkgconfig:${PKG_CONFIG_PATH}
export PKG_CONFIG_PATH=${LIGO_LIBS}/gsl-2.0/lib/pkgconfig:${PKG_CONFIG_PATH}
export PKG_CONFIG_PATH=${lalsuite}/lib/pkgconfig:${PKG_CONFIG_PATH}
export PKG_CONFIG_PATH=${LIGO_LIBS}/metaio-8.3.0/lib/pkgconfig:${PKG_CONFIG_PATH}
#export PKG_CONFIG_PATH=${LIGO_LIBS}/FFTW3/fftw3.3.4/lib/pkgconfig:${PKG_CONFIG_PATH}

export PYTHONPATH=${LIGO_LIBS}/h5py/lib64/python2.6/site-packages:${PYTHONPATH}
export PYTHONPATH=${PYTHONPATH}:${HOME}/src/lscsoft/bayeswave/trunk/utils
export PYTHONPATH=${PYTHONPATH}:${HOME}/Projects/osg_tools/bwb/bin

export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/src
export PATH=${PATH}:${HOME}/Projects/osg_tools/bwb/bin

source ${HOME}/opt/lscsoft/glue/etc/glue-user-env.sh
source ${HOME}/opt/lscsoft/pylal/etc/pylal-user-env.sh

export LAL_DATA_PATH=${HOME}/data/ROM_data

source /nv/hp11/jclark308/.local/etc/dqsegdb-user-env.sh

