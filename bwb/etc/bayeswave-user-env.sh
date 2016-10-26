#!/bin/bash

LIGO_LIBS="/cvmfs/oasis.opensciencegrid.org/ligo/pipeline/bayeswave"

lscsoft="${LIGO_LIBS}/lscsoft/"
lalsuite="${lscsoft}lalsuite-6.38"

source ${lalsuite}/etc/lalsuiterc
#source ${lalsuite}/pylal/etc/pylal-user-env.sh
#source ${lalsuite}/glue/etc/glue-user-env.sh

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

# -------------------------------------------
#               PYTHON & FRIENDS            #
# -------------------------------------------

# Python
export PATH=${LIGO_LIBS}/non-lsc/Python-2.7.5/bin:${PATH}
export LD_LIBRARY_PATH=${LIGO_LIBS}/non-lsc/Python-2.7.5/lib:${LD_LIBRARY_PATH}
export PKG_CONFIG_PATH=${LIGO_LIBS}/non-lsc/Python-2.7.5/lib/pkgconfig:${PKG_CONFIG_PATH}

# Numpy
export PATH=${LIGO_LIBS}/non-lsc/numpy-1.9.1/bin:${PATH}
export PYTHONPATH=${LIGO_LIBS}/non-lsc/numpy-1.9.1/lib/python2.7/site-packages:${PYTHONPATH}

# Scipy
export PATH=${LIGO_LIBS}/non-lsc/scipy-0.12.1/bin:${PATH}
export PYTHONPATH=${LIGO_LIBS}/non-lsc/scipy-0.12.1/lib/python2.7/site-packages:${PYTHONPATH}

# Matplotlib
export PATH=${LIGO_LIBS}/non-lsc/matplotlib-1.2.0/bin:${PATH}
export PYTHONPATH=${LIGO_LIBS}/non-lsc/matplotlib-1.2.0/lib/python2.7/site-packages:${PYTHONPATH}
export GEOS_DIR=${LIGO_LIBS}/non-lsc/geos-3.5.0
export PYTHONPATH=${LIGO_LIBS}/non-lsc/basemap-1.0.7/lib/python2.7/site-packages:${PYTHONPATH}

# Misc dependencies
export PYTHONPATH=${LIGO_LIBS}/non-lsc/setuptools-28.6.1/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${LIGO_LIBS}/non-lsc/healpy-1.9.1/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${LIGO_LIBS}/non-lsc/acor-1.1.1/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${LIGO_LIBS}/non-lsc/h5py-2.6.0/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${LIGO_LIBS}/non-lsc/skyarea-0.1/lib/python2.7/site-packages:${PYTHONPATH}

# -------------------------------------------
#           LALSUITE DEPENDENCIES           #
# -------------------------------------------
# Note: not all of these are required at run-time but all are required for
# building (e.g., swig)

# SWIGPATH (only necessary for building)
export SWIGPATH=${LIGO_LIBS}/non-lsc/swig-3.0.7
export PATH=${SWIGPATH}/bin:${PATH}

# FFTW
export FFTW=${LIGO_LIBS}/non-lsc/fftw-3.3.5
export PATH=${FFTW}/bin:${PATH}
export PKG_CONFIG_PATH=${FFTW}/lib/pkgconfig:${PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=${FFTW}/lib:${LD_LIBRARY_PATH}

# GSL
export GSL=${LIGO_LIBS}/non-lsc/gsl-1.15
export PATH=${GSL}/bin:${PATH}
export PKG_CONFIG_PATH=${GSL}/lib/pkgconfig:${PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=${GSL}/lib:${LD_LIBRARY_PATH}

# LIBFRAME
export LIBFRAME=${LIGO_LIBS}/lscsoft/libframe-8.30
export PATH=${LIBFRAME}/bin:${PATH}
export PKG_CONFIG_PATH=${LIBFRAME}/lib/pkgconfig:${PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=${LIBFRAME}/lib:${LD_LIBRARY_PATH}

# LIBMETAIO
export LIBMETAIO=${LIGO_LIBS}/lscsoft/libmetaio-8.4.0
export PATH=${LIBMETAIO}/bin:${PATH}
export PKG_CONFIG_PATH=${LIBMETAIO}/lib/pkgconfig:${PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=${LIBMETAIO}/lib:${LD_LIBRARY_PATH}

# HDF5
export HDF5=$LIGO_LIBS/non-lsc/hdf5-1.8.16
export PATH=${HDF5}/bin:${PATH}
export LD_LIBRARY_PATH=${HDF5}/lib:${LD_LIBRARY_PATH}
export C_INCLUDE_PATH=${HDF5}/include:${C_INCLUDE_PATH} # for building

# CFITSIO
export CFITSIO=${LIGO_LIBS}/non-lsc/cfitsio-3.39
export LD_LIBRARY_PATH=${CFITSIO}/lib:${LD_LIBRARY_PATH}
export PKG_CONFIG_PATH=${CFITSIO}/lib/pkgconfig:${PKG_CONFIG_PATH}

# HEALPIX
export HEALPIX=${LIGO_LIBS}/non-lsc/Healpix_3.30
export LD_LIBRARY_PATH=${HEALPIX}/lib:${LD_LIBRARY_PATH}
export PKG_CONFIG_PATH=${HEALPIX}/lib:${PKG_CONFIG_PATH}


