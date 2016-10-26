#!/bin/bash

LIGO_LIBS="/cvmfs/oasis.opensciencegrid.org/ligo/pipeline/bayeswave"

lscsoft="${LIGO_LIBS}/lscsoft/"
lalsuite="${lscsoft}lalsuite-v6.38"

export LAL_DATA_PATH=${HOME}/data/ROM_data

export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/src
export PATH=${PATH}:${HOME}/Projects/osg_tools/bwb/bin
export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess
export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess/dist
export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess/skymap
export PATH=${PATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess/skymap/dist
export PATH=${LIGO_LIBS}/Python-2.7.5/bin:${PATH}
export PATH=${LIGO_LIBS}/numpy-1.9.1/bin:${PATH}
export PATH=${LIGO_LIBS}/scipy-0.12.1/bin:${PATH}
export PATH=${LIGO_LIBS}/matplotlib-1.2.0/bin:${PATH}

unset PYTHONPATH
#export PYTHONPATH=${PYTHONPATH}:${LIGO_LIBS}/ligo-gracedb-1.20/lib/python2.6/site-packages
export PYTHONPATH=${PYTHONPATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess
export PYTHONPATH=${PYTHONPATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/postprocess/skymap
export PYTHONPATH=${PYTHONPATH}:${HOME}/data/src/lscsoft/bayeswave/trunk/utils
export PYTHONPATH=${PYTHONPATH}:${HOME}/Projects/osg_tools/bwb/bin

#export PYTHONPATH=${PYTHONPATH}:${HOME}/opt/python-cjson-1.1.0/lib64/python2.6/site-packages
export PYTHONPATH=${LIGO_LIBS}/numpy-1.9.1/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${LIGO_LIBS}/scipy-0.12.1/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${LIGO_LIBS}/matplotlib-1.2.0/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${LIGO_LIBS}/setuptools-28.6.1/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${LIGO_LIBS}/healpy-1.9.1/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${LIGO_LIBS}/acor-1.1.1/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${LIGO_LIBS}/h5py-2.6.0/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${LIGO_LIBS}/basemap-1.0.7/lib/python2.7/site-packages:${PYTHONPATH}
#export PYTHONPATH=${LIGO_LIBS}/skyarea-0.2.1/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=${LIGO_LIBS}/skyarea-0.1/lib/python2.7/site-packages:${PYTHONPATH}


export LD_LIBRARY_PATH=$LIGO_LIBS/Python-2.7.5/lib:${LD_LIBRARY_PATH}

# HDF5
export HDF5=$LIGO_LIBS/hdf5-1.8.16
export PATH=${HDF5}/bin:${PATH}
export LD_LIBRARY_PATH=${HDF5}/lib:${LD_LIBRARY_PATH}
export LIBRARY_PATH=${HDF5}/lib:${LIBRARY_PATH}

# CFITSIO
export CFITSIO=${LIGO_LIBS}/cfitsio-3.39
export LD_LIBRARY_PATH=${CFITSIO}/lib:${LD_LIBRARY_PATH}
export PKG_CONFIG_PATH=${CFITSIO}/lib/pkgconfig:${PKG_CONFIG_PATH}

# HEALPIX
export HEALPIX=${LIGO_LIBS}/Healpix_3.31
export LD_LIBRARY_PATH=${HEALPIX}/lib:${LD_LIBRARY_PATH}
export LIBRARY_PATH=${HEALPIX}/lib:${LIBRARY_PATH}
export PKG_CONFIG_PATH=${HEALPIX}/lib:${PKG_CONFIG_PATH}
export C_INCLUDE_PATH=${HEALPIX}/include:${C_INCLUDE_PATH}

export GEOS_DIR=${LIGO_LIBS}/geos-3.5.0

source ${lscsoft}/lscsoft-user-env.sh
source ${lalsuite}/etc/lalsuiterc
source ${lalsuite}/pylal/etc/pylal-user-env.sh
source ${lalsuite}/glue/etc/glue-user-env.sh

