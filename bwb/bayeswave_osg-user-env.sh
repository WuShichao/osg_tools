#!/bin/bash

export C_INCLUDE_PATH=${LIGO_LIBS}/lalsuite/master/include:/usr/include:/usr/include/gsl:$C_INCLUDE_PATH
export LIBRARY_PATH=${LIGO_LIBS}/lalsuite/master/lib:/usr/lib64:/usr/lib:$LIBRARY_PATH
export LIBRARY_PATH=${LIBRARY_PATH}:/cvmfs/oasis.opensciencegrid.org/ligo/framelib_v8r19p2_root-v5-34-08.patched/Linux-x86_64
export LIBRARY_PATH=${LIBRARY_PATH}:/cvmfs/oasis.opensciencegrid.org/ligo/metaio-8.3.0/lib
export LIBRARY_PATH=${LIBRARY_PATH}:/cvmfs/oasis.opensciencegrid.org/osg/modules/gsl/1.16/lib/
export LIBRARY_PATH=${LIBRARY_PATH}:/cvmfs/oasis.opensciencegrid.org/ligo/FFTW3/fftw3.3.4/lib
export LD_LIBRARY_PATH=${LIGO_LIBS}/lalsuite/master/lib:/usr/lib64:/usr/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cvmfs/oasis.opensciencegrid.org/ligo/framelib_v8r19p2_root-v5-34-08.patched/Linux-x86_64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cvmfs/oasis.opensciencegrid.org/ligo/metaio-8.3.0/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cvmfs/oasis.opensciencegrid.org/osg/modules/gsl/1.16/lib/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cvmfs/oasis.opensciencegrid.org/ligo/FFTW3/fftw3.3.4/lib
