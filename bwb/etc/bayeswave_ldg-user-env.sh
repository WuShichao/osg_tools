#!/bin/bash

export PATH=${PATH}:${HOME}/src/lscsoft/bayeswave/trunk/src
export PATH=${PATH}:${HOME}/Projects/osg_tools/bwb/bin
export LIBRARY_PATH=${HOME}/opt/lscsoft/lalsuite_master/lib:/usr/lib64:/usr/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=${HOME}/opt/lalsuite_master/lib:/usr/lib64:/usr/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=${HOME}/opt/lscsoft/lalsuite_master/include:$C_INCLUDE_PATH
export PYTHONPATH=${PYTHONPATH}:${HOME}/src/lscsoft/bayeswave/trunk/utils

export C_INCLUDE_PATH=${HOME}/opt/lscsoft/lalsuite_master/include:/usr/include:/usr/include/gsl:$C_INCLUDE_PATH
export LIBRARY_PATH=${HOME}/opt/lscsoft/lalsuite_master/lib:/usr/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=${HOME}/opt/lscsoft/lalsuite_master/lib:/usr/lib64:$LD_LIBRARY_PATH

