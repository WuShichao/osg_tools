#!/bin/bash

source ${HOME}/Projects/osg_tools/bwb/etc/bayeswave_ldg-user-env.sh

workdir=${1}    
: ${workdir:?"Need to specify input directory (outputDir from bayeswave)"}

pushd ${workdir}
# Make skymap
python ${HOME}/src/lscsoft/bayeswave/trunk/postprocess/skymap/megasky.py

# Make the output page
python ${HOME}/src/lscsoft/bayeswave/trunk/postprocess/megaplot.py
popd

