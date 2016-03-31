#!/bin/bash

source ~/.local/etc/pycbc-user-env.sh

catalogdir=/data/lvc_nr/GaTech

h5files=`ls ${catalogdir}/GT0901.h5`
echo ${h5files}

pycbc_make_nr_hdf_catalog \
    --output-file GaTechCatalog.xml.gz \
    --input-files ${h5files}
