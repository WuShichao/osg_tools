#!/bin/bash

source ../etc/bayeswave_ldg-user-env.sh

which bayeswave
which bayeswave_post

bayeswave \
    --ifo H1 --H1-flow 16 --H1-cache LALSimAdLIGO \
    --H1-channel LALSimAdLIGO  \
    --inj GaTechIMBBH.xml --event 0 \
    --srate 512 --seglen 4 \
    --trigtime 1126621184 \
    --PSDstart 1126621184 --PSDlength 1024 \
    --NCmin 2 --NCmax 2 --dataseed 1234 \
    --inj-numreldata GaTechCatalog.xml.gz
