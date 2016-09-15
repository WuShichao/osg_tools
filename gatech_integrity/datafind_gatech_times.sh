#!/bin/bash


# --- ER8 
# H1 ER8 h(t) H1_HOFT_C00 1123856384 1126621184

echo "************* ER8 H1_HOFT_C00 *****************"
gw_data_find --observatory H --type H1_HOFT_C00 -s 1123856384 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

echo "************* ER8 L1_HOFT_C00 *****************"
# L1 ER8 h(t) L1_HOFT_C00 1123856384 1126621184
gw_data_find --observatory L --type L1_HOFT_C00 -s 1123856384 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

echo "************* ER8 H1_HOFT_C01 *****************"
# H1 ER8 h(t) H1_HOFT_C01 1125969920 1126621184
gw_data_find --observatory H --type H1_HOFT_C01 -s 1125969920 -e 1126621184 \
    --server=ligo-gftp.pace.gatech.edu:80 --show-times

echo "************* ER8 L1_HOFT_C01 *****************"
# L1 ER8 h(t) L1_HOFT_C01 1126031360 1126621184
gw_data_find --observatory L --type L1_HOFT_C01 -s 1126031360 -e 1126621184 \
    --server=ligo-gftp.pace.gatech.edu:80 --show-times
