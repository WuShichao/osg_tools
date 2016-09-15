#!/bin/bash


# --- ER8 
# H1 ER8 h(t) H1_HOFT_C00 1123856384 1126621184

echo "************* H *****************"
gw_data_find --observatory H --type H1_HOFT_C00 -s 1123856384 -e 1126621184 \
    --show-times

echo "************* L *****************"
# L1 ER8 h(t) L1_HOFT_C00 1123856384 1126621184
gw_data_find --observatory L --type L1_HOFT_C00 -s 1123856384 -e 1126621184 \
    --show-times

