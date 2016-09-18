#!/bin/bash


# --- ER8 
# H1 ER8 h(t) H1_HOFT_C00 1123856384 1126621184

echo "************* ER8 H1_HOFT_C00 *****************"
gw_data_find --observatory H --type H1_HOFT_C00 -s 1123856384 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80 

gw_data_find --observatory H --type H1_HOFT_C00 -s 1123856384 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80  | awk '{ sum+=$4} END {print sum}'



echo "************* ER8 L1_HOFT_C00 *****************"
# L1 ER8 h(t) L1_HOFT_C00 1123856384 1126621184
gw_data_find --observatory L --type L1_HOFT_C00 -s 1123856384 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

gw_data_find --observatory L --type L1_HOFT_C00 -s 1123856384 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80  | awk '{ sum+=$4} END {print sum}'

echo "************* ER8 H1_HOFT_C01 *****************"
# H1 ER8 h(t) H1_HOFT_C01 1125969920 1126621184
gw_data_find --observatory H --type H1_HOFT_C01 -s 1125969920 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

gw_data_find --observatory H --type H1_HOFT_C01 -s 1125969920 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80  | awk '{ sum+=$4} END {print sum}'

echo "************* ER8 L1_HOFT_C01 *****************"
# L1 ER8 h(t) L1_HOFT_C01 1126031360 1126621184
gw_data_find --observatory L --type L1_HOFT_C01 -s 1126031360 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

gw_data_find --observatory L --type L1_HOFT_C01 -s 1126031360 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80  | awk '{ sum+=$4} END {print sum}'

echo "************* ER8 H1_HOFT_C02 *****************"
# H1 ER8 h(t) H1_HOFT_C01 1125969920 1126621184
gw_data_find --observatory H --type H1_HOFT_C02 -s 1125969920 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

gw_data_find --observatory H --type H1_HOFT_C02 -s 1125969920 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80  | awk '{ sum+=$4} END {print sum}'

echo "************* ER8 L1_HOFT_C02 *****************"
# L1 ER8 h(t) L1_HOFT_C01 1126031360 1126621184
gw_data_find --observatory L --type L1_HOFT_C02 -s 1126031360 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

gw_data_find --observatory L --type L1_HOFT_C02 -s 1126031360 -e 1126621184 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80  | awk '{ sum+=$4} END {print sum}'

echo "************* O1 H1_HOFT_C00 *****************"
# --- O1
# H1 01 h(t) H1_HOFT_C00 1126621184 1137258496
gw_data_find --observatory H --type H1_HOFT_C00 -s 1126621184 -e 1137258496 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

gw_data_find --observatory H --type H1_HOFT_C00 -s 1126621184 -e 1137258496 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80  | awk '{ sum+=$4} END {print sum}'

echo "************* O1 L1_HOFT_C00 *****************"
# L1 01 h(t) L1_HOFT_C00 1126621184 1137254400
gw_data_find --observatory L --type L1_HOFT_C00 -s 1126621184 -e 1137254400 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

gw_data_find --observatory L --type L1_HOFT_C00 -s 1126621184 -e 1137254400 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80  | awk '{ sum+=$4} END {print sum}'

echo "************* O1 H1_HOFT_C01 *****************"
# H1 01 h(t) H1_HOFT_C01 1126621184 1137258496
gw_data_find --observatory H --type H1_HOFT_C01 -s 1126621184 -e 1137258496 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

gw_data_find --observatory H --type H1_HOFT_C01 -s 1126621184 -e 1137258496 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80  | awk '{ sum+=$4} END {print sum}'

echo "************* O1 L1_HOFT_C01 *****************"
# L1 01 h(t) L1_HOFT_C01 1126621184 1137258496
gw_data_find --observatory L --type L1_HOFT_C01 -s 1126621184 -e 1137258496 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

gw_data_find --observatory L --type L1_HOFT_C01 -s 1126621184 -e 1137258496 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80  | awk '{ sum+=$4} END {print sum}'

echo "************* O1 H1_HOFT_C02 *****************"
# H1 01 h(t) H1_HOFT_C02 1126621184 1137258496
gw_data_find --observatory H --type H1_HOFT_C02 -s 1126621184 -e 1137258496 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

gw_data_find --observatory H --type H1_HOFT_C02 -s 1126621184 -e 1137258496 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80  | awk '{ sum+=$4} END {print sum}'

echo "************* O1 L1_HOFT_C02 *****************"
# L1 01 h(t) L1_HOFT_C02 1126621184 1137258496
gw_data_find --observatory L --type L1_HOFT_C02 -s 1126621184 -e 1137258496 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

gw_data_find --observatory L --type L1_HOFT_C02 -s 1126621184 -e 1137258496 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80  | awk '{ sum+=$4} END {print sum}'

echo "************* ER9 H1_HOFT_C00 *****************"
# --- ER9
# H1 ER9 h(t) H1_HOFT_C00 1151848448 1152167936
gw_data_find --observatory H --type H1_HOFT_C00 -s 1151848448 -e 1152167936 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80

gw_data_find --observatory H --type H1_HOFT_C00 -s 1151848448 -e 1152167936 \
    --show-times --server=ligo-gftp.pace.gatech.edu:80 | awk '{ sum+=$4} END {print sum}'
