#!/bin/bash

# Data sets
#
# Experiment: produce cache files on GATech & a reference LDG cluster (CIT),
# compare frame listing

# --- ER8 
# H1 ER8 h(t) H1_HOFT_C00 1123856384 1126621184
gw_data_find --observatory H --type H1_HOFT_C00 -s 1123856384 -e 1126621184 \
    --server=ligo-gftp.pace.gatech.edu:80 --url-type=file --lal-cache > ER8_H1_HOFT_C00_GATech.lcf

# L1 ER8 h(t) L1_HOFT_C00 1123856384 1126621184
gw_data_find --observatory L --type L1_HOFT_C00 -s 1123856384 -e 1126621184 \
    --server=ligo-gftp.pace.gatech.edu:80 --url-type=file --lal-cache > ER8_L1_HOFT_C00_GATech.lcf


# H1 ER8 h(t) H1_HOFT_C01 1125969920 1126621184
gw_data_find --observatory H --type H1_HOFT_C01 -s 1125969920 -e 1126621184 \
    --server=ligo-gftp.pace.gatech.edu:80 --url-type=file --lal-cache > ER8_H1_HOFT_C01_GATech.lcf

# L1 ER8 h(t) L1_HOFT_C01 1126031360 1126621184
gw_data_find --observatory L --type L1_HOFT_C01 -s 1126031360 -e 1126621184 \
    --server=ligo-gftp.pace.gatech.edu:80 --url-type=file --lal-cache > ER8_L1_HOFT_C01_GATech.lcf

# H1 ER8 h(t) H1_HOFT_C02 1125969920 1126621184
gw_data_find --observatory H --type H1_HOFT_C02 -s 1125969920 -e 1126621184 \
    --server=ligo-gftp.pace.gatech.edu:80 --url-type=file --lal-cache > ER8_H1_HOFT_C02_GATech.lcf

# L1 ER8 h(t) L1_HOFT_C02 1126031360 1126621184
gw_data_find --observatory L --type L1_HOFT_C02 -s 1126031360 -e 1126621184 \
    --server=ligo-gftp.pace.gatech.edu:80 --url-type=file --lal-cache > ER8_L1_HOFT_C02_GATech.lcf

# --- O1
# H1 01 h(t) H1_HOFT_C00 1126621184 1137258496
gw_data_find --observatory H --type H1_HOFT_C00 -s 1126621184 -e 1137258496 \
    --server=ligo-gftp.pace.gatech.edu:80 --url-type=file --lal-cache > O1_H1_HOFT_C00_GATech.lcf

# L1 01 h(t) L1_HOFT_C00 1126621184 1137254400
gw_data_find --observatory L --type L1_HOFT_C00 -s 1126621184 -e 1137254400 \
    --server=ligo-gftp.pace.gatech.edu:80 --url-type=file --lal-cache > O1_L1_HOFT_C00_GATech.lcf

# H1 01 h(t) H1_HOFT_C01 1126621184 1137258496
gw_data_find --observatory H --type H1_HOFT_C01 -s 1126621184 -e 1137258496 \
    --server=ligo-gftp.pace.gatech.edu:80 --url-type=file --lal-cache > O1_H1_HOFT_C01_GATech.lcf

# L1 01 h(t) L1_HOFT_C01 1126621184 1137258496
gw_data_find --observatory L --type L1_HOFT_C01 -s 1126621184 -e 1137258496 \
    --server=ligo-gftp.pace.gatech.edu:80 --url-type=file --lal-cache > O1_L1_HOFT_C01_GATech.lcf

# H1 01 h(t) H1_HOFT_C02 1126621184 1137258496
gw_data_find --observatory H --type H1_HOFT_C02 -s 1126621184 -e 1137258496 \
    --server=ligo-gftp.pace.gatech.edu:80 --url-type=file --lal-cache > O1_H1_HOFT_C02_GATech.lcf

# L1 01 h(t) L1_HOFT_C02 1126621184 1137258496
gw_data_find --observatory L --type L1_HOFT_C02 -s 1126621184 -e 1137258496 \
    --server=ligo-gftp.pace.gatech.edu:80 --url-type=file --lal-cache > O1_L1_HOFT_C02_GATech.lcf


# --- ER9
# H1 ER9 h(t) H1_HOFT_C00 1151848448 1152167936
gw_data_find --observatory H --type H1_HOFT_C00 -s 1151848448 -e 1152167936 \
    --server=ligo-gftp.pace.gatech.edu:80 --url-type=file --lal-cache > ER9_H1_HOFT_C00_GATech.lcf



