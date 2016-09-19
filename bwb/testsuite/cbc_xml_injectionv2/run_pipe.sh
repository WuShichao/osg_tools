#!/bin/bash
 
source ~/etc/bayeswave_ldg-user-env.sh

rm -r biginj_followup
# Supply a trigger time (can leave this off and use default 1126259462.392)
bwb_pipe.py burst_mdc_injection_BBH_style.ini \
    -I HL-INJECTIONS_78901_BBH1_cWB-1126051217-11206883.xml \
    --workdir ./biginj_followup \
    --followup-injections trigger_list_BBH_set_1.txt \
    --sim-data
