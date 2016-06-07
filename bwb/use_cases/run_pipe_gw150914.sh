#!/bin/bash
#
# Example call to the pipeline, using trigger time from offline CWB CED:
# https://www.atlas.aei.uni-hannover.de/~vedovato/LSC/reports/G184098_ceds_run1/ced/ced_1126259272_800_G184098_ceds_run1_slag0_lag0_1_rMRA/L1H1_1126259461.500_1126259461.500/

bwb_pipe.py usecase_GW150914-C02.ini \
    --workdir usecase_GW150914-C02 \
    --trigger-time 1126259462.392
