#!/bin/bash
 
###############################################################
#                                                             #
# Example BayesWave run sampling from prior                   #   
#                                                             # 
# RUN TIME ~20 minutes                                        #  
#  - Processor: 2 GHz Intel Core i7                           #
#  - Memory: 4 GB 1600 Mhz                                    #
#                                                             #
###############################################################

bwb_pipe.py sample_prior.ini \
    --workdir sample_prior \
    --trigger-time 900000000
