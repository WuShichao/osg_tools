#!/bin/bash
 

# Supply a trigger time (can leave this off and use default 1126259462.392)
bwb_pipe.py internal_injection.ini \
    --workdir internal_injection \
    --trigger-time 900000000
