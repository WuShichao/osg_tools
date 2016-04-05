BayesWavePost \
        --ifo H1 --H1-flow 16 \
        --ifo L1 --L1-flow 16 \
        --H1-cache H1.cache --H1-channel H1:LDAS-STRAIN \
        --L1-cache L1.cache --L1-channel L1:LDAS-STRAIN \
        --srate 256 --seglen 4 \
        --trigtime 1126259462 --PSDstart 1126259462 --PSDlength 8 \
        --0noise
