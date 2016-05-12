# Job da38b60a8658e41f219f03f457fd77b4
bayeswave --L1-flow 32 --PSDlength 128.0 --verbose  --seglen 4.0 --Niter 2000000 --trigtime 1126259462.39 --noClean  --outputDir bayeswave_1126259462_a717f0c5-9a34-44b0-89a3-a602806fd704 --H1-channel H1:DCS-CALIB_STRAIN_C01 --checkpoint  --L1-channel L1:DCS-CALIB_STRAIN_C01 --PSDstart 1126259462.39 --H1-cache datafind/1126259462.39_H1.cache --L1-cache datafind/1126259462.39_L1.cache --ifo H1 --ifo L1 --srate 1024 --H1-flow 32 

# Job fc4b4e11a22475256c7e0c9acdb41230
bayeswave_post --L1-flow 32 --PSDlength 128.0 --seglen 4.0 --trigtime 1126259462.39 --outputDir bayeswave_1126259462_a717f0c5-9a34-44b0-89a3-a602806fd704 --H1-channel H1:DCS-CALIB_STRAIN_C01 --0noise  --L1-channel L1:DCS-CALIB_STRAIN_C01 --PSDstart 1126259462.39 --H1-cache datafind/1126259462.39_H1.cache --L1-cache datafind/1126259462.39_L1.cache --ifo H1 --ifo L1 --srate 1024 --H1-flow 32 

