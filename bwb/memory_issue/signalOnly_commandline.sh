# Job aa505433cf11dee5c52e3d083d43d554
bayeswave --outputDir bayeswave_1126259462_1b7c7e7d-b391-4e86-bac8-c3603eb5283a --L1-flow 32 --PSDlength 128.0 --verbose  --seglen 4.0 --Niter 2000000 --signalOnly  --noClean  --trigtime 1126259462.39 --H1-channel H1:DCS-CALIB_STRAIN_C01 --checkpoint  --L1-channel L1:DCS-CALIB_STRAIN_C01 --PSDstart 1126259462.39 --H1-cache datafind/1126259462.39_H1.cache --L1-cache datafind/1126259462.39_L1.cache --ifo H1 --ifo L1 --srate 1024 --H1-flow 32 

# Job 97471045897e3cfbc7600de637606596
bayeswave_post --L1-flow 32 --PSDlength 128.0 --seglen 4.0 --trigtime 1126259462.39 --outputDir bayeswave_1126259462_1b7c7e7d-b391-4e86-bac8-c3603eb5283a --H1-channel H1:DCS-CALIB_STRAIN_C01 --0noise  --L1-channel L1:DCS-CALIB_STRAIN_C01 --PSDstart 1126259462.39 --H1-cache datafind/1126259462.39_H1.cache --L1-cache datafind/1126259462.39_L1.cache --ifo H1 --ifo L1 --srate 1024 --H1-flow 32 

