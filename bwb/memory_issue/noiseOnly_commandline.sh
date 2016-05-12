# Job 4238ea1f5dafa3ef91deca499c48b56e
bayeswave --L1-flow 32 --PSDlength 128.0 --noiseOnly  --seglen 4.0 --Niter 2000000 --trigtime 1126259462.39 --noClean  --outputDir bayeswave_1126259462_6c846858-7897-4fdb-a314-fe34526e8c99 --H1-channel H1:DCS-CALIB_STRAIN_C01 --checkpoint  --L1-channel L1:DCS-CALIB_STRAIN_C01 --PSDstart 1126259462.39 --H1-cache datafind/1126259462.39_H1.cache --L1-cache datafind/1126259462.39_L1.cache --ifo H1 --ifo L1 --srate 1024 --H1-flow 32 --verbose  

# Job 10e7760242e35fbc6439a85a517bc976
bayeswave_post --L1-flow 32 --PSDlength 128.0 --seglen 4.0 --trigtime 1126259462.39 --outputDir bayeswave_1126259462_6c846858-7897-4fdb-a314-fe34526e8c99 --H1-channel H1:DCS-CALIB_STRAIN_C01 --0noise  --L1-channel L1:DCS-CALIB_STRAIN_C01 --PSDstart 1126259462.39 --H1-cache datafind/1126259462.39_H1.cache --L1-cache datafind/1126259462.39_L1.cache --ifo H1 --ifo L1 --srate 1024 --H1-flow 32 

