# Job ddc983ddf04638a7cbfe5ade0d5891bf
bayeswave --outputDir bayeswave_1126259462_1e0beeb9-deb7-41d9-8a51-34594ed81d57 --L1-flow 32 --PSDlength 128.0 --verbose  --seglen 4.0 --Niter 2000000 --glitchOnly  --noClean  --trigtime 1126259462.39 --H1-channel H1:DCS-CALIB_STRAIN_C01 --checkpoint  --L1-channel L1:DCS-CALIB_STRAIN_C01 --PSDstart 1126259462.39 --H1-cache datafind/1126259462.39_H1.cache --L1-cache datafind/1126259462.39_L1.cache --ifo H1 --ifo L1 --srate 1024 --H1-flow 32 

# Job e17ae1352eb67e0cb56f5eb43790aaf7
bayeswave_post --L1-flow 32 --PSDlength 128.0 --seglen 4.0 --trigtime 1126259462.39 --outputDir bayeswave_1126259462_1e0beeb9-deb7-41d9-8a51-34594ed81d57 --H1-channel H1:DCS-CALIB_STRAIN_C01 --0noise  --L1-channel L1:DCS-CALIB_STRAIN_C01 --PSDstart 1126259462.39 --H1-cache datafind/1126259462.39_H1.cache --L1-cache datafind/1126259462.39_L1.cache --ifo H1 --ifo L1 --srate 1024 --H1-flow 32 

