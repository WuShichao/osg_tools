#!/bin/bash

find /cvmfs/oasis.opensciencegrid.org/ligo/frames/ER8/ -name *H1*C01*.gwf  | lalapps_path2cache > ER8-H1_HOFT_C01-osg.lcf
find /cvmfs/oasis.opensciencegrid.org/ligo/frames/ER8/ -name *L1*C01*.gwf  | lalapps_path2cache > ER8-L1_HOFT_C01-osg.lcf

find /cvmfs/oasis.opensciencegrid.org/ligo/frames/O1/ -name *H1*C01*.gwf  | lalapps_path2cache > O1-H1_HOFT_C01-osg.lcf
find /cvmfs/oasis.opensciencegrid.org/ligo/frames/O1/ -name *L1*C01*.gwf  | lalapps_path2cache > O1-L1_HOFT_C01-osg.lcf
