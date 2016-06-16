#!/bin/bash
#
# Show most recent Done line in dagman.out
workdir=$1
grep -A2 Done $workdir/*man.out | tail -n 3
