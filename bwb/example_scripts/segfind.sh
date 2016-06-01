#!/bin/bash

#   ligolw_segment_query_dqsegdb -q -t https://segments.ligo.org -s <gpsstart> -e
#   <gpsend> -a <x1:james-is-drunk:1>
#
#   -s - start time
#   -e - end time
#   -a - name of flag (with version number, ideally)
#
#   In python you can use the dqsegdb client directly, but you have to parse the
#   json responses manually, which gwpy has already got set up for you:
#
#   >>> from gwpy.segments import DataQualityFlag
#   >>> flag = DataQualityFlag.query('X1:FLAG:1', start, end)
#
#   see https://gwpy.github.io/docs/latest/segments/ for more details.

segdb="https://segments.ligo.org"
gpsstart=1126249362
gpsend=1126269562
flag="H1:DMT-SCIENCE:1"

ligolw_segment_query_dqsegdb -q -t ${segdb} -s ${gpsstart} -e ${gpsend} -a ${flag}






