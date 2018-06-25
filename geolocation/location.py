#!/usr/bin/env python

import re
import os
import socket

machine_ad_file_name = os.getenv('_CONDOR_MACHINE_AD')
try:
    machine_ad_file = open(machine_ad_file_name, 'r')
    machine_ad = machine_ad_file.read()
    machine_ad_file.close()
except TypeError:
    print socket.getfqdn()
    exit(1)

try:
    print re.search(r'GLIDEIN_Gatekeeper = "(.*):\d*/jobmanager-\w*"', machine_ad, re.MULTILINE).group(1)
except AttributeError:
    try:
        print re.search(r'GLIDEIN_Gatekeeper = "(\S+) \S+:9619"', machine_ad, re.MULTILINE).group(1)
    except AttributeError:
        exit(1)
