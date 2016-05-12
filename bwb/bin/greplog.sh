#!/bin/bash

grep updated $1 > tmp
awk '{print $4, $10}' tmp > tmp2
mv tmp2 $1
rm tmp
