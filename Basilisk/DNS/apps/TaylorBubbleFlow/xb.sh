#!/bin/bash

# Extract time and xb from log file and plot

LOGFILE=log

cat $LOGFILE | grep xb | tr -d , | awk -F " " '{ print $6 }' > t.txt
cat $LOGFILE | grep xb | tr -d , | awk -F " " '{ print $12 }' > xb.txt

python3 xb.py
