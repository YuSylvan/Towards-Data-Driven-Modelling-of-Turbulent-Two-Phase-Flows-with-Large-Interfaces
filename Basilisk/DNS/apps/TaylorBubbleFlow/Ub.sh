#!/bin/bash

# Extract time and Ub from log file and plot

LOGFILE=log

cat $LOGFILE | grep Ub | tr -d , | awk -F " " '{ print $6 }' > t.txt
cat $LOGFILE | grep Ub | tr -d , | awk -F " " '{ print $12 }' > Ub.txt

python3 Ub.py
