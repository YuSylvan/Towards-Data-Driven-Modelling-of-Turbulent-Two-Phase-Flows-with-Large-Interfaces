#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

. $WM_PROJECT_DIR/bin/tools/CleanFunctions

cleanCase
rm -fr 0 __pycache__ *.pdf
