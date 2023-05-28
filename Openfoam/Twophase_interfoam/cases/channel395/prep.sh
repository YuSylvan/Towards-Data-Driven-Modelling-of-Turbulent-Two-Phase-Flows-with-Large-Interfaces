#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

. $WM_PROJECT_DIR/bin/tools/RunFunctions

cp -r 0.orig 0

runApplication blockMesh
