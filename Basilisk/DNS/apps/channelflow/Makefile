# Basilisk user source directory (can also be set using the environment)

BASILISK_USER_SOURCE := $(if $(BASILISK_USER_SOURCE),$(BASILISK_USER_SOURCE),$(HOME)/twophase/2phase/src)

# Project configuration

SOURCE := channelFlow
GRID := multigrid3D
COMPILEMPI := 1

# Compilation options

CFLAGS := -Wall -O2 -std=c99
LIBS := -lm

# Cleanup additional files

BASILISK_CLEAN_APP := __pycache__ *.pdf *.npz

# Include the main Makefile

include $(BASILISK_USER_SOURCE)/Makefile.inc
