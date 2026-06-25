#!/bin/bash

# Configuration for Figaro environment

# Binary directory for local builds
export EXECDIR="bin"

# Solver paths
export MOSEDIR="$HYDRADIR/src/lib/MOSE"
export MOSEMASTER="$MOSEDIR/bin/MOSE"
export MOSELOCAL="./$EXECDIR/MOSE"

export FUSSDIR="$HYDRADIR/src/lib/FUSS"
export FUSSMASTER="$FUSSDIR/bin/FUSS"
export FUSSLOCAL="./$EXECDIR/FUSS"

export GROOTDIR="$HYDRADIR/src/lib/GROOT"
export GROOTMASTER="$GROOTDIR/bin/GROOT"
export GROOTLOCAL="./$EXECDIR/GROOT"

export ARESDIR="$HYDRADIR/src/lib/ARES"
export ARESMASTER="$ARESDIR/bin/ARES"
export ARESLOCAL="./$EXECDIR/ARES"