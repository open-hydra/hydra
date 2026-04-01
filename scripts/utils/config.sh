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

export IGLOODIR="$HYDRADIR/src/lib/IGLOO"
export IGLOOMASTER="$IGLOODIR/bin/IGLOO"
export IGLOOLOCAL="./$EXECDIR/IGLOO"

export ICEDIR="$HYDRADIR/src/lib/ICE"
export ICEMASTER="$ICEDIR/bin/ICE"
export ICELOCAL="./$EXECDIR/ICE"

export MOSKADIR="$HYDRADIR/src/lib/MOSKA"
export MOSKAMASTER="$MOSKADIR/bin/MOSKA"
export MOSKALOCAL="./$EXECDIR/MOSKA"

export Q2DDIR="$HYDRADIR/src/lib/Q2D"
export Q2DMASTER="$Q2DDIR/bin/Q2D"
export Q2DLOCAL="./$EXECDIR/Q2D"
