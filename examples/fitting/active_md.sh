#!/bin/bash

# Input Prameters
CORES=$(nproc)
MPI="prun"                                 # MPI
LMP="lmp_intel_cpu_intelmpi -in"           # LAMMPS
VASP='vasp_std'                            # VASP
WRKDIR=~/packages/active_learning/scripts  # Active Learning Scripts

# Start active learning
source $WRKDIR/run.sh
