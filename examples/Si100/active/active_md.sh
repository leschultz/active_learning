#!/bin/bash

# Input Prameters
MPI="prun"                                     # MPI
LMP="lmp_intel_cpu_intelmpi -in"               # LAMMPS
VASP='vasp_std'                                # VASP
WRKDIR=~/packages/active_learning/scripts      # Active Learning Scripts
PREFIT=../train.cfg                            # Initial AIMD to fit to
ACTIVE_PARITY=true                             # Whether to produce parity plots for each loop
MAX_ITER=1000                                  # Maximum iterations when training (default 1000)
CONV_TOL=0.001                                 # Error tolerance when training (default 0.001)

# Start active learning
source $WRKDIR/active_learn.sh                 # Fit the potential
