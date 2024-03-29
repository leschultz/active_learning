#!/bin/bash

# Input Prameters
MPI="mpirun"                                         # MPI
LMP="lmp_intel_cpu_intelmpi -in"                     # LAMMPS
VASP='vasp_std'                                      # VASP
WRKDIR="$(pwd)/../../../../active_learning/scripts"  # Active Learning Scripts
PREFIT=./train.cfg                                   # Initial AIMD to fit to
ACTIVE_PARITY=false                                  # Whether to produce parity plots for each loop
MAX_ITER=1000                                        # Maximum iterations when training (default 1000)
CONV_TOL=0.001                                       # Error tolerance when training (default 0.001)
EARLY_STOP=true                                      # Stop the active learning if MD is stable (not necessarily accurate)

# Start active learning
source $WRKDIR/active_learn.sh                       # Fit the potential
