#!/bin/bash

# Input Prameters
MPI="prun"                                     # MPI
LMP="lmp_intel_cpu_intelmpi -in"               # LAMMPS
VASP='vasp_std'                                # VASP
WRKDIR=~/packages/active_learning/scripts      # Active Learning Scripts
PREFIT=../aimd                                 # Initial AIMD to fit to
ACTIVE_PARITY=true                             # Whether to produce parity plots for each loop

# Start active learning
source $WRKDIR/fit.sh    # Fit the potential
