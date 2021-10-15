#!/bin/bash

# Input Prameters
MPI="prun"                                     # MPI
LMP="lmp_intel_cpu_intelmpi -in"               # LAMMPS
VASP='vasp_std'                                # VASP
WRKDIR=~/packages/active_learning/scripts      # Active Learning Scripts
PREFIT=../train.cfg                            # Initial fit. Every OUTCAR in directory or mlip2 training file.
ACTIVE_PARITY=true                             # Whether to produce parity plots for each loop
FILT_SIGMA=3                                   # Remove energies/atom outside FILT_SIGMA*sigma

# Start active learning
source $WRKDIR/fit.sh    # Fit the potential
