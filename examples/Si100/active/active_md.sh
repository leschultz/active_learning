#!/bin/bash

# Input Prameters
MPI="prun"                                     # MPI
LMP="lmp_intel_cpu_intelmpi -in"               # LAMMPS
VASP='vasp_std'                                # VASP
WRKDIR=~/packages/active_learning/scripts      # Active Learning Scripts
PREFIT=./train.cfg                             # Initial AIMD to fit to
ACTIVE_PARITY=true                             # Whether to produce parity plots for each loop
FILT_SIGMA=3                                   # Remove energies/atom outside FILT_SIGMA*sigma
FILT_SAMPLE=100                                # Only start with 100 random configurations

# Start active learning
source $WRKDIR/fit.sh    # Fit the potential
