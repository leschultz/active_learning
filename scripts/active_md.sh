#!/bin/bash

# Input Prameters
COMP=Si64
CORES=$(nproc)
JOBDIR=~/runs/active_test

# Programs
MPI="prun"                           # MPI
LMP="lmp_intel_cpu_intelmpi -in"     # LAMMPS
VASP='vasp_std'                      # VASP
WRKDIR=~/packages/vasp_runs/scripts  # Active Learning Scripts

# VASP inputs
RECDIR=~/potentials/vasp/
RECPOTS=$WRKDIR/vasp_pots.csv
TYPE=pbe

# MLIP-2 inputs
POTS=~/potentials/mlip-2/
POT=08.mtp

# Start active learning
source $WRKDIR/run.sh

# Make AIMD starting run
mkdir aimd
cd aimd
AIMDDIR=$(pwd)

# Generate VASP run files and get element masses
MASSES=$(aimd_job $WRKDIR $COMP $RECDIR $RECPOTS $TYPE)

# Run AIMD
$MPI $VASP

cd ../
mkdir potential
cd potential
POTDIR=$(pwd)

# Create training file
mlp convert-cfg --input-format=vasp-outcar ../aimd/OUTCAR train.cfg

# Prepare potential
NELS=$(echo $COMP | grep -Eo '[[:alpha:]]+' | wc -w)
MINDIST=$(mlp mindist train.cfg | awk -F ': ' '{print $2}')
cp "${POTS}${POT}" ./curr.mtp
sed -i "s/species_count.*/species_count = $NELS/g" curr.mtp
sed -i "s/min_dist.*/min_dist = $MINDIST/g" curr.mtp

# Create initial potential
$MPI mlp train curr.mtp train.cfg --trained-pot-name=curr.mtp

# Initialize active learning state
mlp calc-grade curr.mtp train.cfg train.cfg out.cfg --als-filename=state.als

cd ..

mkdir md_dft
cd md_dft

# Start infinite loop
ITERS=0
while [ 1 -gt 0 ]
do

mkdir $ITERS
cd $ITERS

mkdir md
cd md

# Copy starting postions to LAMMPS run
$WRKDIR/convert/poscar2lammps.awk $AIMDDIR/POSCAR > input.pos
touch preselected.cfg
md_job $WRKDIR "$MASSES"  # Preapare MD job

# Add potential files
mv $POTDIR/curr.mtp .  # Needed for MD
mv $POTDIR/state.als . # Needed for MD
mv $POTDIR/train.cfg . # Needed for adding frames to training

$LMP md.in  # Has to run in serial because of active learning

# The number of MD steps with extrapolation grade above a threshold in mlip.ini
n_preselected=$(grep "BEGIN" preselected.cfg | wc -l)

# If n_preseleted is greater than zero.
if [ $n_preselected -gt 0 ]; then

    # Add configurations to the training set from preselected
    mlp select-add curr.mtp train.cfg preselected.cfg diff.cfg --als-filename=state.als
    mkdir ../dft
    mv diff.cfg ../dft
    mv curr.mtp ../dft
    mv state.als ../dft
    mv train.cfg ../dft

    # Calculate energies and forces and convert LAMMPS to MLIP format.
    cd ../dft
    dft_job $WRKDIR $COMP $RECDIR $RECPOTS $TYPE $MPI $VASP ../train.cfg
    mv diff.cfg ../
    mv curr.mtp ../
    mv state.als ../
    mv train.cfg ../
    cd ../

    # Re-train the current potential
    $MPI mlp train curr.mtp train.cfg --trained-pot-name=curr.mtp --update-mindist
    
    # Update the active learning state
    mlp calc-grade curr.mtp train.cfg diff.cfg out.cfg --als-filename=state.als
    
    # Move back to potential folder
    mv curr.mtp $POTDIR
    mv state.als $POTDIR
    mv train.cfg $POTDIR
    cd ../

    # Increment counter
    ITERS=$((ITERS+1))

elif  [ $n_preselected -eq 0 ]; then
    exit
fi

done
