#!/bin/bash

# Scripts
source ./run_types/gen_aimd.sh
source ./run_types/gen_dft.sh
source ./run_types/gen_md.sh

# Input Prameters
COMP=Cu5Zr5
TEBEG=2000
TEND=600
DT=100
SEED=$RANDOM
CORES=$(nproc)

# Programs
MPI="prun"
LMP="lmp_intel_cpu_intelmpi -in"
VASP='vasp_std'
WRKDIR=$(pwd)

# VASP inputs
RECDIR=~/potentials/vasp/
RECPOTS=$WRKDIR/vasp_pots.csv
TYPE=pbe

# MLIP-2 inputs
POTS=~/potentials/mlip-2/
POT=08.mtp

# Clean working space
rm -rf ../runs

# Directory where run data lives
cd ..
mkdir runs
cd runs

mkdir $TEBEG
cd $TEBEG

mkdir aimd
cd aimd

# Generate VASP run files and get element masses
MASSES=$(aimd_job $WRKDIR $COMP $TEBEG $TEBEG $RECDIR $RECPOTS $TYPE)

# Run AIMD
$MPI $VASP

cd ../../
mkdir potential
cd potential

# Create training file
mlp convert-cfg --input-format=vasp-outcar ../$TEBEG/aimd/OUTCAR train.cfg

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

# Copy starting postions to LAMMPS run
$WRKDIR/convert/poscar2lammps.awk ../$TEBEG/aimd/POSCAR > input.pos

# Start infinite loop
while [ 1 -gt 0 ]
do

# Run MD
touch preselected.cfg

# Define the masses for classical MD
md $WRKDIR "$MASSES"

$LMP md.in  # Has to run in serial because of lmp potential

# The number of MD steps with extrapolation grade above a threshold in mlip.ini
n_preselected=$(grep "BEGIN" preselected.cfg | wc -l)

# If n_preseleted is greater than zero.
if [ $n_preselected -gt 0 ]; then

    # Add configurations to the training set from preselected
    mlp select-add curr.mtp train.cfg preselected.cfg diff.cfg --als-filename=state.als
    mkdir calculate_ab_initio_ef
    cp diff.cfg calculate_ab_initio_ef

    # Clean files
    rm -f preselected.cfg
    rm -f selected.cfg

    # Calculate energies and forces and convert LAMMPS to MLIP format.
    cd calculate_ab_initio_ef
    dft $WRKDIR $COMP $RECDIR $RECPOTS $TYPE $MPI $VASP
    cd -

    # Re-train the current potential
    $MPI mlp train curr.mtp train.cfg --trained-pot-name=curr.mtp --update-mindist
    
    # Update the active learning state
    mlp calc-grade curr.mtp train.cfg diff.cfg out.cfg --als-filename=state.als
    
    # Clean files
    rm -f diff.cfg
    rm -f out.cfg
    
elif  [ $n_preselected -eq 0 ]; then
    exit
fi

done
