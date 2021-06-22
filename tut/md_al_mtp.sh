#!/bin/bash

# Working variables
MPI="mpirun -n $(nproc)"
LMP="lmp_mpi -in"

# Clean working space
rm -f train.cfg
rm -f curr.mtp
rm -f preselected.cfg
rm -f diff.cfg
rm -f selected.cfg
rm -f out.cfg

# Initialize potential and training data
cp init.mtp curr.mtp
cp train_init.cfg train.cfg

# Initialize active learning state
mlp calc-grade curr.mtp train.cfg train.cfg out.cfg --als-filename=state.als

# Start infinite loop
while [ 1 -gt 0 ]
do

# Run MD
touch preselected.cfg
$LMP in.nb_md

# The number of MD steps with extrapolation grade above a threshold in mlip.ini
n_preselected=$(grep "BEGIN" preselected.cfg | wc -l)

# If n_preseleted is greater than zero.
if [ $n_preselected -gt 0 ]; then

    # Add configurations to the training set from preselected
    mlp select-add curr.mtp train.cfg preselected.cfg diff.cfg --als-filename=state.als
    cp diff.cfg calculate_ab_initio_ef/

    # Clean files
    rm -f preselected.cfg
    rm -f selected.cfg

    # Calculate energies and forces and convert LAMMPS to MLIP format.
    cd calculate_ab_initio_ef/
    ./ab_initio_calculations.sh
    cd ../

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
