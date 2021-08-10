#!/bin/bash

# Clean
rm -rf run

# Scripts
source $WRKDIR/run_types/gen_dft.sh
source $WRKDIR/run_types/gen_md.sh
source $WRKDIR/run_types/gen_parity.sh

TOPDIR=$(pwd)

# Get elements and masses from POTCAR
MASSES=$(cat POTCAR | grep MASS | awk -F ' ' '{print $3}' | sed 's/\;//g')
ELEMS=$(cat POTCAR | grep 'VRHFIN =' | grep -oP '(?<==).*?(?=:)')
NELS=$(echo $ELEMS | grep -Eo '[[:alpha:]]+' | wc -w)

# Search for any OUTCARS that may have ben previously fit
PREFIT=$(find ./$PREFIT -type f -name 'OUTCAR')

# Make run directory
mkdir run
cd run

if (( ${#PREFIT[@]} != 0 )); then

	mkdir potential
	cd potential
        POTDIR=$(pwd)

	touch train.cfg

	for i in $PREFIT; do
		mlp convert-cfg --input-format=vasp-outcar $TOPDIR/$i tr.cfg
		cat tr.cfg >> train.cfg
	done
	rm tr.cfg

else
	echo "Need to provide initial training OUTCARS"
	exit 1
fi

# Prepare potential
MINDIST=$(mlp mindist train.cfg | awk -F ': ' '{print $2}')
cp $TOPDIR/curr.mtp ./curr.mtp
sed -i "s/species_count.*/species_count = $NELS/g" curr.mtp
sed -i "s/min_dist.*/min_dist = $MINDIST/g" curr.mtp

# Create initial potential
$MPI mlp train curr.mtp train.cfg --trained-pot-name=curr.mtp

# Initialize active learning state
mlp calc-grade curr.mtp train.cfg train.cfg out.cfg --als-filename=state.als

cd ..

# Make directory for dft and md
mkdir md_dft
cd md_dft

# Start infinite loop
ITERS=0
while true
do

	# Make a directory for each iteration
	mkdir $ITERS
	cd $ITERS

	mkdir md
	cd md

	# Copy starting postions to LAMMPS run
	$WRKDIR/convert/poscar2lammps.awk $TOPDIR/POSCAR > input.pos
	touch preselected.cfg
	md_job $TOPDIR "$MASSES" "$ELEMS"  # Preapare MD job

	# Add potential files
	mv $POTDIR/curr.mtp .  # Needed for MD
	mv $POTDIR/state.als . # Needed for MD
	mv $POTDIR/train.cfg . # Needed for adding frames to training
	cp $TOPDIR/mlip.ini . # Neede for potential parameters
	mv $POTDIR/out.cfg . # Not needed but copied for completness

	$LMP md.in  # Has to run in serial because of active learning

	# The number of MD steps with extrapolation grade above a threshold in mlip.ini
	n_preselected=$(grep "BEGIN" preselected.cfg | wc -l)

	# If n_preseleted is greater than zero.
	if [ $n_preselected -gt 0 ]; then

	    # Add configurations to the training set from preselected
	    mlp select-add curr.mtp train.cfg preselected.cfg diff.cfg --als-filename=state.als
	    mkdir ../dft
	    cp diff.cfg ../dft
	    cp curr.mtp ../dft
	    cp state.als ../dft
	    cp train.cfg ../dft
	    cp mlip.ini ../dft
	    cp out.cfg ../dft

	    # Calculate energies and forces and convert LAMMPS to MLIP format.
	    cd ../dft
	    dft_job $TOPDIR $MPI $VASP ../train.cfg

	    # Prepare retraining
	    mkdir ../retrain
	    cp diff.cfg ../retrain
	    cp curr.mtp ../retrain
	    cp state.als ../retrain
	    cp train.cfg ../retrain
	    cp mlip.ini ../retrain
	    cd ../retrain

	    # Re-train the current potential
	    $MPI mlp train curr.mtp train.cfg --trained-pot-name=curr.mtp --update-mindist
	    
	    # Update the active learning state
	    mlp calc-grade curr.mtp train.cfg diff.cfg out.cfg --als-filename=state.als
	    
	    # Move back to potential folder
	    cp curr.mtp $POTDIR
	    cp state.als $POTDIR
	    cp train.cfg $POTDIR
	    cp out.cfg $POTDIR
	    cp mlip.ini $POTDIR

	    # Whether to produce on the fly parity plots
	    if [ $ACTIVE_PARITY = true ]; then
		    test  # Genreate test set and parity plots
	    fi

	    cd ../../

	    # Increment counter
	    ITERS=$((ITERS+1))

	elif  [ $n_preselected -eq 0 ]; then

	    # Move potential back to original folder
	    cp curr.mtp $POTDIR
	    cp state.als $POTDIR
	    cp train.cfg $POTDIR
	    cp mlip.ini $POTDIR
	    cp out.cfg $POTDIR

	    test  # Genreate test set and parity plots

	    exit
	fi

done
