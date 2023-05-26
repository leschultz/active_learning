#!/bin/bash

# Turn on for special remove
shopt -s extglob

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

# Make run directory
mkdir run
cd run

mkdir potential
cd potential
POTDIR=$(pwd)
cp "${TOPDIR}/${PREFIT}" .

# Prepare potential
cp $TOPDIR/curr.mtp ./curr.mtp
sed -i "s/species_count.*/species_count = $NELS/g" curr.mtp

"$MPI" mlp train curr.mtp train.cfg --trained-pot-name=curr.mtp --update-mindist --max-iter=${MAX_ITER} --bfgs-conv-tol=${CONV_TOL}

# Initialize active learning state
mlp calc-grade curr.mtp train.cfg train.cfg out.cfg --als-filename=state.als
rm out.cfg

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
	cp $TOPDIR/mlip.ini . # Needed for potential parameters

	$LMP md.in  # Has to run in serial because of active learning

	# The number of MD steps with extrapolation grade above a threshold in mlip.ini
	n_preselected=$(grep "BEGIN" preselected.cfg | wc -l)

	# If n_preseleted is greater than zero.
	if [ $n_preselected -gt 0 ]; then

	    # Add configurations to the training set from preselected
	    mlp select-add curr.mtp train.cfg preselected.cfg diff.cfg --als-filename=state.als
	    rm preselected.cfg selected.cfg
	    mkdir ../dft
	    mv diff.cfg ../dft
	    mv curr.mtp ../dft
	    mv state.als ../dft
	    mv train.cfg ../dft
	    mv mlip.ini ../dft

	    # Calculate energies and forces and convert LAMMPS to MLIP format.
	    cd ../dft
	    dft_job $TOPDIR "$MPI" $VASP ../train.cfg

	    # Prepare retraining
	    mv diff.cfg $POTDIR
	    mv curr.mtp $POTDIR
	    mv state.als $POTDIR
	    mv train.cfg $POTDIR
	    mv mlip.ini $POTDIR
	    cd $POTDIR

	    "$MPI" mlp train curr.mtp train.cfg --trained-pot-name=curr.mtp --update-mindist --max-iter=${MAX_ITER} --bfgs-conv-tol=${CONV_TOL}
	    
	    # Update the active learning state
	    mlp calc-grade curr.mtp train.cfg diff.cfg out.cfg --als-filename=state.als
	    rm out.cfg diff.cfg
	    
	    # Whether to produce on the fly parity plots
	    if [ $ACTIVE_PARITY = true ]; then
		    test  # Genreate test set and parity plots
	    fi

	    # If MD is stable then stop
	    if [ $EARLY_STOP = true ]; then

	    	if grep -q "Total wall time" ../md/log.lammps; then
			test
			rm ../md/log.lammps
			exit
		fi
			
	    fi

	    cd ../md_dft/

	    # Increment counter
	    ITERS=$((ITERS+1))

	# If n_preselected is equal to zero and finish
	elif [ $n_preselected -eq 0 ]; then

	    rm preselected.cfg log.lammps selected.cfg
	    test  # Genreate test set and parity plots
	    exit
	fi

done
