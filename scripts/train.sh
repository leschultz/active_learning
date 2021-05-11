#!/bin/bash

# Input parameters
POTS=/home/nerve/packages/mlip-2/untrained_mtps/
POT=18.mtp
JOB='../aimd'
CORES=$(nproc)

cd $JOB

# The number of elements
COMP=$(cat run.sh | grep COMP= | awk -F '=' '{print $2}')
NELS=$(echo $COMP | grep -Eo '[[:alpha:]]+' | wc -w)
echo $NELS

# Make file to store values
rm -f -- train.cfg
touch train.cfg

# Create training file
for i in $(find ./ -type f -name OUTCAR); do
	mlp convert-cfg --input-format=vasp-outcar $i out.cfg
	cat out.cfg >> train.cfg
	rm out.cfg
done

# Get the minimum distances
MINDIST=$(mlp mindist train.cfg | awk -F ': ' '{print $2}')
echo $MINDIST

# Specify potential parameters
cp "${POTS}${POT}" ./init.mtp
sed -i "s/species_count.*/species_count = $NELS/g" init.mtp
sed -i "s/min_dist.*/min_dist = $MINDIST/g" init.mtp

# Train the potential
mpirun -n $CORES mlp train init.mtp train.cfg --trained-pot-name=potential.mtp
