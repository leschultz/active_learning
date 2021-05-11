#!/bin/bash

POTS=/home/nerve/packages/mlip-2/untrained_mtps/
POT=18.mtp
JOB='./aimd'

cd $JOB

# Make file to store values
rm -f -- train.cfg
touch train.cfg

# The potential choice
cp "${POTS}${POT}" ./init.mtp

for i in $(find ./ -type f -name OUTCAR); do
	mlp convert-cfg --input-format=vasp-outcar $i out.cfg
	cat out.cfg >> train.cfg
	rm out.cfg
done

# Train the potential
mlp train init.mtp train_init.cfg --trained-pot-name=init.mtp
