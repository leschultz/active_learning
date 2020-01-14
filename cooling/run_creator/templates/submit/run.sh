#!/bin/bash

mpirun -n $NSLOTS vasp_std > out.txt

TEMPERATURE=1000
MIN=500
DT=100
until [ $TEMPERATURE -lt $MIN ]
do
	mkdir ../$TEMPERATURE
	python3 ../../../run_creator/scripts/step_down.py "${TEMPERATURE}"
	cd ../$TEMPERATURE
	mpirun -n $NSLOTS vasp_std > out.txt
	let TEMPERATURE-=$DT
done
