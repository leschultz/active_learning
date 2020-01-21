#!/bin/bash

mpirun -n $NSLOTS vasp_std > out.txt

TEMPERATURE=1150
MIN=100
DT=150
until [ $TEMPERATURE -lt $MIN ]
do
	mkdir ../$TEMPERATURE
	python3 ../../../run_creator/scripts/step_down.py "${TEMPERATURE}"
	cd ../$TEMPERATURE
	mpirun -n $NSLOTS vasp_std > out.txt
	let TEMPERATURE-=$DT
done
