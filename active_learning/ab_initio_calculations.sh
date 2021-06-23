#!/bin/bash

MPI=$1
VASP=$2

n_selected=$(grep "BEGIN" diff.cfg | wc -l)
mlp convert-cfg diff.cfg POSCAR --output-format=vasp-poscar

for ((i=0; i<n_selected; i++))
do
    if [ $n_selected -eq 1 ]; then 
	mkdir 0
        mv POSCAR 0/POSCAR
    elif [ $n_selected -gt 1 ]; then
	mkdir $i
        mv POSCAR$i $i/POSCAR
    fi
    cp ../../aimd/POTCAR $i/
    cp ../../aimd/KPOINTS $i/
    cp ../../../../dft_INCAR $i/INCAR

    # Do DFT to get energies and forces
    cd $i
    $MPI $VASP

    mlp convert-cfg OUTCAR calculated.cfg --input-format=vasp-outcar
    cat calculated.cfg >> ../../train.cfg
    cd -

    rm -rf $i

done

