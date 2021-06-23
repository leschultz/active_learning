#!/bin/bash

n_selected=$(grep "BEGIN" diff.cfg | wc -l)
mlp convert-cfg diff.cfg VASP/POSCAR --output-format=vasp-poscar
mkdir -p VASP

for ((i=0; i<n_selected; i++))
do
    if [ $n_selected -eq 1 ]; then 
        mv VASP/POSCAR VASP/0/POSCAR
    elif [ $n_selected -gt 1 ]; then
        mv VASP/POSCAR"$i" VASP/"$i"/POSCAR
    fi
    cp VASP/POSCAR VASP/"$i"/
    cp VASP/INCAR VASP/"$i"/
done

rm diff.cfg 
rm output_ef*
rm lammps_input*
rm log.lammps
rm input.pos
