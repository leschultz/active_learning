#!/bin/bash

n_selected=$(grep "BEGIN" diff.cfg | wc -l)
mlp convert-cfg diff.cfg lammps_input --output-format=lammps-datafile

for ((i=0; i<n_selected; i++))
do
    if [ $n_selected -eq 1 ]; then 
        cp lammps_input input.pos
    elif [ $n_selected -gt 1 ]; then
        cp lammps_input"$i" input.pos
    fi
    lmp_mpi -in calc_ef.in
    python3 convert_lammps_dump_to_cfg.py
    cat output_ef.cfg >> ../train.cfg
done

rm diff.cfg 
rm output_ef*
rm lammps_input*
rm log.lammps
rm input.pos
