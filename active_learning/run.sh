#!/bin/bash

# Input Prameters
COMP=Si64
TEMPER=2000
SEED=$RANDOM
CORES=$(nproc)

# Programs
MPI="prun"
LMP="lmp_intel_cpu_intelmpi -in"
VASP='vasp_std'

# VASP inputs
RECDIR=~/potentials/vasp/
RECPOTS="${RECDIR}vasp_pots.csv"
TYPE=pbe

# MLIP-2 inputs
POTS=~/potentials/mlip-2/
POT=08.mtp

# Clean working space
rm -rf runs

# Directory where run data lives
mkdir runs
cd runs

# Initial VASP run
mkdir $TEMPER
cd $TEMPER

mkdir aimd
cd aimd

# Create initial features
python3 ../../../gen_poscar.py $COMP  # Creates random initial positions
bash ../../../gen_potcar.sh $COMP $TYPE $RECDIR $RECPOTS  # Creates potential

# Get the masses for the elements
MASSES=$(cat POTCAR | grep MASS | awk -F ' ' '{print $3}' | sed 's/\;//g')

touch KPOINTS
cat > KPOINTS <<!
simple
0 0 0
Gamma
 1 1 1
 0 0 0
!

touch INCAR
cat > INCAR <<!
SYSTEM =  $COMP_$TEMPER

# electronic degrees
LREAL = A                      # real space projection
PREC  = Normal                 # chose Low only after tests
EDIFF = 1E-5                   # do not use default (too large drift)
ISMEAR = -1                    # determine how partial occupancies are set for each orbital
SIGMA = 0.130                  # specifies the width of the smearing in eV
ALGO = Very Fast               # recommended for MD (fall back ALGO = Fast)
MAXMIX = 40                    # reuse mixer from one MD step to next
ISYM = 0                       # no symmetry
NELMIN = 4                     # minimum of steps per time step
NELM = 60                      # maximum of steps per time step

# MD
IBRION = 0                     # main molecular dynamics tag
NSW = 300                      # number of MD steps
POTIM = 3                      # time step of MD [fs]
LCHARG = .FALSE.               # no charge density written out
LWAVE = .FALSE.                # no wave function coefficients written out
TEBEG = $TEMPER                # starting temperature for MD
TEEND = $TEMPER                # end temperature for MD

# canonic (Nose) MD
MDALGO = 2                     # switch to select thermostat
SMASS =  3                     # Nose mass
ISIF = 2                       # this tag selects the ensemble in combination with the thermostat
!

# Run AIMD
$MPI $VASP

MASSES=$(cat POTCAR | grep MASS | awk -F ' ' '{print $3}' | sed 's/\;//g') # Delete
cd ..
mkdir potential
cd potential

# Create training file
mlp convert-cfg --input-format=vasp-outcar ../aimd/OUTCAR train.cfg

# Prepare potential
NELS=$(echo $COMP | grep -Eo '[[:alpha:]]+' | wc -w)
MINDIST=$(mlp mindist train.cfg | awk -F ': ' '{print $2}')
cp "${POTS}${POT}" ./curr.mtp
sed -i "s/species_count.*/species_count = $NELS/g" curr.mtp
sed -i "s/min_dist.*/min_dist = $MINDIST/g" curr.mtp

# Create Initial potential
$MPI mlp train curr.mtp train.cfg --trained-pot-name=curr.mtp

# Initialize active learning state
mlp calc-grade curr.mtp train.cfg train.cfg out.cfg --als-filename=state.als

# Copy starting postions to LAMMPS run
../../../poscar2lammps.awk ../aimd/POSCAR > input.pos

# Start infinite loop
while [ 1 -gt 0 ]
do

# Run MD
touch preselected.cfg

# Define the masses for classical MD
COUNTER=1
rm -f masses.txt
touch masses.txt
for i in $MASSES;
do
	echo "mass $COUNTER $i" >> masses.txt
	COUNTER=$((COUNTER + 1))
done

touch in.nb_md
cat > in.nb_md <<!
units           metal
atom_style      atomic
boundary        p p p

read_data       input.pos

include         masses.txt

pair_style mlip mlip.ini
pair_coeff * *

neighbor        0.5 bin
neigh_modify    every 1 delay 5 check yes

timestep        0.001

fix             1 all npt temp ${TEMPER} ${TEMPER} 100.0 iso 0.0 0.0 1000.0

thermo_style    custom step temp press vol pe ke etotal
thermo 1

dump   1 all custom 1 dump.nb id type x y z fx fy fz

run             100000
reset_timestep  0
!

touch mlip.ini
cat > mlip.ini <<!
mtp-filename       curr.mtp
calculate-efs      TRUE
select             TRUE
                   select:threshold        2.1
                   select:threshold-break  10.0
                   select:save-selected    preselected.cfg
                   select:load-state       state.als

select TRUE
select:threshold 2.1
select:threshold-break 10.0
select:save-selected preselected.cfg
select:load-state state.als
!

$LMP in.nb_md  # Has to run in serial because of lmp potential

# The number of MD steps with extrapolation grade above a threshold in mlip.ini
n_preselected=$(grep "BEGIN" preselected.cfg | wc -l)

# If n_preseleted is greater than zero.
if [ $n_preselected -gt 0 ]; then

    # Add configurations to the training set from preselected
    mlp select-add curr.mtp train.cfg preselected.cfg diff.cfg --als-filename=state.als
    mkdir calculate_ab_initio_ef
    cp diff.cfg calculate_ab_initio_ef

    # Clean files
    rm -f preselected.cfg
    rm -f selected.cfg

    # Calculate energies and forces and convert LAMMPS to MLIP format.
    cd calculate_ab_initio_ef
    ../../../../ab_initio_calculations.sh $MPI $VASP
    cd -

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
