#!/bin/bash

# Programs
CORES=$(nproc)
CORES=2  # Delete this
MPI="mpirun -n ${CORES}"
LMP="lmp_mpi -in"
VASP='vasp_std'

# VASP inputs
RECDIR=~/potentials/vasp/
RECPOTS="${RECDIR}vasp_pots.csv"
TYPE=pbe
COMP=Zr5Cu5
TEMPER=2000

# MLIP-2 inputs
POTS=~/packages/mlip-2/untrained_mtps/
POT=08.mtp

# Clean working space
rm -rf initial
rm -f train.cfg
rm -f curr.mtp
rm -f preselected.cfg
rm -f diff.cfg
rm -f selected.cfg
rm -f out.cfg

# Initial VASP run
mkdir $TEMPER
cd $TEMPER

mkdir vasp
cd vasp

# Create initial features
python3 ../../gen_poscar.py $COMP  # Creates random initial positions
bash ../../gen_potcar.sh $COMP $TYPE $RECDIR $RECPOTS  # Creates potential

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
NELM = 10                      # maximum of steps per time step

# MD (do little writing to save disc space)
IBRION = 0                     # main molecular dynamics tag
NSW = 2                        # number of MD steps
POTIM = 3                      # time step of MD [fs]
NWRITE = 0                     # controls output
LCHARG = .FALSE.               # no charge density written out
LWAVE = .FALSE.                # no wave function coefficients written out
TEBEG = $TEMPER                # starting temperature for MD
TEEND = $TEMPER                # end temperature for MD

# canonic (Nose) MD with XDATCAR updated every 10 steps
MDALGO = 2                     # switch to select thermostat
SMASS =  3                     # Nose mass
ISIF = 2                       # this tag selects the ensemble in combination with the thermostat
!

# Run AIMD
$MPI $VASP

cd ..
mkdir pot
cd pot

# Create training file
mlp convert-cfg --input-format=vasp-outcar ../vasp/OUTCAR train.cfg

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

# Start infinite loop
while [ 1 -gt 0 ]
do

# Run MD
touch preselected.cfg
../../pos2lmp.awk ../vasp/POSCAR > input.pos

touch in.nb_md
cat > in.nb_md <<!
variable        Tequil index $TEMPER
variable        seed index 826626413

units           metal
atom_style      atomic
boundary        p p p

read_data        input.pos

mass            * 92.90638

pair_style mlip mlip.ini
pair_coeff * *

neighbor        0.5 bin
neigh_modify    every 1 delay 5 check yes

timestep        0.001

fix             1 all nve
fix             2 all langevin ${Tequil} ${Tequil} 0.1 ${seed} zero yes

thermo_style    custom step temp
thermo 1000

dump   1 all custom 1000 dump.nb id type x y z fx fy fz

run             100000
reset_timestep  0
!

$MPI $LMP in.nb_md

# The number of MD steps with extrapolation grade above a threshold in mlip.ini
n_preselected=$(grep "BEGIN" preselected.cfg | wc -l)

# If n_preseleted is greater than zero.
if [ $n_preselected -gt 0 ]; then

    # Add configurations to the training set from preselected
    mlp select-add curr.mtp train.cfg preselected.cfg diff.cfg --als-filename=state.als
    cp diff.cfg calculate_ab_initio_ef/

    # Clean files
    rm -f preselected.cfg
    rm -f selected.cfg

    # Calculate energies and forces and convert LAMMPS to MLIP format.
    cd ../calculate_ab_initio_ef/
    ./ab_initio_calculations.sh "${MPI} ${LMP}"
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
