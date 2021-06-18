#!/bin/bash

# VASP inputs
RECDIR=~/potentials/vasp/
RECPOTS="${RECDIR}vasp_pots.csv"
TYPE=pbe
COMP=Zr5Cu5
CORES=$(nproc)

# MLIP-2 inputs
POTS=~/packages/mlip-2/untrained_mtps/
POT=08.mtp

python3 gen_poscar.py $COMP
bash gen_potcar.sh $COMP $TYPE $RECDIR $RECPOTS

cp POSCAR POSCAR_start

# Train the initial potential

# Create initial training file
mkdir initial
cd initial
touch INCAR
cat > INCAR <<!

SYSTEM =  $COMP_$i
ISMEAR = 0  ! Gaussian smearing
!

cp ../POSCAR .
cp ../KPOINTS .
cp ../POTCAR .

mpirun -np $CORES vasp_std

cp ./XDATCAR ../POSCAR

cd -

# Create training file
mlp convert-cfg --input-format=vasp-outcar initial/OUTCAR train_init.cfg

# Prepare potential
NELS=$(echo $COMP | grep -Eo '[[:alpha:]]+' | wc -w)
MINDIST=$(mlp mindist train_init.cfg | awk -F ': ' '{print $2}')
cp "${POTS}${POT}" ./init.mtp
sed -i "s/species_count.*/species_count = $NELS/g" init.mtp
sed -i "s/min_dist.*/min_dist = $MINDIST/g" init.mtp

# Create Initial potential
mpirun -n $CORES mlp train init.mtp train_init.cfg --trained-pot-name=init.mtp

# Train the potential
#mpirun -n $CORES mlp train init.mtp train.cfg --trained-pot-name=potential.mtp

: <<'END'

# Run quench
for i in 2000 1800
do

mkdir $i
cd $i

# Do isothermal run
touch INCAR
cat > INCAR <<!

SYSTEM =  $COMP_$i

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
NSW = 10                       # number of MD steps
POTIM = 3                      # time step of MD
NWRITE = 0                     # controls output
LCHARG = .FALSE.               # no charge density written out
LWAVE = .FALSE.                # no wave function coefficients written out
TEBEG = $i                     # starting temperature for MD
TEEND = $i                     # end temperature for MD

# canonic (Nose) MD with XDATCAR updated every 10 steps
MDALGO = 2                     # switch to select thermostat
SMASS =  3                     # Nose mass
ISIF = 2                       # this tag selects the ensemble in combination with the thermostat
!

cp ../POSCAR .
cp ../KPOINTS .
cp ../POTCAR .

mpirun -np $CORES vasp_std

cp ./XDATCAR ../POSCAR

cd -

done

mv POSCAR_start POSCAR
END
