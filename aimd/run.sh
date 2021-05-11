RECPOTS=~/potentials/vasp/vasp_pots.csv
COMP=Zr5Cu5
CORES=$(nproc)

python3 gen_poscar.py $COMP
bash gen_potcar.sh $COMP

POT=$(cat $RECPOTS | grep yes)

echo $POT

cp POSCAR POSCAR_start

for i in 2000 1300
do

mkdir $i
cd $i

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
NSW = 2                        # number of MD steps
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
