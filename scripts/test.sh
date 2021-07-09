# Create potential test
TOPDIR=$(pwd)
POTDIR=$TOPDIR/run/potential
MASSES=$(cat POTCAR | grep MASS | awk -F ' ' '{print $3}' | sed 's/\;//g')
ELEMS=$(cat POTCAR | grep 'VRHFIN =' | grep -oP '(?<==).*?(?=:)')
ATOMS=$(awk 'NR==7{print $0}' POSCAR)

# Get the composition
for i in "${!MASSES[@]}";
do
    COMP="${ELEMS[i]}${ATOMS[i]}"
done

cd $POTDIR
cd ..
rm -rf test
mkdir test
cd test

mkdir md
cd md
cp -r $POTDIR/* .

# Run LAMMPS that selects test trajectories
cp $TOPDIR/md_test.in .
$WRKDIR/convert/poscar2lammps.awk $TOPDIR/POSCAR > input.pos
sed -i '/select/Q' mlip.ini

# Define the masses for classical MD
COUNTER=1
rm -f masses.txt
touch masses.txt
for i in "${!MASSES[@]}";
do
    echo "mass $COUNTER ${MASSES[i]}  # ${ELEMS[i]}" >> masses.txt
    COUNTER=$((COUNTER + 1))
done

$MPI $LMP md_test.in
cd ..

# Convert to POSCARs for VASP
mkdir aimd
cd aimd

python3 $WRKDIR/funcs/convert.py ../md/dump.atom $COMP

# Prepare test configurations
touch test.cfg
for i in $(find ./ -type f -name POSCAR);
do
    cd $(dirname $i)

    cp $TOPDIR/KPOINTS .
    cp $TOPDIR/TEST_INCAR ./INCAR
    cp $TOPDIR/POTCAR .

    $MPI $VASP

    mlp convert-cfg --input-format=vasp-outcar OUTCAR test.cfg 
    cat test.cfg >> ../test.cfg
    cd -
done

cd ..
mkdir ml
cd ml

cp ../md/curr.mtp .
cp ../aimd/test.cfg .
mlp calc-efs curr.mtp test.cfg test_pred.cfg

python3 $WRKDIR/funcs/parity.py
