test ()
{
# Create potential test
ATOMS=$(awk 'NR==7{print $0}' $TOPDIR/POSCAR)

# Convert to array
ELEARR=($ELEMS)
ATARR=($ATOMS)
MASARR=($MASSES)

# Get the composition
COMP=""
for i in "${!ELEARR[@]}";
do
    COMP="$COMP${ELEARR[i]}${ATARR[i]}"
done

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
for i in "${!MASARR[@]}";
do
    echo "mass $COUNTER ${MASARR[i]}  # ${ELEARR[i]}" >> masses.txt
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
    cp $TOPDIR/INCAR .
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
cd ..

}
