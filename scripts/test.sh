# Create potential test
TOPDIR=$(pwd)
POTDIR=$TOPDIR/run/potential
MASSES=$(cat POTCAR | grep MASS | awk -F ' ' '{print $3}' | sed 's/\;//g')
ELEMS=$(cat POTCAR | grep 'VRHFIN =' | grep -oP '(?<==).*?(?=:)')

cd $POTDIR
cd ..
cp -r $POTDIR potential_test
cd potential_test

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

# Convert to POSCARs for VASP
python3 $WRKDIR/funcs/convert.py dump.atom $COMP

touch test.cfg

for i in $(find ./ -type f -name POSCAR);
do
    cd $(dirname $i)

    cp ../KPOINTS .
    cp ../TEST_INCAR ./INCAR
    cp ../POTCAR .

    $MPI $VASP

    mlp convert-cfg --input-format=vasp-outcar OUTCAR test.cfg 
    cat test.cfg >> ../potential/test.cfg
    cd -
done

mlp calc-errors curr.mtp test.cfg > test_info.txt
mlp calc-efs curr.mtp test.cfg test_pred.cfg

python3 ../../parity.py
