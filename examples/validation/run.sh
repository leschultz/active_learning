COMP=Si50
VASP='vasp_std'
LMP='lmp_intel_cpu_intelmpi'
MPI='prun'

$MPI $LMP < md.in

mkdir runs
cd runs

python3 ../convert.py ../dump.atom $COMP

mkdir potential
cp ../curr.mtp ./potential
touch ./potential/test.cfg

for i in $(find ./ -type f -name POSCAR);
do
    cd $(dirname $i)

    cp ../../KPOINTS .
    cp ../../INCAR .
    cp ../../POTCAR .

    $MPI $VASP

    mlp convert-cfg --input-format=vasp-outcar OUTCAR test.cfg 
    cat test.cfg >> ../potential/test.cfg
    cd -
done

cd potential
mlp calc-errors curr.mtp test.cfg > test_info.txt
mlp calc-efs curr.mtp test.cfg test_pred.cfg

python3 ../../parity.py
