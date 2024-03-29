dft_job ()
{

TOPDIR=$1   # Location of scripts
MPI=$2      # MPI
VASP=$3     # VASP
TRAIN=$4    # The training file

n_selected=$(grep "BEGIN" diff.cfg | wc -l)
mlp convert-cfg diff.cfg POSCAR --output-format=vasp-poscar

for ((n=0; n<n_selected; n++))
do

    mkdir $n  # Create the directory for DFT
    if [ $n_selected -eq 1 ]; then 
        mv POSCAR $n/POSCAR
    elif [ $n_selected -gt 1 ]; then
        mv POSCAR$n $n/POSCAR
    fi

    cd $n

    # Create input files
    cp $TOPDIR/KPOINTS .
    cp $TOPDIR/POTCAR .
    cp $TOPDIR/INCAR .

    # Do DFT to get energies and forces
    $MPI $VASP

    # Clean stuff
    rm -v !("OUTCAR"|"POSCAR"|"KPOINTS"|"INCAR"|"POTCAR")

    mlp convert-cfg OUTCAR calculated.cfg --input-format=vasp-outcar --last > log.txt

    # Skip if warning happen
    if ! grep --quiet WARNING log.txt
    then
        cat calculated.cfg >> $TRAIN
    fi

    rm log.txt calculated.cfg

    cd ../

done

}
