dft_job ()
{

WRKDIR=$1   # Location of scripts
COMP=$2     # The composition
RECDIR=$3   # The directory with VASP potentials
RECPOTS=$4  # The the file location for recomended VASP potentials
TYPE=$5     # The type of potential to consider
MPI=$6      # MPI
VASP=$7     # VASP
TRAIN=$8    # The training file

# Load needed functions
source $WRKDIR/funcs/gen_potcar.sh

n_selected=$(grep "BEGIN" diff.cfg | wc -l)
mlp convert-cfg diff.cfg POSCAR --output-format=vasp-poscar

for ((i=0; i<n_selected; i++))
do
    if [ $n_selected -eq 1 ]; then 
	mkdir 0
        mv POSCAR 0/POSCAR
    elif [ $n_selected -gt 1 ]; then
	mkdir $i
        mv POSCAR$i $i/POSCAR
    fi

    cd $i
    cp $WRKDIR/templates/vasp/KPOINTS .
    cp $WRKDIR/templates/vasp/dft_INCAR ./INCAR
    pot $COMP $TYPE $RECDIR $RECPOTS

    # Do DFT to get energies and forces
    $MPI $VASP

    mlp convert-cfg OUTCAR calculated.cfg --input-format=vasp-outcar
    cat calculated.cfg >> $TRAIN
    cd -

    rm -rf $i

done

}
