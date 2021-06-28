aimd_job ()
{

WRKDIR=$1   # Location of scripts
COMP=$2     # The composition
RECDIR=$3   # The directory with VASP potentials
RECPOTS=$4  # The the file location for recomended VASP potentials
TYPE=$5     # The type of potential to consider

# Load needed functions
source $WRKDIR/funcs/gen_potcar.sh

# Create the poscar
python3 $WRKDIR/funcs/gen_poscar.py $COMP

# Create poteantial
pot $COMP $TYPE $RECDIR $RECPOTS

# Create KPOINTS
cp $WRKDIR/templates/vasp/KPOINTS .

# Create INCAR
cp $WRKDIR/templates/vasp/INCAR .

# Get the masses for the elements
MASSES=$(cat POTCAR | grep MASS | awk -F ' ' '{print $3}' | sed 's/\;//g')
echo $MASSES

}
