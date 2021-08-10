md_job ()
{

TOPDIR=$1   # Location of job
MASSES=$2   # The masses for each element
ELEMS=$3    # The actual elements

# Convert to arrays
MASARR=($MASSES)
ELEARR=($ELEMS)

# Define the masses for classical MD
COUNTER=1
rm -f masses.txt
touch masses.txt
for i in "${!MASARR[@]}";
do
    echo "mass $COUNTER ${MASARR[i]}  # ${ELEARR[i]}" >> masses.txt
    COUNTER=$((COUNTER + 1))
done

# The potential to use
cp $TOPDIR/mlip.ini .

# The LAMMPS input file
cp $TOPDIR/md.in .
}
