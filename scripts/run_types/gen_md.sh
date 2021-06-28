md_job ()
{

WRKDIR=$1   # Location of scripts
MASSES=$2   # The masses for each element

# Define the masses for classical MD
COUNTER=1
rm -f masses.txt
touch masses.txt
for i in $MASSES;
do
    echo "mass $COUNTER $i" >> masses.txt
    COUNTER=$((COUNTER + 1))
done

# The potential to use
cp $WRKDIR/templates/lammps/mlip.ini .

# The LAMMPS input file
cp $WRKDIR/templates/lammps/md.in .
}
