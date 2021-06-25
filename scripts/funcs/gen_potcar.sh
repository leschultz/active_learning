pot ()
{

COMP=$1     # Composition as element and number (e.g. Zr5Cu50Al45)
TYPE=$2     # Type of VASP potential 
RECDIR=$3   # Directory containing VASP potentials
RECPOTS=$4  # File containing recomended potentials

POT=$(echo $COMP | grep -Eo '[[:alpha:]]+')

touch POTCAR
for i in $POT; do
    REC=$(cat $RECPOTS | grep $i | grep yes | awk -F ',' '{print $2}')
    REC=$(find $RECDIR -type d -name "${REC}" | grep $TYPE)
    REC="${REC}/POTCAR"
    cat $REC >> POTCAR
done

}
