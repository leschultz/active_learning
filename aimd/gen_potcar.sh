RECDIR=~/potentials/vasp/
RECPOTS="${RECDIR}vasp_pots.csv"
COMP=$1

POT=$(echo $COMP | grep -Eo '[[:alpha:]]+')

touch POTCAR
for i in $POT; do
	REC=$(cat $RECPOTS | grep $i | grep yes | awk -F ',' '{print $2}')
        echo $REC
	REC=$(find $RECDIR -type d -name "${REC}" | grep pbe)
	REC="${REC}/POTCAR"
	cat $REC >> POTCAR
done
