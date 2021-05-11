COMP=$1
TYPE=$2
RECDIR=$3
RECPOTS=$4

POT=$(echo $COMP | grep -Eo '[[:alpha:]]+')

touch POTCAR
for i in $POT; do
	REC=$(cat $RECPOTS | grep $i | grep yes | awk -F ',' '{print $2}')
	REC=$(find $RECDIR -type d -name "${REC}" | grep $TYPE)
	REC="${REC}/POTCAR"
	cat $REC >> POTCAR
done
