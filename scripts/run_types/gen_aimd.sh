aimd_job ()
{

TOPDIR=$1   # Location of scripts

# Create input files
cp $TOPDIR/POTCAR .          # POTCAR
cp $TOPDIR/KPOINTS .         # KPOINTS
cp $TOPDIR/POSCAR .          # POSCAR
cp $TOPDIR/MD_INCAR ./INCAR  # INCAR
}
