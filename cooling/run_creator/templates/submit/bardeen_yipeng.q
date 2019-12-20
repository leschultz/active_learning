#!/bin/bash

# Declare a name for this job
#$ -N parallel

# Request the queue for this job (yipeng.q, morgan.q)
#$ -q yipeng.q

# Resources
# yipeng.q: 2.2 GHz/core (up to 20 cores)
# morgan.q: 2.5 GHz/core (up to 16 cores) 

# Request up to 48 hours (yipeng.q) or 168 hours (morgan.q) of wall time (hh:mm:ss)
#$ -l h_rt=336:00:00

# Run the job from the directory of submission. Uncomment only if you don't want the defults.
#$ -cwd

# Request processors saved as variable NSLOTS (replace <num> with up to 20 (yipeng.q) or 16 (morgan.q))
#$ -pe intel-mpi 20

# combine SGE standard output and error files
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID

# Transfer all your environment variables.
#$ -V

# Load intel64 MPI environment variables
source /share/apps/intel/parallel_studio_xe_2016.4.072/psxevars.sh intel64

# The executable for parallel jobs
mpirun -n $NSLOTS vasp_std > out.txt

DT="150"  # Change in temperature
MIN="200"  # Minimum temperature
python3 ../../../scripts/step_down.py "${DT}" "${MIN}"
