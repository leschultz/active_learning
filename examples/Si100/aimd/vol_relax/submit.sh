#PBS -S /bin/bash
#PBS -q cpu16m64
#PBS -l select=1:ncpus=16:mpiprocs=16
#PBS -l walltime=72:00:00
#PBS -N job

# available queues: cpu16m64, cpu20m128
# replace the 1 in select= with the number of nodes (but use only 1 - no working IB)
# replace the 16 in ncpus and mpiprocs with 20 depending on the queue
# specify walltime in the HH:MM:SS format
# replace JOB_NAME with your desired job name

cd $PBS_O_WORKDIR 

module load vasp

# use prun with no arguments, no matter what MPI environment is used
date
time prun vasp_std
date
