#!/bin/sh
#SBATCH --partition=shared
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=4000
#SBATCH --error=job.e.%J
#SBATCH --output=job.o.%J

module load lammps vasp/5.4.4 mlp

./active_md.sh
