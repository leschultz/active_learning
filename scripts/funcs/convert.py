from ase.io.lammpsrun import read_lammps_dump_text
from ase.io import read, write
import sys
import os

if __name__ == '__main__':

    a = read(sys.argv[1], index=':')

    count = 0
    for i in a:

        i.symbols = sys.argv[2]
        os.system('mkdir {}'.format(count))
        write('{}/POSCAR'.format(count), i)

        count += 1
