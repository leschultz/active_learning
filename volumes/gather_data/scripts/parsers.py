import pandas as pd
import numpy as np


def incar(path):
    '''
    Parse INCAR file.

    inputs:
        path = The file path.

    outputs:
        param = The input parameters.
    '''

    param = {}
    with open(path) as f:
        for line in f:
            line = line.strip().split(' ')

            if '=' in line:
                param[line[0]] = line[2]

    return param


def poscar(path):
    '''
    Parse POSCAR file.

    inputs:
        path = The file path.

    outputs:
        lattice = The lattice parameters.
        coords = The atom coordinates.
    '''

    lattice = []
    coords = []

    begin_coords = 0
    count = 0
    with open(path) as f:
        for line in f:
            line = line.strip().split(' ')

            if (count > 1) and (count < 5):
                lattice.append(list(map(lambda i: float(i), line)))

            if begin_coords == 1:
                coords.append(line)

            if line[0] == 'direct':
                begin_coords = 1

            count += 1

    lattice = np.array(lattice)
    coords = pd.DataFrame(coords, columns=['x', 'y', 'z', 'element'])

    return lattice, coords


def outcar(path):
    '''
    Parse OUTCAR file.

    inputs:
        path = The file path.

    outputs:
        volumes = The volume data.
        pressures =  The pressure data.
        temperatures = The temperature data.
        etotal = The total energy of the system.
    '''

    volumes = []
    pressures = []
    temperatures = []
    etotal = []

    iteration = 0
    with open(path) as f:
        for line in f:

            line = line.strip().split(' ')
            line = [i for i in line if i != '']

            if 'POSCAR' in line:
                composition = ''.join(line[2:])

            if 'Iteration' in line:
                iteration += 1

            if iteration > 0:

                if ('volume' in line) and ('cell' in line):
                    volumes.append(float(line[-1]))

                if ('external' in line) and ('pressure' in line):
                    pressures.append(float(line[3]))

                if '(temperature' in line:
                    temperatures.append(float(line[5]))

                if ('free' in line) and ('TOTEN' in line):
                    etotal.append(float(line[4]))

    volumes = np.array(volumes)
    pressures = np.array(pressures)
    temperatures = np.array(temperatures)
    etotal = np.array(etotal)

    return composition, volumes, pressures, temperatures, etotal
