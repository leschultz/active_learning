import pandas as pd
import numpy as np


def output(path):
    '''
    Parse the VASP print to screen file.

    inputs:
        path = The file path.

    outputs:
        error = Whether or not there was an error.
    '''

    # Check if the run has any errors
    error = False
    with open(path) as f:
        for line in f:
            if 'ERROR' in line:
                error = True

    return error


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
    '''

    df = {}  # Avoid unusual printing occasionally seen

    finished = False
    iteration = 0
    with open(path) as f:
        for line in f:
            line = line.strip().split(' ')
            line = [i for i in line if i != '']

            if 'Voluntary' in line:
                finished = True

            if 'POSCAR' in line:
                composition = ''.join(line[2:])

            if 'Iteration' in line:
                index = line.index('Iteration')+1
                iteration = int(line[index].split('(')[0])

                df[iteration] = {
                                 'volume': np.nan,
                                 'pressure': np.nan,
                                 'temperature': np.nan,
                                 }

            if iteration > 0:

                if ('volume' in line) and ('cell' in line):
                    df[iteration]['volume'] = float(line[-1])

                if ('external' in line) and ('pressure' in line):
                    df[iteration]['pressure'] = float(line[3])

                if '(temperature' in line:
                    df[iteration]['temperature'] = float(line[5])

                if 'energy(sigma->0)' in line:
                    df[iteration]['total_energy'] = float(line[-1])

    # Create dataframe
    df = pd.DataFrame(df).T

    # Prepare values for export
    volumes = df['volume'].values
    pressures = df['pressure'].values
    temperatures = df['temperature'].values
    total_energy = df['total_energy'].values

    # Return nan for not finished jobs
    if not finished:
        volumes[:] = np.nan
        pressures[:] = np.nan
        temperatures[:] = np.nan
        total_energy[:] = np.nan

    return composition, volumes, pressures, temperatures, total_energy


def oszicar(path):
    '''
    Parse OSZICAR file.

    inputs:
        path = The file path.

    outputs:
        T = The temperature.
        E = The total free energy.
        F = The total free energy (no kinetic contribution).
        E0 = The energy for sigma -> 0.
        EK = The kinetic energy.
        SP = The potential energy of the Nose thermostat.
        SK = The kinetic energy of the Nose thermostat.
    '''

    T = []
    E = []
    F = []
    E0 = []
    EK = []
    SP = []
    SK = []

    with open(path) as f:
        for line in f:
            line = line.strip().split(' ')
            line = [i for i in line if i != '']

            if 'T=' in line:

                T.append(float(line[2]))
                E.append(float(line[4]))
                F.append(float(line[6]))
                E0.append(float(line[8]))
                EK.append(float(line[10]))
                SP.append(float(line[12]))
                SK.append(float(line[14]))

    return T, E, F, E0, EK, SP, SK
