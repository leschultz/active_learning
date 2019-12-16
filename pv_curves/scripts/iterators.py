from matplotlib import pyplot as pl

from os.path import join
from math import floor
import pandas as pd
import numpy as np
import parsers


def analysis(path, incar, poscar, outcar, fraction):
    '''
    The analysis for a run.

    inputs:
        path = Path to run.
        incar = The INCAR file.
        poscar = The POSCAR file.
        outcar = The OUTCAR file.
        fraction = The amount of data to average.

    outputs:
        volume = The average system volume
        pressure = The average system pressure
    '''

    # INCAR parameters
    params = parsers.incar(join(path, incar))
    
    # POSCAR paramters
    lattice, coords = parsers.poscar(join(path, poscar))

    # OUTCAR paramters
    iterations, volumes, pressures = parsers.outcar(join(path, outcar))

    # Average the last percent amount of holds data
    start = floor(fraction*len(pressures))
    volume = np.mean(volumes[start:])
    pressure = np.mean(pressures[start:])

    return volume, pressure


def iterate(paths, *args, **kwargs):
    '''
    Analyze every run.

    inputs:
        paths = Paths to all runs.

    outputs:
        df = Data frame containing all pertinent information.
    '''

    runs = []
    volumes = []
    pressures = []

    counts = str(len(paths))
    count = 1
    for i in paths:
        print('Run '+'('+str(count)+'/'+counts+')'+': '+i)
        volume, pressure = analysis(i, *args)

        runs.append(i)
        volumes.append(volume)
        pressures.append(pressure)

        count += 1

    df = {
          'run': runs,
          'volume': volumes,
          'pressure': pressures,
            }

    df = pd.DataFrame(df)

    return df
