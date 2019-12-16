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
        volume = The average system volume.
        pressure = The average system pressure.
        temprature = The average system temperature.
        start_temp = The starting temperature defined in INCAR.
        end_temp = The ending temperature defined in INCAR.
    '''

    # INCAR parameters
    params = parsers.incar(join(path, incar))
    start_temp = params['TEBEG']  # Starting temperature
    end_temp = params['TEEND']  # Ending temperature

    # POSCAR paramters
    lattice, coords = parsers.poscar(join(path, poscar))

    # OUTCAR paramters
    composition, volumes, pressures, temperatures = parsers.outcar(join(
                                                                        path,
                                                                        outcar
                                                                        ))

    # Average the last percent amount of holds data
    start = floor(fraction*len(pressures))
    volume = np.mean(volumes[start:])
    pressure = np.mean(pressures[start:])
    temperature = np.mean(temperatures[start:])

    return composition, volume, pressure, temperature, start_temp, end_temp


def iterate(paths, *args, **kwargs):
    '''
    Analyze every run.

    inputs:
        paths = Paths to all runs.

    outputs:
        df = Data frame containing all pertinent information.
    '''

    runs = []
    compositions = []
    volumes = []
    pressures = []
    temperatures = []
    temp_is = []
    temp_fs = []

    counts = str(len(paths))
    count = 1
    for i in paths:

        # Status
        print('Run '+'('+str(count)+'/'+counts+')'+': '+i)

        # Run data
        j, k, l, m, n, o = analysis(i, *args)

        runs.append(i)
        compositions.append(j)
        volumes.append(k)
        pressures.append(l)
        temperatures.append(m)
        temp_is.append(n)
        temp_fs.append(o)

        count += 1

    df = {
          'run': runs,
          'composition': compositions,
          'volume': volumes,
          'pressure': pressures,
          'temperature': temperatures,
          'start_temperature': temp_is,
          'end_temperature': temp_fs,
          }

    df = pd.DataFrame(df)

    return df
