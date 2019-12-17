from matplotlib import pyplot as pl
from os.path import join
from math import floor
import pandas as pd
import numpy as np
import parsers


def analysis(path, incar, poscar, outcar, fraction, show_plots):
    '''
    The analysis for a run.

    inputs:
        path = Path to run.
        incar = The INCAR file.
        poscar = The POSCAR file.
        outcar = The OUTCAR file.
        fraction = The amount of data to average.
        show_plots = Wheter to display analysis plots.

    outputs:
        comp = The composition.
        vol = The average system volume.
        press = The average system pressure.
        temp = The average system temperature.
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
    comp, vol, press, temp = parsers.outcar(join(path, outcar))

    # Find minium data length because of runs stopping abruptly
    cut = min(map(len, [vol, press, temp]))
    vol = vol[:cut]
    press = press[:cut]
    temp = temp[:cut]

    # Average the last percent amount of holds data
    start = floor(fraction*len(press))
    volume = np.mean(vol[start:])
    pressure = np.mean(press[start:])
    temperature = np.mean(temp[start:])

    if show_plots:
        n = 3  # Number of plots
        fig, ax = pl.subplots(n)

        x = np.array(range(len(vol)))*float(params['POTIM'])  # Time

        ax[0].plot(x, vol, color='b', label=r'Data')
        ax[1].plot(x, press, color='b', label=r'Data')
        ax[2].plot(x, temp, color='b', label=r'Data')

        ax[0].axhline(volume, xmin=fraction, color='g', label=r'Mean Data')
        ax[1].axhline(pressure, xmin=fraction, color='g', label=r'Mean Data')
        ax[2].axhline(temperature, xmin=fraction, color='g', label=r'Mean Data')

        ax[0].set_ylabel(r'Volume $[\AA^{3}]$')
        ax[1].set_ylabel(r'Pressure $[\AA^{3}]$')
        ax[2].set_ylabel(r'Temperature $[\AA^{3}]$')
        
        for i in range(n):
            ax[i].axvline(
                          start*float(params['POTIM']),
                          linestyle=':',
                          color='r',
                          label='Averaging start'
                          )
            ax[i].legend(loc='best')

        ax[-1].set_xlabel(r'Time $[fs]$')
        fig.tight_layout()
        pl.show()
        pl.close('all')

    return comp, vol, press, temp, start_temp, end_temp


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
