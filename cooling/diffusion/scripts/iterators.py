from ovito.modifiers import CalculateDisplacementsModifier
from ovito.io import import_file

from msd import msd, self_diffusion
import functions
import parsers

from matplotlib import pyplot as pl
from os.path import join
from math import floor
import pandas as pd
import numpy as np

pl.rcParams["figure.figsize"] = [8, 6]  # Define plot size

# For calculating displacements
msd_modifier = CalculateDisplacementsModifier()
msd_modifier.assume_unwrapped_coordinates = True


def analysis(
             path,
             out_print,
             incar,
             poscar,
             outcar,
             xdatcar,
             fraction,
             save_plots=False
             ):
    '''
    The analysis for a run.

    inputs:
        path = Path to run.
        out_print = The VASP print to screen file.
        incar = The INCAR file.
        poscar = The POSCAR file.
        outcar = The OUTCAR file.
        xdatcar = The XDATCAR file.
        fraction = The amount of data to average.
        save_plots = Whether to save analysis plots.

    outputs:
        comp = The composition.
        vol = The average system volume.
        press = The average system pressure.
        temp = The average system temperature.
        start_temp = The starting temperature defined in INCAR.
        end_temp = The ending temperature defined in INCAR.
        diff = The self diffusion.
    '''

    # Check if the run has any errors
    error = parsers.output(join(path, out_print))

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

    # Compute diffusion
    pipeline = import_file(join(path, xdatcar))
    pipeline.modifiers.append(msd_modifier)
    pipeline.modifiers.append(msd)

    time_modifier = float(params['POTIM'])
    frames = np.array(range(pipeline.source.num_frames))

    # Remove fraction data from beginning of frames
    frames = frames[start:]

    pipeline.modifiers[0].reference_frame = frames[0]  # Origin

    msd_calc = []
    for frame in frames:

        out = pipeline.compute(frame)
        msd_calc.append(out.attributes['msd'])

    times = frames*time_modifier
    times -= np.min(times)

    diff, m, b = self_diffusion(times, msd_calc, 6)

    if save_plots:

        # Save thermodynamic plot
        n = 3  # Number of plots
        fig, ax = pl.subplots(n)

        x = np.array(range(len(vol)))*time_modifier  # Time
        xmin = x[start]/(max(x)-min(x))  # Fraction not exact due to rounding

        ax[0].plot(x, vol, color='b', label=r'Data')
        ax[1].plot(x, press, color='b', label=r'Data')
        ax[2].plot(x, temp, color='b', label=r'Data')

        ax[0].axhline(
                      volume,
                      xmin=xmin,
                      color='g',
                      label=r'Mean Data'
                      )
        ax[1].axhline(
                      pressure,
                      xmin=xmin,
                      color='g',
                      label=r'Mean Data'
                      )
        ax[2].axhline(
                      temperature,
                      xmin=xmin,
                      color='g',
                      label=r'Mean Data'
                      )

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
            ax[i].legend(loc='upper left')

        ax[-1].set_xlabel(r'Time $[fs]$')
        fig.tight_layout()

        name = path.strip('../')
        name = name.split('/')
        name = '_'.join(name)
        therm_name = name+'_thermodynamic.png'

        # Save figure
        fig.savefig(join(save_plots, therm_name))

        pl.close('all')

        # Save MSD plot

        times_fit = np.linspace(min(times), max(times))
        msd_fit = m*times_fit+b

        fig, ax = pl.subplots()

        ax.plot(
                times,
                msd_calc,
                linestyle='none',
                marker='.',
                label='Data'
                )

        ax.plot(
                times_fit,
                msd_fit,
                linestyle='none',
                marker='.',
                label='Linear Fit (m='+str(m)+', b='+str(b)+')'
                )

        ax.set_xlabel(r'Time $[fs]$')
        ax.set_ylabel(r'MSD $[\AA^{2}]$')

        ax.legend()

        fig.tight_layout()

        msd_name = name+'_msd.png'

        # Save figure
        fig.savefig(join(save_plots, msd_name))

        pl.close('all')

    # If there is an error, use none of the paresd values
    if error:
        volume = np.nan
        pressure = np.nan
        temperature = np.nan
        diff = np.nan

    return comp, volume, pressure, temperature, start_temp, end_temp, diff


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
    diffs = []

    counts = str(len(paths))
    count = 1
    for i in paths:

        # Status
        print('Run '+'('+str(count)+'/'+counts+')'+': '+i)

        # Run data
        j, k, l, m, n, o, p = analysis(i, *args)

        runs.append(i)
        compositions.append(j)
        volumes.append(k)
        pressures.append(l)
        temperatures.append(m)
        temp_is.append(n)
        temp_fs.append(o)
        diffs.append(p)

        count += 1

    df = {
          'run': runs,
          'composition': compositions,
          'volume': volumes,
          'pressure': pressures,
          'temperature': temperatures,
          'start_temperature': temp_is,
          'end_temperature': temp_fs,
          'self_diffusion': diffs
          }

    df = pd.DataFrame(df)

    return df
