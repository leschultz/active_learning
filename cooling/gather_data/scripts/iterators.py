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
import re

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
             oszicar,
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
        oszicar = The OSZICAR file.
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
    comp, vol, press, temp, etot = parsers.outcar(join(path, outcar))

    # Count atoms from composition
    atoms = re.split('(\d+)', comp)
    atoms = sum([int(i) for i in atoms if i.isdigit()])

    # OSZICAR parameters
    T, E, F, E0, EK, SP, SK = parsers.oszicar(join(path, oszicar))

    etot = E  # Use the energy values from OSZICAR instead

    # Find minium data length because of runs stopping abruptly
    cut = min(map(len, [vol, press, temp, etot]))
    vol = vol[:cut]
    press = press[:cut]
    temp = temp[:cut]
    etot = etot[:cut]

    # Average the last percent amount of holds data
    start = floor(fraction*len(press))
    volume = np.mean(vol[start:])
    pressure = np.mean(press[start:])
    temperature = np.mean(temp[start:])
    etotal = np.mean(etot[start:])

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
        n = 4  # Number of plots
        fig, ax = pl.subplots(n)

        x = np.array(range(len(vol)))*time_modifier  # Time
        xmin = x[start]/(max(x)-min(x))  # Fraction not exact due to rounding

        ax[0].plot(x, vol, color='b', label=r'Data')
        ax[1].plot(x, press, color='b', label=r'Data')
        ax[2].plot(x, temp, color='b', label=r'Data')
        ax[3].plot(x, etot, color='b', label=r'Data')

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
        ax[3].axhline(
                      etotal,
                      xmin=xmin,
                      color='g',
                      label=r'Mean Data'
                      )

        ax[0].set_ylabel(r'Volume $[\AA^{3}]$')
        ax[1].set_ylabel(r'Pressure $[kB]$')
        ax[2].set_ylabel(r'Temperature $[K]$')
        ax[3].set_ylabel(r'Total Energy $[eV]$')

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
        etotal = np.nan
        diff = np.nan

    return comp, volume, pressure, temperature, etotal, start_temp, end_temp, diff, atoms


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
    etots = []
    temp_is = []
    temp_fs = []
    diffs = []
    atoms = []

    counts = str(len(paths))
    count = 1
    for i in paths:

        # Status
        print('Run '+'('+str(count)+'/'+counts+')'+': '+i)

        # Run data
        out = analysis(i, *args)

        runs.append(i)
        compositions.append(out[0])
        volumes.append(out[1])
        pressures.append(out[2])
        temperatures.append(out[3])
        etots.append(out[4])
        temp_is.append(out[5])
        temp_fs.append(out[6])
        diffs.append(out[7])
        atoms.append(out[8])

        count += 1

    df = {
          'run': runs,
          'composition': compositions,
          'volume': volumes,
          'pressure': pressures,
          'temperature': temperatures,
          'total_energy': etots,
          'start_temperature': temp_is,
          'end_temperature': temp_fs,
          'self_diffusion': diffs,
          'atoms': atoms,
          }

    df = pd.DataFrame(df)

    return df
