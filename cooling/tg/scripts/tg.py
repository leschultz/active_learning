from scipy.constants import physical_constants
from scipy.interpolate import interp1d
from matplotlib import pyplot as pl
from os.path import join
from knee import opt
import pandas as pd
import numpy as np
import functions
import re

# Input paramters
df = '../../gather_data/data/data.csv'  # Path to starting data
order = 1  # The order of the interpolation
density = 10000  # The number of points to use in iterpolation
data_save = '../data'  # The data save location
save_plots = '../figures'  # The figures save folder

df = pd.read_csv(df)

# Only include isothermal holds
df = df.loc[df['start_temperature'] == df['end_temperature']]

# Sort by temperature
df = df.sort_values(by=['end_temperature'])

# Energy normalized by atoms
df['total_energy/atoms'] = df['total_energy']/df['atoms']

# Calculate E-3KT (not showing expected behavior)
kb = physical_constants['Boltzmann constant in eV/K'][0]
df['E-3kT'] = df['total_energy/atoms']-3*kb*df['end_temperature']

groups = df.groupby(['composition', 'atoms'])

compositions = []
atoms = []
tgs = []

for group, values in groups:

    print(group)

    x = values['end_temperature']
    y = values['E-3kT']

    f = interp1d(x, y, kind=order)  # Linear interpolation

    xfit = np.linspace(min(x), max(x), density)
    yfit = f(xfit)

    xcut, endpoints, middle_rmse = opt(xfit, yfit)

    compositions.append(group[0])
    atoms.append(group[1])
    tgs.append(xcut)

    if save_plots:

        # Create saving directory for figures
        functions.create_dir(save_plots)

        save_name = join(save_plots, str(group))

        # Plot glass transition curve
        fig, ax = pl.subplots()

        ax.plot(
                x,
                y,
                linestyle='none',
                marker='.',
                label='Data',
                color='b'
                )
        ax.plot(
                xfit,
                yfit,
                label='Interpolation (order='+str(order)+')',
                color='g'
                )
        ax.axvline(
                   x=xcut,
                   label='Tg='+str(xcut)+' [K]',
                   color='r'
                   )

        ax.set_xlabel(r'Temperature $[K]$')
        ax.set_ylabel(r'E-3kT $[eV/atom]$')
        ax.legend()

        fig.tight_layout()
        fig.savefig(save_name+'_tg.png')

        # Plot error curve
        fig, ax = pl.subplots()

        ax.plot(endpoints, middle_rmse, label='Data', color='g')
        ax.axvline(x=xcut, label='Minimum', color='r')

        ax.set_xlabel(r'Temperature $[K]$')
        ax.set_ylabel(r'RMSE $[eV]$')
        ax.legend()

        fig.tight_layout()
        fig.savefig(save_name+'_rmse.png')

    # Create data saving directory
    functions.create_dir(data_save)

    data = {
            'composition': compositions,
            'atoms': atoms,
            'tg': tgs,
            }

    data = pd.DataFrame(data)
    data.to_csv(join(data_save, 'data.csv'), index=False)
