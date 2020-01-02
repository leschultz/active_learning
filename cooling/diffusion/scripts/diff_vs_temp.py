from scipy.constants import physical_constants
from scipy.interpolate import interp1d
from matplotlib import pyplot as pl
from os.path import join
import pandas as pd
import numpy as np
import functions
import re

# Input paramters
df = '../../gather_data/data/data.csv'  # Path to starting data
order = 1  # The order of the interpolation
density = 10000  # The number of points to use in iterpolation
save_plots = '../figures'  # The figures save folder

df = pd.read_csv(df)

# Only include isothermal holds
df = df.loc[df['start_temperature'] == df['end_temperature']]

# Sort by temperature
df = df.sort_values(by=['end_temperature'])

groups = df.groupby(['composition', 'atoms'])

for group, values in groups:

    print(group)

    x = 1/values['end_temperature']
    y = df['self_diffusion']

    f = interp1d(x, y, kind=order)  # Linear interpolation

    xfit = np.linspace(min(x), max(x), density)
    yfit = f(xfit)

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

        ax.set_xlabel(r'1/Temperature $[1/K]$')
        ax.set_ylabel(r'Diffusion $[\AA^{2}/fs]$')
        ax.legend()

        fig.tight_layout()
        fig.savefig(save_name+'_diffusion.png')
