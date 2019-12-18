from matplotlib import pyplot as pl
from scipy.stats import linregress
from os.path import join
import pandas as pd
import numpy as np
import functions

# Input parameters
df = '../../lv_curves/data/data.csv'  # Path to data file
data_save_dir = '../data'  # The data save folder
data_save_name = 'data.csv'  # The data save name
save_plots = '../figures'  # The figures save folder

# Load Data
df = pd.read_csv(df)

# Group data
groups = df.groupby(['composition'])

# Create figures directory
functions.create_dir(save_plots)

compositions = []
slopes = []
intercepts = []

for group, values in groups:
    name = str(group)
    print('Grouped by: '+name)

    x = values['temperature'].values
    y = values['cube_length'].values

    m, b, r, p, std = linregress(x, y)

    compositions.append(np.unique(values['composition'])[0])
    slopes.append(m)
    intercepts.append(b)

    if save_plots:

        equation = '('+str(m)+')x+'+'('+str(b)+')'
        xfit = np.linspace(min(x), max(x))
        yfit = m*xfit+b

        fig, ax = pl.subplots()

        ax.plot(xfit, yfit, label='Data')
        ax.plot(x, y, linestyle='none', marker='.', label=equation)

        ax.set_ylabel(r'Cube Length $[\AA]$')
        ax.set_xlabel(r'Temperature $[K]$')

        ax.legend()

        fig.tight_layout()

        fig.savefig(join(save_plots, name))
        pl.close('all')

df = {
      'composition': compositions,
      'slope': slopes,
      'intercept': intercepts,
      }

df = pd.DataFrame(df)

# Create directory and save data
functions.create_dir(data_save_dir)
data_save = join(data_save_dir, data_save_name)
df.to_csv(data_save, index=False)

print('Saved: '+str(data_save))
