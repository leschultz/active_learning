from matplotlib import pyplot as pl
from scipy.optimize import curve_fit
from scipy.stats import linregress
from os.path import join
import pandas as pd
import numpy as np
import functions

# Input parameters
df = '../../gather_data/data/data.csv'  # Path to data file
density = 10000  # The number of points to use in fit
data_save_dir = '../data'  # The data save folder
data_save_name = 'data.csv'  # The data save name
save_plots = '../figures'  # The figures save folder

# Load Data
df = pd.read_csv(df)
df = df.sort_values(by=['volume', 'start_temperature', 'end_temperature'])

# Group data
groups = df.groupby(['composition', 'start_temperature', 'end_temperature'])

# Create figures directory
functions.create_dir(save_plots)

compositions = []
temperatures = []
lengths = []
length_errors = []
pressures = []

for group, values in groups:
    name = str(group)
    print('Grouped by: '+name)

    y = values['pressure'].values
    x = values['volume'].values
    x **= 1/3

    for i in range(y.shape[0]):
        if y[i] <= 0:
            break

    pos = i-1  # Positive pressure
    neg = i  # Negative pressure

    xnew = [x[neg], x[pos]]
    ynew = [y[neg], y[pos]]

    m, b, r, p, std = linregress(xnew, ynew)

    xfit = np.linspace(min(xnew), max(xnew), density)
    yfit = np.array(list(map(lambda x: m*x+b, xfit)))

    index = functions.nearest(0, yfit)

    compositions.append(np.unique(values['composition'])[0])
    temperatures.append(np.unique(values['end_temperature'])[0])
    lengths.append(xfit[index])
    length_errors.append(x[neg]-x[pos])
    pressures.append(yfit[index])

    if save_plots:
        fig, ax = pl.subplots()

        ax.plot(
                xfit,
                yfit,
                label='Extrapolation'
                )
        ax.plot(
                xfit[index],
                yfit[index],
                marker='8',
                linestyle='none',
                label='Zero Pressure',
                )
        ax.plot(
                x,
                y,
                linestyle='none',
                marker='*',
                label='Data'
                )

        ax.set_xlabel(r'Cube Length $[\AA]$')
        ax.set_ylabel(r'Pressure $[kB]$')

        ax.legend()

        fig.savefig(join(save_plots, name))
        pl.close('all')

df = {
      'composition': compositions,
      'temperature': temperatures,
      'cube_length': lengths,
      'cube_length_error': length_errors,
      'pressure': pressures,
      }

df = pd.DataFrame(df)

# Create directory and save data
functions.create_dir(data_save_dir)
data_save = join(data_save_dir, data_save_name)
df.to_csv(data_save, index=False)

print('Saved: '+str(data_save))
