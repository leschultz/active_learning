from scipy.interpolate import interp1d
from matplotlib import pyplot as pl
from scipy.stats import linregress
from os.path import join
import pandas as pd
import numpy as np
import functions

# Input parameters
df = '../../gather_data/data/data.csv'  # Path to data file
density = 100000  # The number of points to use in fit
data_save_dir = '../data'  # The data save folder
data_save_name = 'data.csv'  # The data save name
save_plots = '../figures'  # The figures save folder
order = 3  # The order of the spline fit

# Load Data
df = pd.read_csv(df)
df = df.dropna()  # Drop bad jobs
df = df.sort_values(by=['volume', 'start_temperature', 'end_temperature'])

# Group data
groups = df.groupby(['composition', 'start_temperature', 'end_temperature'])

# Create figures directory
functions.create_dir(save_plots)

compositions = []
temperatures = []
lengths = []
length_errors = []
pressure_errors = []

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

    xfit = np.linspace(min(x), max(x), density)

    # Attempt spline
    spline_fail = False
    try:

        f2 = interp1d(x, y, kind=order)
        yfit = f2(xfit)

        index = functions.nearest(0, yfit)
        lengths.append(xfit[index])

    except Exception:
        spline_fail = True

        m, b, r, p, std = linregress(xnew, ynew)
        yfit = np.array(list(map(lambda x: m*x+b, xfit)))

        index = functions.nearest(0, yfit)
        lengths.append(xfit[index])

    length_errors.append(min(abs(xnew-xfit[index])))
    pressure_errors.append(min(abs(ynew-xfit[index])))

    compositions.append(np.unique(values['composition'])[0])
    temperatures.append(np.unique(values['end_temperature'])[0])

    if save_plots:

        if spline_fail:
            label = '('+str(m)+')x+('+str(b)+')'

        else:
            label = 'Spline Order '+str(order)

        fig, ax = pl.subplots()

        ax.plot(
                xfit,
                yfit,
                label=label
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
      'pressure_error': pressure_errors,
      }

df = pd.DataFrame(df)

# Create directory and save data
functions.create_dir(data_save_dir)
data_save = join(data_save_dir, data_save_name)
df.to_csv(data_save, index=False)

print('Saved: '+str(data_save))
