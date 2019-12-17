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
show_plots = False  # Display plots or not
save_dir = '../data'  # The save folder
save_name = 'data.csv'  # The save name

df = pd.read_csv(df)
df = df.sort_values(by=['volume'])
groups = df.groupby(['composition'])

compositions = []
temperatures = []
volumes = []

for group, values in groups:
    print('Grouped by: '+str(group))

    y = values['pressure'].values
    x = values['volume'].values

    for i in range(y.shape[0]):
        if y[i] <= 0:
            break

    pos = i-1
    neg = i

    xnew = [x[neg], x[pos]]
    ynew = [y[neg], y[pos]]

    m, b, r, p, std = linregress(xnew, ynew)
    line = lambda x: m*x+b

    xfit = np.linspace(min(xnew), max(xnew), density)
    yfit = np.array(list(map(line, xfit)))

    index = functions.nearest(0, yfit)

    compositions.append(np.unique(values['composition'])[0])
    temperatures.append(np.unique(values['end_temperature'])[0])
    volumes.append(xfit[index])

    if show_plots:
        fig, ax = pl.subplots()
        ax.plot(xfit, yfit, linestyle='none', marker='.')
        ax.plot(xfit[index], yfit[index], marker='8')
        ax.plot(x, y, linestyle='none', marker='*')
        ax.set_xlabel(r'Volumes $[\AA^{3}]$')
        ax.set_ylabel(r'Pressure $[kB]$')
        pl.show()
        pl.close('all')

df = {
      'composition': compositions,
      'temperature': temperatures,
      'volume': volumes,
      }

df = pd.DataFrame(df)

# Create directory and save data
functions.create_dir(save_dir)
save = join(save_dir, save_name)
df.to_csv(save, index=False)

print(df)
