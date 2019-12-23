from matplotlib import pyplot as pl
from pymatgen import Composition
from scipy import stats

import pandas as pd
import numpy as np
import sys
import re
import os

figdir = '../figures'
datadir = '../data'
df = '../data/dmax.csv'

df = pd.read_csv(df)
df = df.sort_values(by=['dmax'])

species = []
for i in df['composition']:
    specie = Composition(i)
    specie  = set(str(specie).split(' '))
    species.append(specie)

df['species'] = species
df['log10(dmax^2)'] = np.log10(df['dmax']**2)

if not os.path.exists(datadir):
    os.makedirs(datadir)

df.to_csv('../data/dmax.csv', index=False)

x = df['log10(dmax^2)']

fig, ax = pl.subplots()

for i in [10, 25, 50]:
    ax.hist(
            x,
            bins=i,
            edgecolor='black',
            label='Data Counts: Bins='+str(i)
            )

ax.set_ylabel('Number of Data Points [-]')
ax.set_xlabel(r'$log_{10}(Dmax^{2})$ $[log_{10}(mm^{2})]$')

ax.legend()

fig.tight_layout()

if not os.path.exists(figdir):
    os.makedirs(figdir)

fig.savefig(os.path.join(figdir, 'counts_with_tl_with_tg'))

print('Number of Compositions: '+str(df.shape[0]))
print('Max Tl: '+str(max(df['tl']))+' [K]')
print('Min Tl: '+str(min(df['tl']))+' [K]')
print('Max Tg: '+str(max(df['tg']))+' [K]')
print('Min Tg: '+str(min(df['tg']))+' [K]')
