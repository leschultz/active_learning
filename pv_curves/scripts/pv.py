from os.path import join
import iterators
import functions

# Input paramters
incar = 'INCAR'  # INCAR file
poscar = 'POSCAR'  # POSCAR file
outcar = 'OUTCAR'  # OUTCAR file
vol_dir = '../../volumes'  # Run directory
fraction = 0.5  # The fraction of hold data to average
save_dir = '../data'  # The save folder
save_name = 'data.csv'  # The save name

# Paths for run files
incar_paths = functions.finder(incar, vol_dir)
poscar_paths = functions.finder(poscar, vol_dir)

# Paths containing all relevant files
paths = incar_paths.intersection(poscar_paths)

# Iterate for every run
df = iterators.iterate(
                       paths,
                       incar,
                       poscar,
                       outcar,
                       fraction,
                       )

functions.create_dir(save_dir)
df.to_csv(join(save_dir, save_name), index=False)
