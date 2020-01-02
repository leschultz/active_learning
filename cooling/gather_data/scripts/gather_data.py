from os.path import join
import iterators
import functions

# Input paramters
out_print = 'out.txt'  # The VASP print to screen file
incar = 'INCAR'  # INCAR file
poscar = 'POSCAR'  # POSCAR file
outcar = 'OUTCAR'  # OUTCAR file
xdatcar = 'XDATCAR'  # XDATCAR file
oszicar = 'OSZICAR'  # OSZICAR file
vol_dir = '../../run_sets'  # Run directory
fraction = 0.5  # The fraction of hold data to average
data_save_dir = '../data'  # The data save folder
data_save_name = 'data.csv'  # The data save name
save_plots = '../figures'  # The figures save folder

# Paths for run files
incar_paths = functions.finder(incar, vol_dir)
poscar_paths = functions.finder(poscar, vol_dir)

# Paths containing all relevant files
paths = incar_paths.intersection(poscar_paths)

# Create directory to save plots
functions.create_dir(save_plots)

# Iterate for every run
df = iterators.iterate(
                       paths,
                       out_print,
                       incar,
                       poscar,
                       outcar,
                       xdatcar,
                       oszicar,
                       fraction,
                       save_plots,
                       )

# Create directory and save data
functions.create_dir(data_save_dir)
data_save = join(data_save_dir, data_save_name)
df.to_csv(data_save, index=False)

print('Saved: '+data_save)
