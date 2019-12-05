from ast import literal_eval
import numpy as np
import functions
import sys
import os

sys.argv[1:] = [literal_eval(i) for i in sys.argv[1:]]  # Convert types

# Input paramters
a = sys.argv[1]  # Cubic lattice constant
elements = sys.argv[2]  # Elements
numbers = sys.argv[3]  # Number corresponding to each element
coords = sys.argv[4]  # Starting coordinates
save = sys.argv[5]  # Location to save POSCAR
runs = sys.argv[6]  # The number of runs to generate

coords = np.loadtxt(coords)  # Load starting coordinates
coords = coords/coords.max(axis=0)  # Make fractional

# Generates structure
structure = functions.gen_cubic(
                                elements,
                                numbers,
                                coords,
                                a,
                                )

# Random seeds
seeds = np.random.randint(low=100000, high=999999, size=runs)

for seed in seeds:
    run = os.path.join(save, str(seed))
    functions.create_dir(run)
    structure.to(fmt='poscar', filename=os.path.join(run, 'POSCAR'))  # Save POSCAR

    print('Generated: '+run)
