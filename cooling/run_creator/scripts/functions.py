from pymatgen import Lattice, Structure
from math import ceil
import numpy as np
import os


def create_dir(path):
    '''
    Create a directory if it does not already exist

    inputs:
        path = The directory to create.
    '''

    if not os.path.exists(path):
        os.makedirs(path)


def finder(name, source):
    '''
    Find the diretories with a file.

    inputs:
        name = The generic name for files to get path of.
        source = The parent directory of all files.

    outputs:
        paths = The matching paths.
    '''

    # Count all mathching paths
    paths = []
    for item in os.walk(source):

        if name not in item[2]:
            continue

        paths.append(item[0])

    paths = set(paths)  # Because ordoer should not matter

    return paths


def gen_cubic(elements, numbers, l):
    '''
    Create a cubic structure.

    inputs:
        elements = A list of elements.
        numbers = The number of species for each element.

    outputs:
        l = The length for the sides of the cubic volume.
    '''

    atoms = sum(numbers)  # The number of atoms

    # Element type for each coordinate
    e = []
    for i, j in zip(elements, numbers):
        e += [i]*j

    # Make starting structure
    n = ceil(atoms**(1/3))
    lattice = Lattice.cubic(l)
    structure = Structure(lattice, [elements[0]], [[0, 0, 0]])
    structure.make_supercell([n, n, n])

    # Truncate and shuffle coordinates
    coords = structure.frac_coords[:atoms]
    np.random.shuffle(coords)

    # Generate and save structure
    structure = Structure(lattice, e, coords)

    return structure
