import pymatgen as mg
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


def gen_cubic(elements, numbers, coords, a):
    '''
    Generate a cubic supercell.

    inputs:
        elements = List of elements.
        numbers = The number of atoms for each element.
        coords = Starting coordinates.
        a = Cubic lattice paramter.

    outputs:
        structure = Pymatgen structure object.
    '''

    # Assign species for each coordinate
    species = []
    for i in range(len(elements)):
        species += [elements[i]]*numbers[i]

    coords = coords[:len(species)]  # Remove left over coordinates

    lattice = mg.Lattice.cubic(a)  # Generate cubic lattice
    structure = mg.Structure(lattice, species, coords)  # Generate structure

    return structure
