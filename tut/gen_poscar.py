from pymatgen.core import Lattice, Structure, Element
from math import ceil

import numpy as np
import sys
import re


def gen_cubic(elements, numbers):
    '''
    Create a cubic structure.

    inputs:
        elements = A list of elements.
        numbers = The number of species for each element.

    outputs:
        structure = The structure generated
    '''

    atoms = sum(numbers)  # The number of atoms

    # Get maximum radii
    rads = []
    for el in elements:
        el = Element(el)

        rad = []

        # Gather applicable radii
        try:
            rad.append(el.atomic_radius)
        except Exception:
            pass
        try:
            rad.append(el.average_anionic_radius)
        except Exception:
            pass
        try:
            rad.append(el.average_cationic_radius)
        except Exception:
            pass
        try:
            rad.append(el.average_ionic_radius)
        except Exception:
            pass
        try:
            rad.append(el.metallic_radius)
        except Exception:
            pass

        rads.append(max(rad))  # Maximum radius

    side = 2.0*max(rads)*atoms**(1/3)

    # Element type for each coordinate
    e = []
    for i, j in zip(elements, numbers):
        e += [i]*j

    # Make starting structure
    n = ceil(atoms**(1/3))
    lattice = Lattice.cubic(side)
    structure = Structure(lattice, [elements[0]], [[0, 0, 0]])
    structure.make_supercell([n, n, n])

    # Truncate and shuffle coordinates
    coords = structure.frac_coords
    np.random.shuffle(coords)
    coords = coords[:atoms]

    # Generate and save structure
    structure = Structure(lattice, e, coords)

    return structure


if __name__ == '__main__':

    comp = sys.argv[1]
    comp = re.split('(\d+)', comp)
    comp = [i for i in comp if i != '']

    elements = comp[0::2]
    numbers = list(map(lambda i: int(i), comp[1::2]))

    structure = gen_cubic(elements, numbers)

    structure.to(filename='POSCAR')
