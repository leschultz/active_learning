from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core import Structure
from maml.apps.pes import MTPotential


if __name__ == '__main__':

    data = 'aimd/2000/vasprun.xml'

    data = Vasprun(data)
    print(data.structures)
    print(dir(data))

    print(data.final_energy)
    train_structures = data.structures
    train_energies = [d['output']['energy'] for d in data]
    train_forces = [d['output']['forces'] for d in data]
    train_stresses = [d['output']['stress'] for d in data]
