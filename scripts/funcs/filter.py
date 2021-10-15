from decimal import Decimal
import numpy as np
import sys


def read_cfg(name):
    '''
    Read the traning cfg file from MLIP.

    inputs:
        name = The name of the cfg file.

    outputs:
        frames = The parsed frames.
    '''

    frames = []  # The fame containing all info

    # Conditions data acquisition
    frame_switch = False
    etot_switch = False
    size_switch = False
    stress_switch = False
    atom_switch = False
    cell_switch = False

    # Read data
    with open(name, 'r') as handle:
        for line in handle:

            if 'BEGIN_CFG' in line:
                frame = {}
                frame_switch = True
            elif 'END_CFG' in line:
                frames.append(frame)
                frame_switch = False

            elif 'Energy' in line:
                etot_switch = True
            elif etot_switch is True:
                frame['Energy'] = line.strip()
                etot_switch = False
                atom_switch = False

            elif 'Size' in line:
                size_switch = True
            elif size_switch is True:
                frame['Size'] = line.strip()
                size_switch = False

            elif 'PlusStress' in line:
                frame['PlusStress'] = {}
                stress_switch = True
            elif stress_switch is True:
                vals = line.strip().split()
                frame['PlusStress']['xx'] = vals[0]
                frame['PlusStress']['yy'] = vals[1]
                frame['PlusStress']['zz'] = vals[2]
                frame['PlusStress']['yz'] = vals[3]
                frame['PlusStress']['xz'] = vals[4]
                frame['PlusStress']['xy'] = vals[5]
                stress_switch = False

            elif 'Supercell' in line:
                frame['Supercell'] = {}
                cell_switch = 1
            elif cell_switch == 1:
                vals = line.strip().split()
                frame['Supercell']['xx'] = vals[0]
                frame['Supercell']['xy'] = vals[1]
                frame['Supercell']['xz'] = vals[2]
                cell_switch += 1
            elif cell_switch == 2:
                vals = line.strip().split()
                frame['Supercell']['yx'] = vals[0]
                frame['Supercell']['yy'] = vals[1]
                frame['Supercell']['yz'] = vals[2]
                cell_switch += 1
            elif cell_switch == 3:
                vals = line.strip().split()
                frame['Supercell']['zx'] = vals[0]
                frame['Supercell']['zy'] = vals[1]
                frame['Supercell']['zz'] = vals[2]
                cell_switch = False

            elif 'AtomData' in line:
                line = line.strip().split()[1:]
                frame['atoms'] = {i: [] for i in line}
                atom_switch = True
            elif atom_switch is True:
                vals = line.strip().split()
                frame['atoms']['id'].append(vals[0])
                frame['atoms']['type'].append(vals[1])
                frame['atoms']['cartes_x'].append(vals[2])
                frame['atoms']['cartes_y'].append(vals[3])
                frame['atoms']['cartes_z'].append(vals[4])
                frame['atoms']['fx'].append(vals[5])
                frame['atoms']['fy'].append(vals[6])
                frame['atoms']['fz'].append(vals[7])

            elif 'mindist' in line:
                vals = line.strip().split()[-1]
                frame['mindist'] = vals

    return frames


def remove(frames, n):
    '''
    Flags a configuration if the normalized energy is outside of n*sigma.

    inputs:
        frames = The parsed frames from MLIP.
        n = The multiple for sigma filtering.

    outputs:
        frames = The subset of selected frames.
    '''

    energies = [i['Energy'] for i in frames]
    energies = list(map(Decimal, energies))  # Keep precision
    energies = np.array(energies)

    sizes = [i['Size'] for i in frames]
    sizes = list(map(int, sizes))

    energies = energies/sizes

    mean = np.mean(energies)
    std = np.std(energies)
    cut = n*std

    upper = mean+cut
    lower = mean-cut

    cond = (energies >= lower) & (energies <= upper)

    frames = [i for i, j in zip(frames, cond) if j == True]

    return frames


def write(frames, name):
    '''
    Write the filtered cfg file.

    inputs:
        frames = The frames to write.
        name = The output file name.
    '''

    with open(name, 'w') as handle:

        for i in frames:

            handle.write('BEGIN_CFG\n')

            handle.write(' Size\n')
            handle.write('    {}\n'.format(i['Size']))

            handle.write(' Supercell\n')
            xx = i['Supercell']['xx']
            xy = i['Supercell']['xy']
            xz = i['Supercell']['xz']
            yx = i['Supercell']['yx']
            yy = i['Supercell']['yy']
            yz = i['Supercell']['yz']
            zx = i['Supercell']['zx']
            zy = i['Supercell']['zy']
            zz = i['Supercell']['zz']

            handle.write('        {}      {}      {}\n'.format(xx, xy, xz))
            handle.write('        {}      {}      {}\n'.format(yx, yy, yz))
            handle.write('        {}      {}      {}\n'.format(zx, zy, zz))

            handle.write(' AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n')
            a = i['atoms']['id']
            b = i['atoms']['type']
            c = i['atoms']['cartes_x']
            d = i['atoms']['cartes_y']
            e = i['atoms']['cartes_z']
            f = i['atoms']['fx']
            g = i['atoms']['fy']
            h = i['atoms']['fz']
            for a in zip(a, b, c, d, e, f, g, h):
                handle.write('            {}    {}       {}      {}      {}    {}   {}    {}\n'.format(*a))

            handle.write(' Energy\n')
            handle.write('       {}\n'.format(i['Energy']))

            handle.write(' PlusStress:  xx          yy          zz          yz          xz          xy\n')
            xx = i['PlusStress']['xx']
            yy = i['PlusStress']['yy']
            zz = i['PlusStress']['zz']
            yz = i['PlusStress']['yz']
            xz = i['PlusStress']['xz']
            xy = i['PlusStress']['xy']
            handle.write('       {}   {}   {}    {}     {}    {}\n'.format(xx, yy, zz, yz, xz, xy))

            handle.write(' Feature   EFS_by       VASP\n')

            if 'mindist' in i.keys():
                handle.write(' Feature   mindist      {}\n'.format(i['mindist']))

            handle.write('END_CFG\n')


if __name__ == '__main__':

    print('*'*79)
    print('Filtering')
    frames = read_cfg(sys.argv[1])
    start = len(frames)
    frames = remove(frames, int(sys.argv[2]))
    end = len(frames)
    write(frames, sys.argv[3])

    print('Kept {} from {}'.format(end, start))
