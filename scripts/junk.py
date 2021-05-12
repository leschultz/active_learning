from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core import Structure

from maml.apps.pes import MTPotential

from multiprocessing import Pool
from functools import partial
from tqdm import tqdm

import sys
import os


def parallel(func, x, *args, **kwargs):
    '''
    Run some function in parallel.
    '''

    pool = Pool(os.cpu_count())
    part_func = partial(func, *args, **kwargs)

    data = list(tqdm(pool.imap(part_func, x), total=len(x), file=sys.stdout))
    pool.close()
    pool.join()

    return data


def find(loc, match):
    '''
    Find all xml files.
    '''

    paths = []
    for i in os.walk(loc):
        for j in i[2]:
            if j.split('.')[-1] == match:
                paths.append(os.path.join(i[0], j))

    return paths


def gather_step(step):
    '''
    Gather the energies, forces, stresses, and structure for each MD step.
    '''

    # The stress as xx, yy, zz, yz, xz, and xy
    stress = step['stress']
    stress = [
              stress[0][0],
              stress[1][1],
              stress[2][2],
              stress[1][2],
              stress[0][2],
              stress[0][1]
              ]

    data = {}
    data['structure'] = step['structure']
    data['outputs'] = {}
    data['outputs']['energy'] = step['e_wo_entrp']
    data['outputs']['forces'] = step['forces']
    data['outputs']['stress'] = stress

    return data


def gather_job(job):
    '''
    Gather data for each AIMD job.
    '''

    job = Vasprun(job)
    data = []
    job = (job.ionic_steps)
    data = parallel(gather_step, job[:2])

    return data


def train(data):

    train_structures = [i['structure'] for i in data]
    train_energies = [i['outputs']['energy'] for i in data]
    train_forces = [i['outputs']['forces'] for i in data]
    train_stresses = [i['outputs']['stress'] for i in data]

    model = MTPotential()
    model.train(
                train_structures=train_structures,
                train_energies=train_energies,
                train_forces=train_forces,
                train_stresses=None,
                #unfitted_mtp='10.mtp',
                energy_weight=1,
                force_weight=0.01,
                stress_weight=0.0,
                max_dist=5.0,
                max_iter=1000
                )

    

    return model


if __name__ == '__main__':

    runs = '/mnt/c/Users/Nerve/Desktop/jun_md'
    match = 'xml'

    jobs = find(runs, match)[:2]
    njobs = len(jobs)
    data = []
    count = 1
    for i in jobs:
        print('{} out of {}: {}'.format(count, njobs, i))
        i = os.path.join(runs, i)
        data += gather_job(i)
        count += 1


    train(data)
