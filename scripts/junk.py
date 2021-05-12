from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core import Structure
from maml.apps.pes import MTPotential
from copy import deepcopy

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


def gather_step(step):
    '''
    Gather the energies, forces, stresses, and structure for each MD step.
    '''

    data = {}
    data['structure'] = step['structure'].as_dict()
    data['output'] = {}
    data['output']['energy'] = step['e_wo_entrp']
    data['output']['forces'] = step['forces']
    data['output']['stress'] = step['stress']

    return data


def gather_job(job):
    '''
    Gather data for each AIMD job.
    '''

    job = Vasprun(job)
    data = []
    job = (job.ionic_steps)
    data = parallel(gather_step, job)

    return data


if __name__ == '__main__':

    runs = '/mnt/c/Users/Nerve/Desktop/jun_md'
    data  = []
    for i in os.listdir(runs):
        i = os.path.join(runs, i)
        print('-'*79)
        data += gather_job(i)


    print(data[0])
