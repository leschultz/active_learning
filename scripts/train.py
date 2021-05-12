from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core import Structure

from maml.apps.pes import MTPotential
from sklearn import metrics

from multiprocessing import Pool
from functools import partial
from tqdm import tqdm

import pandas as pd
import numpy as np
import random
import plots
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


def eval_metrics(y, y_pred):
    '''
    Evaluate standard prediction metrics.
    '''

    rmse = metrics.mean_squared_error(y, y_pred)**0.5
    rmse_sig = rmse/np.std(y)
    mae = metrics.mean_absolute_error(y, y_pred)
    r2 = metrics.r2_score(y, y_pred)

    results = {}
    results[r'$RMSE$'] = rmse
    results[r'$RMSE/\sigma$'] = rmse_sig
    results[r'$MAE$'] = mae
    results[r'$R^{2}$'] = r2

    return results


def eval_metrics_prop(y, y_pred):
    '''
    Gather metrics per property.
    '''

    df = pd.DataFrame()
    y_norm = y['y_orig']/y['n']
    y_norm_pred = y_pred['y_orig']/y_pred['n']

    df['y'] = y_norm
    df['y_pred'] = y_norm_pred
    df['dtype'] = y['dtype']

    data = {}
    for i, j in df.groupby(['dtype']):

        z = j['y'].values
        z_pred = j['y_pred'].values

        res = eval_metrics(z, z_pred)

        data[i] = eval_metrics(z, z_pred)
        data[i][i] = list(z)
        data[i][i+'_pred'] = list(z_pred)

    return data


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
    data = parallel(gather_step, job)

    return data


def train(data, train_split):

    train_split = int(len(data)*train_split)
    random.shuffle(data)

    structures = [i['structure'] for i in data]
    energies = [i['outputs']['energy'] for i in data]
    forces = [i['outputs']['forces'] for i in data]
    stresses = [i['outputs']['stress'] for i in data]

    # Train split
    train_structures = structures[:train_split]
    train_energies = energies[:train_split]
    train_forces = forces[:train_split]
    train_stresses = stresses[:train_split]

    # Test split
    test_structures = structures[train_split:]
    test_energies = energies[train_split:]
    test_forces = forces[train_split:]
    test_stresses = stresses[train_split:]

    # Train a model for evaluation
    train_model = MTPotential()
    train_model.train(
                      train_structures=train_structures,
                      train_energies=train_energies,
                      train_forces=train_forces,
                      train_stresses=None,
                      #unfitted_mtp='10.mtp',
                      energy_weight=1,
                      force_weight=0.01,
                      stress_weight=0.0,
                      max_dist=5.0,
                      max_iter=10
                      )

    # Model evaluation
    eval_train = train_model.evaluate(
                                      test_structures=train_structures,
                                      test_energies=train_energies,
                                      test_forces=train_forces
                                      )

    y_train, y_train_pred = eval_train

    eval_test = train_model.evaluate(
                                     test_structures=test_structures,
                                     test_energies=test_energies,
                                     test_forces=test_forces
                                     )

    y_test, y_test_pred = eval_test

    evaluation = {}
    evaluation['train'] = eval_metrics_prop(y_train, y_train_pred)
    evaluation['test'] = eval_metrics_prop(y_test, y_test_pred)

    for ikey, ivalue in evaluation.items():
        for jkey, jvalue in ivalue.items():
            name = '{}_{}'.format(ikey, jkey)
            plots.parity(jvalue, jkey, name)

    model = MTPotential()
    model.train(
                train_structures=structures,
                train_energies=energies,
                train_forces=forces,
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

    runs = '/mnt/c/Users/Nerve/Desktop/md_jun'
    match = 'xml'
    train_split = 0.8

    jobs = find(runs, match)
    njobs = len(jobs)
    data = []
    count = 1
    for i in jobs:
        print('{} out of {}: {}'.format(count, njobs, i))
        i = os.path.join(runs, i)
        data += gather_job(i)
        count += 1

    train(data, train_split)
