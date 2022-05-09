from matplotlib import pyplot as pl
from sklearn import metrics

import pandas as pd
import numpy as np
import json


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


def parity(mets, y, y_pred, name, units, save):
    '''
    Make a paroody plot.
    '''

    # Parody plot
    label = r''
    for key, value in mets.items():
        label += r'{}={:.2}'.format(key, value)
        label += '\n'
    label = label[:-2]  # Compensate for last break

    fig, ax = pl.subplots()
    ax.plot(
            y,
            y_pred,
            linestyle='none',
            marker='.',
            zorder=0,
            label='Data'
            )

    ax.text(
            0.65,
            0.05,
            label,
            transform=ax.transAxes,
            bbox=dict(facecolor='none', edgecolor='black')
            )

    limits = []
    min_range = min(min(y), min(y_pred))
    max_range = max(max(y), max(y_pred))
    span = max_range-min_range
    limits.append(min_range-0.1*span)
    limits.append(max_range+0.1*span)

    # Line of best fit
    ax.plot(
            limits,
            limits,
            label=r'$y=\hat{y}$ Line',
            color='k',
            linestyle=':',
            zorder=1
            )

    ax.set_aspect('equal')
    ax.set_ylim(limits)
    ax.set_xlim(limits)
    ax.legend()

    ax.set_ylabel('Predicted {} {}'.format(name, units))
    ax.set_xlabel('True {} {}'.format(name, units))

    fig.tight_layout()
    fig.savefig(save)
    pl.close(fig)

    data = {}
    data['y_pred'] = y_pred
    data['y'] = y
    data['metrics'] = mets

    jsonfile = '{}.json'.format(save)
    with open(jsonfile, 'w') as handle:
        json.dump(data, handle)


def load(true, pred):
    f, e, s = parse(true)
    f_pred, e_pred, s_pred = parse(pred)
    return f, f_pred, e, e_pred, s, s_pred


def parse(data):

    forces = []
    energies = []
    stresses = []

    econd = False
    scond = False
    with open(data, 'r') as handle:
        for i in handle:
            i = i.strip().split(' ')
            i = [j for j in i if j != '']

            if len(i) == 9:
                header = i[-3:]
                atoms = 0
            elif len(i) == 8:
                points = list(map(float, i[-3:]))
                element = int(i[1])
                forces.append({
                               'element': element,
                               'fx': points[0],
                               'fy': points[1],
                               'fz': points[2],
                               })
                atoms += 1

            elif 'Energy' in i:
                econd = True
            elif econd:
                energies.append(float(i[0])/atoms)
                econd = False
            elif 'PlusStress:' in i:
                scond = True
            elif scond:
                points = list(map(float, i))
                for i in points:
                    stresses.append(i)
                scond = False

    return forces, energies, stresses


def main(true, pred, save):

    f, f_pred, e, e_pred, s, s_pred = load(true, pred)

    # Prepare data for per_elemnt analysis
    dff = pd.DataFrame(f)
    dff_pred = pd.DataFrame(f_pred)
    dff['set'] = 'test'
    dff_pred['set'] = 'pred'
    df = pd.concat([dff, dff_pred])

    for i, j in df.groupby('element'):
        k = j.loc[j['set'] == 'test']
        l = j.loc[j['set'] == 'pred']
        elfx = k['fx'].tolist()
        elfy = k['fy'].tolist()
        elfz = k['fz'].tolist()
        elf = elfx+elfy+elfz

        elfx_pred = l['fx'].tolist()
        elfy_pred = l['fy'].tolist()
        elfz_pred = l['fz'].tolist()
        elf_pred = elfx_pred+elfy_pred+elfz_pred

        fmets = eval_metrics(elf, elf_pred)
        fxmets = eval_metrics(elfx, elfx_pred)
        fymets = eval_metrics(elfy, elfy_pred)
        fzmets = eval_metrics(elfz, elfz_pred)

        parity(fmets, elf, elf_pred, 'Force', r'[eV/$\AA$]', save+'_Force'+str(i))
        parity(fxmets, elfx, elfx_pred, 'Force', r'[eV/$\AA$]', save+'_xForce'+str(i))
        parity(fymets, elfy, elfy_pred, 'Force', r'[eV/$\AA$]', save+'_yForce'+str(i))
        parity(fzmets, elfz, elfz_pred, 'Force', r'[eV/$\AA$]', save+'_zForce'+str(i))

    fx = [i['fx'] for i in f]
    fy = [i['fy'] for i in f]
    fz = [i['fz'] for i in f]
    f = fx+fy+fz

    fx_pred = [i['fx'] for i in f_pred]
    fy_pred = [i['fy'] for i in f_pred]
    fz_pred = [i['fz'] for i in f_pred]
    f_pred = fx_pred+fy_pred+fz_pred

    fmets = eval_metrics(f, f_pred)
    fxmets = eval_metrics(fx, fx_pred)
    fymets = eval_metrics(fy, fy_pred)
    fzmets = eval_metrics(fz, fz_pred)
    emets = eval_metrics(e, e_pred)
    smets = eval_metrics(s, s_pred)

    parity(fxmets, fx, fx_pred, 'Force', r'[eV/$\AA$]', save+'_xForce')
    parity(fymets, fy, fy_pred, 'Force', r'[eV/$\AA$]', save+'_yForce')
    parity(fzmets, fz, fz_pred, 'Force', r'[eV/$\AA$]', save+'_zForce')
    parity(fmets, f, f_pred, 'Force', r'[eV/$\AA$]', save+'_Force')
    parity(emets, e, e_pred, 'Energy', '[eV/atom]', save+'_Energy')
    parity(
           smets,
           s,
           s_pred,
           'Volume Times Virial Stress',
           '[eV]',
           save+'_Stress'
           )


if __name__ == '__main__':
    main('test.cfg', 'test_pred.cfg', 'test')
    main('train.cfg', 'train_pred.cfg', 'train')
