from matplotlib import pyplot as pl
from sklearn import metrics

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


def parity(mets, y, y_pred, name, units):
    '''
    Make a paroody plot.
    '''

    # Parody plot
    label = r'$RMSE/\sigma=$'
    label += '{:.2}'.format(mets[r'$RMSE/\sigma$'])
    label += '\n'
    label += r'$RMSE=$'
    label += '{:.2}'.format(mets[r'$RMSE$'])
    label += '\n'
    label += r'$MAE=$'
    label += '{:.2}'.format(mets[r'$MAE$'])
    label += '\n'
    label += r'$R^{2}=$'
    label += '{:.2}'.format(mets[r'$R^{2}$'])

    fig, ax = pl.subplots()
    ax.plot(
            y_pred,
            y,
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
    limits.append(min(min(y), min(y_pred))-0.25)
    limits.append(max(max(y), max(y_pred))+0.25)

    # Line of best fit
    ax.plot(
            limits,
            limits,
            label=r'$45^{\circ}$ Line',
            color='k',
            linestyle=':',
            zorder=1
            )

    ax.set_aspect('equal')
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.legend()

    ax.set_xlabel('Predicted {} {}'.format(name, units))
    ax.set_ylabel('Actual {} {}'.format(name, units))

    fig.tight_layout()
    fig.savefig(name)
    pl.close(fig)

    data = {}
    data['y_pred'] = y_pred
    data['y'] = y
    data['metrics'] = mets

    jsonfile = '{}.json'.format(name)
    with open(jsonfile, 'w') as handle:
        json.dump(data, handle)


def load(true, pred):
    f, e = parse(true)
    f_pred, e_pred = parse(pred)
    return f, f_pred, e, e_pred


def parse(data):

    coords = []
    energies = []

    cond = False
    with open(data, 'r') as handle:
        for i in handle:
            i = i.strip().split(' ')
            i = [j for j in i if j != '']

            if len(i) == 9:
                header = i[-3:]
            elif len(i) == 8:
                points = list(map(float, i[-3:]))
                for i in points:
                    coords.append(i)
            elif 'Energy' in i:
                cond = True
            elif cond:
                energies.append(float(i[0]))
                cond = False

    return coords, energies


def main():

    true = 'test.cfg'
    pred = 'test_pred.cfg'

    f, f_pred, e, e_pred = load(true, pred)

    fmets = eval_metrics(f, f_pred)
    emets = eval_metrics(e, e_pred)

    parity(fmets, f, f_pred, 'Force', r'[eV/$AA$]')
    parity(emets, e, e_pred, 'Energy', '[eV]')


if __name__ == '__main__':
    main()
