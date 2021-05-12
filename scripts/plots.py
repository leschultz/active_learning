from matplotlib import pyplot as pl
import jsonlines


def parity(data, prop, name):
    '''
    Make a paroody plot.
    '''

    # Parody plot
    label = r'$RMSE/\sigma=$'
    label += '{:.2}'.format(data[r'$RMSE/\sigma$'])
    label += '\n'
    label += r'$RMSE=$'
    label += '{:.2}'.format(data[r'$RMSE$'])
    label += '\n'
    label += r'$MAE=$'
    label += '{:.2}'.format(data[r'$MAE$'])
    label += '\n'
    label += r'$R^{2}=$'
    label += '{:.2}'.format(data[r'$R^{2}$'])

    x = data[prop+'_pred']
    y = data[prop]

    fig, ax = pl.subplots()
    ax.plot(
            x,
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
    limits.append(min(min(y), min(x))-0.25)
    limits.append(max(max(y), max(x))+0.25)

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

    ax.set_xlabel('Predicted '+prop)
    ax.set_ylabel('Actual '+prop)

    fig.tight_layout()
    fig.savefig(name)
    pl.close(fig)

    with jsonlines.open(name+'.jsonl', mode='w') as writer:
        writer.write(data)
