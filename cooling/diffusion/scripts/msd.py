from scipy.stats import linregress
import numpy as np


def msd(frame, data):
    '''
    MSD modifier from ovito.
    '''

    displacement_magnitudes = data.particles['Displacement Magnitude']
    msd = np.sum(displacement_magnitudes**2)/len(displacement_magnitudes)
    data.attributes['msd'] = msd


def self_diffusion(x, y, ddof):
    '''
    Calculate self diffusion from MSD curve.
    inputs:
        x = The time.
        y = The MSD.
        ddof = Degrees of freedom from dimension of space.
    outputs:
        d = Diffusion coefficient.
        m = The slope of the fit.
    '''

    m, b, r, p, err = linregress(x, y)  # Fit linear line
    d = m/ddof  # Divide by degrees of freedom

    return d, m, b
