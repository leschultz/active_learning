from scipy.interpolate import UnivariateSpline as interpolate
from matplotlib import pyplot as pl

from scipy.stats import linregress
from os.path import join

import numpy as np
import os


def rmse(act, pred):
    '''
    Calculate the root mean squared error (RMSE).

    inputs:
        act = The actual values.
        pred = The predicted values.

    outputs:
        e = The RMSE.
    '''

    act = np.array(act)
    pred = np.array(pred)

    e = (pred-act)**2
    e = e.mean()
    e = np.sqrt(e)

    return e


def opt(x, y):
    '''
    Linearely fit data from both ends of supplied data and calculate RMSE.

    Left fit:
        Start at a pivot (x[0]) and then fit multiple lines from x[n] until
        x[1]. Calculate RMSE for each fit.

    Right fit:
        Start at a pivot (x[n]) and then fit multiple lines from x[0] until
        x[n-1]. Calculate RMSE for each fit.

    After the fits are completed, their errors are averaged for each
    x point included in a fit excluding the pivots.

    The minimum averaged RMSE is considered to be the optimal point for
    fitting two lines.

    inputs:
        x = Horizontal axis data.
        y = Vertical axis data.

    outputs:
        xcut = The chosen x point.
        endpoints = The non-pivot x values.
        middlermse = The average RMSE.
    '''

    # Make data into numpy arrays
    x = np.array(x)
    y = np.array(y)

    n = len(x)  # The length of data

    ldata = []  # Left fits non-pivot and RMSE
    rdata = []  # Right fits non-pivot and RMSE

    for i in list(range(n-1)):

        # Left fit
        xl = x[:n-i]
        yl = y[:n-i]

        ml, il, _, _, _ = linregress(xl, yl)
        yfitl = ml*xl+il
        rmsel = rmse(yl, yfitl)  # Left RMSE

        # Right fit
        xr = x[i:]
        yr = y[i:]

        mr, ir, _, _, _ = linregress(xr, yr)
        yfitr = mr*xr+ir
        rmser = rmse(yr, yfitr)  # Right RMSE

        # Exclude pivot points from data
        if i > 0:
            ldata.append((xl[-1], rmsel))
            rdata.append((xr[0], rmser))

    # Align data based on pivot point
    ldata = np.array(ldata)
    rdata = np.array(rdata[::-1])

    middle_rmse = (ldata[:, 1]+rdata[:, 1])/2  # Mean RMSE
    mcut = np.argmin(middle_rmse)  # Minimum RMSE index
    xcut = ldata[mcut, 0]  # x data with minimum RMSE

    endpoints = ldata[:, 0]  # Non-pivot points for fits

    return xcut, endpoints, middle_rmse
