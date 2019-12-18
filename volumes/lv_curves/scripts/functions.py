import numpy as np
import os


def create_dir(path):
    '''
    Create a directory if it does not already exist

    inputs:
        path = The directory to create.

    outputs:
        NA
    '''

    if not os.path.exists(path):
        os.makedirs(path)


def nearest(value, x):
    '''
    Find the nearest index in list to a desired value.

    inputs:
        value = Value of interest.
        x = List of data.

    outputs:
        index = The index of the value closet to the value of interest.
    '''

    index = np.abs(x-value).argmin()

    return index
