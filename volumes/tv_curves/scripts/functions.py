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
