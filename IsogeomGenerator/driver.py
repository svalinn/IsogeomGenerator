"""Driver script for the two overarching steps in an isosurface
geometry creation.
"""

import sys
import os
import shutil
import warnings
import numpy as np
import math as m
import meshio

import visit as v


def generate_levels(N, minN, maxN, mode='lin'):
    """Auto-generate evenly-spaced level values between the min and max
    value.

    Input:
    ------
        N: int or float, number of levels (int) to generate
            (lin or log mode); or the ratio (float) to use to separate
            levels (ratio mode).
        minN: float, minimum level value
        maxN: float, maximum level value
        mode: str, options are 'lin' (default), 'log', or 'ratio'.
            lin: N linearly spaced values between minN and maxN
            log: N logarithmically spaced values between minN and maxN
            ratio: levels that are spaced by a constant ratio N.
                minN will be used as minimum level value and the maximum
                level value is less than or equal to maxN.
    """
    if mode == 'lin':
        levels = list(np.linspace(minN, maxN,
                                  num=N, endpoint=True))
        return levels

    if mode == 'log':
        base = 10.
        start = m.log(minN, base)
        stop = m.log(maxN, base)
        levels = list(np.logspace(start, stop, num=N,
                                  endpoint=True, base=base))
        return levels

    if mode == 'ratio':
        # set minN as the minimum and get all other values until maxN
        tmpmax = 0.
        levels = [minN]
        while tmpmax < maxN:
            next_val = levels[-1] * float(N)
            if next_val <= maxN:
                levels.append(next_val)
                tmpmax = next_val
            else:
                break
        return levels

    raise RuntimeError("Level generation mode {} not " +
                       "recognized.".format(mode))


def generate_volumes(ivdb, filename, data,
                     dbname=os.getcwd() + "/tmp",
                     levelinfo=None):
    """Creates an STL file for each isovolume. N+1 files are
    generated and stored in the dbname folder.

    Input:
    ------
        ivdb: IsoVolDatabase.IvDb object to use with database generation
        filename: string, path to vtk file with the mesh
        data: string, name of the data whose values exist on the
            mesh (will be used to generate isocontours and
            isovolumes)
        dbname: (optional), string, name of folder to store created
            surface files. Must be absolute path!
            default: a folder called 'tmp' in the current directory
        levelinfo: string or list of floats, level value information.
            string: relative path to file with level information. Each
                line of the file should have exactly one float to be
                used as a level value.
            list: list of user-defined values to use for contour levels
    """
    # initialize attributes
    if data is not None:
        # overwrite previously set variable
        ivdb.data = data
    if ivdb.data is None:
        raise RuntimeError("Name of 'data' not provided.")

    if ivdb.db is None:
        ivdb.db = dbname

    # gather level information
    if ivdb.levels is None:
        # level information has not been previously assigned
        # check type of input
        if type(levelinfo) is list:
            ivdb.assign_levels(levelinfo)
        elif type(levelinfo) is str:
            ivdb.read_levels(levelinfo)
        elif type(levelinfo) is None:
            raise RuntimeError("Information for assigning levels " +
                               "must be provided (levelinfo).")
        else:
            raise RuntimeError("Type of levelinfo is not recognized.")

    # Generate isovolumes using VisIT
    try:
        v.LaunchNowin()
    except:
        pass

    # create volumes
    v.OpenDatabase(filename)
    print("Generating isovolumes...")
    ivdb.generate_vols()
    print("...Isovolumes files generated!")
    v.CloseComputeEngine()

    # write levels to file in database
    ivdb.write_levels()

    # set flag for indicating completion
    ivdb.completed = True
