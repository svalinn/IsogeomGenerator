import sys
import os
import shutil
import warnings
import numpy as np
import math as m


class IsoGeomGen(object):

    def __init__(self, levels=None, data=None, db=None):

        # check type of level input
        if type(levels) is list:
            self.assign_levels(levels)
        elif levels is not None:
            raise RuntimeWarning("Type of levels provided not allowed. " +
                                 "Provide list of floats.")
        else:
            self.levels = None
        self.data = data
        self.db = db

    def assign_levels(self, levels):
        """User defines the contour levels to be used as the isosurfaces.

        Input:
        ------
            levels: list of floats, list of user-defined values to use
                for contour levels
        """
        # make sure values are floats
        levels = [float(i) for i in levels]
        self.levels = sorted(levels)

    def read_levels(self, levelfile):
        """Read level values from a file. One value per line only.

        Input:
        ------
            levelfile: str, relative path to file with level information.
        """
        levels = []
        try:
            f = open(levelfile, 'r')
        except IOError:
            raise RuntimeError("Level file {} does not " +
                               "exist.".format(levelfile))

        levels = []
        lines = f.readlines()
        for line in lines:
            levels.append(line)

        self.assign_levels(levels)
