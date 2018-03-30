import visit as v
import sys
import os
import numpy as np
import math as m

class IsoVolumes(object):

    def __init__(self, filename, data='ww_n'):
        pass

    def assign_levels(self, levels):
        self.levels = levels

    def generate_levels(self, N, minN, maxN, log):
        # eventually determine the max/min/N/log automatically from
        # data provided
        if log:
            base = 10.
            start = m.log(minN, base)
            stop = m.log(maxN, base)
            self.levels = np.logspace(start, stop, num=N, endpoint=True, base=base)
        else:
            self.levels = np.linspace(minN, maxN, num=N, endpoint=True)

    def generate_volumes(self):
        pass
