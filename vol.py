import visit as v
import sys
import os
import numpy as np
import math as m

class IsoVolumes(object):
    """This class contains methods to create a set of isovolumes given
    a cartesian mesh with data.
    Parameters:
    -----------
        filename: string, path to vtk file with the mesh
        data: string, name of the data to use from the file
        log: bool (optional), True (default) to use log scaling, False
            to use linear (applies to entire problem)
        dbname: string (optional), folder name to store volume data
            files to be generated. Default "volumes" in current working
            directory.
    """

    def __init__(self, filename, data, log=True, dbname="volumes"):
        self.data = data
        self.f = filename
        self.log = log
        self.db = dbname

    def assign_levels(self, levels):
        self.levels = sorted(levels)
        self.minN = min(self.levels)
        self.maxN = max(self.levels)

    def generate_levels(self, N, minN, maxN):
        # eventually determine the max/min/N/log automatically from
        # data provided

        self.minN = minN
        self.maxN = maxN

        if self.log:
            base = 10.
            start = m.log(minN, base)
            stop = m.log(maxN, base)
            self.levels = np.logspace(start, stop, num=N, endpoint=True, base=base)
        else:
            self.levels = np.linspace(minN, maxN, num=N, endpoint=True)

    def _plot_pseudocolor(self):
        # add the pseudocolor plot to contour
        v.AddPlot("Pseudocolor", self.data)
        att = v.PseudocolorAttributes()
        if self.log:
            att.SetScaling(1)
        else:
            att.SetScaling(0)

        # min/max for the pseudocolor plot
        att.minFlag = True
        att.min = self.minN
        att.maxFlag = True
        att.max = self.maxN

        # plot
        v.SetPlotOptions(att)
        v.DrawPlots()

    def _get_isovol(self, lbound, ubound, i):

        # generate
        v.AddOperator("Isovolume")
        att = v.IsovolumeAttributes()
        att.lbound = lbound
        att.ubound = ubound
        v.SetOperatorOptions(att)
        v.DrawPlots()

        # export to dbname folder
        e = v.ExportDBAttributes()
        cwd = os.getcwd()
        e.dirname = cwd + "/" + self.db
        e.db_type = "STL"
        e.filename = str(i)
        e.variables = self.data
        v.ExportDatabase(e)

        # delete the operator
        v.RemoveLastOperator()

    def generate_volumes(self):

        # create folder to store data

        if not os.path.isdir(self.db):
            os.mkdir(self.db)

        # launch database in VisIt
        v.LaunchNowin()
        v.OpenDatabase(self.f)

        # plot the pseudocolor plot
        self._plot_pseudocolor()

        # get the minimum level
        lbound = 0.0
        ubound = self.levels[0]
        self._get_isovol(lbound, ubound, 0)

        for l in self.levels[1:]:

            # get index of current level
            i = self.levels.index(l)

            # assign
            lbound = self.levels[i-1]
            ubound = l

            self._get_isovol(lbound, ubound, i)

        # get maximum level
        lbound = self.levels[-1]
        ubound = 1000.
        self._get_isovol(lbound, ubound, i+1)

        # close everything
        v.DeleteAllPlots()
        v.Close()

