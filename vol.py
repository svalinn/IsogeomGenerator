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
    """

    def __init__(self, filename, data):
        self.data = data
        self.f = filename


    def assign_levels(self, levels):
        """User defines the contour levels to be used in the isovolumes.

        Input:
        ------
            levels: list of floats, list of user-defined values to use
                for contour levels
        """
        self.levels = sorted(levels)
        self.minN = min(self.levels)
        self.maxN = max(self.levels)


    def generate_levels(self, N, minN, maxN, log=True):
        """Auto-generate evenly-spaced level values to use given a
        user-defined maximum, minimum, and number of levels to use.

        Input:
        ------
            N: int, number of levels generate
            minN: float, minimum level value
            maxN: float, maximum level value
            log: bool (optional), True to generate evenly spaced levels
                on a log scale (default), False to use linear scale.
        """

        # set min/max values
        self.minN = minN
        self.maxN = maxN

        # generate evenly spaced values
        if log:
            base = 10.
            start = m.log(minN, base)
            stop = m.log(maxN, base)
            self.levels = np.logspace(start, stop, num=N,
                                      endpoint=True, base=base)
        else:
            self.levels = np.linspace(minN, maxN, num=N, endpoint=True)


    def _plot_pseudocolor(self):
        """Plots the data on a pseudocolor plot to use."""

        # add the pseudocolor plot to contour
        v.AddPlot("Pseudocolor", self.data)
        att = v.PseudocolorAttributes()

        # min/max for the pseudocolor plot
        att.minFlag = True
        att.min = self.minN
        att.maxFlag = True
        att.max = self.maxN

        # plot
        v.SetPlotOptions(att)
        v.DrawPlots()


    def _get_isovol(self, lbound, ubound, i):
        """Set selection for isovolume and export as STL."""

        # generate isovolume
        v.AddOperator("Isovolume")
        att = v.IsovolumeAttributes()
        att.lbound = lbound
        att.ubound = ubound
        v.SetOperatorOptions(att)

        # set operator setting to only get surfaces meshes
        v.AddOperator("ExternalSurface")

        # draw plot
        v.DrawPlots()

        # export current volume to folder
        e = v.ExportDBAttributes()
        cwd = os.getcwd()
        e.dirname = cwd + "/" + self.db
        e.db_type = "VTK"
        e.filename = str(i)
        e.variables = self.data
        v.ExportDatabase(e)

        # delete the operators (external surface/isovolume selection)
        v.RemoveAllOperators()


    def generate_volumes(self, dbname="volumes"):
        """Generates the isosurfaces volumes inbetween the contour levels.
        Data files are exported as STLs and saved in the folder dbname.
        Files will be named based on their index corresponding to their
        level values (0.stl is lowest values).

        Input:
        ------
        dbname: string (optional), folder name to store volume data
            files to be generated. Default "volumes" in current working
            directory.
        """

        # set folder to save exported data
        self.db = dbname

        # make sure levels have been set before proceding
        if not self.levels:
            print("No contour levels have been set!")
            print("Please use assign_levels or generate_levels to set levels.")
            return

        # create folder to store data if it does not already exist
        if not os.path.isdir(self.db):
            os.mkdir(self.db)

        # launch VisIt and open database
        v.LaunchNowin()
        v.OpenDatabase(self.f)

        # plot the pseudocolor data inorder to get volumes
        self._plot_pseudocolor()

        # get the minimum isovolume level
        lbound = 0.0
        ubound = self.levels[0]
        self._get_isovol(lbound, ubound, 0)

        # iterate over all isovolume levels
        for l in self.levels[1:]:

            # get index of current level
            i = self.levels.index(l)

            # assign bounds
            lbound = self.levels[i-1]
            ubound = l

            # get volume
            self._get_isovol(lbound, ubound, i)

        # get maximum isovolume level
        lbound = self.levels[-1]
        ubound = 1000.
        self._get_isovol(lbound, ubound, i+1)

        # close everything
        v.DeleteAllPlots()
        v.Close()
