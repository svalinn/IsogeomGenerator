import sys
import os
import shutil
import warnings
import numpy as np
import math as m
import meshio

from isg_gen import IsoGeomGen

import visit as v


class IvDb(IsoGeomGen):
    """VisIt Step
    """

    def __init__(self, levels=None, data=None, db=None):
        # initialize attributes
        super(IvDb, self).__init__(levels, data, db)
        self.completed = False

        # set db:
        if self.db is None:
            self.db = os.getcwd() + "/tmp"

    def generate_vols(self, filename):
        """Generates the isosurface volumes between the level values.
        Data files are exported as STLs and saved in the folder db.
        Files will be named based on their index corresponding to their
        level values (0.stl is lowest values).

        Input:
        ------
            filename: string, path to vtk file with the mesh
        """
        # create folder to store data if it does not already exist
        i = 0
        while os.path.isdir(self.db):
            i += 1
            new_dir = self.db.rstrip("/").rstrip("-{}".format(i - 1)) +\
                "-{}/".format(i)
            warnings.warn("Database {} exists. Using {} " +
                          "instead.".format(self.db, new_dir))
            self.db = new_dir
        os.makedirs(self.db + "/vols/")

        # launch VisIt
        try:
            v.LaunchNowin()
        except:
            pass

        # open file
        v.OpenDatabase(filename)

        # read data using meshio to get min and max and make sure levels are
        # within data bounds:
        mf = meshio.read(filename)
        mindata = min(mf.cell_data['hexahedron'][self.data])
        maxdata = max(mf.cell_data['hexahedron'][self.data])
        arbmin = mindata - 10  # lower than lowest data
        arbmax = maxdata + 10  # higher than highest data
        all_levels = list(self.levels)
        for level in all_levels:
            if (level <= mindata) or (level >= maxdata):
                warnings.warn("Level {} is out of data bounds.".format(level))
                self.levels.remove(level)
        if len(self.levels) == 0:
            raise RuntimeError("No data exists within provided levels.")

        # plot the pseudocolor data in order to get volumes
        v.AddPlot("Pseudocolor", self.data)
        v.DrawPlots()

        # iterate over all isovolume levels
        for i, l in enumerate(self.levels):
            res = 0
            while res == 0:

                # lower bound
                if i == 0:
                    lbound = arbmin
                else:
                    lbound = self.levels[i - 1]

                # upper bound
                ubound = self.levels[i]

                # get volume
                # res = 0 if no level found (should update to next level)
                res, ubound = self.__get_isovol(lbound, ubound, i)

        # get maximum isovolume level
        lbound = self.levels[-1]
        ubound = arbmax
        self.__get_isovol(lbound, ubound, i + 1)

        # delete plots
        v.DeleteAllPlots()

        # close visit
        v.CloseComputeEngine()

    def __get_isovol(self, lbound, ubound, i):
        """Gets the volume selection for isovolume and export just the
        outer surface of the volume as STL.

        Input:
        ------
            lbound: float, lower boundary value for the isovolume
            ubound: float, upper boundary value for the isovolume
            i: int, surface number
        """
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
        e.dirname = self.db + "/vols/"
        e.db_type = "STL"
        e.filename = str(i)
        e.variables = self.data
        export_res = v.ExportDatabase(e)

        # check if exporting was successful or not and adjust values
        if export_res == 0:
            # export not successful because there was no data
            # get new upper bound
            warn_message = "Warning: no data to export between " \
                + "{} and {}.\n".format(lbound, ubound) \
                + "Increasing upper bound to next selected level."
            warnings.warn(warn_message)
            index = self.levels.index(ubound)
            ubound_old = ubound
            ubound = self.levels[index + 1]
            self.levels.remove(ubound_old)
        # delete the operators
        v.RemoveAllOperators()

        return export_res, ubound

    def write_levels(self):
        """Write the final level values used to a file that can be used by
        read_levels().
        """
        # store levels as a string
        level_str = ""
        for level in self.levels:
            level_str += str(level)
            level_str += "\n"

        # write to a file
        filepath = self.db + '/levelfile'
        with open(filepath, "w") as f:
            f.write(level_str)
