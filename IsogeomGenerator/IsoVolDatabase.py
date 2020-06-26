import sys
import os
import shutil
import warnings
import numpy as np
import math as m
import meshio

import visit as v


class IvDb(object):
    """VisIt Step
    """

    def __init__(self, levels=None):
        # initialize levels list to be populated in various ways
        self.levels = levels
        self.completed = False

    def generate_volumes(self, filename, data,
                         dbname=os.getcwd() + "/tmp",
                         levels=None, levelfile=None, genmode=None, N=None,
                         minN=None, maxN=None):
        """Creates an STL file for each isovolume. N+1 files are
        generated and stored in the dbname folder.

        Input:
        ------
            filename: string, path to vtk file with the mesh
            data: string, name of the data whose values exist on the
                mesh (will be used to generate isocontours and
                isovolumes)
            dbname: (optional), string, name of folder to store created
                surface files. Must be absolute path!
                default: a folder called 'tmp' in the current directory

            level information: One of the following 3 sets of inputs must be
                provided if one of the level assignment methods was not called
                previously (assign_levels(), read_levels(), generate_levels()):
                Option 1: Assign levels
                    levels: list of floats, list of user-defined values to use
                        for contour levels
                Option 2: Read levels from file
                    levelfile: str, relative path to file with level
                        information. Each line of the file should have exactly
                        one float to be used as a level value.
                Option 3: Generate Levels
                    genmode: str, options are 'lin', 'log', or 'ratio'.
                        lin: N linearly spaced values between minN and maxN
                        log: N logarithmically spaced values between minN and
                            maxN
                        ratio: levels that are spaced by a constant ratio N.
                            minN will be used as minimum level value and the
                            maximum level value is less than or equal to maxN.
                    N: int or float, number of levels (int) to generate (lin
                        or log mode); or the ratio (float) to use to separate
                        levels (ratio mode).
                    minN: float, minimum level value
                    maxN: float, maximum level value
        """
        # initialize attributes
        self.data = data
        self.db = dbname

        # gather level information
        if self.levels is None:
            # level information has not been previously assigned
            if levels is not None:
                self.assign_levels(levels)
            elif levelfile is not None:
                self.read_levels(levelfile)
            elif genmode is not None:
                if (N is None) or (minN is None) or (maxN is None):
                    raise RuntimeError("Generating levels requires input " +
                                       "values for N, minN, and maxN")
                else:
                    self.generate_levels(N, minN, maxN, mode=genmode)
            else:
                raise RuntimeError("Information for assigning levels " +
                                   "must be provided.")

        # read data using meshio to get min and max and make sure levels are
        # within data bounds:
        mf = meshio.read(filename)
        mindata = min(mf.cell_data['hexahedron'][self.data])
        maxdata = max(mf.cell_data['hexahedron'][self.data])
        self.arbmin = mindata - 10  # lower than lowest data
        self.arbmax = maxdata + 10  # higher than highest data
        all_levels = list(self.levels)
        for level in all_levels:
            if (level <= mindata) or (level >= maxdata):
                warnings.warn("Level {} is out of data bounds.".format(level))
                self.levels.remove(level)

        if len(self.levels) == 0:
            raise RuntimeError("No data exists within provided levels.")

        # Generate isovolumes using VisIT
        try:
            v.LaunchNowin()
        except:
            pass

        # create volumes
        v.OpenDatabase(filename)
        print("Generating isovolumes...")
        self.__generate_vols()
        print("...Isovolumes files generated!")
        v.CloseComputeEngine()

        # write levels to file in database
        self.__write_levels()

        # set flag for indicating completion
        self.completed = True

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
            levels.append(float(line))

        self.levels = sorted(levels)

    def generate_levels(self, N, minN, maxN, mode='lin'):
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
            self.levels = list(np.linspace(minN, maxN,
                                           num=N, endpoint=True))
        elif mode == 'log':
            base = 10.
            start = m.log(minN, base)
            stop = m.log(maxN, base)
            self.levels = list(np.logspace(start, stop, num=N,
                                           endpoint=True, base=base))
        elif mode == 'ratio':
            # set minN as the minimum value and get all other values until maxN
            tmpmax = 0.
            self.levels = [minN]
            while tmpmax < maxN:
                next_val = self.levels[-1] * float(N)
                if next_val <= maxN:
                    self.levels.append(next_val)
                    tmpmax = next_val
                else:
                    break
        else:
            raise RuntimeError("Level generation mode {} not " +
                               "recognized.".format(mode))

    def __generate_vols(self):
        """Generates the isosurface volumes between the contour levels.
        Data files are exported as STLs and saved in the folder dbname.
        Files will be named based on their index corresponding to their
        level values (0.stl is lowest values).
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

        # plot the pseudocolor data in order to get volumes
        v.AddPlot("Pseudocolor", self.data)
        v.DrawPlots()

        # iterate over all isovolume levels
        for i, l in enumerate(self.levels):
            res = 0
            while res == 0:

                # lower bound
                if i == 0:
                    lbound = self.arbmin
                else:
                    lbound = self.levels[i - 1]

                # upper bound
                ubound = self.levels[i]

                # get volume
                # res = 0 if no level found (should update to next level)
                res, ubound = self.__get_isovol(lbound, ubound, i)

        # get maximum isovolume level
        lbound = self.levels[-1]
        ubound = self.arbmax
        self.__get_isovol(lbound, ubound, i + 1)

        # delete plots
        v.DeleteAllPlots()

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

    def __write_levels(self):
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
