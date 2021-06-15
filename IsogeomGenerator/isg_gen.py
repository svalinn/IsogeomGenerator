import os


class IsoGeomGen(object):
    """Parent class for common member variables and methods used in
    ivdb.IsoVolDatabase() and isg.IsoSurfGen() classes.

    Attributes:
    -----------
        levels: list of floats, values used for isosurface values
        data: string, name of data on mesh
        db: string, path to database folder with isovolume files

    Methods:
    --------
        read_levels(): read level values from a file or list and set
            attribute
    """

    def __init__(self, levels=None, data=None, db=None, extents=None):
        """Create IsoGoemGen object

        Input:
        ------
            levels: (optional), list of floats or string, values to use
                for isosurface values. If string, then it is the path
                to levelfile where one value is given per line.
            data: (optional), string, name of data on the mesh
            db: (optional), string, path to database folder with
                isovolume files
            extents: (optional) list of list of floats, minimum and
                maximum values for x, y, and z in mesh. Must be
                structured like [[xmin, ymin, zmin], [xmax, ymax, zmax]]
        """
        self.levels = levels
        self.data = data
        self.db = db

        # set extents
        if extents is not None:
            self.xmin = extents[0][0]
            self.ymin = extents[0][1]
            self.zmin = extents[0][2]
            self.xmax = extents[1][0]
            self.ymax = extents[1][1]
            self.zmax = extents[1][2]
        else:
            self.xmin = None
            self.ymin = None
            self.zmin = None
            self.xmax = None
            self.ymax = None
            self.zmax = None

        # set levels
        if self.levels is not None:
            self.read_levels(levels)

        # set db default
        if self.db is None:
            self.db = os.getcwd() + "/tmp"

    def read_levels(self, levels):
        """Read list of levels or file to assign as levels attribute

        Input:
        ------
            levels: list of floats or string, values to use for levels
                attribute. If string, then it is the path to levelfile
                where one value is given per line.
        """
        if type(levels) is list:
            self.__assign_levels(levels)
            return
        elif type(levels) is str:
            try:
                f = open(levels, 'r')
            except IOError:
                raise RuntimeError("Level file {} does not " +
                                   "exist.".format(levels))
            levels_list = []
            lines = f.readlines()
            for line in lines:
                levels_list.append(line)
            self.__assign_levels(levels_list)
            return
        else:
            raise RuntimeError("Type of levels provided not allowed. " +
                               "Provide list of floats or file.")

    def __assign_levels(self, levels):
        """assign list as levels attribute

        Input:
        ------
            levels: list of floats, values to use for levels
        """
        # make sure values are floats
        levels = [float(i) for i in levels]
        self.levels = sorted(levels)
