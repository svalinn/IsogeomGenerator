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
        assign_levels(): assigns level values as an attribute
        read_levels(): read level values from a file
    """

    def __init__(self, levels=None, data=None, db=None):
        """Create IsoGoemGen object

        Input:
        ------
            levels: (optional), list of floats, values to use for
                isosurface values
            data: (optional), string, name of data on the mesh
            db: (optional), string, path to database folder with
                isovolume files
        """
        # check type of level input
        if type(levels) is list:
            self.assign_levels(levels)
        elif levels is not None:
            raise RuntimeError("Type of levels provided not allowed. " +
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
                for isosurface values
        """
        # make sure values are floats
        levels = [float(i) for i in levels]
        self.levels = sorted(levels)

    def read_levels(self, levelfile):
        """Read level values from a file. One value per line only.

        Input:
        ------
            levelfile: str, path to file with level information.
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
