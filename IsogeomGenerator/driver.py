"""Driver script for the two overarching steps in an isosurface
geometry creation.
"""

import os
import warnings
import numpy as np
import math as m


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

    Returns:
    --------
        levels: list of floats, evenly spaced values according to mode
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


def generate_volumes(ivdb, filename, data=None, db=os.getcwd() + "/tmp",
                     levelinfo=None):
    """Creates an STL file for each isovolume. N+1 files are
    generated and stored in the dbname folder.

    Input:
    ------
        ivdb: IsoVolDatabase.IvDb object to use with database generation
        filename: string, path to vtk file with the mesh
        data: (optional) string, name of the data whose values exist on
            the mesh (will be used to generate isosurfaces and
            isovolumes)
        db: (optional), string, name of folder to store created
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
        ivdb.db = db

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

    # create volumes
    print("Generating isovolumes...")
    ivdb.generate_vols(filename)
    print("...Isovolumes files generated!")

    # write levels to file in database
    ivdb.write_levels()

    # set flag for indicating completion
    ivdb.completed = True


def create_geometry(isogeom, ivdb=None, data=None, dbname=None,
                    levelfile=None, tag_for_viz=False, norm=1.0,
                    merge_tol=1e-5, tags=None, sname=None, sdir=None):
    """Over-arching function to do all steps to create a single
    isosurface geometry for DAGMC using pyMOAB.

    Input:
    ------
        isogeom: IsoGeom object, initialized object to use to create
            full isosurface geometry
        ivdb: IsoVolDatabase object (optional), pass a completed
            IsoVolDatabase object with already set levels and database
            information. If passed, will attempt to get necessary info
            from the object instead of dbname and levelfile.
        data: string (optional), name of the data. Required if no
            ivdb was provided.
        tag_for_viz: bool (optional), True to tag each triangle on
            every surface with the data value. Needed to visualize
            values in VisIt. Default=False.
        norm: float (optional), default=1. All data values will be
            multiplied by the normalization factor.
        merge_tol: float (optional), default=1e-5 cm. Merge tolerance for
            mesh based merge of coincident surfaces. Recommended to be
            1/10th the mesh voxel size.
        dbname: (optional), string, name of folder to store created
            surface files. Must be absolute path.
            default: a folder called 'tmp' in the current directory
        levelfile: str (optional), relative path to file with level
            information. Each line of the file should have exactly
            one float to be used as a level value. If not provided, it
            will be looked for in the database folder (dbname). If levels
            already exist on the ivdb or isogeom object, levelfile will be
            ignored.
        tags: (optional), dict, set of names and values to tag on the
            geometry root set. Dictionary should be structured with each
            key as a tag name (str) and with a single value (float) for the
            tag. Example: {'NAME1':1.0, 'NAME2': 2.0}
        sname: (optional), str, name of file (including extension) for the
            written geometry file. Acceptable file types are VTK and H5M.
            Default name: isogeom.h5m
    """
    if ivdb is not None:
        isogeom.read_isovol(ivdb)

    # set data name
    if data is not None:
        # check if data was already set in the init
        if isogeom.data is not None:
            warnings.warn("New variable data will overwrite " +
                          "previously provided data.")
        isogeom.data = data
    # if still not set, raise error
    if isogeom.data is None:
        raise RuntimeError("Variable 'data' must be provided.")

    # set database info
    if dbname is not None:
        # check if database was already set in init
        if isogeom.db is not None:
            warnings.warn("New variable dbname will overwrite " +
                          "previously provided dbname.")
        isogeom.db = dbname
    # if still not set, use default location
    if isogeom.db is None:
        isogeom.db = os.getcwd() + "/tmp"

    # check that the database folder exists:
    if not os.path.exists(isogeom.db + '/vols/'):
        raise RuntimeError('Database {} does not '.format(isogeom.db) +
                           'contain an isovolume database.')

    # get level information from file
    if levelfile is None:
        levelfile = isogeom.db + '/levelfile'
    if isogeom.levels is None:
        isogeom.read_levels(levelfile)
    else:
        warnings.warn("levels already set, ignoring levelfile.")

    print("Reading database...")
    isogeom.read_database()
    print("... Reading complete!")

    # Step 1: Separate Isovolume Surfaces
    print("Separating isovolumes...")
    isogeom.separate_isovols()
    print("...Separation complete!")

    # Step 2: Merge Coincident Surfaces
    print("Merging surfaces...")
    isogeom.imprint_merge(norm, merge_tol)
    print("...Merging complete!")

    # Step 3: Assign Parent-Child Relationship
    isogeom.make_family()

    if tag_for_viz:
        print('Tagging triangles with data...')
        isogeom.tag_for_viz()
        print('... tags complete')

    if tags is not None:
        isogeom.set_tags(tags)

    if sdir is None:
        sdir = isogeom.db
    if sname is None:
        sname = 'isogeom.h5m'

    isogeom.write_geometry(sname, sdir)
