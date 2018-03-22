import visit as v
import sys
import os
import numpy as np

def _unpack_params(params):
    """Unpack the parameters to define the 3D contour for isosurfaces.
    If a certain parameter is not provided, the default value will be
    set.

    Parameters:
    -----------
        params: dictionary, dict with zero or more of the following keys
            and values:
            'N': (optional) int, number of contours to generate;
                default=10
            'log': (optional) boolean, true to use log scale, false to
                use linear; default is true (log scale)
            'minval': (optional) float, minimum value to use for
                isosurfaces; if not provided, minval=1.e-10 when
                log=True, else minval=0.
            'maxval': (optional) float, maximum value to use for
                isosurfaces; if not provided, to be chosen by VisIt

            example: {'N': 6, 'log': True, 'minval': 1e-5, 'maxval': .5}
            If using all default values, params={}.

    Returns:
    --------
        log: bool, true to use log scale, false for linear; default=True
        minval: float, minimum value to use for contours; default=0
            if log=False, or default=1.e-10 if log=True.
        maxval: float or None, if given in params, maxval is float; if
            not provided by params, set to None to be chosen by VisIt
        N: int, number of contours to use; default=10
    """

    # determine if log scale or linear
    if 'log' in params.keys():
        log = params['log']
    else:
        log = False

    # minimum and maximum isosurface values to use
    if 'minval' not in params.keys():
        if log:
            minval = 1.e-10
        else:
            minval = 0.
    else:
        minval = params['minval']

    # set maximum isosurface value if any
    if 'maxval' in params.keys():
        maxval = params['maxval']
    else:
        maxval = None

    # set number of contour levels
    if 'N' in params.keys():
        N = params['N']
    else:
        N = 10

    return log, minval, maxval, N


def _plot_isosurfs(data, N, log, minval, maxval):
    """Get the isosurface contours according to desired parameters.

    Paramters:
    ----------
        data: string, name of data to use to make contours
            example: "ww_n" for neutron weight windows
        N: int, number of contours to generate, default=10
        log: boolean, true to use log scale, false to use
            inear; default is true (log scale)
        minval: float, minimum value to use for isosurfaces;
            if not provided, to be chosen by VisIt
        maxval: float or None, maximum value to use for isosurfaces;
            if None, to be chosen by VisIt
    """

    v.AddPlot("Contour", data)

    # adjust settings
    att = v.ContourAttributes()

    # log or linear scale
    if log:
        att.SetScaling(1)
    else:
        att.SetScaling(0)

    # minimum isosurface values to use
    att.minFlag = True
    att.min = minval

    # set maximum isosurface value if given
    if maxval:
        att.maxFlag = True
        att.max = maxval

    # number of evenly spaced contours
    att.contourNLevels = N

    # generate contours
    v.SetPlotOptions(att)
    v.DrawPlots()


def _export_isosurf_db(dbname, data):
    """Export all 3D contours as a database. This will save a folder
    with each contour as a different .vtk file.

    Parameters:
    -----------
        dbname: string, name of the database to use for folder and files
        data: string, variable name that is contoured to be exported
    """

    e = v.ExportDBAttributes()
    e.db_type = "VTK"
    e.filename = dbname
    e.variables = data
    v.ExportDatabase(e)


def _get_boundaries(filename):
    """Get boundaries of the domain by opening one file

    Parameters:
    -----------
        filename: string, path to single .vtk file

    Returns:
    --------
        bounds: list of floats, list of bounds in this order:
            [xmin, xmax, ymin, ymax, zmin, zmax]
    """

    f = open(filename)
    end = False
    while not end:
        line = f.readline()
        if line.split(" ")[0] == 'avtOriginalBounds':
            bline = f.readline()
            bound_list = bline.split(" ")
            bound_list.remove("\n")
            bounds = [float(x) for x in bound_list]
            end = True
    f.close()

    return bounds


def _generate_slices(dbname):
    """Exports the vtk files for the intersection of countours with the
    domain bounds.

    Parameters:
    -----------
        dbname: string, name of the database (must correspond to folder
            name and .vtk files within the folder)
    """

    # determine number of contour files in the database
    path, dirs, files = os.walk(dbname).next()
    N = len(files)

    # get domain boundaries
    bounds = _get_boundaries(dbname + "/" + dbname + ".0.vtk")
    xb = bounds[0:2]
    yb = bounds[2:4]
    zb = bounds[4:6]
    bounds = {'x': xb, 'y': yb, 'z': zb}

    # open database
    v.OpenDatabase(dbname + "/" + dbname + ".*.vtk database")

    # load and draw meshes
    v.AddPlot("Mesh", "mesh")
    v.DrawPlots()

    # suppress upcoming error and warning messages from visit (we
    # know that some slices will yield no data)
    v.SuppressMessages(1)

    # create directory to hold slices
    dirname = dbname + '-slices'
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    # for each domain boundary, check if there is an intersection with
    # each contour and export if there is.
    for axis in bounds.keys():

        # slice for each bound on the current axis
        for b in bounds[axis]:

            # create slice
            v.AddOperator("Slice")
            atts = v.SliceAttributes()

            # do not project to 2D to preserve 3D coordinates on slices
            atts.project2d = 0

            # set as arbitrary axis and set normal direction depending
            # on axis (x, y, or z)
            atts.axisType = 3
            if axis == 'x':
                norm = np.array((1, 0, 0))
            elif axis == 'y':
                norm = np.array((0, 1, 0))
            elif axis == 'z':
                norm = np.array((0, 0, 1))

            # determine if min or max bound
            # set surface normal accordingly
            if b == min(bounds[axis]):
                slice_type = 'min'
                atts.normal = tuple(-1*norm)
            else:
                slice_type = 'max'
                atts.normal = tuple(norm)

            # set slice location and draw
            atts.originType = atts.Intercept
            atts.originIntercept = b
            v.SetOperatorOptions(atts)
            v.DrawPlots()

            # try to export for each time step
            for t in range(0, N):

                # set time step
                v.SetTimeSliderState(t)

                # try export
                e = v.ExportDBAttributes()
                e.db_type = "VTK"
                e.dirname = dirname
                e.filename = 's{}-{}{}'.format(t, axis, slice_type)
                e.variables = "ww_n"
                v.ExportDatabase(e)

            # remove slice to be able to re-slice
            v.RemoveLastOperator()

    # unsuppress messages
    v.SuppressMessages(4)


def GenerateIsosurfaceContours(f, data, dbname, params={}):
    """Generates a set of .vtk files that correspond to the following:
        VTK Database called <dbname>: contains one vtk file for each
            3D generated isosurface (# of files <= N).
        VTK Files in folder "<dbname>-slices": if isosurface intersects
            with a domain boundary, a single vtk file is generated of
            the 2D countour of the instersection. Files have the
            following naming scheme:
                s{c}-{a}{m}.vtk, where
                c = countour # corresponding to original 3D isosurface
                    in the <dbname> database
                a = domain boundary axis (either x, y, or z)
                m = "min" or "max" if slice is at the minimum or maximum
                    boundary, respectively.
            There is the possibilty of up to 6 vtk files for each 3D
            isosurface.

    Paramters:
    ----------
        f: string, path to the data file to contour (.vtk)
        data: string, name of data to use to make contours
            example: "ww_n" for neutron weight windows
        dbname: string, name of the database to use for folder and files
        params: (optional) dictionary, with the following optional keys
            and values:
            'N': (optional) int, number of contours to generate;
                default=10
            'log': (optional) boolean, true to use log scale, false to
                use linear; default is true (log scale)
            'minval': (optional) float, minimum value to use for
                isosurfaces; if not provided, minval=1.e-10 when
                log=True, else minval=0.
            'maxval': (optional) float, maximum value to use for
                isosurfaces; if not provided, to be chosen by VisIt

            example: {'N': 6, 'log': True, 'minval': 1e-5, 'maxval': .5}
            If using all default values, no params argument is necessary.
    """

    log, minval, maxval, N = _unpack_params(params)

    v.LaunchNowin()

    v.OpenDatabase(f)
    _plot_isosurfs(data, N, log, minval, maxval)
    _export_isosurf_db(dbname, data)

    v.DeleteAllPlots()

    _generate_slices(dbname)

    v.Close()


def SavePlot3d(f, data, params, name="plot3d"):
    """Save a PNG image of the 3D contour plot.

    Parameters:
    -----------
        f: string, path to the data file to contour (.vtk)
        data: string, name of data to use to make contours
            example: "ww_n" for neutron weight windows
        params: (optional) dictionary, with the following optional keys
            and values:
            'N': (optional) int, number of contours to generate;
                default=10
            'log': (optional) boolean, true to use log scale, false to
                use linear; default is true (log scale)
            'minval': (optional) float, minimum value to use for
                isosurfaces; if not provided, minval=1.e-10 when
                log=True, else minval=0.
            'maxval': (optional) float, maximum value to use for
                isosurfaces; if not provided, to be chosen by VisIt

            example: {'N': 6, 'log': True, 'minval': 1e-5, 'maxval': .5}
            If using all default values, no params argument is necessary.
        name: (optional), string, file name to save as name.png.
            Default="plot3d".
    """

    # unpack parameters
    log, minval, maxval, N = _unpack_params(params)

    v.LaunchNowin()

    v.OpenDatabase(f)

    # plot the 3D contours
    _plot_isosurfs(data, N, log, minval, maxval)

    # save the plot
    s = v.SaveWindowAttributes()
    s.format = s.PNG
    s.fileName = name
    s.width, s.height = 1024,768
    s.screenCapture = 0
    v.SetSaveWindowAttributes(s)
    v.SaveWindow()

    v.Close()


def SavePlot2d(f, data, params, axis, val, name="plot2d"):
    """Saves a PNG image of the 2D slice plot along axis at location
    val.

    Parameters:
    -----------
        f: string, path to the data file to contour (.vtk)
        params: (optional) dictionary, with the following optional keys
            and values:
            'N': (optional) int, number of contours to generate;
                default=10
            'log': (optional) boolean, true to use log scale, false to
                use linear; default is true (log scale)
            'minval': (optional) float, minimum value to use for
                isosurfaces; if not provided, minval=1.e-10 when
                log=True, else minval=0.
            'maxval': (optional) float, maximum value to use for
                isosurfaces; if not provided, to be chosen by VisIt

            example: {'N': 6, 'log': True, 'minval': 1e-5, 'maxval': .5}
            If using all default values, no params argument is necessary.
        axis: string, axis to slice along ('x', 'y', or 'z')
        val: float, value on axis to slice at
        name: (optional), string, file name to save as name.png.
            Default="plot2d".
    """

    # unpack parameters
    log, minval, maxval, N = _unpack_params(params)

    v.LaunchNowin()

    v.OpenDatabase(f)

    # plot 3D in order to slice
    _plot_isosurfs(data, N, log, minval, maxval)

    # add slice
    v.AddOperator("Slice")
    atts = v.SliceAttributes()

    # set axis
    if axis in ['x', 'X']:
        ax = 0
    elif axis in ['y', 'Y']:
        ax = 1
    elif axis in ['y', 'Y']:
        ax = 2
    else:
        print('Axis {} not recognized.'.format(axis))
        return

    atts.axisType = ax
    atts.originIntercept = val

    v.SetOperatorOptions(atts)

    # redraw to include slice operator
    v.DrawPlots()

    # save plot
    s = v.SaveWindowAttributes()
    s.format = s.PNG
    s.fileName = name
    s.width, s.height = 1024,768
    s.screenCapture = 0
    v.SetSaveWindowAttributes(s)
    v.SaveWindow()

    v.RemoveAllOperators()

    v.Close()
