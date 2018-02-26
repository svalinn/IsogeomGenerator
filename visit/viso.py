import visit as v
import sys
import os



def _load_vtk(f):
    """Load a .vtk file into visit that contains weight window data.

    Paramaters:
    -----------
        f: string, path to .vtk file to open
    """
    v.OpenDatabase(f)


def _get_isosurfs(data, N, log, minval, maxval):
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
        print(maxval)
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


def _export_slices(dbname):
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

            # determine axis
            if axis == 'x':
                atts.axisType = 0
            elif axis == 'y':
                atts.axisType = 1
            elif axis == 'z':
                atts.axisType = 2

            # determine if min or max bound
            # need to adjust bounds by 1e-14 to ensure inclusivity
            # within visit slicing program
            if b == min(bounds[axis]):
                slice_type = 'min'
                b += 1.e-14
            else:
                slice_type = 'max'
                b -= 1.e-14

            # set slice location and draw
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


def SavePlot2d():
    pass


def SavePlot3d(plot3D=True, plot2D=True, basename="plot", **kwargs):
    """Save a PNG image of the plot if desired. 3D contour and 2D slice
    of contour plot possible.

    Parameters:
    -----------
        plot3D: (optional) bool, True to save as a 3D image of the contours;
            default=True
        plot2D: (optional) bool, True to save as a 2D image of the contours;
            default=True
        basename: (optional), string, base name to use for file names.
            Default="plot". Files saved as plot-2d.png and plot-3d.png
            for 2D and 3D plots, respectively.
        axis: string, required if 2D plot, axis to slice along
        val: float, required if 2D plot, value to slice at
    """

    # plot 3D first
    if plot3D:
        v.DrawPlots()
        s = v.SaveWindowAttributes()
        s.format = s.PNG
        s.fileName = basename + "-3d"
        s.width, s.height = 1024,768
        s.screenCapture = 0
        v.SetSaveWindowAttributes(s)
        v.SaveWindow()

    # plot 2D next
    if plot2D:
        # check that appropriate kwargs supplied
        # if not supplied, return no plot
        if ('axis' and 'val') not in kwargs:
            print("Must supply slicing axis and location for 2D plot.")
            return
        else:
            axis = kwargs['axis']
            # check that "axis" is either x, y, or z
            if axis not in ['x', 'y', 'z', 'X', 'Y', 'Z']:
                print("Axis {} is invalid.".format(axis))
                return
            else:
                v.AddOperator("Slice")
                atts = v.SliceAttributes()

                # set axis
                if axis in ['x', 'X']:
                    ax = 0
                elif axis in ['y', 'Y']:
                    ax = 1
                elif axis in ['y', 'Y']:
                    ax = 2

                atts.axisType = ax
                val = kwargs['val']
                atts.originIntercept = val

                v.SetOperatorOptions(atts)

                # redraw to include slice operator
                v.DrawPlots()

                # save plot
                s = v.SaveWindowAttributes()
                s.format = s.PNG
                s.fileName = basename + "-2d"
                s.width, s.height = 1024,768
                s.screenCapture = 0
                v.SetSaveWindowAttributes(s)
                v.SaveWindow()

                # remove the slice to continue working
                v.RemoveLastOperator()

                # redraw plot to get back to 3D
                v.DrawPlots()


def GenerateIsosurfaceContours(f, data, dbname, N=10, log=True, **kwargs):
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
        N: (optional) int, number of contours to generate, default=10
        log: (optional) boolean, true to use log scale, false to use
            linear; default is true (log scale)
        minval: (optional) float, minimum value to use for isosurfaces;
            if not provided, minval=1.e-10 when log=True, else minval=0.
        maxval: (optional) float, maximum value to use for isosurfaces;
            if not provided, to be chosen by VisIt
    """

    # minimum and maximum isosurface values to use
    if 'minval' not in kwargs:
        if log:
            minval = 1.e-10
        else:
            minval = 0.
    else:
        minval = kwargs['minval']

    # set maximum isosurface value if any
    if 'maxval' in kwargs:
        maxval = kwargs['maxval']
    else:
        maxval = None

    v.LaunchNowin()

    _load_vtk(f)
    _get_isosurfs(data, N, log, minval, maxval=maxval)
    _export_isosurf_db(dbname, data)

    v.DeleteAllPlots()

    _export_slices(dbname)


