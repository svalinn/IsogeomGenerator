import visit as v
import sys
import os

v.LaunchNowin()

def load_vtk(f):
    """Load a .vtk file into visit that contains weight window data.

    Paramaters:
    -----------
        f: string, path to .vtk file to open
    """
    v.OpenDatabase(f)


def get_contours(data, N=10, log=True, **kwargs):
    """Get the isosurface contours according to desired parameters.

    Paramters:
    ----------
        data: string, name of data to use to make contours
            example: "ww_n" for neutron weight windows
        N: (optional) int, number of contours to generate, default=10
        log: (optional) boolean, true to use log scale, false to use
            inear; default is true (log scale)
        minval: (optional) float, minimum value to use for isosurfaces;
            if not provided, to be chosen by VisIt
        maxval: (optional) float, maximum value to use for isosurfaces;
            if not provided, to be chosen by VisIt
    """

    v.AddPlot("Contour", data)

    # adjust settings
    att = v.ContourAttributes()

    # log or linear scale
    if log:
        att.SetScaling(1)
    else:
        att.SetScaling(0)

    # minimum and maximum surface values to use
    if 'minval' in kwargs:
        att.minFlag = True
        att.min = kwargs['minval']
    if 'maxval' in kwargs:
        att.maxFlag = True
        att.max = kwargs['maxval']

    # number of evenly spaced contours
    att.contourNLevels = N

    v.SetPlotOptions(att)
    v.DrawPlots()


def save_plots(plot3D=True, plot2D=True, basename="plot", **kwargs):
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
                print(atts)

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


def export_complete_db(dbname, data):
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


def export_slices(dbname):
    """Exports the vtk files for the intersection of countours with the
    domain bounds.

    Parameters:
    -----------
        dbname: string, name of the database (must correspond to folder
            name and .vtk files within the folder)
    """

    # open database
    v.OpenDatabase(dbname + "/" + dbname + ".*.vtk database")

    # load meshes
    v.AddPlot("Mesh", "mesh")

    # determine number of contour files in the database
    path, dirs, files = os.walk(dbname).next()
    N = len(files)

    # get domain boundaries


    # for each contour, get intersections with boundaries


def main():

    # load the file
    # this should be "expanded_tags.vtk"
    f = sys.argv[1]
    data = "ww_n"
    dbname = "contours"

    #load_vtk(f)

    # to capture the full problem
    # get_contours(data, N=10, log=True, minval=5.e-9, maxval=0.5)

    # to get clean lines, use these values
    #get_contours(data, N=6, log=True, minval=1.e-6, maxval=0.5)

    #save_plots(plot3D=True, plot2D=True, basename="plot", axis='y', val=0.0)

    #export_complete_db(dbname, data)

    # delete all active plots before moving on to slices
    v.DeleteAllPlots()

    export_slices(dbname)

if __name__ == "__main__":
    main()

