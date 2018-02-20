import visit as v
import sys

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
    print(att)


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

def export_db():
    # save the database as OBJ

    e = v.ExportDBAttributes()
    e.db_type="WavefrontOBJ"
    e.filename="test_export_obj"
    e.variables="ww_n"
    v.ExportDatabase(e)


def main():

    # load the file
    # this should be "expanded_tags.vtk"
    f = sys.argv[1]

    load_vtk(f)
    get_contours("ww_n", N=10, log=True, minval=5.e-9, maxval=0.5)

    save_plots(plot3D=True, plot2D=True, basename="plot", axis='y', val=0.0)
    #save_plots()

    #export_db()


if __name__ == "__main__":
    main()

