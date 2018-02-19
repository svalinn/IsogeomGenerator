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
        att.max = 0.5

    v.SetPlotOptions(att)

    # number of evenly spaced contours - TO DO

    # draw the plot
    #v.DrawPlots()

    print(att)


def save_plot():
    #save the plot
    s = v.SaveWindowAttributes()
    s.format = s.PNG
    s.fileName = "test"
    s.width, s.height = 1024,768
    s.screenCapture = 0
    v.SetSaveWindowAttributes(s)
    name = v.SaveWindow()
    print(name)


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

    #save_plot()

    #export_db()


if __name__ == "__main__":
    main()

