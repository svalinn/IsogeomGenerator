import visit as v
import sys

v.LaunchNowin()


def get_contours():
    # get the contours from the file
    # to add: ability to choose number of contour levels
    # add plot
    v.AddPlot("Contour", "ww_n")

    # adjust settings
    att = v.ContourAttributes()
    att.SetScaling(1)
    att.minFlag = True
    att.maxFlag = True
    att.min = 5e-9
    att.max = 0.5
    v.SetPlotOptions(att)
    #print(att)

    # draw the plot
    v.DrawPlots()

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
    v.OpenDatabase(f)

    get_contours()

    #save_plot()

    export_db()


if __name__ == "__main__":
    main()

