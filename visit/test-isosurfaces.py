import viso as v

def main():
    f = "../magic/expanded_tags.vtk"
    data = "ww_n"
    dbname = "isosurfaces"
    params = {'N': 6, 'log': True, 'minval':1.e-6, 'maxval':0.5}

    v.GenerateIsosurfaceContours(f, data, dbname, params)
    v.SavePlot3d(f, data, params, name="plot3d")
    axis='x'
    val=0
    v.SavePlot2d(f, data, params, axis, val, name="plot2d")


if __name__ == '__main__':
    main()
