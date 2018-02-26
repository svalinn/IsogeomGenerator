import viso as v

def main():
    f = "../magic/expanded_tags.vtk"
    data = "ww_n"
    dbname = "isosurfaces"
    v.GenerateIsosurfaceContours(f, data, dbname, N=6, log=True, minval=1.e-6, maxval=0.5)


if __name__ == '__main__':
    main()
