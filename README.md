# Isosurface Geometry Generator

This tool can be used to generate an isosurface geometry from any 3D mesh
tagged with scalar values (eg, a Cartesian weight window mesh used for
Monte Carlo particle transport).
The tool will create a meshed surface geometry (DAGMC-compliant) using specified isosurface
values.

## Dependencies:

* Python 2.7
* [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/)
* [MOAB](https://sigma.mcs.anl.gov/moab-library/) v5.1+ with PyMOAB enabled

## Installation:

From the source directory, run the following to install:

      pip install . --user

---

## Python Module Usage

To use, import the Python module with `from IsogeomGenerator import vol`.

### Steps:

To generate a full geometry from a starting Cartesian mesh file, run the following steps in order:

1. **Set contour levels values:** This will assign the values that will be used for the isosurface selections
    in the mesh file. The values can be set using one of three methods:

    * `assign_levels(levels)`: User defines the contour levels to be used as the isosurfaces.

        Input: `levels`, list of floats, list of user-defined values

    * `read_levels(levelfile)`: Read level values from a file. One value per line only.

        Input: `levelfile`, string, relative path to file with level information.

    * `generate_levels(N, minN, maxN, mode='lin')`: Auto-generate evenly-spaced level values between the min and max value.

        Input:
        * `N`: int or float, number of levels (int) to generate
            (`'lin'` or `'log'` mode); or the ratio (float) to use to separate
            levels (`'ratio'` mode).
        * `minN`: float, minimum level value
        * `maxN`: float, maximum level value
        * `mode` (optional): string, options are `'lin'` (default), `'log'`, or `'ratio'`.
            * `'lin'`: `N` linearly spaced values between `minN` and `maxN`
            * `'log'`: `N` logarithmically spaced values between `minN` and `maxN`
            * `ratio`: levels that are spaced by a constant ratio `N`.
                `minN` will be used as minimum level value and the maximum
                level value is less than or equal to `maxN`.

2. **Generate the isovolume database:** This step will use the level values assigned
    in step 1 to create isovolumes in the mesh file using VisIt.
    Two consecutive level values are used to generate an isovolume in the meshfile
    and the surface of that volume is exported as a meshed surface. Each isovolume
    is exported as one file in the database.

    * `generate_volumes(self, filename, data, dbname='/tmp/')`: Creates an STL file for each isovolume. Files are generated in the `dbname` folder.

        Input:
        * `filename`: string, path to vtk file with the mesh file
        * `data`: string, name of the data whose values exist on the mesh
        (will be used to generate the isovolumes)
        * `dbname`: (optional), string, Absolute path to the folder to store created
        surface files. Default: a folder called `tmp/` in the current directory.

3. **Create the DAGMC isosurface geometry:**

    * `create_geometry(tag_for_viz=False, norm=1.0, merge_tol=1e-5, dbname='/tmp/', tags=None, sname=None, sdir=None)`:
    Creates a DAGMC-compliant isosurface geometry from the isovolume files in the database created in step 2 using MOAB.

        Input:
        * `tag_for_viz`: bool (optional), True to tag each triangle on
            every surface with the data value. Needed to visualize
            values in VisIt. Default=False.
        * `norm`: float (optional), default=1. All data values will be
            multiplied by the normalization factor.
        * `merge_tol`: float (option), default=1e-5 cm. Merge tolerance for
            mesh based merge of coincident surfaces.
        * `dbname`: string (optional), name of folder to store created
            surface files. Must be absolute path.
            default: a folder called 'tmp' in the current directory
        * `tags`: dict (optional), set of names and values to tag on the
            geometry root set. Dictionary should be structured with each
            key as a tag name (str) and with a single value (float) for the
            tag. Example: `{'NAME1':1.0, 'NAME2': 2.0}`
        * `sname`: string (optional), name of file (including extension) for the
            written geometry file. Acceptable file types are VTK and H5M.
            Default name: `isogeom.h5m`

----

## Command Line Tool

The steps for creating an isosurface geometry can be done on the command line
with the `generate_isogeom` command. This tool can be run in three different
modes:

* `full`: this will run both the visit step then the moab step (described below). Command:

      generate_isogeom full <meshfile> <dataname> [options]

* `visit`: starting from a Cartesian mesh file, this will generate only a database of
separate isosurface files from VisIt (step 2 above). Command:

      generate_isogeom visit <meshfile> <dataname> [options]

* `moab`: this will start from the database (generated in the visit step) to create
a DAGMC-compliant isosurface geometry using PyMOAB (step 3 above).
Command:

      generate_isogeom moab [options]

To view the three different modes, run `generate_isogeom --help`.

### Options
Each mode has different required or optional arguments needed at run time.
Below is a table summary of those options. This information is also available
via terminal help message: (`generate_isogeom [mode] --help`).

* `X` = required
* `O` = optional
* `-` = Not allowed

|_**Command Option**_ | | | | _**Run Mode**_ | | |
|-|-|-|-|-|----------|-|
| | **Option** | **Description** | **Default value** | `full` | `visit` | `moab` |
| *Mesh file information* | | | | | | |
| Cartesian Mesh File |`meshfile` | Relative path to the Cartesian mesh file that will be used to generate isosurfaces. | | `X` | `X` | `-` |
| Data Name |`dataname` | The name of the scalar data on the Cartesian mesh file to use for the isosurfaces. | | `X` | `X` | `-` |
| *Level value information* | _One of the following options is required: `-lf`, `-lv`, `-gl`_ | _These options set the values that will be used for the isosurfaces in the mesh file._ | | `X` | `X` | `X` |
| Level File | `-lf`/`--levelfile` `LEVELFILE` | Relative path to file containing values to use for isosurface levels. File should be structured to have one value per line. | | `O` | `O` | `O` |
| Level Values | `-lv`/`--levelvalues` `VAL [VAL VAL]` | List of values used to generate isosurfaces in VisIt. | | `O` | `O` | `O` |
| Generate Levels | `-gl`/`--generatelevels` `ratio`/`log`/`lin` |  Specifies the mode for generating level values to be used for the isosurfaces. If used, values for the minimum and maximum levels (`-lx`) and the ratio or number of levels (`-N`) are also required. Options are: (1) `ratio`: `N` is the ratio between levels ranging from the min value up to, but not exceeding, the max value. (2) `log`: `N` is the number of levels to be evenly spaced logarithmically between the min and max values. (3) `lin`: `N` is the number of levels to be evenly spaced linearly between the min and max values. | | `O` | `O` | `-` |
| Number of Levels | `-N`/`--numlevels` `N` | (Required if using `-gl`, otherwise not allowed) | | `X`/`-` | `X`/`-` | `-` |
| Level Min/Max| `-lx`/`--levelextrema` `MIN_VAL MAX_VAL` | (Required if using `-gl`, otherwise not allowed) | | `X`/`-` | `X`/`-` | `-` |
| *Geometry information* | | _These options specify information needed to generate a DAGMC geometry from the isosurfaces produced from VisIt._ | | | | |
| Database Path | `-db`/`--database` `DATABASE_PATH` | Relative path to folder where isosurface meshfiles from VisIt are to be saved (`visit` step) or read from (`moab` step). | A folder called `tmp/` in the current directory: `${PWD}/tmp/` | `O` | `O` | `O` |
| Mesh-based Merge Tolerance | `-m`/`--mergetol` `TOL` | Merge tolerance for mesh based merging of coincident surfaces. | `1e-5` | `O` | `-` | `O` |
| Surface Normalization Factor | `-n`/`--norm` `NORM_FACTOR` |  All level values will be multiplied by this normalization factor when the geometry is generated. | `1.0` | `O` | `-` | `O` |
| Visualization Tag | `-v`/`--viz` | If set, surfaces generated will be tagged with data for visualization purposes in VisIt. Note: this will increase the file size as each individual facet will be tagged with data. | | `O` | `-` | `O` |
| Isosurface Geometry Filename | `-g`/`--geomfile` `GEOM_FILENAME` | Filename to write generated isosurface geometry file. Must be either a .h5m or .vtk file name. | `isogeom.h5m` | `O` | `-` | `O` |
| Save Location | `-sp`/`--savepath` `PATH` | Absolue path to folder to write generated geometry file. | Database Path | `O` | `-` | `O` |
| Extra Tag Information | `-t`/`--tag` `TAGNAME TAGVAL` | Information to tag on the whole geometry. First entry must be the name for the tag (string). Second entry must be the value for the tag (will be tagged as float). Option can be set more than once to set more tags. | | `O` | `-` | `O` |

### Example Usage

All of these examples will assume a starting Cartesian mesh file called `cw_mesh` with the desired data called `wwn`.

* Run all the steps start to finish, defining 3 values for the level information at runtime, and tag for visualization:

      generate_isogeom full cw_mesh wwn -lv 0.1 5.2 12.3 --viz

* Generate a geometry start to finish, with 5 levels logarithmically spaced from 1e-5 to 1e+3. Also tag the geometry two metadata tags called E1 and E2 with values of 1.0 and 10.0, respectively:

      generate_isogeom full meshfile my_data -gl log -lx 1e-5 1e+3 -N 5 -t E1 1.0 -t E2 10.0

* Generate only the isovolume database using levels that logarithmically spaced between 0.1 and 1e+14 and specifying where to write the generated database:

      generate_isogeom visit cw_mesh wwn -gl log -lx 0.1 1e+14 -db my_database/

* Generate an isovolume database in the default location using levels between 1.0 2e+4 that are spaced with a ratio of 20:

      generate_isogeom visit cw_mesh wwn -gl ratio -lx 1.0 2.e4 -N 20

* Generate an isosurface geometry using the levelfile and database located in my_database/, specifying a file name for file produced:

      generate_isogeom moab -lf my_database/levelfile -db my_database -g geom1.h5m

* Generate a geometry from a database located in 'my_isogeom/', read the level info from a file called 'levelinfo', mutliply all data by a factor of 2e4, and save the file as 'my_isogeom.vtk' in a new folder called 'output_folder/':

      generate_isogeom moab -db my_isogeom/ -lf levelinfo -n 2e4 -g my_isogeom.vtk -sp output_folder/
