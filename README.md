# Isosurface Geometry Generator

This tool can be used to generate an isosurface geometry from any 3D mesh
tagged with scalar values (eg, a Cartesian weight window mesh used for
Monte Carlo particle transport).
The tool will create a meshed surface geometry using specified isosurface
values.

## Dependencies:

* Python 2.7
* [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/)
* [MOAB](https://sigma.mcs.anl.gov/moab-library/) v5.1+ with PyMOAB enabled

## Installation:

From the source directory, run `pip install . --user`.

## Python Module Usage

To use, import the Python module with `from IsogeomGenerator import vol`.

### Steps:

1. Set contour levels with `assign_levels()`, `read_levels()` or `generate_levels()`
2. Generate the isovolume files using `generate_volumes()`
3. Create the MOAB geometry with `create_geometry()`

## Command Line Tool

The steps for creating an isosurface geometry can be done on the command line
with the `generate_isogeom` command. This tool can be run in three different
modes:

* full: this will run both the visit step then the moab step.
(command: `generate_isogeom full ...`)
* visit: starting from a Cartesian mesh file, this will generate only a database of
separate isosurface files from VisIt (step 2 above). (command: `generate_isogeom visit ...`)
* moab: this will start from the database (generated in the visit step) to create
a DAGMC-compliant isosurface geometry using PyMOAB (step 3 above).
(command: `generate_isogeom moab ...`)

To view the three different modes, run `generate_isogeom --help`.

### Use
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
