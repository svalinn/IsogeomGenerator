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

Each mode has different required or optional arguments needed at run time.
Below is a table summary of those options. This information is also available
via terminal help message: (`generate_isogeom [mode] --help`).

* `X` = required
* `O` = optional
* `-` = Not allowed

|||||| Mode ||
|------|
| | Option | Description | Default value | full | visit | moab |
| Mesh file information |||||||
| |`meshfile` ||| `X` | `X` | `-` |
| |`dataname` ||| `X` | `X` | `-` |
| Level value information | One of the following options is required: `-lf`, `-lv`, `-gl` | | | `X` | `X` | `X` |
| | `-lf/--levelfile LEVELFILE` | | | `O` | `O` | `O` |
| | `-lv/--levelvalues VAL [VAL VAL]` | | | `O` | `O` | `O` |
| | `-gl/--generatelevels ratio/log/lin` | | | `O` | `O` | `-` |
| | `-N/--numlevels N` | (Required if using `-gl`, otherwise not allowed) | | `X`/`-` | `X`/`-` | `-` |
| | `-lx/-levelextrema MIN_VAL MAX_VAL` | (Required if using `-gl`, otherwise not allowed) | | `X`/`-` | `X`/`-` | `-` |
| Geometry Information | | | | | | |
| | `-db/--database DATABASE_PATH` | | | `O` | `O` | `O` |
| | `-m/--mergetol TOL` | | | `O` | `-` | `O` |
| | `-n/--norm NORM_FACTOR` | | | `O` | `-` | `O` |
| | `-v/--viz` | | | `O` | `-` | `O` |
| | `-g/--geomfile GEOM_FILENAME` | | | `O` | `-` | `O` |
| | `-sp/--savepath PATH` | | | `O` | `-` | `O` |
