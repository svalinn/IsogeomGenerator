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

From the source directory, run `pip install -e . --user`.

## Usage

To use, import the Python module with `import IsogeomGenerator`.

Required input file: a 3D mesh file (such as a Cartesian mesh) with scalar
values defined on each mesh voxel.

### Steps:

    1. Set contour levels with `assign_levels()` or `generate_levels()`
    2. Generate the isovolume files using `generate_volumes()`
    3. Create the MOAB geometry with `create_geometry()`
    4. Write to a file with `write_geometry()`
