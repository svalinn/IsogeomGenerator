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

From the source directory, run `pip install -e . --user`
