"""For a given isosurface geometry file, apply mesh refinement techniques.

Mesh Refinements techniques:
----------------------------

    * Decimate: use fewer triangles in the mesh surfaces
    * smooth: remove roughness from the surfaces
"""

import argparse
from pymoab import core

formatter = argparse.RawDescriptionHelpFormatter


def decimate_triangles(factor):
    """Reduce the number of triangles in the surfaces by user specified
    amount
    """
    pass


def smooth_surfaces(factor):
    """Smooth surfaces by a user-specified factor"""
    pass


def refine_all():
    """smooth and decimate surfaces. Smoothing will be applied first"""
    pass


def parse_arguments():
    description = "temp description"
    usage = "temp usage"
    examples = "temp examples"
    parser = argparse.ArgumentParser(description=description,
                                     usage=usage,
                                     epilog=examples,
                                     formatter_class=formatter)
    parser.add_argument('geomfile',
                        action='store',
                        nargs=1,
                        type=str,
                        help='Relative path an isosurface geometry file to ' +
                        'be refined.'
                        )
    parser.add_argument('-d', '--decimate',
                        action='store_true',
                        required=False,
                        dest='decimate',
                        help='If set, decimate surface triangles by the ' +
                        'decimation factor df.')
    parser.add_argument('-df', '--decimatefactor',
                        action='store',
                        required=False,
                        default=[0.5],
                        metavar='DECIMATION_FACTOR',
                        dest='deci_factor',
                        help='Decimation factor for triangles (0 < df < 1). ' +
                        'Higher value means higher reduction in triangles. ' +
                        'Default value is 0.5.',
                        type=float
                        )
    parser.add_argument('-s', '--smooth',
                        action='store_true',
                        required=False,
                        dest='smooth',
                        help='If set, smooth surfaces by the ' +
                        'smoothing factor')
    parser.add_argument('-sf', '--smoothfactor',
                        action='store',
                        required=False,
                        default=[0.5],
                        metavar='SMOOTH_FACTOR',
                        dest='smooth_factor',
                        help='Smoothing relaxation factor for surfaces (0 < sf < 1). ' +
                        'Higher value means more smoothing. ' +
                        'Default value is 0.5.',
                        type=float
                        )
    args = parser.parse_args()
    return args


def main():
    args = parse_arguments()

    if args.smooth:
        smooth_surfaces(args.smooth_factor[0])
    if args.decimate:
        decimate_triangles(args.decimate_factor[0])


if __name__ == "__main__":
    main()
