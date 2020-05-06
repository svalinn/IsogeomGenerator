import argparse
import os
import vol as v


def set_level_options(parser, moab):
    """Sets options for specifying level values. If moab mode, do not
    give options for generating levels.
    """
    level_group = parser.add_mutually_exclusive_group(required=True)
    level_group.add_argument('-lf', '--levelfile',
                             action='store',
                             nargs=1,
                             type=str,
                             help='Relative path to file containing values ' +
                             'to use for isosurface levels. '
                             'File should be structured to have one value ' +
                             'per line.'
                             )
    level_group.add_argument('-lv', '--levelvalues',
                             action='store',
                             nargs='+',
                             default=None,
                             metavar='VAL',
                             type=float,
                             help='List of values used to generate ' +
                             'isosurfaces in VisIt.'
                             )

    # only have option to generate levels if not moab mode
    if not moab:
        level_group.add_argument('-gl', '--generatelevels',
                                 action='store',
                                 nargs=1,
                                 choices=['ratio', 'log', 'lin'],
                                 default=None,
                                 metavar='ratio/log/lin',
                                 type=str,
                                 help='Specifies the mode for generating ' +
                                 'level values to be used for the ' +
                                 'isosurfaces. '
                                 'If used, values for the minimum and ' +
                                 'maximum levels (-lex) and the ratio or ' +
                                 'number of levels (-N) are also required. ' +
                                 'Options are: ' +
                                 '(1) ratio: N is the ratio between levels ' +
                                 'ranging from the min value upto, but not ' +
                                 'exceeding, the max value. ' +
                                 '(2) log: N is the number of levels to be ' +
                                 'evenly spaced logarithmically between the ' +
                                 'min and max values. ' +
                                 '(3) lin: N is the number of levels to be ' +
                                 'evenly spaced linearly between the min ' +
                                 'and max values.'
                                 )
        parser.add_argument('-lx', '--levelextrema',
                            action='store',
                            nargs=2,
                            required=False,
                            default=None,
                            metavar='MIN_VAL MAX_VAL',
                            dest='extN',
                            type=float,
                            help='float, minimum and maximum values to use ' +
                            'for generating the set of level values to be ' +
                            'used for the isosurfaces.'
                            )
        parser.add_argument('-N', '--numlevels',
                            action='store',
                            nargs=1,
                            required=False,
                            default=None,
                            metavar='N',
                            dest='N',
                            type=float,
                            help='If generating levels (-gl), it is either ' +
                            'the ratio between adjacent level values (ratio ' +
                            'mode), or the number of levels to generate ' +
                            '(log or lin mode).'
                            )


def set_visit_only_options(parser):
    """set options specific to the Visit step.
    """
    parser.add_argument('meshfile',
                        action='store',
                        nargs=1,
                        type=str,
                        help='Relative path to the Cartesian mesh file ' +
                        '(vtk format) that will be used to generate ' +
                        'isosurfaces.'
                        )
    parser.add_argument('dataname',
                        action='store',
                        nargs=1,
                        type=str,
                        help='String representing the name of the scalar ' +
                        'data on the Cartesian mesh file to use for the ' +
                        'isosurfaces.'
                        )


def set_moab_only_options(parser):
    """Set options specific to the MOAB step
    """
    parser.add_argument('-m', '--mergetol',
                        action='store',
                        nargs=1,
                        required=False,
                        default=1e-5,
                        metavar='TOL',
                        dest='mergetol',
                        type=float,
                        help='Merge tolerance for mesh based merging of ' +
                        'coincident surfaces. Default=1e-5.'
                        )
    parser.add_argument('-n', '--norm',
                        action='store',
                        nargs=1,
                        required=False,
                        default=1,
                        metavar='NORM_FACTOR',
                        dest='norm',
                        type=float,
                        help='All level values will be multiplied by this ' +
                        'normalization factor when the geometry is ' +
                        'generated. ' +
                        'Default=1'
                        )
    parser.add_argument('-v', '--viz',
                        action='store_true',
                        required=False,
                        dest='tagviz',
                        help='If set, surfaces generated will be tagged ' +
                        'with data for visualization purposes.'
                        )
    parser.add_argument('-g', '--geomfile',
                        action='store',
                        nargs=1,
                        required=False,
                        default=None,
                        metavar='GEOM_FILENAME',
                        dest='geomfile',
                        type=str,
                        help='Filename to write generated isosurface ' +
                        'geometry file. ' +
                        'Must be either a .h5m or .vtk file name. ' +
                        'Default name is isogeom.h5m.'
                        )
    parser.add_argument('-sp', '--savepath',
                        action='store',
                        nargs=1,
                        required=False,
                        default=None,
                        metavar='PATH',
                        dest='savepath',
                        type=str,
                        help='Absolue path to folder to write generated ' +
                        'geometry file. If not set, file will be saved in ' +
                        'the database folder (-db).'
                        )
    parser.add_argument('-t', '--tag',
                        action='append',
                        nargs=2,
                        required=False,
                        default=None,
                        metavar='TAGNAME TAGVAL',
                        dest='tags',
                        help='Information to tag on the whole geometry. ' +
                        'First entry must be the name for the tag (string). ' +
                        'Second entry must be the value for the tag (will ' +
                        'be tagged as float). ' +
                        'Option can be set more than once to set more tags.'
                        )


def set_shared_options(parser, moab=False):
    """set options that are in both the moab and visit steps.
    """
    set_level_options(parser, moab)
    parser.add_argument('-db', '--database',
                        action='store',
                        nargs=1,
                        required=False,
                        default="/tmp",
                        metavar='DATABASE_PATH',
                        dest='db',
                        type=str,
                        help='Relative path to folder where isosurface ' +
                        'meshfiles from VisIt are located or to be saved. ' +
                        'Default is a folder named tmp/ in the current ' +
                        'directory.'
                        )


def parse_arguments():
    """parse user args
    """

    parser = argparse.ArgumentParser(description='Generate isosurface ' +
                                     'geometry from a Cartesian mesh file ' +
                                     'with VisIt and MOAB')
    subparsers = parser.add_subparsers(title='Run Steps',
                                       help='Select which steps to run for ' +
                                       'generating the geometry.')

    # set full mode options
    full_parser = subparsers.add_parser('full', help='Start-to-finish ' +
                                        'generation from a Cartesian mesh ' +
                                        'file to a DAGMC-compliant geometry.')
    set_visit_only_options(full_parser)
    set_shared_options(full_parser)
    set_moab_only_options(full_parser)
    full_parser.set_defaults(which='full')

    # set visit only mode options
    visit_parser = subparsers.add_parser('visit', help='Only generate the ' +
                                         'isosurface mesh files using VisIt.')
    set_visit_only_options(visit_parser)
    set_shared_options(visit_parser)
    visit_parser.set_defaults(which='visit')

    # set moab only mode options
    moab_parser = subparsers.add_parser('moab', help='Only generate the ' +
                                        'DAGMC-compliant geometry with MOAB ' +
                                        'starting from the VisIt mesh files.')
    set_shared_options(moab_parser, moab=True)
    set_moab_only_options(moab_parser)
    moab_parser.set_defaults(which='moab')

    args = parser.parse_args()
    return args


def check_level_gen(args):
    """Check that correct args were supplied if generating level
    """
    if (args.extN is None) or (args.N is None):
        raise RuntimeError("Min/Max level values (-lx) and number of levels " +
                           "(-N) must be set to use --generatelevels option.")


def get_levels(args, g):
    """Get level information depending on supplied args.

    Input:
    ------
        args: argparse ArgumentParser
        g: IsoVolume instance
    """
    # collect level information:
    if args.levelfile is not None:
        # option 1: read from file
        g.read_levels(args.levelfile[0])
    elif args.levelvalues is not None:
        # option 2: set values at run time
        g.assign_levels(args.levelvalues)
    elif args.generatelevels is not None:
        # option 3: generate levels
        check_level_gen(args)
        minN = min(args.extN)
        maxN = max(args.extN)
        N = args.N[0]
        g.generate_levels(N, minN, maxN, mode=args.generatelevels[0])
    else:
        raise RuntimeError("Mode for setting level information is not " +
                           "recognized")


def process_tags(tags):
    """Convert the provided tag information to a dictionary to pass to the
    geometry creation stepself.

    Input:
    -----
        tags: list of lists, [[tagname1, tagval1], [tagname2, tagval2], ...]

    Return:
    -------
        tagdict: dict, key=TAGNAME, value=TAGVALUE
    """
    tagdict = {}
    for tagset in tags:
        key = tagset[0]
        val = float(tagset[1])
        tagdict[key] = val

    return tagdict


def main():

    # get args
    args = parse_arguments()
    print(args)

    # create instance
    g = v.IsoVolume()

    # get level info regardless of mode
    get_levels(args, g)

    # get database information
    db = os.getcwd() + '/' + args.db[0]

    # run steps depending on mode
    mode = args.which
    if mode == 'full' or mode == 'visit':
        # generate isosurfaces
        meshfile = args.meshfile[0]
        dataname = args.dataname[0]
        g.generate_volumes(meshfile, dataname, dbname=db)

    if mode == 'full' or mode == 'moab':
        if args.tags is not None:
            tags = process_tags(args.tags)
        else:
            tags = None

        # create/write geometry
        g.create_geometry(tag_for_viz=args.tagviz[0],
                          norm=args.norm[0],
                          merg_tol=args.merg_tol[0],
                          tags=tags,
                          dbname=db,
                          sname=args.geomfile[0],
                          sdir=args.savepath[0])


if __name__ == "__main__":
    main()
