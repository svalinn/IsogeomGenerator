import argparse
import os
from IsogeomGenerator import driver, isg, ivdb

"""This is a script that can be installed for a user to easily run all the
steps to create an isosurface geometry from a mesh file with scalar data from
the command line. This script will parse the command line options and call the
necessary methods from the vol.py file. There are three modes for running this
script: full, visit, and moab. Each are documented below.
"""
formatter = argparse.RawDescriptionHelpFormatter


def set_level_options(parser, moab):
    """Sets options for specifying level values.

    Three options are available (two available for MOAB mode). Exactly one
    option must be supplied. If in MOAB mode (moab=True), do not give options
    for generating levels.

    Input:
    ------
        parser: ArgumentParser object to attach options to
        moab: bool, True if setting MOAB args
    """
    level_group = parser.add_mutually_exclusive_group(required=True)
    level_group.add_argument('-lf', '--levelfile',
                             action='store',
                             nargs=1,
                             default=[None],
                             type=str,
                             help='Relative path to file containing values ' +
                             'to use for isosurface levels. '
                             'File should be structured to have one value ' +
                             'per line.'
                             )
    level_group.add_argument('-lv', '--levelvalues',
                             action='store',
                             nargs='+',
                             default=[None],
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
                                 default=[None],
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
                            default=[None],
                            metavar=('MIN_VAL', 'MAX_VAL'),
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
                            default=[None],
                            metavar='N',
                            dest='N',
                            type=float,
                            help='If generating levels (-gl), it is either ' +
                            'the ratio between adjacent level values (ratio ' +
                            'mode), or the number of levels to generate ' +
                            '(log or lin mode).'
                            )


def set_visit_only_options(parser):
    """Set options specific to the VisIt step.

    Input:
    ------
        parser: ArgumentParser object to attach options to
    """
    parser.add_argument('meshfile',
                        action='store',
                        nargs=1,
                        type=str,
                        help='Relative path to the Cartesian mesh file ' +
                        '(vtk format) that will be used to generate ' +
                        'isosurfaces.'
                        )


def set_moab_only_options(parser):
    """Set options specific to the MOAB step.

    Input:
    ------
        parser: ArgumentParser object to attach options to
    """
    parser.add_argument('-m', '--mergetol',
                        action='store',
                        nargs=1,
                        required=False,
                        default=[1e-5],
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
                        default=[1.0],
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
                        default=[None],
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
                        default=[None],
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
                        default=[],
                        metavar=('TAGNAME', 'TAGVAL'),
                        dest='tags',
                        help='Information to tag on the whole geometry. ' +
                        'First entry must be the name for the tag (string). ' +
                        'Second entry must be the value for the tag (will ' +
                        'be tagged as float). ' +
                        'Option can be set more than once to set more tags.'
                        )


def set_shared_options(parser, moab=False):
    """Set options that are for both the MOAB and VisIt steps.

    Input:
    ------
        parser: ArgumentParser object to attach options to
        moab: bool, True if setting MOAB args (default=False)
    """
    set_level_options(parser, moab)
    parser.add_argument('-db', '--database',
                        action='store',
                        nargs=1,
                        required=False,
                        default=["/tmp"],
                        metavar='DATABASE_PATH',
                        dest='db',
                        type=str,
                        help='Relative path to folder where isosurface ' +
                        'meshfiles from VisIt are located or to be saved. ' +
                        'Default is a folder named tmp/ in the current ' +
                        'directory.'
                        )
    parser.add_argument('dataname',
                        action='store',
                        nargs=1,
                        type=str,
                        help='String representing the name of the scalar ' +
                        'data on the Cartesian mesh file to use for the ' +
                        'isosurfaces.'
                        )


def parse_arguments():
    """Parse user args

    There are three subparsers, one for each mode: full, visit, and moab.
    Full mode runs both the visit and moab steps. Each parser should have a
    full help message, simplified usage statement, and examples.
    """
    mode_examples = """
To view full options for each mode, use 'generate_isogeom MODE -h'.

Example usage:
    (1) Run all the steps start to finish (full mode) starting with meshfile
        'cw_mesh', scalar data 'wwn', and defining 3 values for the level
        information at runtime:

        generate_isogeom full cw_mesh wwn -lv 0.1 5.2 12.3

    (2) Run just the first step (visit mode), generating logarithmically spaced
        levels between 0.1 and 1e+14 and specifying where to write the
        generated database:

        generate_isogeom visit cw_mesh wwn -gl log -lx 0.1 1e14 -db my_database

    (3) Run only the second step (moab mode), using the levelfile and database
        from the MOAB step, and specifying a file name for file produced:

        generate_isogeom moab -lf my_database/levelfile -db my_database
            -g geom1.h5m
"""
    mode_description = """
Use this to generate a full isosurface geometry from a starting Cartesian mesh
file containing scalar data using VisIt and MOAB. This tool can be run in three
different modes:
    full: run both steps starting from the Cartesian mesh file to produce
        a full DAGMC-compliant isosurface geom. This step first runs the visit
        step then the moab step.
    visit: run only the first step using VisIt. This will generate a database
        of individual mesh isosurfaces from the Cartesian mesh fileself.
    moab: run only the second step using MOAB. This will generate a full DAGMC-
        compliant isosurface geometry starting from the database generated from
        the visit step.
"""
    parser = argparse.ArgumentParser(description=mode_description,
                                     usage='generate_isogeom MODE [OPTIONS]',
                                     epilog=mode_examples,
                                     formatter_class=formatter)
    subparsers = parser.add_subparsers(title='Modes',
                                       help='Select which steps to run for ' +
                                       'generating the geometry.')

    # set full mode options
    full_description = """
Start-to-finish generation from a Cartesian mesh file to a DAGMC-compliant
geometry.

Levels information must be provided with either the -lf, -lv, or -gl option.
If using the -gl option (generate levels), then options -lx and -N must also be
provided.
"""
    full_usage = \
        'generate_isogeom full meshfile dataname [-lf/-lv/-gl] [OPTIONS]'
    full_examples = """
Example Usage:
    (1) Create an isosurface geometry called 'my_isogeom.h5m' with assigned
        level values of 0.1 0.4 and 1.0, and tag the surfaces with data for
        vizualization:

        generate_isogeom full meshfile my_data -lv 0.1 0.4 1.0
            -g my_isogeom.h5m --viz

    (2) Generate a geometry with 5 levels lograthmically spaced from 1e-5 and
        1e+3. Also tag the geometry two metadata tags called E1 and E2 with
        values of 1.0 and 10.0, respectively:

        generate_isogeom full meshfile my_data -gl log -lx 1e-5 1e+3 -N 5
            -t E1 1.0 -t E2 10.0

    (3) Store the generated database in a different folder called 'my_isogeom/'
        and read level information from a file called 'levelfile' located in
        the current directory:

        generate_isogeom full meshfile my_data -lf levelfile -db my_isogeom/
    """
    full_parser = subparsers.add_parser('full',
                                        description=full_description,
                                        usage=full_usage,
                                        epilog=full_examples,
                                        formatter_class=formatter)
    set_visit_only_options(full_parser)
    set_shared_options(full_parser)
    set_moab_only_options(full_parser)
    full_parser.set_defaults(which='full')

    # set visit only mode options
    visit_description = """
Only generate the isosurface mesh file database using VisIt.

Levels information must be provided with either the -lf, -lv, or -gl option.
If using the -gl option (generate levels), then options -lx and -N must also be
provided.
"""
    visit_usage = \
        'generate_isogeom visit meshfile dataname [-lf/-lv/-gl] [OPTIONS]'
    visit_examples = """
Example Usage:
    (1) Generate a database located at 'my_database/' with assigned
        level values of 0.1 0.4 and 1.0:

        generate_isogeom visit meshfile my_data -lv 0.1 0.4 1.0
            -db my_isogeom/

    (2) Generate a database in the default location using levels between 1.0
        2e+4 that are spaced with a ratio of 20:

        generate_isogeom visit meshfile my_data -gl ratio -lx 1.0 2.e4 -N 20

    (3) Generate a database in the default location using 15 levels between 1.0
        2e+4 that are spaced logarithmically:

        generate_isogeom visit meshfile my_data -gl log -lx 1.0 2.e4 -N 15

    (4) Generate a database in a folder called 'my_isogeom/' and read the level
        information from a file in the current directory called 'levelfile':

        generate_isogeom visit meshfile my_data -lf levelfile -db my_isogeom/
    """
    visit_parser = subparsers.add_parser('visit',
                                         description=visit_description,
                                         usage=visit_usage,
                                         epilog=visit_examples,
                                         formatter_class=formatter)
    set_visit_only_options(visit_parser)
    set_shared_options(visit_parser)
    visit_parser.set_defaults(which='visit')

    # set moab only mode options
    moab_description = """
Only generate the DAGMC-compliant geometry with MOAB starting from the VisIt
mesh file database.

Levels information must be provided with either the -lf or -lv option.
"""
    moab_usage = 'generate_isogeom moab [-lf/-lv] [OPTIONS]'
    moab_examples = """
Example Usage:
    (1) Create an isosurface geometry called 'my_isogeom.h5m' with assigned
        level values of 0.1 0.4 and 1.0, and tag the surfaces with data for
        vizualization (assume default database location):

        generate_isogeom moab -lv 0.1 0.4 1.0 -g my_isogeom.h5m --viz

    (2) Generate a geometry from a database located in 'my_isogeom/', read the
        level info from a file called 'levelinfo', mutliply all data by a
        factor of 2e4, and save the file as 'my_isogeom.vtk' in a new folder
        called 'output_folder/':

        generate_isogeom moab -db my_isogeom/ -lf levelinfo -n 2e4
            -g my_isogeom.vtk -sp output_folder/

    (3) Generate a geometry from a database in the default location, read
        levels from a file called 'levelfile' located in the database, tag the
        geometry two metadata tags called E1 and E2 with values of 1.0 and
        10.0, respectively, and tag the geometry with the level information for
        vizualization:

        generate_isogeom moab -lf tmp/levelfile -t E1 1.0 -t E2 10.0 -v
    """
    moab_parser = subparsers.add_parser('moab',
                                        description=moab_description,
                                        usage=moab_usage,
                                        epilog=moab_examples,
                                        formatter_class=formatter)
    set_shared_options(moab_parser, moab=True)
    set_moab_only_options(moab_parser)
    moab_parser.set_defaults(which='moab')

    args = parser.parse_args()
    return args


def check_level_gen(args):
    """Check that correct args were supplied for generating levels.

    If the min/max values or N are not supplied, then raise error.

    Input:
    ------
        args: set of ArgumentParser args
    """
    if (args.extN[0] is None) or (args.N[0] is None):
        raise RuntimeError("Min/Max level values (-lx) and number of levels " +
                           "(-N) must be set to use --generatelevels option.")


def get_levels(args):
    """Get level information depending on supplied args.

    Input:
    ------
        args: set of ArgumentParser args
    """
    # collect level information:
    if args.levelfile[0] is not None:
        # option 1: read from file, set levels to file name
        levels = args.levelfile[0]
    elif args.levelvalues[0] is not None:
        # option 2: set list of values
        levels = args.levelvalues
    elif args.generatelevels is not None:
        # option 3: generate levels
        check_level_gen(args)
        minN = min(args.extN[0])
        maxN = max(args.extN[0])
        levels = driver.generate_levels(args.N[0], minN, maxN,
                                        mode=args.generatelevels[0])
    else:
        raise RuntimeError("Mode for setting level information is not " +
                           "recognized")
    return levels


def process_tags(tags):
    """Process the provided tag information to correct format.

    Converts the list of information to a dictionary to be able to pass to the
    geometry creation step.

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

    # generate level info if necessary
    # levels is either a list of values or path to file
    levels = get_levels(args)

    # get database information
    db = os.getcwd() + '/' + args.db[0]

    # get dataname
    data = args.dataname[0]

    # run steps depending on mode
    mode = args.which
    visit_modes = ["full", "visit"]
    moab_modes = ["full", "moab"]
    iv = None  # initialize

    if mode in visit_modes:
        iv = ivdb.IvDb(levels=levels, data=data, db=db)
        driver.generate_volumes(iv, args.meshfile[0])

    if mode in moab_modes:
        if args.tags != []:
            tags = process_tags(args.tags)
        else:
            tags = None

        # pass IvDb info if object exists from previous step
        if iv is not None:
            ig = isg.IsGm(ivdb=iv)
        else:
            ig = isg.IsGm(levels=levels, data=data, db=db)

        # create geom from info
        driver.create_geometry(ig,
                               tag_for_viz=args.tagviz,
                               norm=args.norm[0],
                               merge_tol=args.mergetol[0],
                               tags=tags,
                               sname=args.geomfile[0],
                               sdir=args.savepath[0])


if __name__ == "__main__":
    main()
