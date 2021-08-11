"""For a given isosurface geometry file, apply mesh refinement techniques
using the VTK library. Decimation (decrease number of mesh triangles)
and surface smoothing can be applied. If both are applied, smoothing is
applied first as it can introduce additional triangles in some cases.

Mesh Refinements techniques:
----------------------------
    * Decimate: use fewer triangles on the mesh surfaces (vtkDecimatePro)
    * Smooth: decrease surface roughness (vtkWindowedSincPolyDataFilter)
"""

import os
import argparse
import warnings
import numpy as np
from pymoab import core, types, rng
import vtk


def parse_arguments():
    """Parser for user input arguments for the commandline

    Inputs:
    -------
        None

    Return:
    -------
        args: dictionary of parsed arguments
    """
    formatter = argparse.RawDescriptionHelpFormatter

    description = """
This tool will refine surfaces of an isosurface mesh geometry. Options
for refinement include smoothing and decimating. Smoothing with smooth mesh
surfaces (decrease surface roughness). Decimation will decrease the number
of triangles on the surface. If both smoothing and decimation are applied,
smoothing will be done first as it can cause an increase in triangles. If
neither decimation or smoothing is specified, both will be applied by
default using 0.5 factors. Smoothing is only applied to interior surfaces.
Both refinement methods are built on VTK."""

    usage = "refine_isogeom geomfile [OPTIONS]"

    examples = """
Examples:

    1) using default refinement factors (0.5 decimation, 0.5 smoothing)
        refine_isogeom my_geom.h5m
    2) applying only decimation with a 0.7 target reduction
        refine_isogeom my_geom.h5m -d 0.7
    3) applying only smoothing with a 0.4 smoothing factor
        refine_isogeom my_geom.h5m -s 0.4
    4) applying both decimation and smoothing with different factors
        refine_isogeom my_geom.h5m -s 0.4 -d 0.9
"""

    parser = argparse.ArgumentParser(description=description,
                                     usage=usage,
                                     epilog=examples,
                                     formatter_class=formatter)
    parser.add_argument('geomfile',
                        action='store',
                        nargs=1,
                        type=str,
                        help='Relative path an isosurface geometry file to ' +
                        'be refined. Must be an HDF5 (.h5m) file.'
                        )
    parser.add_argument('-d', '--decimate',
                        action='store',
                        required=False,
                        default=None,
                        metavar='DECIMATION_FACTOR',
                        dest='deci_factor',
                        help='Decimation factor for triangles (0 < df < 1). ' +
                        'Higher value means higher reduction in triangles.',
                        type=float
                        )
    parser.add_argument('-s', '--smooth',
                        action='store',
                        required=False,
                        default=None,
                        metavar='SMOOTH_FACTOR',
                        dest='smooth_factor',
                        help='Smoothing relaxation factor for surfaces (0 < sf < 1). ' +
                        'Higher value means more smoothing. ' +
                        'Default value is 0.5.',
                        type=float
                        )
    parser.add_argument('-o', '--output',
                        action='store',
                        required=False,
                        type=str,
                        nargs=1,
                        default=['refined_geom.h5m'],
                        metavar='OUTPUT_FILE',
                        dest='output',
                        help='Name to be used for the refined output file. ' +
                        'Default is refined_geom.h5m')
    parser.add_argument('-k', '--keep',
                        action='store_true',
                        required=False,
                        dest='keep',
                        help='If set, temporary surfaces files that are ' +
                        'generated will not be deleted. (Use for debugging).'
                        )
    args = parser.parse_args()
    return args


def apply_filters(fname, df, sf, surf_type):
    """Apply the decimation and smoothing filters using the VTK library
    and the user specified refinement factors. Decimation uses the
    vtkDecimatePro method. Smoothing uses the vtkWindowSincPolyDataFilter
    method. Both refinement processes maintain boundary vertices and
    topology. If both filters are applied, smoothing will be applied first.

    Input:
    ------
        fname: string, name of the surface vtk file that will be refined
        df: float or None, factor for target reduction in the decimation
            process. Must be between 0 and 1. If set to None, no
            decimation will be applied. A higher value indicates a
            greater reduction in the number of surface triangles.
        sf: float or None, smoothing factor for the smoothing process.
            Must be between 0 and 1. If set to None, no smoothing will
            be applied. A higher value indicates more smoothing.
        surf_type: str, 'interior' or 'exterior', used to indicate if
            the surface being refined is an interior or exterior surface.
            Only interior surfaces will have smoothing applied if smoothing is requested.

    Return:
    -------
        outfile: string, name of the output file from VTK after refinement
    """
    # read in surface w/ vtk as an unstructured grid
    usg = vtk.vtkUnstructuredGridReader()
    usg.SetFileName(fname)
    usg.Update()
    usgport = usg.GetOutputPort()

    # convert USG to polydata
    gf = vtk.vtkGeometryFilter()
    gf.SetInputConnection(usgport)
    gf.Update()
    pdport = gf.GetOutputPort()

    # only apply smoothing filter to interior surfaces to avoid smoothing
    # the corners of the geometry
    if sf is not None and surf_type == 'interior':
        smoothfilter = vtk.vtkWindowedSincPolyDataFilter()
        smoothfilter.SetInputConnection(pdport)
        smoothfilter.SetNumberOfIterations(20)
        smoothfilter.BoundarySmoothingOff  # don't apply to boundaries
        smoothfilter.NonManifoldSmoothingOff()  # don't collapse topology/vols
        smoothfilter.Update()
        pdport = smoothfilter.GetOutputPort()

    if df is not None:
        # setup decimate operator
        decifilter = vtk.vtkDecimatePro()
        decifilter.SetInputConnection(pdport)
        # set decimate options
        decifilter.SetTargetReduction(df)  # target reduction
        # preserve topology (splitting or hole elimination not allowed)
        decifilter.SetPreserveTopology(True)
        decifilter.SetSplitting(False)  # no mesh splitting allowed
        # no boundary vertex (eddge/curve) deletion allowed
        decifilter.SetBoundaryVertexDeletion(False)
        decifilter.Update()
        pdport = decifilter.GetOutputPort()

    # convert polydata back to USG
    appendfilter = vtk.vtkAppendFilter()
    appendfilter.SetInputConnection(pdport)
    appendfilter.Update()
    usg_new = vtk.vtkUnstructuredGrid()
    usg_new.ShallowCopy(appendfilter.GetOutput())

    outfile = fname.split('.')[0] + '_refined.vtk'
    writer1 = vtk.vtkGenericDataObjectWriter()
    writer1.SetFileName(outfile)
    writer1.SetInputData(usg_new)
    writer1.Update()
    writer1.Write()

    return outfile


def get_viz_info(mb, surf):
    """Use MOAB to determine if the triangles of the surface were
    previously tagged individually with a data value for visualization.
    If so, retrieve the data name and tag.

    Input:
    ------
        mb: moab instance with geometry file loaded
        surf: entity handle for a surface to be examined

    Return:
    -------
        data_tag: entity handle or None, if tag is present, the entity
            handle for the data tag.
    """
    # check if data tag is on triangles for viz (only do this one time)
    all_tags = mb.tag_get_tags_on_entity(
        mb.get_entities_by_type(surf, types.MBTRI)[0])
    non_names = ['GEOM_DIMENSION', 'GLOBAL_ID',
                 'CATEGORY', 'GEOM_SENSE_2', 'SURF_TYPE']
    data_name = None
    for tag in all_tags:
        tag_name = tag.get_name()
        if tag_name not in non_names:
            data_name = tag_name
            break

    data_tag = None
    if data_name is not None:
        data_tag = mb.tag_get_handle(data_name, size=1,
                                     tag_type=types.MB_TYPE_DOUBLE,
                                     storage_type=types.MB_TAG_SPARSE,
                                     create_if_missing=False)

    return data_tag


def refine_surfaces(filename, df, sf, output, keep):
    """For a given isogeom file, iterate through each surface applying
    the refinement factors specified. A new geometry file will be
    written out that has gone through the specified refinement.

    Inputs:
    -------
        filename: string, name of the HDF5 (.h5m) isogeometry file to
            to be refined.
        df: float or None, factor for target reduction in the decimation
            process. Must be between 0 and 1. If set to None, no
            decimation will be applied. A higher value indicates a
            greater reduction in the number of surface triangles.
        sf: float or None, smoothing factor for the smoothing process.
            Must be between 0 and 1. If set to None, no smoothing will
            be applied. A higher value indicates more smoothing.
        output: string, name of output file to be written after the
            refinement process
        keep: bool, if set to True, the temporary surface files written
            to disk during the process will not be deleted.

    Return:
    -------
        None
    """
    # load as a moab instance
    mb = core.Core()
    mb.load_file(filename)
    rs = mb.get_root_set()
    # get necessary tag information
    dim_tag = mb.tag_get_handle('GEOM_DIMENSION', size=1,
                                tag_type=types.MB_TYPE_INTEGER,
                                storage_type=types.MB_TAG_SPARSE,
                                create_if_missing=False)
    surf_tag = mb.tag_get_handle('SURF_TYPE', size=32,
                                 tag_type=types.MB_TYPE_OPAQUE,
                                 storage_type=types.MB_TAG_SPARSE,
                                 create_if_missing=False)
    # get all surfaces
    all_surfs = mb.get_entities_by_type_and_tag(rs, types.MBENTITYSET,
                                                dim_tag, [2])
    # get data tag if present on the triangles (for rewriting viz data later)
    data_tag = get_viz_info(mb, all_surfs[0])

    # iterate over all surfaces refining one at a time
    for surf in all_surfs:
        print('Refining surface {}'.format(surf))

        # get surface type (interior or exterior)
        surf_type = mb.tag_get_data(surf_tag, surf)

        # get all triangles and vertices (these will be replaced)
        tris = mb.get_entities_by_type(surf, types.MBTRI)
        verts = mb.get_entities_by_type(surf, types.MBVERTEX)

        # write surface to .vtk file to use in vtk
        fname = 'tmp_{}.vtk'.format(surf)
        mb.write_file(fname, [surf])

        # apply vtk filters
        outfile = apply_filters(fname, df, sf, surf_type)

        # read in vtk file again
        new_surf = mb.create_meshset()
        mb.load_file(outfile, new_surf)

        # get all triangles and vertices on new surf
        tris_new = mb.get_entities_by_type(new_surf, types.MBTRI)
        verts_new = mb.get_entities_by_type(new_surf, types.MBVERTEX)

        #  remove old tris/verts from surfs and vols
        mb.remove_entities(surf, tris)
        mb.remove_entities(surf, verts)
        vols = mb.get_parent_meshsets(surf)
        for vol in vols:
            mb.remove_entities(vol, tris)
            mb.remove_entities(vol, verts)
        mb.add_entities(surf, tris_new)
        mb.add_entities(surf, verts_new)

        # get value of data tag on surface & tag new triangles
        if data_tag is not None:
            data_val = mb.tag_get_data(data_tag, surf)[0][0]
            vals = np.full(len(tris_new), data_val)
            mb.tag_set_data(data_tag, tris_new, vals)

        # delete temporary files
        if not keep:
            os.remove(fname)
            os.remove(outfile)

    # write full geometry - only the necessary entities
    print("Writing refined geometry")
    all_vols = mb.get_entities_by_type_and_tag(
        rs, types.MBENTITYSET, dim_tag, [3])
    all_sets = rng.unite(all_vols, all_surfs)

    # check file extension of save name:
    ext = output.split(".")[-1]
    if ext.lower() not in ['h5m', 'vtk']:
        warnings.warn("File extension {} ".format(ext) +
                      " not recognized. File will be saved as type .h5m.")
        output = output.split(".")[0] + ".h5m"
    mb.write_file(output, all_sets)


def main():
    args = parse_arguments()

    # use defaults for both if not set
    if (args.smooth_factor is None) and (args.deci_factor is None):
        warnings.warn('No refinement factors provided. Applying default ' +
                      'decimation factor (0.5) and smooth factor (0.5).')
        args.smooth_factor = 0.5
        args.deci_factor = 0.5

    refine_surfaces(args.geomfile[0], args.deci_factor, args.smooth_factor,
                    args.output[0], args.keep)


if __name__ == "__main__":
    main()
