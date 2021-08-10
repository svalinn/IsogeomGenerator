"""For a given isosurface geometry file, apply mesh refinement techniques
using the VTK library. Decimation (decrease number of mesh triangles)
and surface smoothing can be applied. If both are applied, smoothing is
applied first as it can introduce additional triangles in some cases.

Mesh Refinements techniques:
----------------------------
    * Decimate: use fewer triangles on the mesh surfaces (vtkDecimatePro)
    * Smooth: decrease surface roughness (vtkWindowedSincPolyDataFilter)
"""

import argparse
import os
import numpy as np
from pymoab import core, types, rng
import vtk


def parse_arguments():
    formatter = argparse.RawDescriptionHelpFormatter
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
                        action='store',
                        required=False,
                        default=None,
                        metavar='DECIMATION_FACTOR',
                        dest='deci_factor',
                        help='Decimation factor for triangles (0 < df < 1). ' +
                        'Higher value means higher reduction in triangles. ' +
                        'Default value is 0.5.',
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
    args = parser.parse_args()
    return args


def apply_filters(fname, df, sf, surf_type):
    # read in suface w/ vtk
    # read unstructured grid
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
    print('wrote {}'.format(outfile))

    return outfile


def get_viz_info(mb, surf):
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

    return data_tag, data_name


def refine_surfaces(filename, df, sf):

    # load as a moab instance
    mb = core.Core()
    mb.load_file(filename)
    rs = mb.get_root_set()
    dim_tag = mb.tag_get_handle('GEOM_DIMENSION', size=1,
                                tag_type=types.MB_TYPE_INTEGER,
                                storage_type=types.MB_TAG_SPARSE,
                                create_if_missing=False)
    surf_tag = mb.tag_get_handle('SURF_TYPE', size=32,
                                 tag_type=types.MB_TYPE_OPAQUE,
                                 storage_type=types.MB_TAG_SPARSE,
                                 create_if_missing=False)

    all_surfs = mb.get_entities_by_type_and_tag(rs, types.MBENTITYSET,
                                                dim_tag, [2])

    data_tag, data_name = get_viz_info(mb, all_surfs[0])

    # iterate over all surfaces refining one at a time
    for surf in all_surfs:
        print('refining surface {}'.format(surf))

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

        # get value of data tag on surface & tag tris again if already tagged
        if data_name is not None:
            data_val = mb.tag_get_data(data_tag, surf)[0][0]
            vals = np.full(len(tris_new), data_val)
            mb.tag_set_data(data_tag, tris_new, vals)

        # delete temporary files
        os.remove(fname)
        os.remove(outfile)

    # wrtite full geometry - only the necessary entities
    print("Writing refined geometry")
    all_vols = mb.get_entities_by_type_and_tag(
        rs, types.MBENTITYSET, dim_tag, [3])
    all_sets = rng.unite(all_vols, all_surfs)
    mb.write_file('refined_geom.h5m', all_sets)


def main():
    args = parse_arguments()
    # use defaults for both if not set
    if (args.smooth_factor is None) and (args.deci_factor is None):
        print('Applying default decimation factor (0.5) and smooth factor (0.5).')
        args.smooth_factor = 0.5
        args.deci_factor = 0.5

    refine_surfaces(args.geomfile[0], args.deci_factor, args.smooth_factor)


if __name__ == "__main__":
    main()
