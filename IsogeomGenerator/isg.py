import sys
import os
import shutil
import warnings
import numpy as np
import math as m

from isg_gen import IsoGeomGen

from pymoab import core, types
from pymoab.rng import Range
from pymoab.skinner import Skinner


class IsGm(IsoGeomGen):
    """Class containing necessary methods for generating a full mesh
    geometry from a database of isosurface volumes. This class uses the
    python interface for MOAB (pyMOAB).

    Attributes:
    -----------
        levels: list of floats, values used for isosurface values
        data: string, name of data on mesh
        db: string, path to database folder with isovolume files
        mb: MOAB core instance, contains all information of the mesh
            geometry
        isovol_meshsets: dictionary, information relating curve,
            surface, and volume entity handles to each other.
        val_tag: MOAB tag entity handle, tag for surface value
        sense_tag: MOAB tag entity handle, tag for surface sense

    Methods:
    --------
    """

    def __init__(self, ivdb=None, levels=None, data=None, db=None,
                 extents=None):
        """Create IsGm object. Information provided by an ivdb object
        will overwrite other data provided.

        Input:
        ------
            ivdb: (object), IvDb object, must be a completed object
            levels: (optional), list of floats or string, values to use
                for isosurface values. If string, then it is the path
                to levelfile where one value is given per line.
            data: (optional), string, name of data on the mesh
            db: (optional), string, path to database folder with
                isovolume files
            extents: (optional) list of list of floats, minimum and
                maximum values for x, y, and z in mesh. Must be
                structured like [[xmin, ymin, zmin], [xmax, ymax, zmax]]
        """
        # initialize variables
        super(IsGm, self).__init__(levels, data, db, extents)

        # if ivdb object is provided, overwrite with that info
        if ivdb is not None:
            self.read_ivdb(ivdb)

        # set MOAB related attributes
        self.mb = core.Core()
        self.isovol_meshsets = {}
        self.val_tag = \
            self.mb.tag_get_handle(self.data, size=1,
                                   tag_type=types.MB_TYPE_DOUBLE,
                                   storage_type=types.MB_TAG_SPARSE,
                                   create_if_missing=True)
        self.sense_tag = \
            self.mb.tag_get_handle('GEOM_SENSE_2', size=2,
                                   tag_type=types.MB_TYPE_HANDLE,
                                   storage_type=types.MB_TAG_SPARSE,
                                   create_if_missing=True)

    def read_ivdb(self, ivdb):
        """read information from IvDb object.

        Input:
        ------
            ivdb: completed IvDb object
        """
        # check that object is completed
        if not ivdb.completed:
            raise RuntimeError("Incomplete IvDb object was " +
                               "provided. Please run 'generate_volumes()'.")
        # set values
        self.levels = ivdb.levels
        self.db = ivdb.db
        self.data = ivdb.data
        self.xmin = ivdb.xmin
        self.xmax = ivdb.xmax
        self.ymin = ivdb.ymin
        self.ymax = ivdb.ymax
        self.zmin = ivdb.zmin
        self.zmax = ivdb.zmax

    def read_database(self):
        """Read the files from the database and initialize the meshset info.
        """
        # check that levels exist:
        if self.levels is None:
            raise RuntimeError("Object must have levels defined.")

        # check there are correct number of files:
        file_list = sorted(os.listdir(self.db + "/vols/"))
        if len(self.levels) != len(file_list):
            raise RuntimeError("Number of levels does not match number of " +
                               "isovolume files in the database.")

        # read files
        for f in file_list:
            # get file name
            fpath = self.db + "/vols/" + f
            i = int(f.strip(".stl"))  # must be an integer

            # load file and create EH for file-set
            fs = self.mb.create_meshset()
            self.mb.load_file(fpath, file_set=fs)

            # initiate dictionary
            iv_info = (i, fs)
            self.isovol_meshsets[iv_info] = {}

            # add value min/max info (min, max)
            if i == 0:
                self.isovol_meshsets[iv_info]['bounds'] =\
                    (None, self.levels[i])
            elif i == len(self.levels):
                self.isovol_meshsets[iv_info]['bounds'] =\
                    (self.levels[i - 1], None)
            else:
                self.isovol_meshsets[iv_info]['bounds'] =\
                    (self.levels[i - 1], self.levels[i])

    def separate_isovols(self):
        """Split isosurfaces into different surfaces for exterior vs
        interior surfaces. Exterior surfaces are those in which full
        triangles are on the planes defining the bounding box of the
        geometry. For each isovolume in the database, separate any disjoint
        surfaces into unique single surfaces.
        """
        for iv_info in self.isovol_meshsets.keys():
            # extract isovolume information
            iso_id = iv_info[0]
            fs = iv_info[1]

            # get set of all vertices for the isosurface
            all_verts = self.mb.get_entities_by_type(fs, types.MBVERTEX)

            # separate the verts into interior vs exterior
            verts_interior = []
            verts_exterior = []

            for vert in all_verts:
                c = self.mb.get_coords(vert)
                if (c[0] == self.xmin) or (c[0] == self.xmax) or \
                    (c[1] == self.ymin) or (c[1] == self.ymax) or \
                    (c[2] == self.zmin) or (c[2] == self.zmax):
                    verts_exterior.append(vert)
                else:
                    verts_interior.append(vert)

            # get sets of interior vs exterior tris
            # all tris connected to interior verts are automatically
            # considered interior tris
            tris_interior = self.mb.get_adjacencies(
                verts_interior, 2, op_type=1)

            # exterior tris are only those with all vertices on an
            # exterior surface
            tris_exterior = self.__get_surf_triangles(verts_exterior)

            # create interior and exterior surface meshsets
            surf_exterior = self.mb.create_meshset()
            self.mb.add_entities(surf_exterior, tris_exterior)
            self.mb.add_entities(surf_exterior, verts_exterior)

            surf_interior = self.mb.create_meshset()
            self.mb.add_entities(surf_interior, tris_interior)
            self.mb.add_entities(surf_interior, verts_interior)

            # separate meshset into disjoint surfaces based on their
            # connectedness
            ext_surfs = self.__separate(surf_exterior)
            int_surfs = self.__separate(surf_interior)

            # tag the surfaces with whether they are interior or exterior
            surf_type_tag = \
                self.mb.tag_get_handle('SURF_TYPE', size=32,
                                       tag_type=types.MB_TYPE_OPAQUE,
                                       storage_type=types.MB_TAG_SPARSE,
                                       create_if_missing=True)
            for s in ext_surfs:
                self.mb.tag_set_data(surf_type_tag, s, 'exterior')
            for s in int_surfs:
                self.mb.tag_set_data(surf_type_tag, s, 'interior')

            # store separate surface entity handles
            self.isovol_meshsets[iv_info]['surfs_EH'] = []
            self.isovol_meshsets[iv_info]['surfs_EH'].extend(ext_surfs)
            self.isovol_meshsets[iv_info]['surfs_EH'].extend(int_surfs)

    def imprint_merge(self, norm):
        """Uses PyMOAB to check if surfaces are coincident. Creates a
        single surface where surfaces are coincident values are tagged
        on each surface. Surface senses are also determined and tagged.

        Input:
        ------
            norm: float, All data values will be multiplied by this factor.
        """
        # get list of all original isovolumes
        all_vols = sorted(self.isovol_meshsets.keys())
        for i, isovol in enumerate(all_vols):
            if i != len(self.levels) - 1:
                # do not need to check the last isovolume because it
                # will be checked against its neighbor already
                self.__compare_surfs(isovol, all_vols[i + 1], norm)

        # if a surface doesn't have a value tagged after merging
        # give it a value of 0 and tag forward sense
        for isovol in all_vols:
            for surf in self.isovol_meshsets[isovol]['surfs_EH']:

                # tag val=0
                try:
                    val = self.mb.tag_get_data(self.val_tag, surf)
                except:
                    val = 0.0
                    self.mb.tag_set_data(self.val_tag, surf, val)
                    verts = \
                        self.mb.get_entities_by_type(surf,
                                                     types.MBVERTEX)
                    tris = self.__get_surf_triangles(verts)
                    self.mb.add_entities(surf, tris)

                # tag fwd sense
                try:
                    sense = self.mb.tag_get_data(self.sense_tag, surf)
                except:
                    fwd = isovol[1]
                    bwd = np.uint64(0)
                    self.mb.tag_set_data(self.sense_tag,
                                         surf, [fwd, bwd])

    def make_family(self):
        """Makes the correct parent-child relationships with volumes
        and surfaces. Tags geometry type, category, and ID on surfaces
        and volumes.
        """
        # create geometry dimension, category, and global id tags
        geom_dim = \
            self.mb.tag_get_handle('GEOM_DIMENSION', size=1,
                                   tag_type=types.MB_TYPE_INTEGER,
                                   storage_type=types.MB_TAG_SPARSE,
                                   create_if_missing=True)
        category = \
            self.mb.tag_get_handle('CATEGORY', size=32,
                                   tag_type=types.MB_TYPE_OPAQUE,
                                   storage_type=types.MB_TAG_SPARSE,
                                   create_if_missing=True)
        global_id = \
            self.mb.tag_get_handle('GLOBAL_ID', size=1,
                                   tag_type=types.MB_TYPE_INTEGER,
                                   storage_type=types.MB_TAG_SPARSE,
                                   create_if_missing=True)

        vol_id = 0
        surf_id = 0
        completed_surf_list = []
        for v in self.isovol_meshsets.keys():
            vol_eh = v[1]

            # tag volume
            self.mb.tag_set_data(geom_dim, vol_eh, 3)
            self.mb.tag_set_data(category, vol_eh, 'Volume')
            vol_id += 1
            self.mb.tag_set_data(global_id, vol_eh, vol_id)

            for surf_eh in self.isovol_meshsets[v]['surfs_EH']:
                # create relationship
                self.mb.add_parent_child(vol_eh, surf_eh)

                # tag surfaces (if not already done)
                if surf_eh not in completed_surf_list:
                    self.mb.tag_set_data(geom_dim, surf_eh, 2)
                    self.mb.tag_set_data(category, surf_eh, 'Surface')
                    surf_id += 1
                    self.mb.tag_set_data(global_id, surf_eh, surf_id)
                    completed_surf_list.append(surf_eh)

    def tag_for_viz(self):
        """Tags all triangles on all surfaces with the data value for
        that surface. This is for vizualization purposes.
        """
        for isovol in self.isovol_meshsets.keys():
            for surf in self.isovol_meshsets[isovol]['surfs_EH']:
                # get the tagged data
                val = self.mb.tag_get_data(self.val_tag, surf)

                # get the triangles
                tris = self.mb.get_entities_by_type(surf,
                                                    types.MBTRI)

                # create data array
                num = len(tris)
                data = np.full((num), val)

                # tag the data
                self.mb.tag_set_data(self.val_tag, tris, data)

    def set_tags(self, tags):
        """Set provided tag values on the root set.

        Tag size and type will be determined by the TAGVALUE. Any
        TAGNAME that is not a string will be converted to a string.

        Input:
        ------
            tags: dict, key=TAGNAME, value=TAGVALUE, TAGNAMEs should be
                strings and TAGVALUEs can be single values or lists of
                any type.
        """
        rs = self.mb.get_root_set()
        for tagname, tagval in tags.items():
            # get tag size
            try:
                taglength = len(tagval)  # strings or lists of ints/floats
            except:
                taglength = 1  # can't get length on a single int/float values

            # get datatype:
            if taglength > 1:
                mbtype = types.pymoab_data_type(type(tagval[0]))
            else:
                mbtype = types.pymoab_data_type(type(tagval))

            # create tag
            tag = self.mb.tag_get_handle(str(tagname), size=taglength,
                                         tag_type=mbtype,
                                         storage_type=types.MB_TAG_SPARSE,
                                         create_if_missing=True)
            self.mb.tag_set_data(tag, rs, tagval)

    def write_geometry(self, sname, sdir):
        """Writes out the geometry stored in memory.

        Input:
        ------
            sname: string, name of file to save written file
            sdir: string, absolute path for writing file
        """
        # check file extension of save name:
        ext = sname.split(".")[-1]
        if ext.lower() not in ['h5m', 'vtk']:
            warnings.warn("File extension {} ".format(ext) +
                          " not recognized. File will be saved as type .h5m.")
            sname = sname.split(".")[0] + ".h5m"
        # save the file
        save_location = sdir + "/" + sname
        self.mb.write_file(save_location)
        print("Geometry file written to {}.".format(save_location))

    def __separate(self, ms):
        """For a given surface meshset, separate meshset into unique and
        disjoint surfaces based on their connectedness.

        Input:
        ------
            ms: meshset entity handle

        Returns:
        --------
            surf_list: list of entity handles for the set of unique and
                disjoint surfaces in the meshset
        """
        surf_list = []

        # get set of all vertices for the isosurface
        all_verts = self.mb.get_entities_by_type(ms, types.MBVERTEX)

        while len(all_verts) > 0:
            # get full set of connected verts starting from a seed
            verts = [all_verts[0]]
            verts_check = [all_verts[0]]
            vtmp_all = set(verts[:])

            # gather set of all vertices that are connected to the seed
            while True:
                # check adjancency and connectedness of new vertices
                vtmp = self.mb.get_adjacencies(self.mb.get_adjacencies(
                                               verts_check, 2, op_type=1),
                                               0, op_type=1)

                # add newly found verts to all list
                vtmp_all.update(set(vtmp))

                # check if different from already found verts
                if len(list(vtmp_all)) == len(verts):
                    # no more vertices are connected, so full surface
                    # has been found
                    break
                else:
                    # update vertices list to check only newly found
                    # vertices
                    verts_check = vtmp_all.difference(verts)
                    verts = list(vtmp_all)

            # get the connected set of triangles that make the single
            # surface and store into a unique meshset
            tris = self.__get_surf_triangles(verts)
            surf = self.mb.create_meshset()
            self.mb.add_entities(surf, tris)
            self.mb.add_entities(surf, verts)

            # store surfaces in completed list
            surfs_list.append(surf)

            # remove surface from original meshset
            self.mb.remove_entities(ms, tris)
            self.mb.remove_entities(ms, verts)

            # resassign vertices that remain
            all_verts = self.mb.get_entities_by_type(ms, types.MBVERTEX)

        return surf_list

    def __get_surf_triangles(self, verts_good):
        """This function will take a set of vertex entity handles and
        return the set of triangles for which all vertices of all
        triangles are in the set of vertices.

        Input:
        ------
            verts_good: list of entity handles, list of vertices to
                compare against. Only triangles will be returned whose
                complete set of vertices are in this list.

        Returns:
        --------
            tris: list of entity handles, EHs for the triangles for
                which all three vertices are in the verts list.
        """
        tris_all = self.mb.get_adjacencies(verts_good, 2, op_type=1)
        verts_all = self.mb.get_connectivity(tris_all)
        verts_bad = set(verts_all) - set(verts_good)

        if verts_bad:
            # not an empty set
            tris_bad = self.mb.get_adjacencies(list(verts_bad), 2, op_type=1)
            tris_good = set(tris_all) - set(tris_bad)
            return list(tris_good)
        else:
            # empty set so all tris are good
            return tris_all

    def __list_coords(self, eh):
        """Gets list of all coords as a list of tuples for an entity
        handle eh.

        Input:
        ------
            eh: MOAB entity handle for meshset to retrieve coordinates

        Returns:
        --------
            coords: dictionary, key is the MOAB entity handle for the
                vertice and the value is a tuple of the coordinate
                (x, y, z).
        """
        # list of all entity handles for all vertices
        all_verts_eh = self.mb.get_entities_by_type(eh, types.MBVERTEX)
        coords = {}
        for v in all_verts_eh:
            coord = tuple(self.mb.get_coords(v))
            coords[v] = coord
        return coords

    def __compare_surfs(self, v1, v2, norm):
        """finds coincident surfaces between two isovolumes.

        Input:
        ------
            v1/2: tuple, corresponds to the dictionary keys for two
                isovolumes in self.isovol_meshsets that will be compared
            norm: float, All data values will be multiplied by this factor.
        """
        print("comparing surfaces in isovolumes {} and {}.".format(
            v1[0], v2[0]))

        # tag to indicate if interior or exterior surface
        surf_type_tag = \
            self.mb.tag_get_handle('SURF_TYPE', size=32,
                                   tag_type=types.MB_TYPE_OPAQUE,
                                   storage_type=types.MB_TAG_SPARSE,
                                   create_if_missing=False)

        # store dict of matched surfaces
        #   key: surf to remove in v2
        #   value: surf in v1 to replace it with
        surfs_to_remove = {}

        # compare all interior surfaces in v1 (s1) to all
        # interior surfaces in v2 (s2)
        for s1 in self.isovol_meshsets[v1]['surfs_EH']:

            # check if s1 is interior surf
            surf_type = self.mb.tag_get_data(surf_type_tag, s1)
            if surf_type != 'interior':
                continue

            # get list of all coordinates in s1
            verts1 = self.__list_coords(s1)

            for s2 in self.isovol_meshsets[v2]['surfs_EH']:
                # check surf type
                surf_type = self.mb.tag_get_data(surf_type_tag, s2)
                if surf_type != 'interior':
                    continue

                # get surface 2 vertices
                verts2 = self.__list_coords(s2)

                # check if any vertex matches between the two surfaces
                match = False
                for coord1 in verts1.values():
                    if coord2 in verts2.values():
                        # exact match, don't need to check anymore
                        match = True
                        break

                if match:
                    # match was found so s1 and s2 are coincident
                    # delete s2 (all triangles and all vertices that are
                    # not on the outer curve)
                    tris2 = self.__get_surf_triangles(verts2.keys())
                    self.mb.delete_entities(tris2)
                    verts_save = sk.find_skin(s2, tris2, True, False)
                    verts_delete = list(set(verts2.keys()) - set(verts_save))
                    self.mb.delete_entities(verts_delete)
                    surfs_to_remove[s2] = s1

        # remove the matched surfaces from volume 2
        # assign sense and value tags
        for s2, s1 in surfs_to_remove.items():
            self.isovol_meshsets[v2]['surfs_EH'].remove(s2)
            self.isovol_meshsets[v2]['surfs_EH'].append(s1)

            # assign sense tag to surface
            # [forward=v1, backward=v2]
            fwd = v1[1]
            bwd = v2[1]
            self.mb.tag_set_data(self.sense_tag, s1, [fwd, bwd])

            # tag the new surface with the shared value
            shared = \
                list(set(self.isovol_meshsets[v1]['bounds']) &
                     set(self.isovol_meshsets[v2]['bounds']))
            if not(bool(shared)):
                warnings.warn("No matching value for volumes " +
                              "{} and {}".format(v1, v2))
                val = 0.0
            else:
                val = shared[0] * norm
            self.mb.tag_set_data(self.val_tag, s1, val)
