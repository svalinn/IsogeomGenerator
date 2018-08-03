import sys
import os
import shutil
import numpy as np
import math as m

import visit as v
from pymoab import core, types


class IsoVolume(object):
    """This class contains methods to create a DAGMC geometry of
    isovolumes given any Cartesian mesh with tagged data.

    Users should follow the following sequence to completely build a
    geometry file:
        (1) Set contour levels with assign_levels() or generate_levels()
        (2) Generate the isovolume files using generate_volumes()
        (3) Create the MOAB geometry with create_geometry()
        (4) Write to a file with write_geometry()
    """

    def __init__(self):
        pass


    def assign_levels(self, levels):
        """User defines the contour levels to be used in the isovolumes.

        Input:
        ------
            levels: list of floats, list of user-defined values to use
                for contour levels
        """
        # make sure values are floats
        levels = [float(i) for i in levels]
        self.levels = sorted(levels)
        self.minN = min(self.levels)
        self.maxN = max(self.levels)
        self.N = len(self.levels)


    def generate_levels(self, N, minN, maxN, log=True):
        """Auto-generate evenly-spaced level values to use given a
        user-defined maximum, minimum, and number of levels to use.

        Input:
        ------
            N: int, number of levels generate
            minN: float, minimum level value
            maxN: float, maximum level value
            log: bool (optional), True to generate evenly spaced levels
                on a log scale (default), False to use linear scale.
        """

        # set min/max values
        self.minN = minN
        self.maxN = maxN
        self.N = N

        # generate evenly spaced values
        if log:
            base = 10.
            start = m.log(self.minN, base)
            stop = m.log(self.maxN, base)
            self.levels = list(np.logspace(start, stop, num=self.N,
                                      endpoint=True, base=base))
        else:
            self.levels = list(np.linspace(self.minN, self.maxN,
                                      num=self.N, endpoint=True))


    def _plot_pseudocolor(self):
        """Plots the data on a pseudocolor plot to use."""

        # add the pseudocolor plot to contour
        v.AddPlot("Pseudocolor", self.data)
        att = v.PseudocolorAttributes()

        # min/max for the pseudocolor plot
        att.minFlag = True
        att.min = self.minN
        att.maxFlag = True
        att.max = self.maxN

        # plot
        v.SetPlotOptions(att)
        v.DrawPlots()


    def _get_isovol(self, lbound, ubound, i):
        """Gets the volume selection for isovolume and export just the
        outer surface of the volume as STL.

        Input:
        ------
            lbound: float, lower boundary value for the isovolume
            ubound: float, upper boundary value for the isovolume
            i: int, surface number
        """

        # generate isovolume
        v.AddOperator("Isovolume")
        att = v.IsovolumeAttributes()
        att.lbound = lbound
        att.ubound = ubound
        v.SetOperatorOptions(att)

        # set operator setting to only get surfaces meshes
        v.AddOperator("ExternalSurface")

        # draw plot
        v.DrawPlots()

        # export current volume to folder
        e = v.ExportDBAttributes()
        e.dirname = self.db + "/vols/"
        e.db_type = "STL"
        e.filename = str(i)
        e.variables = self.data
        v.ExportDatabase(e)

        # delete the operators
        v.RemoveAllOperators()


    def _generate_vols(self):
        """Generates the isosurface volumes between the contour levels.
        Data files are exported as STLs and saved in the folder dbname.
        Files will be named based on their index corresponding to their
        level values (0.stl is lowest values).
        """

        # create folder to store data if it does not already exist
        if not os.path.isdir(self.db):
            os.mkdir(self.db)
        if os.path.isdir(self.db + "/vols/"):
            # make sure folder is empty by removing it first
            shutil.rmtree(self.db + "/vols/")
        os.mkdir(self.db + "/vols/")

        # plot the pseudocolor data inorder to get volumes
        self._plot_pseudocolor()

        # get the minimum isovolume level
        lbound = 0.0
        ubound = self.levels[0]
        self._get_isovol(lbound, ubound, 0)

        # iterate over all isovolume levels
        for l in self.levels[1:]:

            # get index of current level
            i = self.levels.index(l)

            # assign bounds
            lbound = self.levels[i-1]
            ubound = l

            # get volume
            self._get_isovol(lbound, ubound, i)

        # get maximum isovolume level
        lbound = self.levels[-1]
        ubound = 1.e50
        self._get_isovol(lbound, ubound, i+1)

        # delete plots
        v.DeleteAllPlots()


    def generate_volumes(self, filename, data,
                            dbname=os.getcwd()+"/tmp"):
        """Creates an STL file for each isovolume. N+1 files are
        generated and stored in the dbname folder.

        Input:
        ------
            filename: string, path to vtk file with the mesh
            data: string, name of the data whose values exist on the
                mesh (will be used to generate isocontours and
                isovolumes)
            dbname: (optional), string, name of folder to store created
                surface files. Must be absolute path!
                default: a folder called 'tmp' in the current directory
        """

        self.data = data
        self.db = dbname

        # make sure levels have been set before proceding
        try:
            self.levels
        except:
            print("ERROR: No contour levels have been set. " +\
                "Please use assign_levels or generate_levels to set.")
            sys.exit()

        # Generate isovolumes using VisIT
        v.LaunchNowin()
        v.OpenDatabase(filename)
        print("Generating isovolumes...")
        self._generate_vols()
        print("...Isovolumes files generated!")
        v.Close()


    def _separate(self, iv_info):
        """Separates a given surface into separate surfaces. All
        resulting surfaces are disjoint surfaces that made up the
        original surface.

        Input:
        ------
            iv_info: tuple, (iso_id, fs) where iso_id is the name of the
                loaded isovolume file (without the extension) and fs is
                a MOAB EntityHandle corresponding to the file_set for
                the loaded isovolume file.
        """

        # extract isovolume information
        iso_id = iv_info[0]
        fs = iv_info[1]

        # get set of all vertices for the isosurface
        all_verts = self.mb.get_entities_by_type(fs, types.MBVERTEX)

        # initiate list to store separate surface entity handles
        self.isovol_meshsets[iv_info]['surfs_EH'] = []

        # separate the surfaces
        print("separating isovolume {}".format(iso_id))
        while len(all_verts) > 0:
            # get full set of connected verts starting from a seed
            verts = [all_verts[0]]

            # gather set of all vertices that are connected to the seed
            while True:
                # this step takes too long for large surfaces
                # check adjancency and connectedness of vertices
                vtmp = self.mb.get_adjacencies(self.mb.get_adjacencies(
                                               verts, 2, op_type=1),
                                               0, op_type=1)
                if len(vtmp) == len(verts):
                    # no more vertices are connected, so full surface
                    # has been found
                    break
                else:
                    # update vertices list to include newly found
                    # connected vertices
                    verts = vtmp

            # get the connected set of triangles that make the single
            # surface and store into a unique meshset
            tris = self.mb.get_adjacencies(verts, 2, op_type=1)
            surf = self.mb.create_meshset()
            self.mb.add_entities(surf, tris)
            self.mb.add_entities(surf, verts)

            self.isovol_meshsets[iv_info]['surfs_EH'].append(surf)

            # remove surface from original meshset
            self.mb.remove_entities(fs, tris)
            self.mb.remove_entities(fs, verts)

            # resassign vertices that remain
            all_verts = self.mb.get_entities_by_type(fs, types.MBVERTEX)


    def _separate_isovols(self):
        """For each isovolume in the database, separate any disjoint
        surfaces into unique single surfaces.
        """

        for f in sorted(os.listdir(self.db + "/vols/")):
            # get file name
            fpath = self.db + "/vols/" + f
            i = int(f.strip(".stl")) # must be an integer

            # load file and create EH for file-set
            fs = self.mb.create_meshset()
            self.mb.load_file(fpath, file_set=fs)

            # initiate dictionary
            iv_info = (i, fs)
            self.isovol_meshsets[iv_info] = {}

            # separate
            self._separate(iv_info)

            # add value min/max info (min, max)
            if i == 0:
                self.isovol_meshsets[iv_info]['bounds'] =\
                    (None, self.levels[i])
            elif i == self.N:
                self.isovol_meshsets[iv_info]['bounds'] =\
                    (self.levels[i-1], None)
            else:
                self.isovol_meshsets[iv_info]['bounds'] =\
                    (self.levels[i-1], self.levels[i])


    def _list_coords(self, eh):
        """Gets list of all coords as a list of tuples for an entity
        handle eh.

        Input:
        ------
            eh: MOAB entity handle for meshset to retrieve coordinates

        Returns:
        --------
            coords: dictionary, key is the MOAB entity handle for the
                vertice and the value is a tuple of the coordinate
                (x, y, z)
        """

        # list of all entity handles for all vertices
        all_verts_eh = self.mb.get_entities_by_type(eh, types.MBVERTEX)
        coords = {}
        for v in all_verts_eh:
            coords[v] = tuple(self.mb.get_coords(v))

        return coords


    def _get_matches(self, vertsA, vertsB):
        """Collects the set of entity handles in set of vertsA and their
        coordinates that exist in vertsB.

        Input:
        ------
            vertsA/B: dictionary, key is the MOAB entity handle for the
                vertice and the value is a tuple of the coordinate
                (x, y, z)

        Returns:
        --------
            sA_match_eh: list of MOAB entity handles, the entity handles
                for set vertsA that exist is vertsB
            sA_match_coords: list of tuples, each entry is the
                corresponding coordinate for the EH in sA_match_eh
        """

        sA_match_eh = []
        sA_match_coords = []
        for vert in vertsA.items():
            eh = vert[0]
            coord = vert[1]
            if coord in vertsB.values():
                sA_match_eh.append(eh)
                sA_match_coords.append(coord)

        return sA_match_eh, sA_match_coords


    def _compare_surfs(self, v1, v2):
        """finds coincident surfaces between two isovolumes.

        Input:
        ------
            v1/2: tuple, corresponds to the dictionary keys for two
                isovolumes in self.isovol_meshsets that will be compared
        """

        print("comparing surfaces in isovolumes {} and {}.".format(
            v1[0], v2[0]))

        match_surfs = []

        # compare all surfaces in v1 (s1) to all surfaces in v2 (s2)
        for s1 in self.isovol_meshsets[v1]['surfs_EH']:
            # get list of all coordinates in s1
            verts1 = self._list_coords(s1)

            for s2 in self.isovol_meshsets[v2]['surfs_EH']:

                # get list of all coordinates in s2
                verts2 = self._list_coords(s2)

                # compare vertices and gather sets for s1 and s2
                # that are coincident
                s1_match_eh, s1_match_coords = self._get_matches(verts1,
                                                                 verts2)

                if s1_match_eh != []:
                    # matches were found, so continue

                    # must also collect the corresponding entity handles
                    # for s2 so they can be properly updated
                    s2_match_eh, s2_match_coords = self._get_matches(
                                                    verts2, verts1)

                    # check that the set of coordinates match for each
                    if set(s1_match_coords) != set(s2_match_coords):
                        print("Sets of coincident coords do not match!")

                    # create new coincident surface
                    tris1 = self.mb.get_adjacencies(s1_match_eh, 2,
                                                    op_type=1)
                    surf = self.mb.create_meshset()
                    self.mb.add_entities(surf, tris1)
                    self.mb.add_entities(surf, s1_match_eh)

                    # assign sense tag to surface
                    # [forward=v1, backward=v2]
                    fwd = v1[1]
                    bwd = v2[1]
                    self.mb.tag_set_data(self.sense_tag, surf,
                                            [fwd, bwd])

                    # get s2 tris to delete (no new surface needed)
                    tris2 = self.mb.get_adjacencies(s2_match_eh, 2,
                                                    op_type=1)

                    # delete verts/tris from original surfaces
                    self.mb.remove_entities(s1, tris1)
                    self.mb.remove_entities(s1, s1_match_eh)
                    self.mb.remove_entities(s2, tris1)
                    self.mb.remove_entities(s2, s2_match_eh)

                    # tag the new surface with the shared value
                    shared = \
                        list(set(self.isovol_meshsets[v1]['bounds']) \
                        & set(self.isovol_meshsets[v2]['bounds']))
                    if not(bool(shared)):
                        print('no matching value!', v1, v2)
                        val = 0.0
                    else:
                        val = shared[0]

                    self.mb.tag_set_data(self.val_tag, surf, val)

                    # add new surface to coincident surface list
                    match_surfs.append(surf)

                    # check if original surfaces are empty (no vertices)
                    # if so delete empty meshset and remove from list
                    s2_remaining = self.mb.get_entities_by_type(s2,
                                                        types.MBVERTEX)
                    if len(s2_remaining) == 0:
                        # delete surface from list and mb instance
                        self.isovol_meshsets[v2]['surfs_EH'].remove(s2)

                    s1_remaining = self.mb.get_entities_by_type(s1,
                                                        types.MBVERTEX)
                    if len(s1_remaining) == 0:
                        # delete from list and mb instance and move to
                        # next surf
                        self.isovol_meshsets[v1]['surfs_EH'].remove(s1)
                        break

        # After all comparisons have been made, add surfaces to lists
        self.isovol_meshsets[v1]['surfs_EH'].extend(match_surfs)
        self.isovol_meshsets[v2]['surfs_EH'].extend(match_surfs)


    def _imprint_merge(self):
        """Uses PyMOAB to check if surfaces are coincident. Creates a
        single surface where surfaces are coincident values are tagged
        on each surface. Surface senses are also determined and tagged.
        """

        # set up surface tag information (value and sense)
        self.val_tag = self.mb.tag_get_handle(self.data, size=1,
                        tag_type=types.MB_TYPE_DOUBLE,
                        storage_type=types.MB_TAG_SPARSE,
                        create_if_missing=True)
        self.sense_tag = self.mb.tag_get_handle('GEOM_SENSE_2', size=2,
                        tag_type=types.MB_TYPE_HANDLE,
                        storage_type=types.MB_TAG_SPARSE,
                        create_if_missing=True)

        # get list of all original isovolumes
        all_vols = sorted(self.isovol_meshsets.keys())
        for i, isovol in enumerate(all_vols):

            if i != self.N:
                # do not need to check the last isovolume because it
                # will be checked against its neighbor already
                self._compare_surfs(isovol, all_vols[i+1])

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

                # tag fwd sense
                try:
                    sense = self.mb.tag_get_data(self.sense_tag, surf)
                except:
                    fwd = isovol[1]
                    bwd = np.uint64(0) # this needs to be a valid EH... fix TBD
                    self.mb.tag_set_data(self.sense_tag, surf, [fwd, bwd])


    def _make_family(self):
        """Makes the correct parent-child relationships with volumes
        and surfaces. Tags geometry type and category on surfaces and
        volumes.
        """
        # create geometry dimension and category tags
        geom_dim = self.mb.tag_get_handle('GEOM_DIMENSION', size=1,
                        tag_type=types.MB_TYPE_INTEGER,
                        storage_type=types.MB_TAG_SPARSE,
                        create_if_missing=True)
        category = self.mb.tag_get_handle('CATEGORY', size=32,
                        tag_type=types.MB_TYPE_OPAQUE,
                        storage_type=types.MB_TAG_SPARSE,
                        create_if_missing=True)

        for v in self.isovol_meshsets.keys():
            vol_eh = v[1]

            # tag volume
            self.mb.tag_set_data(geom_dim, vol_eh, 3)
            self.mb.tag_set_data(category, vol_eh, 'Volume')

            for surf_eh in self.isovol_meshsets[v]['surfs_EH']:
                # create relationship
                self.mb.add_parent_child(vol_eh, surf_eh)

                # tag surfaces
                self.mb.tag_set_data(geom_dim, surf_eh, 2)
                self.mb.tag_set_data(category, surf_eh, 'Surface')


    def _tag_groups(self):
        """Tags the surfaces with metadata as a group with value
            {data}_{value}. Example: 'wwn_0.005'. The original data
            name will be stripped of any underscores.
        """

        # create group entity set per data value & tag as a 'Group'

        # category tag
        category = self.mb.tag_get_handle('CATEGORY', size=32,
                        tag_type=types.MB_TYPE_OPAQUE,
                        storage_type=types.MB_TAG_SPARSE,
                        create_if_missing=True)

        # tag for group metadata
        tag_name = self.mb.tag_get_handle('NAME', size=32,
                                tag_type=types.MB_TYPE_OPAQUE,
                                storage_type=types.MB_TAG_SPARSE,
                                create_if_missing=True)

        # strip underscores from base data name
        data = self.data.replace('_', '')

        # create meshsets (starting w/ 0.0)
        data_groups = {}

        data_groups[0.0] = self.mb.create_meshset()
        self.mb.tag_set_data(category, data_group[0.0], 'Group')
        name = '{}_{}'.format(data, 0.0)
        self.mb.tag_set_data(tag_name, data_group[0.0], name)

        for val in self.levels:
            data_groups[val] = self.mb.create_meshset()
            self.mb.tag_set_data(category, data_group[val], 'Group')
            name = '{}_{}'.format(data, val)
            self.mb.tag_set_data(tag_name, data_group[val], name)

        # add surfs to groups
        for isovol in self.isovol_meshsets.keys():
            for surf in self.isovol_meshsets[isovol]:
                # get the tagged data
                val_data = tag_get_data(self.val_tag, surf)

                # add to group with that same data
                self.add_entities(data_groups[val_data], surf)


    def _tag_for_viz(self):
        """Tags all triangles on all surfaces with the data value for
        that surface. This is for vizualization purposes.
        """
        for isovol in self.isovol_meshsets.keys():
            for surf in self.isovol_meshsets[isovol]:
                # get the tagged data
                val = tag_get_data(self.val_tag, surf)

                # get the triangles
                verts = self.mb.get_entities_by_type(surf,
                            types.MBVERTEX, recur=True)
                tris = self.mb.get_adjacencies(verts, 2, op_type=1)

                # create data array
                num = len(tris)
                data = np.full((1, num), val)

                # tag the data
                self.mb.tag_set_data(self.val_tag, tris, data)


    def create_geometry(self, tag_groups=False, tag_for_viz=False):
        """Over-arching function to do all steps to create a single
        isovolume geometry for DAGMC.

        Input:
        ------
            tag_groups: bool (optional), True to tag surfaces in groups
                with NAMES '{data}_{value}' where data is the data name
                and value is the value for that surface. Default=False.
            tag_for_viz: bool (optional), True to tag each triangle on
                every surface with the data value. Needed to visualize
                values in VisIt. Default=False.
        """

        # check that database is identified
        try:
            self.db
        except:
            print("ERROR: Database not found. " +\
                "Please run generate_volumes first.")
            sys.exit()

        # Step 1: Separate Isovolume Surfaces
        self.mb = core.Core()
        self.isovol_meshsets = {}
        print("Separating isovolumes...")
        self._separate_isovols()
        print("...Separation complete!")

        # Step 2: Merge Coincident Surfaces
        print("Merging surfaces...")
        self._imprint_merge()
        print("...Merging complete!")

        # Step 3: Assign Parent-Child Relationship
        self._make_family()

        if tag_groups:
            print('Tagging data groups...')
            self._tag_groups()
            print('... Data groups complete')

        if tag_for_viz:



    def write_geometry(self, sname="", sdir=""):
        """Writes out the geometry stored in memory.

        Input:
        ------
            sname: (optional), string, name of file to save written file
                default: 'geom-DATA.h5m' where DATA is the name of the
                data used to create the geometry (if known).
        """

        # check that there is a geometry in memory
        try:
            self.mb
        except:
            print("ERROR: No geometry in memory to write to file. " + \
                "Please run create_geometry first.")
            sys.exit()

        # create default filename if unassigned.
        if sname == "":
            try:
                self.data
            except:
                print("WARNING: Data name is unknown! " + \
                    "File will be saved as geom-isovols.h5m.")
                sname="geom-isovols.h5m"
            else:
                sname="geom-{}.h5m".format(self.data)

        if sdir == "":
            try:
                self.db
            except:
                print("WARNING: Database location is unknown! " + \
                    "File will be saved in tmp/ folder.")
                sdir = os.getcwd() + "/tmp"
            else:
                sdir = self.db

        # check file extension of save name:
        ext = sname.split(".")[-1]
        if ext.lower() not in ['h5m', 'vtk']:
            print("WARNING: File extension {} ".format(ext) +\
                " not recognized. File will be saved as type .h5m.")
            self.sname = sname.split(".")[0] + ".h5m"

        # save the file
        save_location = sdir + "/" + sname
        self.mb.write_file(save_location)

        print("Geometry file written to {}.".format(save_location))
