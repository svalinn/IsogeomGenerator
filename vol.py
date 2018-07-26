import visit as v
import sys
import os
import shutil
import numpy as np
import math as m
from pymoab import core, types
from pymoab.rng import Range


class IsoVolume(object):
    """This class contains methods to create a DAGMC geometry of
    isovolumes given any cartesian mesh with tagged data.

    Parameters:
    -----------
        filename: string, path to vtk file with the mesh
        data: string, name of the data whose values exist on the mesh
            (will be used to generate isocontours and isovolumes)
        dbname: (optional), string, name of folder to store created
            surface files. Must be absolute path!
    """

    def __init__(self, filename, data, dbname=os.getcwd() + "/tmp"):
        self.data = data
        self.f = filename
        self.db = dbname


    def assign_levels(self, levels):
        """User defines the contour levels to be used in the isovolumes.

        Input:
        ------
            levels: list of floats, list of user-defined values to use
                for contour levels
        """
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


    def _generate_volumes(self):
        """Generates the isosurface volumes between the contour levels.
        Data files are exported as STLs and saved in the folder dbname.
        Files will be named based on their index corresponding to their
        level values (0.stl is lowest values).
        """

        # make sure levels have been set before proceding
        if self.levels is None or []:
            print("No contour levels have been set!")
            print("Please use assign_levels or generate_levels to set levels.")
            return

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

            # add ww min/max info (min, max)
            if i == 0:
                self.isovol_meshsets[iv_info]['ww_bounds'] = (None, self.levels[i])
            elif i == self.N:
                self.isovol_meshsets[iv_info]['ww_bounds'] = (self.levels[i-1], None)
            else:
                self.isovol_meshsets[iv_info]['ww_bounds'] = (self.levels[i-1], self.levels[i])


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

            v1/2: tuple, corresponds to the dictionary keys for two
                isovolumes in self.isovol_meshsets
        """

        print("comparing surfaces in isovolumes {} and {}.".format(v1[0], v2[0]))

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
                s1_match_eh, s1_match_coords = self._get_matches(verts1, verts2)

                if s1_match_eh != []:
                    # matches were found, so continue

                    # must also collect the corresponding entity handles for
                    # s2 so they can be properly updated
                    s2_match_eh, s2_match_coords = self._get_matches(verts2, verts1)

                    # check that the set of coordinates match for each
                    if set(s1_match_coords) != set(s2_match_coords):
                        print("Sets of coincident coords do not match!!")

                    # create new coincident surface
                    tris1 = self.mb.get_adjacencies(s1_match_eh, 2, op_type=1)
                    surf = self.mb.create_meshset()
                    self.mb.add_entities(surf, tris1)
                    self.mb.add_entities(surf, s2_match_eh)

                    # get s2 tris to delete (no new surface needed)
                    tris2 = self.mb.get_adjacencies(s2_match_eh, 2, op_type=1)

                    # delete verts/tris from original surfaces
                    self.mb.remove_entities(s1, tris1)
                    self.mb.remove_entities(s1, s1_match_eh)
                    self.mb.remove_entities(s2, tris1)
                    self.mb.remove_entities(s2, s2_match_eh)

                    # tag the new surface with the shared WW value
                    shared_ww = list(set(self.isovol_meshsets[v1]['ww_bounds']) & set(self.isovol_meshsets[v2]['ww_bounds']))
                    if not(bool(shared_ww)):
                        print('no matching ww value!', v1, v2)
                        ww_val = 0.0
                    else:
                        ww_val = float(shared_ww[0])

                    #ww_tag = self.mb.tag_get_handle('ww', size=1, tag_type=types.MB_TYPE_DOUBLE,
                    #        storage_type=types.MB_TAG_DENSE, create_if_missing=True)
                    self.mb.tag_set_data(self.ww_tag, surf, ww_val)

                    # add new surface to coincident surface list
                    match_surfs.append(surf)

                    # check if original surfaces are empty (no vertices)
                    # if so delete empty surf meshset and remove from list
                    s2_remaining = self.mb.get_entities_by_type(s2, types.MBVERTEX)
                    if len(s2_remaining) == 0:
                        # delete surface from list and mb instance
                        #self.mb.delete_entities(s2)
                        self.isovol_meshsets[v2]['surfs_EH'].remove(s2)

                    s1_remaining = self.mb.get_entities_by_type(s1, types.MBVERTEX)
                    if len(s1_remaining) == 0:
                        # delete from list and mb instance and move to next surf
                        #self.mb.delete_entities(s1)
                        self.isovol_meshsets[v1]['surfs_EH'].remove(s1)
                        break

        # After all comparisons have been made, add surfaces to lists
        self.isovol_meshsets[v1]['surfs_EH'].extend(match_surfs)
        self.isovol_meshsets[v2]['surfs_EH'].extend(match_surfs)


    def _imprint_merge(self):
        """Uses PyMOAB to check if surfaces are coincident. Creates a
        single surface where surfaces are coincident. Weight-Window
        values are tagged on each surface.
        """
        # set up weight window tag information
        self.ww_tag = self.mb.tag_get_handle('ww', size=1, tag_type=types.MB_TYPE_DOUBLE,
                        storage_type=types.MB_TAG_DENSE, create_if_missing=True)

        # get list of all original isovolumes
        all_vols = sorted(self.isovol_meshsets.keys())
        for i, isovol in enumerate(all_vols):

            if i != self.N:
                # do not need to check the last isovolume because it
                # will be checked against its neighbor already
                self._compare_surfs(isovol, all_vols[i+1])

        # if a surface doesn't have a WW tagged value after merging
        # give it a ww value of 0
        for isovol in all_vols:
            for surf in self.isovol_meshsets[isovol]['surfs_EH']:
                # determine if surf already tagged with ww data,
                # if not then tag with 0
                try:
                    ww_val = self.mb.tag_get_data(self.ww_tag, surf)
                except:
                    ww_val = 0.0
                    self.mb.tag_set_data(self.ww_tag, surf, ww_val)


    def create_geometry(self):
        """Over-arching function to do all steps to create a single
        isovolume geometry for DAGMC.
        """

        ###########################################
        # Step 1: Generate Isovolumes using VisIT #
        ###########################################
        v.LaunchNowin()
        v.OpenDatabase(self.f)
        print("Generating isovolumes...")
        self._generate_volumes()
        print("...Isovolumes files generated!")
        v.Close()

        #######################################
        # Step 2: Separate Isovolume Surfaces #
        #######################################
        self.mb = core.Core()
        self.isovol_meshsets = {}
        print("Separating isovolumes...")
        self._separate_isovols()
        print("...Separation complete!")

        #####################################
        # Step 3: Merge Coincident Surfaces #
        #####################################
        print("Merging surfaces...")
        self._imprint_merge()
        print("...Merging complete!")

        # Step 4: create parent-child meshsets for each isovolume

        # Step 5: write out as single h5m
