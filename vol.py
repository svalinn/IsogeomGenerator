import sys
import os
import shutil
import numpy as np
import math as m

import visit as v
from pymoab import core, types
from pymoab import rng
from pymoab.rng import Range
from pymoab.skinner import Skinner

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


    def generate_levels(self, N, minN, maxN, log=True, ratio=False):
        """Auto-generate evenly-spaced level values to use given a
        user-defined maximum, minimum, and number of levels or set
        constant ratio for spacing.

        Input:
        ------
            N: int or float, number of levels (int) to generate
                (ratio=False); or the ratio (float) to use to separate
                levels (ratio=True).
            minN: float, minimum level value
            maxN: float, maximum level value
            log: bool (optional), True to generate evenly spaced levels
                on a log scale (default), False to use linear scale.
            ratio: bool (optional), True to generate levels that are
                spaced by some constant ratio N. If True, log input is
                ignored. If True, minN will be used as minimum level
                value, maximum level value is less than or equal to
                maxN.
        """

        # set min value
        self.minN = minN

        if not ratio:
            # ratio not being used
            # set max value
            self.maxN = maxN
            self.N = int(N)

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

        else:
            # ratio being used
            # set minN as the minimum value and get all other values
            # until maxN
            self.maxN = 0.
            self.levels = [self.minN]

            while self.maxN < maxN:
                next_val = self.levels[-1]*float(N)
                if next_val <= maxN:
                    self.levels.append(next_val)
                    self.maxN = next_val
                else:
                    break

            self.N = len(self.levels)


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
        try:
            v.LaunchNowin()
        except:
            print("VisIt already launched.")
        v.OpenDatabase(filename)
        print("Generating isovolumes...")
        self._generate_vols()
        print("...Isovolumes files generated!")
        v.CloseComputeEngine()


    def create_geometry(self, elower, eupper, tag_groups=False, tag_for_viz=False, norm=1.0,
                        merge_tol=1e-5, facet_tol=1e-3):
        """Over-arching function to do all steps to create a single
        isovolume geometry for DAGMC.

        Input:
        ------
            e_lower: double, energy group lower bound
            e_upper: double, energy group upper bound
            tag_groups: bool (optional), True to tag surfaces in groups
                with NAMES '{data}_{value}' where data is the data name
                and value is the value for that surface. Default=False.
            tag_for_viz: bool (optional), True to tag each triangle on
                every surface with the data value. Needed to visualize
                values in VisIt. Default=False.
            norm: float (optional), default=1. All ww values will be
                multiplied by the normalization factor.
            merge_tol: float (option), default=1e-5 cm. Merge tolerance for mesh
                based merge of coincident surfaces. Recommended to be
                1/10th the mesh voxel size.
        """

        self.norm = norm
        self.merge_tol = merge_tol
        self.facet_tol = facet_tol
        self.el = e_lower
        self.eu = e_upper

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
        #self._get_skins()
        self._make_family()

        if tag_groups:
            print('Tagging data groups...')
            self._tag_groups()
            print('... tags complete')

        if tag_for_viz:
            print('Tagging triangles with data...')
            self._tag_for_viz()
            print('... tags complete')

        # add void materials
        self._add_mats()

        # tag energy bounds on the root set
        el_tag = self.mb.tag_get_handle('E_LOW_BOUND', size=1,
                        tag_type=types.MB_TYPE_DOUBLE,
                        storage_type=types.MB_TAG_SPARSE,
                        create_if_missing=True)
        ep_tag = self.mb.tag_get_handle('E_UP_BOUND', size=1,
                        tag_type=types.MB_TYPE_DOUBLE,
                        storage_type=types.MB_TAG_SPARSE,
                        create_if_missing=True)



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

        # write out only the ranges stored
        self.mb.write_file(save_location)

        print("Geometry file written to {}.".format(save_location))


    ################################################################
    ################# Extra functions for VisIT step ###############
    ################################################################

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


    def _update_levels(self, value):
        """Removes a value from the levels list and resets N.

        Input:
        ------
            value: float, value to remove
        """

        self.levels.remove(value)
        self.N = len(self.levels)


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
        draw_res = v.DrawPlots()
        if draw_res == 0:
            sys.exit("Error with producing isovolume")

        # export current volume to folder
        e = v.ExportDBAttributes()
        e.dirname = self.db + "/vols/"
        e.db_type = "STL"
        e.filename = str(i)
        e.variables = self.data
        export_res = v.ExportDatabase(e)

        if export_res == 0:
            # export not successful because there was no data
            # get new upper bound
            warn_message = "Warning: no data to export between {} and {}.\n".format(lbound, ubound) \
                + "Increasing upper bound to next selected level."
            print(warn_message)
            if ubound in self.levels:
                self._update_levels(ubound)
            else:
                # it is the arbitrary upper level set and is not needed
                self._update_levels(self.levels[-1])

        # delete the operators
        v.RemoveAllOperators()

        return export_res, ubound


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
            res = 0
            while res == 0:
                # get index of current level
                i = self.levels.index(l)

                # assign bounds
                lbound = self.levels[i-1]
                ubound = l

                # get volume
                # res = 0 if no level found (should update to next level)
                res = self._get_isovol(lbound, ubound, i)

        # get maximum isovolume level
        lbound = self.levels[-1]
        ubound = 1.e200
        self._get_isovol(lbound, ubound, i+1)

        # delete plots
        v.DeleteAllPlots()


    ##############################################################
    ############## Functions for PyMOAB step #####################
    ##############################################################

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
            verts_check = [all_verts[0]]
            vtmp_all = set(verts[:])

            # gather set of all vertices that are connected to the seed
            while True:
                # this step takes too long for large surfaces
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
            tris = self._get_surf_triangles(verts)
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

        tol_set = False
        facet_tol_tag = self.mb.tag_get_handle('FACETING_TOL', size=1,
                        tag_type=types.MB_TYPE_DOUBLE,
                        storage_type=types.MB_TAG_SPARSE,
                        create_if_missing=True)
        resabs_tag = self.mb.tag_get_handle('GEOMETRY_RESABS', size=1,
                        tag_type=types.MB_TYPE_DOUBLE,
                        storage_type=types.MB_TAG_SPARSE,
                        create_if_missing=True)

        for f in sorted(os.listdir(self.db + "/vols/")):
            # get file name
            fpath = self.db + "/vols/" + f
            i = int(f.strip(".stl")) # must be an integer

            # load file and create EH for file-set
            fs = self.mb.create_meshset()
            self.mb.load_file(fpath, file_set=fs)

            if not tol_set:
                # only tag on one fileset
                self.mb.tag_set_data(facet_tol_tag, fs, self.facet_tol)
                self.mb.tag_set_data(resabs_tag, fs, self.merge_tol)
                tol_set = True

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


    def _list_coords(self, eh, invert=False):
        """Gets list of all coords as a list of tuples for an entity
        handle eh.

        Input:
        ------
            eh: MOAB entity handle for meshset to retrieve coordinates
            invert: bool, default=False, True to invert keys and values
                in returned coords dict

        Returns:
        --------
            coords: dictionary, key is the MOAB entity handle for the
                vertice and the value is a tuple of the coordinate
                (x, y, z). If invert=True, keys and values are switched.
        """

        # list of all entity handles for all vertices
        all_verts_eh = self.mb.get_entities_by_type(eh, types.MBVERTEX)
        coords = {}
        for v in all_verts_eh:
            coord = tuple(self.mb.get_coords(v))

            if invert:
                # invert is true
                key = coord
                value = v
            else:
                key = v
                value = coord

            coords[key] = value

        return coords


    def _get_matches(self, vertsA, vertsB):
        """Collects the set of entity handles and coordinates in set of
        vertsA and vertsB that match within the specified absolute
        tolerance (self.merge_tol).

        Input:
        ------
            vertsA: dictionary, key is the MOAB entity handle for the
                vertice and the value is a tuple of the coordinate
                (x, y, z)
            vertsB: dictionary, key is a tuple of the coordinate and the
                value is the MOAB entity handle for the coordinate

        Returns:
        --------
            sA_match_eh: list of MOAB entity handles, the entity handles
                for set vertsA that exist is vertsB
            sA_match_coords: list of tuples, each entry is the
                corresponding coordinate for the EH in sA_match_eh
            sB_match_eh: list of MOAB entity handles, the entity handles
                for set vertsB that exist is vertsA
            sB_match_coords: list of tuples, each entry is the
                corresponding coordinate for the EH in sB_match_eh
        """

        sA_match_eh = []
        sA_match_coords = []
        sB_match_eh = []
        sB_match_coords = []

        match_dict = {}

        bcoords = vertsB.keys()

        # get exact matches
        for vert in vertsA.items():
            ehA = vert[0]
            coord = vert[1]
            if coord in bcoords:
                # exact match
                sA_match_eh.append(ehA)
                sA_match_coords.append(coord)
                sB_match_coords.append(coord)
                sB_match_eh.append(vertsB[coord])

                match_dict[vertsB[coord]] = ehA

            else:
                # check approx
                tf = np.isclose(coord, bcoords, rtol=0, atol=self.merge_tol)

                # get index of matches if they exist
                ix = np.where(zip(*tf)[0])[0]
                iy = np.where(zip(*tf)[1])[0]
                iz = np.where(zip(*tf)[2])[0]

                # get index if only x y and z match
                index_set = list(set(ix) & set(iy) & set(iz))

                if index_set != []:
                    index = index_set[0]

                    # get the close match coordinate in the bcoords list
                    bcoord = bcoords[index]
                    ehB = vertsB[bcoord]

                    sA_match_eh.append(ehA)
                    sA_match_coords.append(coord)
                    sB_match_eh.append(ehB)
                    sB_match_coords.append(bcoord)

                    match_dict[ehB] = ehA


        return sA_match_eh, sA_match_coords, sB_match_eh, sB_match_coords, match_dict


    def _get_surf_triangles(self, verts_good):
        """This function will take a set of vertice entity handles and
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
        sk = Skinner(self.mb)

        # compare all surfaces in v1 (s1) to all surfaces in v2 (s2)
        for s1 in self.isovol_meshsets[v1]['surfs_EH']:
            # get list of all coordinates in s1
            verts1 = self._list_coords(s1)

            # initialize list of curves
            if s1 not in self.surf_curve.keys():
                self.surf_curve[s1] = []

            for s2 in self.isovol_meshsets[v2]['surfs_EH']:
                if s2 not in self.surf_curve.keys():
                    self.surf_curve[s2] = []

                # get list of all coordinates in s2 (inverted)
                verts2_inv = self._list_coords(s2, invert=True)

                # compare vertices and gather sets for s1 and s2
                # that are coincident
                s1_match_eh, s1_match_coords, s2_match_eh, s2_match_coords, match_dict =\
                    self._get_matches(verts1, verts2_inv)

                if s1_match_eh != []:
                    # matches were found, so continue

                    # get only tris1 that have all match vertices
                    tris1 = self._get_surf_triangles(s1_match_eh)

                    # get s2 tris to delete (no new surface needed)
                    tris2 = self._get_surf_triangles(s2_match_eh)

                    # create new coincident surface
                    surf = self.mb.create_meshset()
                    self.mb.add_entities(surf, tris1)
                    self.mb.add_entities(surf, s1_match_eh)
                    self.surf_curve[surf] = []

                    # get skin of new merged surf (gets curve)

                    curve_verts = sk.find_skin(surf, tris1, True, False)
                    curve_edges = sk.find_skin(surf, tris1, False, False)

                    # if curve_verts/edges is empty, closed surf is created
                    # so no new curve is needed
                    if len(curve_verts) > 0:
                        # if not empty, make new curve
                        curve = self.mb.create_meshset()
                        self.mb.add_entities(curve, curve_verts)
                        self.mb.add_entities(curve, curve_edges)
                        self.surf_curve[s1].append(curve)
                        self.surf_curve[s2].append(curve)
                        self.surf_curve[surf].append(curve)

                        # remove merged verts and tris from each already existing surf
                        for vert_delete in s2_match_eh:

                            # get all triangles connected to the vert to be deleted
                            tris_adjust = self.mb.get_adjacencies(vert_delete, 2, op_type=1)

                            # get the vert that will replace the deleted vert
                            replacement = match_dict[vert_delete]

                            # for every tri to be deleted, replace vert by setting connectivity
                            for tri in tris_adjust:
                                tri_verts = self.mb.get_connectivity(tri)
                                new_verts = [0, 0, 0]
                                for i, tv in enumerate(tri_verts):
                                    if tv == vert_delete:
                                        new_verts[i] = replacement
                                    else:
                                        new_verts[i] = tv

                                # set connectivity
                                self.mb.set_connectivity(tri, new_verts)

                    # remove from both sets (already in new surface)
                    self.mb.remove_entities(s1, tris1)
                    self.mb.remove_entities(s1, s1_match_eh)
                    self.mb.remove_entities(s2, tris2)
                    self.mb.remove_entities(s2, s2_match_eh)

                    # delete surf 2 (repeats)
                    self.mb.delete_entities(tris2)
                    #self.mb.delete_entities(s2_match_eh)

                    # TAG INFORMATION

                    # assign sense tag to surface
                    # [forward=v1, backward=v2]
                    fwd = v1[1]
                    bwd = v2[1]
                    self.mb.tag_set_data(self.sense_tag, surf,
                                            [fwd, bwd])

                    # tag the new surface with the shared value
                    shared = \
                        list(set(self.isovol_meshsets[v1]['bounds']) \
                        & set(self.isovol_meshsets[v2]['bounds']))
                    if not(bool(shared)):
                        print('no matching value!', v1, v2)
                        val = 0.0
                    else:
                        val = shared[0]*self.norm

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

        # create dictionary of curves to match to surfaces:
        # key = surf eh, value = list of child curve eh
        self.surf_curve = {}

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
                    verts = self.mb.get_entities_by_type(surf,
                            types.MBVERTEX)
                    tris = self._get_surf_triangles(verts)
                    self.mb.add_entities(surf, tris)

                # tag fwd sense
                try:
                    sense = self.mb.tag_get_data(self.sense_tag, surf)
                except:
                    fwd = isovol[1]
                    bwd = np.uint64(0)
                    self.mb.tag_set_data(self.sense_tag,
                                        surf, [fwd, bwd])


    def _get_skins(self):
        """Get the curves from the skins of surfaces and their vertices
        to create as meshsets.
        """
        sk = Skinner(self.mb)

        vol_keys = self.isovol_meshsets.keys()

        # get skins of all surfaces (unmerged)
        for v1 in vol_keys:
            print("*!*!*!*! Volume {} {}".format(v1[0], v1[1]))
            self.isovol_meshsets[v1]['curves'] = {}
            for s1 in self.isovol_meshsets[v1]['surfs_EH']:
                # get skin of surface
                tris = self.mb.get_entities_by_type(s1, types.MBTRI)
                skin_verts = sk.find_skin(s1, tris, True, False)
                skin_edges = sk.find_skin(s1, tris, False, False)



                # create curve from skin
                # create meshset, add verts, add edges, tag w/ geom dimension
                if len(skin_verts) > 0:
                    print("!!!!!! surf eh: {}".format(s1))
                    print("verts: ", list(skin_verts))
                    print("edges: ", list(skin_edges))
                    curve = self.mb.create_meshset()
                    self.mb.add_entities(curve, skin_verts)
                    self.mb.add_entities(curve, skin_edges)

                    # save skin information
                    self.isovol_meshsets[v1]['curves'][s1] = [curve]
                else:
                    # no curves, give default value -1
                    self.isovol_meshsets[v1]['curves'][s1] = []

    #def _imprint_merge_curves(self):
    #
    #    # imprint and merge curves
    #    for iv, v1 in enumerate(vols_keys):
    #        # compare curves of each surface of a volume to all the other
    #        # surfaces in that volume
    #        for iss, s1 in enumerate(self.isovol_meshsets[v1]['surfs']):
    #
    #            for curve in self.isovol_meshsets[v1]['skins'][s1]
    #            curve1 = self.isovol_meshsets[v1]['skins'][s1]
    #            if curve1 != -1:
    #                for jss, s2 in enumerate(self.isovol_meshsets[v1]['surfs'][iss+1:]):
    #                    curve2 = self.isovol_meshsets[v1]['skins'][s2]
    #                    if curve2 != -1:
    #
    #                        # find if any verts/edges have identical handles
    #                        verts1 = self.mb.get_entities_by_type(curve1, types.MBVERTEX)
    #                        edges1 = self.mb.get_entities_by_type(curve1, types.MBEDGE)
    #
    #                        verts2 = self.mb.get_entities_by_type(curve2, types.MBVERTEX)
    #                        edges2 = self.mb.get_entities_by_type(curve2, types.MBEDGE)
    #
    #                        edge_match = rng.intersect(edges1, edges2)
    #                        verts_match = rng.intersect(
    #
    #
    #                        # imprint - get matches and get diffs
    #                        # separate matches -> create new curves
    #                        # add curve to each curve list
    #                        # separate diffs -> create new curves




    def _make_family(self):
        """Makes the correct parent-child relationships with volumes
        and surfaces. Tags geometry type, category, and ID on surfaces
        and volumes.
        """
        # create geometry dimension, category, and global id tags
        geom_dim = self.mb.tag_get_handle('GEOM_DIMENSION', size=1,
                        tag_type=types.MB_TYPE_INTEGER,
                        storage_type=types.MB_TAG_SPARSE,
                        create_if_missing=True)
        category = self.mb.tag_get_handle('CATEGORY', size=32,
                        tag_type=types.MB_TYPE_OPAQUE,
                        storage_type=types.MB_TAG_SPARSE,
                        create_if_missing=True)
        global_id = self.mb.tag_get_handle('GLOBAL_ID', size=1,
                        tag_type=types.MB_TYPE_INTEGER,
                        storage_type=types.MB_TAG_SPARSE,
                        create_if_missing=True)

        vol_id = 0
        surf_id = 0
        curve_id = 0

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

                # tag surfaces
                self.mb.tag_set_data(geom_dim, surf_eh, 2)
                self.mb.tag_set_data(category, surf_eh, 'Surface')
                surf_id += 1
                self.mb.tag_set_data(global_id, surf_eh, surf_id)

                # tag the curve of the surface
                #curve_eh = self.isovol_meshsets[v]['curves'][surf_eh]
                #if curve_eh != []:
                #    self.mb.tag_set_data(geom_dim, curve_eh, 1)
                #    self.mb.tag_set_data(category, curve_eh, 'Curve')
                #    curve_id += 1
                #    self.mb.tag_set_data(global_id, curve_eh, curve_id)
                #    self.mb.add_parent_child(surf_eh, curve_eh[0])

        curve_id = 0
        for s in self.surf_curve.keys():
            for c in self.surf_curve[s]:
                # create relationship
                self.mb.add_parent_child(s, c)

                 # tag curves
                self.mb.tag_set_data(geom_dim, c, 1)
                self.mb.tag_set_data(category, c, 'Curve')
                curve_id += 1
                self.mb.tag_set_data(global_id, c, curve_id)


    def _tag_groups(self):
        """Tags the surfaces with metadata as a group with value
            {data}_{value}. Example: 'wwn_0.005'. The original data
            name will be stripped of any underscores.
        """

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
        self.mb.tag_set_data(category, data_groups[0.0], 'Group')
        name = '{}_{}'.format(data, 0.0)
        self.mb.tag_set_data(tag_name, data_groups[0.0], name)

        for val in self.levels:
            val = val*self.norm
            data_groups[val] = self.mb.create_meshset()
            self.mb.tag_set_data(category, data_groups[val], 'Group')
            name = '{}_{}'.format(data, val)
            self.mb.tag_set_data(tag_name, data_groups[val], name)

        # add surfs to groups
        for isovol in self.isovol_meshsets.keys():
            for surf in self.isovol_meshsets[isovol]['surfs_EH']:
                # get the tagged data (get one value from the array)
                val_data = self.mb.tag_get_data(self.val_tag, surf)
                val = float(val_data[0])

                # add to group with that same data
                self.mb.add_entities(data_groups[val], [surf])


    def _tag_for_viz(self):
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


    def _add_mats(self):
        """Assign material void to all volumes.
        """

        print("Assigning void materials..")

        # create tags for adding materials to groups
        name_tag = self.mb.tag_get_handle('NAME', size=32,
                    tag_type=types.MB_TYPE_OPAQUE,
                    storage_type=types.MB_TAG_SPARSE,
                    create_if_missing=True)
        category = self.mb.tag_get_handle('CATEGORY', size=32,
                    tag_type=types.MB_TYPE_OPAQUE,
                    storage_type=types.MB_TAG_SPARSE,
                    create_if_missing=True)
        mat = 'mat:Vacuum'
        group_ms = self.mb.create_meshset()

        # add all volumes to group
        for isovol in self.isovol_meshsets.keys():
            self.mb.add_entities(group_ms, [isovol[1]])

        # assign as a group and assign material 0
        self.mb.tag_set_data(name_tag, group_ms, mat)
        self.mb.tag_set_data(category, group_ms, 'Group')

        print("... Done assigning materials!")
