import visit as v
import sys
import os
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

        # generate evenly spaced values
        if log:
            base = 10.
            start = m.log(self.minN, base)
            stop = m.log(self.maxN, base)
            self.levels = list(np.logspace(start, stop, num=N,
                                      endpoint=True, base=base))
        else:
            self.levels = list(np.linspace(self.minN, self.maxN,
                                      num=N, endpoint=True))


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
        if not os.path.isdir(self.db + "/vols/"):
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
            rootname = f.strip(".stl")

            # load file and create EH for file-set
            fs = self.mb.create_meshset()
            self.mb.load_file(fpath, file_set=fs)

            # initiate dictionary
            iv_info = (rootname, fs)
            self.isovol_meshsets[iv_info] = {}

            # separate
            self._separate(iv_info)


    def merge_surfs(self):
        """use pymoab to create parent/child meshsets for volumes.
        """
        pass


    def create_geometry(self):
        """Over-arching function to do all steps to create a single
        isovolume geometry for DAGMC.
        """
        # either step 0 or 1a - calculate average WW value to use in
        # each volume and tag
        # create tuples of (vol_id, ww_val)

        ###########################################
        # Step 1: Generate Isovolumes using VisIT #
        ###########################################
        v.LaunchNowin()
        v.OpenDatabase(self.f)
        print("generating isovolumes...")
        self._generate_volumes()
        print("isovolumes complete!")
        v.Close()


        ####################################################
        # Step 2: Separate isovolumes into unique surfaces #
        ####################################################
        self.mb = core.Core()
        self.isovol_meshsets = {}
        print("separating isovolumes...")
        self._separate_isovols()
        print("separation complete!")
        print(self.isovol_meshsets.items())

        # Step 3: Merge Isovolume Surfaces

        # for every isovolume:
        #     for every separate_vol in isovolume:
        #          check if parts match every separate_vol in isvols[i+1]

        # Step 4: create parent-child meshsets for each isovolume

        # Step 5: write out as single h5m
