import visit as v
import sys
import os
import numpy as np
import math as m
from pymoab import core, types
from pymoab.rng import Range


class IsoVolumes(object):
    """This class contains methods to create a set of isovolumes given
    a cartesian mesh with data.

    Parameters:
    -----------
        filename: string, path to vtk file with the mesh
        data: string, name of the data to use from the file
    """

    def __init__(self, filename, data, dbname="tmp"):
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
            self.levels = np.logspace(start, stop, num=N,
                                      endpoint=True, base=base)
        else:
            self.levels = np.linspace(self.minN, self.maxN, num=N, endpoint=True)



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
        """Set selection for isovolume and export as STL."""

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
        cwd = os.getcwd()
        e.dirname = cwd + "/vols/" + self.db
        e.db_type = "STL"
        e.filename = str(i)
        e.variables = self.data
        v.ExportDatabase(e)

        # delete the operators (external surface/isovolume selection)
        v.RemoveAllOperators()


    def _generate_volumes(self):
        """Generates the isosurfaces volumes inbetween the contour levels.
        Data files are exported as STLs and saved in the folder dbname.
        Files will be named based on their index corresponding to their
        level values (0.stl is lowest values).
        """

        # make sure levels have been set before proceding
        if not self.levels:
            print("No contour levels have been set!")
            print("Please use assign_levels or generate_levels to set levels.")
            return

        # create folder to store data if it does not already exist
        if not os.path.isdir(self.db):
            os.mkdir(self.db)

        # launch VisIt and open database
        v.LaunchNowin()
        v.OpenDatabase(self.f)

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
        ubound = 1000.
        self._get_isovol(lbound, ubound, i+1)

        # close everything
        v.DeleteAllPlots()
        v.Close()


    def _generate_contours(self):
        """Generates the single surface contours that match the generated
        isovolumes.
        """

        # draw matching contours
        v.LaunchNowin()
        v.OpenDatabase(self.f)
        v.AddPlot("Contour", self.data)
        att = v.ContourAttributes()
        att.SetContourMethod(1)
        att.SetContourValue(tuple(self.levels))
        v.SetPlotOptions(att)
        v.DrawPlots()

        # export
        e = v.ExportDBAttributes()
        cwd = os.getcwd()
        e.dirname = cwd + "/" + self.db
        e.db_type = "STL"
        e.filename = "c"
        e.variables = self.data
        v.ExportDatabase(e)

        # close everything
        v.DeleteAllPlots()
        v.Close()


    def _separate_isovols(self):
        """for each volume surf in the database, separate into single
        surface files.
        """
        for f in os.listdir(self.db + "/vols/"):

            # get file
            fpath = os.getcwd() + "/" + self.db + "/vols/" + f
            rootname = f.strip(".stl")

            # start pymoab instance
            mb = core.Core()

            # load file and get set of all verts
            mb.load_file(fpath)
            root_set = mb.get_root_set()
            all_verts = mb.get_entities_by_type(root_set, types.MBVERTEX)

            i = 0
            while len(all_verts) > 0:

                # get full set of connected verts starting from a seed
                verts = [all_verts[0]]
                while True:
                    # this step takes too long for complex surfaces
                    vtmp = mb.get_adjacencies(mb.get_adjacencies(verts, 2, op_type=1), 0, op_type=1)
                    if set(list(vtmp)) == set(list(verts)):
                        break
                    else:
                        verts = vtmp

                # get the connected set of triangles that make the single surf
                tris = mb.get_adjacencies(verts, 2, op_type=1)
                surf = mb.create_meshset()
                mb.add_entities(surf, tris)

                # write to file
                r_surf = Range(surf)
                mb.write_file(self.db + "/{}-{}.stl".format(rootname, i), r_surf)

                # remove surface from meshset
                mb.delete_entities(tris)
                mb.delete_entities(verts)

                # resassign vertices
                all_verts = mb.get_entities_by_type(root_set, types.MBVERTEX)
                i += 1

            # add line here to delete the original file f since it is now
            # broken into multiple


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

        # Step 1: Generate Isovolumes using VisIT
        self._generate_volumes()

        # Step 2: Generate corresponding Isocontours using VisIT
        self._generate_contours()

        # Step 3: Separate isovolumes into Surfaces w/ PyMOAB
        self._separate_isovols()

        # Step 4: separate isocontours into separate files

        # Step 5: create parent-child meshsets for each isovolume

        # Step 6: create isocontour meshsets

        # Step 7: merge isocontour/isovol surfaces (see algorithm in notes)

        # Step 8: write out as single h5m