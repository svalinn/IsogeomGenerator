import os
from os import listdir
from os.path import isfile, join
import filecmp
import shutil

import unittest

from IsogeomGenerator import vol

test_dir = os.getcwd() + "/tests/test_files/"
ww_file = test_dir + "cwwm.vtk"

class TestIsogeom(unittest.TestCase):

    def test_assign_levels(self):
        """assign predetermined levels (out of order)"""
        v = vol.IsoVolume()
        levels = [0.1, 0.2, 0.05]
        exp = sorted(levels)
        v.assign_levels(levels)
        assert(v.levels == exp)
        assert(v.minN == min(levels))
        assert(v.maxN == max(levels))
        assert(v.N == len(levels))

    def test_generate_levels_linear(self):
        """generate linearly spaced levels"""
        v = vol.IsoVolume()
        N = 6
        minN = 5
        maxN = 15
        exp = [5., 7., 9., 11., 13., 15.]
        v.generate_levels(N, minN, maxN, log=False)
        assert(v.levels == exp)
        assert(v.minN == min(exp))
        assert(v.maxN == max(exp))
        assert(v.N == len(exp))

    def test_generate_levels_log(self):
        """generate lograthmically spaced levels"""
        v = vol.IsoVolume()
        N = 6
        minN = 1
        maxN = 1e5
        exp = [1., 10., 1.e2, 1.e3, 1.e4, 1.e5]
        v.generate_levels(N, minN, maxN, log=True)
        assert(v.levels == exp)
        assert(v.minN == min(exp))
        assert(v.maxN == max(exp))
        assert(v.N == len(exp))

    def test_generate_levels_ratio_1(self):
        """generate levels by ratio, max included"""
        v = vol.IsoVolume()
        N = 5
        minN = 1
        maxN = 625
        exp = [1., 5., 25., 125., 625.]

        v.generate_levels(N, minN, maxN, ratio=True)
        assert(v.levels == exp)
        assert(v.minN == min(exp))
        assert(v.maxN == max(exp))
        assert(v.N == len(exp))

    def test_generate_levels_ratio_2(self):
        """generate levels by ratio, max not included"""
        v = vol.IsoVolume()
        N = 5
        minN = 1
        maxN = 700
        exp = [1., 5., 25., 125., 625.]

        v.generate_levels(N, minN, maxN, ratio=True)
        assert(v.levels == exp)
        assert(v.minN == min(exp))
        assert(v.maxN == max(exp))
        assert(v.N == len(exp))

    def test_generate_levels_ratio_3(self):
        """generate levels by ratio, override log=True"""
        v = vol.IsoVolume()
        N = 5
        minN = 1.0
        maxN = 125.0
        exp = [1., 5., 25., 125.]

        v.generate_levels(N, minN, maxN, ratio=True, log=True)
        assert(v.levels == exp)
        assert(v.minN == min(exp))
        assert(v.maxN == max(exp))
        assert(v.N == len(exp))

    def test_generate_volumes_in_range(self):
        """Generate all volume files when min and max are both within
        range of the data"""
        g = vol.IsoVolume()
        g.generate_levels(4, 8.e-7, 1.7e-6, log=False)
        db = test_dir + "/test-1"
        g.generate_volumes(ww_file, 'ww_n', dbname=db)

        comp_vols = test_dir + "/vols-1/vols"
        comp_files = [f for f in listdir(comp_vols) if isfile(join(comp_vols, f))]
        gen_vols = db + "/vols"

        # check that files produced are the same
        res = filecmp.cmpfiles(comp_vols, gen_vols, comp_files)
        match_list = res[0]
        non_match = res[1]
        assert(match_list == comp_files)
        assert(non_match == [])

        shutil.rmtree(db)

    def test_generate_volumes_max_out_of_range(self):
        """Generate all volume files when max is out of range of the
        data"""
        g = vol.IsoVolume()
        g.generate_levels(4, 8.e-7, 3.e-6, log=False)
        db = test_dir + "/test-2"
        g.generate_volumes(ww_file, 'ww_n', dbname=db)

        comp_vols = test_dir + "/vols-2/vols"
        comp_files = [f for f in listdir(comp_vols) if isfile(join(comp_vols, f))]
        gen_vols = db + "/vols"

        # check that files produced are the same
        res = filecmp.cmpfiles(comp_vols, gen_vols, comp_files)
        match_list = res[0]
        non_match = res[1]
        assert(match_list == comp_files)
        assert(non_match == [])

        shutil.rmtree(db)


    def test_generate_volumes_min_out_of_range(self):
        """Generate all volume files when min is out of range of the
        data"""
        g = vol.IsoVolume()
        g.generate_levels(4, 5.e-7, 1.7e-6, log=False)
        db = test_dir + "/test-3"
        g.generate_volumes(ww_file, 'ww_n', dbname=db)

        comp_vols = test_dir + "/vols-3/vols"
        comp_files = [f for f in listdir(comp_vols) if isfile(join(comp_vols, f))]
        gen_vols = db + "/vols"

        # check that files produced are the same
        res = filecmp.cmpfiles(comp_vols, gen_vols, comp_files)
        match_list = res[0]
        non_match = res[1]
        assert(match_list == comp_files)
        assert(non_match == [])

        shutil.rmtree(db)


    def test_generate_volumes_not_in_range(self):
        """Generate all volume files when min and max are out of range
        of the data"""
        g = vol.IsoVolume()
        g.generate_levels(4, 5.e-7, 3.e-6, log=False)
        db = test_dir + "/test-4"
        g.generate_volumes(ww_file, 'ww_n', dbname=db)

        comp_vols = test_dir + "/vols-4/vols"
        comp_files = [f for f in listdir(comp_vols) if isfile(join(comp_vols, f))]
        gen_vols = db + "/vols"

        # check that files produced are the same
        res = filecmp.cmpfiles(comp_vols, gen_vols, comp_files)
        match_list = res[0]
        non_match = res[1]
        assert(match_list == comp_files)
        assert(non_match == [])

        shutil.rmtree(db)


    # def test_generate_volumes_single_level():
    #     # test that two levels are generated if a single value is passed
    #     # (code needs work to accomplish this)
    #     pass

    def test_generate_volumes_no_levels(self):
        """Try to generate volumes without assigning levels first and
        catch error."""
        pass

    ## Can these only be run if all the other above tests pass?
    def test_create_geometry(self):
        """Test geom creation from existing STL files."""
        pass

    def test_create_geometry_no_database(self):
        """Try to create geometry without a database from VisIt."""
        pass

    def test_create_geometry_out_of_range(self):
        """Create a geometry when some levels were out of range"""
        pass

    def test_create_geometry_with_tags(self):
        """Create geometry with the tags and viz"""
        pass
