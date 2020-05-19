import os
from os import listdir
from os.path import isfile, join
import filecmp
import shutil

import unittest
import pytest

from IsogeomGenerator import vol


test_dir = os.getcwd() + "/tests/test_files/"
ww_file = test_dir + "cwwm.vtk"


class TestDatabase(unittest.TestCase):

    def test_assign_levels(self):
        """assign predetermined levels (out of order)"""
        v = vol.IsoVolDatabase()
        levels = [0.1, 0.2, 0.05]
        exp = sorted(levels)

        v.assign_levels(levels)
        assert(v.levels == exp)

    def test_generate_levels_linear(self):
        """generate linearly spaced levels"""
        v = vol.IsoVolDatabase()
        N = 6
        minN = 5
        maxN = 15
        exp = [5., 7., 9., 11., 13., 15.]

        v.generate_levels(N, minN, maxN, mode='lin')
        assert(v.levels == exp)

    def test_generate_levels_log(self):
        """generate lograthmically spaced levels"""
        v = vol.IsoVolDatabase()
        N = 6
        minN = 1
        maxN = 1e5
        exp = [1., 10., 1.e2, 1.e3, 1.e4, 1.e5]

        v.generate_levels(N, minN, maxN, mode='log')
        assert(v.levels == exp)

    def test_generate_levels_ratio_1(self):
        """generate levels by ratio, max included"""
        v = vol.IsoVolDatabase()
        N = 5
        minN = 1
        maxN = 625
        exp = [1., 5., 25., 125., 625.]

        v.generate_levels(N, minN, maxN, mode='ratio')
        assert(v.levels == exp)

    def test_generate_levels_ratio_2(self):
        """generate levels by ratio, max not included"""
        v = vol.IsoVolDatabase()
        N = 5
        minN = 1
        maxN = 700
        exp = [1., 5., 25., 125., 625.]

        v.generate_levels(N, minN, maxN, mode='ratio')
        assert(v.levels == exp)


# parametrized tests:
# Min and Max w/in data bounds: (4, 8.e-7, 1.7e-6, 1)
# Max out of data bounds: (4, 8.e-7, 3.e-6, 2)
# Min out of data bounds: (4, 5.e-7, 1.7e-6, 3)
# Min and Max out of data bounds: (4, 5.e-7, 3.e-6, 4)
@pytest.mark.parametrize("N,minN,maxN,id", [(4, 8.e-7, 1.7e-6, 1),
                                            (4, 8.e-7, 3.e-6, 2),
                                            (4, 5.e-7, 1.7e-6, 3),
                                            (4, 5.e-7, 3.e-6, 4)])
def test_generate_volumes(N, minN, maxN, id):
    """Generate all volume files when min and max are both within
    range of the data"""
    # Expected results:
    exp_vols_dir = test_dir + "/vols-{}/vols".format(id)
    common_files = [f for f in listdir(exp_vols_dir)
                    if isfile(join(exp_vols_dir, f))]
    exp_levelfile = test_dir + "/vols-{}/levelfile".format(id)
    # Generate the volumes
    g = vol.IsoVolDatabase()
    g.generate_levels(N, minN, maxN, mode='lin')
    db = test_dir + "/test-{}".format(id)
    g.generate_volumes(ww_file, 'ww_n', dbname=db)
    gen_vols_dir = db + "/vols"
    levelfile = db + "/levelfile"
    # check that files produced are the same
    res = filecmp.cmpfiles(exp_vols_dir, gen_vols_dir, common_files)
    match_list = res[0]
    non_match = res[1]
    assert(match_list == common_files)
    assert(non_match == [])
    # check that the level files are the same
    res = filecmp.cmp(exp_levelfile, levelfile)
    assert(res)
    shutil.rmtree(db)


    # def test_generate_volumes_single_level():
    #     # test that two levels are generated if a single value is passed
    #     # (code needs work to accomplish this)
    #     pass

    def test_generate_volumes_no_levels(self):
        """Try to generate volumes without assigning levels first and
        catch error."""
        pass
