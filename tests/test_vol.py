import os
from os import listdir
from os.path import isfile, join
import filecmp
import shutil
import pytest

from IsogeomGenerator import vol


test_dir = os.getcwd() + "/tests/test_files/"
ww_file = test_dir + "cwwm.vtk"


def test_assign_levels():
    """assign predetermined levels (out of order)"""
    v = vol.IsoVolDatabase()
    levels = [0.1, 0.2, 0.05]
    exp = sorted(levels)

    v.assign_levels(levels)
    assert(v.levels == exp)


# Generate Levels parametrized tests:
# linear: (6, 5, 15, 'lin', [5., 7., 9., 11., 13., 15.]
# log: (6, 1, 1e5, 'log', [1, 10, 1e2, 1e3, 1e4, 1e5])
# ratio, max included: (5, 1, 625, 'ratio', [1., 5., 25., 125., 625.])
# ratio, max not included: (5, 1, 700, 'ratio', [1., 5., 25., 125., 625.])
@pytest.mark.parametrize("N,minN,maxN,mode,exp",
                         [(6, 5, 15, 'lin', [5., 7., 9., 11., 13., 15.]),
                          (6, 1, 1e5, 'log', [1, 10, 1e2, 1e3, 1e4, 1e5]),
                          (5, 1, 625, 'ratio', [1., 5., 25., 125., 625.]),
                          (5, 1, 700, 'ratio', [1., 5., 25., 125., 625.])])
def test_generate_levels(N, minN, maxN, mode, exp):
    """generate levels with different modes"""
    v = vol.IsoVolDatabase()
    v.generate_levels(N, minN, maxN, mode=mode)
    assert(v.levels == exp)


# Generate Volumes parametrized tests:
# Min and Max w/in data bounds: (4, 8.e-7, 1.7e-6, 1)
# Max out of data bounds: (4, 8.e-7, 3.e-6, 2)
# Min out of data bounds: (4, 5.e-7, 1.7e-6, 3)
# Min and Max out of data bounds: (4, 5.e-7, 3.e-6, 4)
@pytest.mark.parametrize("N,minN,maxN,id", [(4, 8.e-7, 1.7e-6, 1),
                                            (4, 8.e-7, 3.e-6, 2),
                                            (4, 5.e-7, 1.7e-6, 3),
                                            (4, 5.e-7, 3.e-6, 4)])
def test_generate_volumes(N, minN, maxN, id):
    """Generate all isovolume files with different bounds"""
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


def test_generate_volumes_assign_levels():
    """Assign levels at time of generating levels"""
    # Expected results:
    exp_vols_dir = test_dir + "/vols-assign/vols"
    common_files = [f for f in listdir(exp_vols_dir)
                    if isfile(join(exp_vols_dir, f))]
    exp_levelfile = test_dir + "/vols-assign/levelfile"
    # Generate the volumes
    g = vol.IsoVolDatabase()
    levels = [8e-07, 1.2e-6, 1.7e-06]
    db = test_dir + "/test-assign"
    g.generate_volumes(ww_file, 'ww_n', dbname=db, levels=levels)
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


def test_generate_volumes_levelfile():
    """Read levels from file at time of generating levels"""
    # Expected results:
    exp_vols_dir = test_dir + "/vols-assign/vols"
    common_files = [f for f in listdir(exp_vols_dir)
                    if isfile(join(exp_vols_dir, f))]
    exp_levelfile = test_dir + "/vols-assign/levelfile"
    # Generate the volumes
    g = vol.IsoVolDatabase()
    levelfile_in = test_dir + "/vols-assign/levelfile"
    db = test_dir + "/test-read"
    g.generate_volumes(ww_file, 'ww_n', dbname=db, levelfile=levelfile_in)
    gen_vols_dir = db + "/vols"
    levelfile_out = db + "/levelfile"
    # check that files produced are the same
    res = filecmp.cmpfiles(exp_vols_dir, gen_vols_dir, common_files)
    match_list = res[0]
    non_match = res[1]
    assert(match_list == common_files)
    assert(non_match == [])
    # check that the level files are the same
    res = filecmp.cmp(exp_levelfile, levelfile_out)
    assert(res)
    shutil.rmtree(db)


def test_generate_volumes_genmode():
    pass


def test_generate_volumes_preset():
    pass


def test_generate_volumes_no_levels(self):
    """Try to generate volumes without assigning levels first and
    catch error."""
    g = vol.IsoVolDatabase()
    g.generate_volumes(ww_file, 'ww_n')
