"""Tests for the ivdb module"""

from os import listdir, remove, getcwd
from os.path import isfile, isdir, join
import filecmp
import shutil
import pytest

from IsogeomGenerator import ivdb


# Set up test files and expected results
test_dir = getcwd() + "/tests/test_files/"
test_mesh = test_dir + "test_mesh.vtk"
data = 'dname'
levels = [15, 5, 25, 35]
exp_db = test_dir + "/exp-test/"
exp_vols_dir = exp_db + "/vols"
common_files = [f for f in listdir(exp_vols_dir)
                if isfile(join(exp_vols_dir, f))]
exp_levelfile = exp_db + "/levelfile"
exp_levels = [5, 15, 25, 35]
exp_geom = test_dir + '/exp-isogeom.h5m'


def __compare_levelfiles(f1, f2):
    """compares the levelfile contents regardless of string format"""
    fs1 = open(f1, 'r')
    fs2 = open(f2, 'r')

    levels1 = []
    lines = fs1.readlines()
    for line in lines:
        levels1.append(float(line))

    levels2 = []
    lines = fs2.readlines()
    for line in lines:
        levels2.append(float(line))

    if levels1 == levels2:
        return True
    else:
        return False


def test_init_none():
    r0 = r1 = r2 = r3 = False
    iv = ivdb.IvDb()
    if iv.levels is None:
        r0 = True
    if iv.data is None:
        r1 = True
    if iv.db == getcwd() + "/tmp":
        r2 = True
    if iv.completed is False:
        r3 = True
    assert(all([r0, r1, r2, r3]))


def test_init_input():
    r0 = r1 = r2 = r3 = False
    iv = ivdb.IvDb(levels=levels, data=data, db=exp_db)
    if iv.levels == exp_levels:
        r0 = True
    if iv.data == data:
        r1 = True
    if iv.db == exp_db:
        r2 = True
    if iv.completed is False:
        r3 = True
    assert(all([r0, r1, r2, r3]))


def test_init_input_file():
    iv = ivdb.IvDb(levels=exp_levelfile)
    assert(iv.levels == exp_levels)


# Generate Volumes parametrized tests:
#   (1) Min and Max w/in data bounds
#   (2) Mid level no data
#   (3) Min and max out of data bounds
@pytest.mark.parametrize("levels,id", [([5, 15, 25, 35], 1),
                                       ([5, 15, 25, 28, 35], 2)])
def test_generate_vols(levels, id):
    """Generate all isovolume files with different bounds.
    Some should raise warnings. Make sure generated files are correct.
    """
    # assert flags
    r0 = r1 = r2 = False
    # test database path
    db = test_dir + "/test-gen-vols-{}".format(id)
    if isdir(db):
        shutil.rmtree(db)
    # init ivdb obj
    iv = ivdb.IvDb(levels=levels, data=data, db=db)
    # data out of range should produce warning
    with pytest.warns(None) as warn_info:
        iv.generate_vols(test_mesh)
    # assert number of warnings raised
    if id == 2:
        # no data found warning
        if len(warn_info) == 1:
            r0 = True
    else:
        if len(warn_info) == 0:
            r0 = True
    # check that files produced are the same
    gen_vols_dir = db + "/vols"
    res = filecmp.cmpfiles(exp_vols_dir, gen_vols_dir, common_files)
    match_list = res[0]
    non_match = res[1]
    if match_list == common_files:
        r1 = True
    if non_match == []:
        r2 = True
    # remove files
    shutil.rmtree(iv.db)
    # check results
    assert(all([r0, r1, r2]))


def test_check_levels():
    # Generate the volumes
    iv = ivdb.IvDb(levels=exp_levels, data=data)
    arbmin, arbmax = iv._IvDb__check_levels(test_mesh)
    exp_min = -10
    exp_max = 50
    assert(arbmin == exp_min)
    assert(arbmax == exp_max)


def test_check_levels_outofbounds():
    """level values provided are not within data bounds, so raise error"""
    # Generate the volumes
    iv = ivdb.IvDb(levels=[-5, 5, 15, 25, 35, 45], data=data)
    # data out of range should produce warning
    with pytest.warns(None) as warn_info:
        iv._IvDb__check_levels(test_mesh)
    assert(len(warn_info) == 2)


@pytest.mark.filterwarnings("ignore:Level")
def test_check_levels_nodata():
    """level values provided are not within data bounds, so raise error"""
    # Generate the volumes
    iv = ivdb.IvDb(levels=[-5, 0, 45], data=data)
    # data out of range should produce warning
    with pytest.raises(RuntimeError) as error_info:
        iv._IvDb__check_levels(test_mesh)
    assert "No data exists" in str(error_info)
