"""Tests for the ivdb module"""

from os import listdir, remove, getcwd, mkdir
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


def test_generate_vols():
    """Generate all isovolume files.
    """
    # assert flags
    r1 = r2 = False
    # test database path
    db = test_dir + "/test-gen-vols"
    if isdir(db):
        shutil.rmtree(db)
    # init ivdb obj
    iv = ivdb.IvDb(levels=levels, data=data, db=db)
    iv.generate_vols(test_mesh)
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
    assert(all([r1, r2]))


def test_make_db_dir():
    db = test_dir + "/test-mkdir"
    if isdir(db):
        shutil.rmtree(db)
    iv = ivdb.IvDb(levels=levels, data=data, db=db)
    iv._IvDb__make_db_dir()
    res = False
    if isdir(db):
        res = True
    shutil.rmtree(iv.db)
    assert(res)


def test_make_db_dir_exists():
    r0 = r1 = r2 = False
    db = test_dir + "/test-direxists"
    if isdir(db):
        shutil.rmtree(db)
    mkdir(db)
    db_exp = test_dir + "/test-direxists-1/"
    iv = ivdb.IvDb(levels=levels, data=data, db=db)
    with pytest.warns(None) as warn_info:
        iv._IvDb__make_db_dir()
    if len(warn_info) == 1:
        r0 = True
    if iv.db == db_exp:
        r1 = True
    if isdir(iv.db):
        r2 = True
    shutil.rmtree(iv.db)
    shutil.rmtree(db)
    assert(all([r0, r1, r2]))


def test_check_levels():
    """test check levels, all data good"""
    r0 = r1 = r2 = False
    iv = ivdb.IvDb(levels=exp_levels, data=data)
    arbmin, arbmax = iv._IvDb__check_levels(test_mesh)
    exp_min = -10
    exp_max = 50
    if arbmin == exp_min:
        r0 = True
    if arbmax == exp_max:
        r1 = True
    if iv.levels == exp_levels:
        r2 = True
    assert(all([r0, r1, r2]))


def test_check_levels_outofbounds():
    """test check levels, data out of bounds"""
    r0 = r1 = False
    iv = ivdb.IvDb(levels=[-5, 5, 15, 25, 35, 45], data=data)
    # data out of range should produce warning
    with pytest.warns(None) as warn_info:
        iv._IvDb__check_levels(test_mesh)
    if len(warn_info) == 2:
        r0 = True
    if iv.levels == exp_levels:
        r1 = True
    assert(all([r0, r1]))


@pytest.mark.filterwarnings("ignore:Level")
def test_check_levels_nodata():
    """check levels, no levels"""
    r0 = r1 = False
    iv = ivdb.IvDb(levels=[-5, 0, 45], data=data)
    # no levels remaining, creates error
    with pytest.raises(RuntimeError) as error_info:
        iv._IvDb__check_levels(test_mesh)
    if "No data exists" in str(error_info):
        r0 = True
    if iv.levels == []:
        r1 = True
    assert(all([r0, r1]))


def test_write_levels():
    db = test_dir + "/test-write-levels/"
    iv = ivdb.IvDb(levels=exp_levels, data=data, db=db)
    if isdir(db):
        shutil.rmtree(db)
    mkdir(db)
    iv.write_levels()
    levelfile_out = db + "/levelfile"
    if isfile(levelfile_out):
        r0 = True
    r1 = __compare_levelfiles(levelfile_out, exp_levelfile)
    shutil.rmtree(db)
    assert(all([r0, r1]))
