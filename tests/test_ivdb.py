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
# there are two different level lists expected: one from the init and
# one after the first step is fully complete (arbmax has been added to
# the level list for the exp_levels_post)
exp_levels_init = [5, 15, 25, 35]
exp_levels_post = [5, 15, 25, 35, 45]


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
    if iv.levels == exp_levels_init:
        r0 = True
    if iv.data == data:
        r1 = True
    if iv.db == exp_db:
        r2 = True
    if iv.completed is False:
        r3 = True
    assert(all([r0, r1, r2, r3]))


def test_init_input_file():
    # since we are reading from a file, we expect it to be the same values
    # as the file. Since the file matches exp_levels_post, the test checks
    # against that list, rather than exp_levels_init even though it is only
    # the init
    iv = ivdb.IvDb(levels=exp_levelfile)
    assert(iv.levels == exp_levels_post)


def test_generate_vols():
    """Generate all isovolume files."""
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


def test_generate_vols_single():
    """Generate all isovolume files with single volume"""
    # assert flags
    r1 = r2 = False
    # test database path
    db = test_dir + "/test-gen-vols-single"
    if isdir(db):
        shutil.rmtree(db)
    # init ivdb obj
    iv = ivdb.IvDb(levels=[25], data=data, db=db)
    iv.generate_vols(test_mesh)
    gen_vols_dir = db + "/vols"
    # expected files
    exp_db = test_dir + "/exp-single/"
    exp_vols_dir = exp_db + "/vols"
    common_files = [f for f in listdir(exp_vols_dir)
                    if isfile(join(exp_vols_dir, f))]
    res = filecmp.cmpfiles(exp_vols_dir, gen_vols_dir, common_files)
    match_list = res[0]
    non_match = res[1]
    # check that files produced are the same
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
    db_exp = test_dir + "/test-direxists_1/"
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
    iv = ivdb.IvDb(levels=levels, data=data)
    arbmin, arbmax = iv._IvDb__check_levels(test_mesh)
    exp_min = -10
    exp_max = 50
    if arbmin == exp_min:
        r0 = True
    if arbmax == exp_max:
        r1 = True
    if iv.levels == exp_levels_init:
        # arbmax has not yet been added to levels list so should match
        # init levels
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
    if iv.levels == exp_levels_init:
        # arbmax has not yet been added to levels list so should match
        # init levels
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
    # this test goes right to the write step and skips that addition
    # of the arbmax value in the process, so we will init with a level
    # list that already includes the arbmax value and expect it to match
    # the level file that would be present after the complete process
    iv = ivdb.IvDb(levels=exp_levels_post, data=data, db=db)
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


def test_get_isovol():
    """test get_isovol"""
    r0 = r1 = r2 = False
    db = test_dir + "/test-isovol/"
    iv = ivdb.IvDb(levels=levels, data=data, db=db)
    if isdir(db):
        shutil.rmtree(db)
    mkdir(db)
    mkdir(db + '/vols/')
    # launch VisIt
    import visit
    try:
        visit.LaunchNowin()
    except:
        pass
    visit.OpenDatabase(test_mesh)
    visit.AddPlot('Pseudocolor', data)
    # run __get_isovol
    lbound = 5
    ubound = 15
    i = 1
    export_res, ubound_out = iv._IvDb__get_isovol(lbound, ubound, i)
    # close VisIt
    visit.CloseComputeEngine()
    # check returned values
    if export_res == 1:
        r0 = True
    if ubound_out == ubound:
        r1 = True
    # check that vol file produced are the same
    gen_vol = db + "/vols/1.stl"
    exp_vol = exp_vols_dir + "/1.stl"
    r2 = filecmp.cmp(gen_vol, exp_vol)
    shutil.rmtree(db)
    assert(all([r0, r1, r2]))


def test_get_isovol_nodata():
    """test get_isovol no data present"""
    r0 = r1 = r2 = rs = r4 = False
    db = test_dir + "/test-isovol-nodata/"
    iv = ivdb.IvDb(levels=[5, 15, 25, 28, 35], data=data, db=db)
    if isdir(db):
        shutil.rmtree(db)
    mkdir(db)
    mkdir(db + '/vols/')
    # launch VisIt
    import visit
    try:
        visit.LaunchNowin()
    except:
        pass
    visit.OpenDatabase(test_mesh)
    visit.AddPlot('Pseudocolor', data)
    # run __get_isovol
    lbound = 25
    ubound = 28
    i = 3
    with pytest.warns(None) as warn_info:
        export_res, ubound_out = iv._IvDb__get_isovol(lbound, ubound, i)
    # close VisIt
    visit.CloseComputeEngine()
    # check returned/changed values
    if export_res == 0:
        r0 = True
    if ubound_out == 35:
        r1 = True
    if iv.levels == [5, 15, 25, 35]:
        r2 = True
    gen_vol = db + "/vols/3.stl"
    if not isfile(gen_vol):
        # no file should be generated
        r3 = True
    if len(warn_info) == 1:
        r4 = True
    shutil.rmtree(db)
    assert(all([r0, r1, r2, r3, r4]))
