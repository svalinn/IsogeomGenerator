"""tests for the IsoSurfGeom module"""
from os import listdir, remove, getcwd, mkdir
from os.path import isfile, isdir, join
import pytest
from pymoab import core
import numpy as np

from IsogeomGenerator import isg, ivdb

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


def __ivdb_obj(completed):
    iv = ivdb.IvDb(levels=exp_levels, data=data, db=exp_db)
    iv.completed = completed
    return iv


def test_init_none():
    r0 = r1 = r2 = r3 = r4 = False
    ig = isg.IsGm()
    if ig.levels is None:
        r0 = True
    if ig.data is None:
        r1 = True
    if ig.db == getcwd() + "/tmp":
        r2 = True
    if isinstance(ig.mb, type(core.Core())):
        r3 = True
    if ig.isovol_meshsets == {}:
        r4 = True
    assert(all([r0, r1, r2, r3, r4]))


def test_init_input():
    r0 = r1 = r2 = False
    ig = isg.IsGm(levels=levels, data=data, db=exp_db)
    if ig.levels == exp_levels:
        r0 = True
    if ig.data == data:
        r1 = True
    if ig.db == exp_db:
        r2 = True
    assert(all([r0, r1, r2]))


def test_init_ivdb():
    """test that info is taken from ivdb"""
    r0 = r1 = r2 = False
    iv = __ivdb_obj(True)
    ig = isg.IsGm(ivdb=iv)
    if ig.levels == exp_levels:
        r0 = True
    if ig.data == data:
        r1 = True
    if ig.db == exp_db:
        r2 = True
    assert(all([r0, r1, r2]))


def test_init_input_ivdb():
    """test that info from ivdb overwrites other input"""
    r0 = r1 = r2 = False
    iv = __ivdb_obj(True)
    ig = isg.IsGm(ivdb=iv, levels=[0, 2], data='nonsense', db='fake_db')
    if ig.levels == exp_levels:
        r0 = True
    if ig.data == data:
        r1 = True
    if ig.db == exp_db:
        r2 = True
    assert(all([r0, r1, r2]))


def test_init_input_file():
    ig = isg.IsGm(levels=exp_levelfile)
    assert(ig.levels == exp_levels)


def test_read_ivdb():
    """read info from ivdb obj"""
    iv = __ivdb_obj(True)
    ig = isg.IsGm()
    ig.read_ivdb(iv)
    r0 = r1 = r2 = False
    if ig.levels == exp_levels:
        r0 = True
    if ig.data == data:
        r1 = True
    if ig.db == exp_db:
        r2 = True
    assert(all([r0, r1, r2]))


def test_read_ivdb_incomplete():
    """raise error if incomplete ivdb obj"""
    iv = __ivdb_obj(False)
    ig = isg.IsGm()
    with pytest.raises(RuntimeError) as error_info:
        ig.read_ivdb(iv)
    assert "Incomplete IvDb object" in str(error_info)


def test_read_database():
    """check that meshsets are properly populated with read_database"""
    # create obj and read database
    ig = isg.IsGm(levels=levels, data=data, db=exp_db)
    ig.read_database()
    # expected meshset entity handles
    ehs = [12682136550675316737,
           12682136550675316738,
           12682136550675316739,
           12682136550675316740,
           12682136550675316741]
    # setup truth array
    res = np.full(len(ehs) + 1, False)
    # check that meshsets exist in the moab instance
    for r, eh in enumerate(ehs):
        try:
            # any moab call that will work if meshet exists, else
            # it will fail
            ig.mb.get_child_meshsets(eh)
        except RuntimeError:
            pass
        else:
            res[r] = True
    # check meshets and bound information are in dictionary
    exp_meshsets = {(0, ehs[0]): {'bounds': (None, 5.0)},
                    (1, ehs[1]): {'bounds': (5.0, 15.0)},
                    (2, ehs[2]): {'bounds': (15.0, 25.0)},
                    (3, ehs[3]): {'bounds': (25.0, 35.0)},
                    (4, ehs[4]): {'bounds': (35.0, None)}}
    if sorted(ig.isovol_meshsets) == sorted(exp_meshsets):
        res[-1] = True
    # assert all pass
    assert(all(res))


def test_read_database_error():
    """read_database throws error if num levels and files mismatch"""
    # create obj and read database
    ig = isg.IsGm(levels=[300], data=data, db=exp_db)
    with pytest.raises(RuntimeError) as error_info:
        ig.read_database()
    assert "does not match number" in str(error_info)
