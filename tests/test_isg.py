"""tests for the IsoSurfGeom module"""
from os import listdir, remove, getcwd, mkdir
from os.path import isfile, isdir, join
import pytest
from pymoab import core

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
