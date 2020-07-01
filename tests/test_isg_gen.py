# test parent class

import os
import pytest

from IsogeomGenerator import isg_gen

test_dir = os.getcwd() + "/tests/test_files/"
db = test_dir + "/exp-test/"
levelfile = db + "/levelfile"
levels = [15, 5, 25, 35]
exp_levels = [5, 15, 25, 35]
data = 'dname'


def test_init_none():
    r0 = r1 = r2 = False
    ig = isg_gen.IsoGeomGen()
    if ig.levels is None:
        r0 = True
    if ig.data is None:
        r1 = True
    if ig.db is None:
        r2 = True
    assert(all([r0, r1, r2]))


def test_init_input():
    r0 = r1 = r2 = False
    ig = isg_gen.IsoGeomGen(levels=levels, data=data, db=db)
    if ig.levels == exp_levels:
        r0 = True
    if ig.data == data:
        r1 = True
    if ig.db == db:
        r2 = True
    assert(all([r0, r1, r2]))


def test_init_levels_error():
    with pytest.raises(RuntimeError) as error_info:
        ig = isg_gen.IsoGeomGen(levels=1)
    assert "Type of levels" in str(error_info)


def test_assign_levels():
    ig = isg_gen.IsoGeomGen()
    ig.assign_levels(levels)
    assert(ig.levels == exp_levels)


def test_read_levels():
    ig = isg_gen.IsoGeomGen()
    ig.read_levels(levelfile)
    assert(ig.levels == exp_levels)


def test_read_levels_error():
    ig = isg_gen.IsoGeomGen()
    with pytest.raises(RuntimeError) as error_info:
        ig.read_levels('fake_file')
    assert "does not exist" in str(error_info)
