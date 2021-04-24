# test parent class

import os
import pytest

from IsogeomGenerator import isg_gen

test_dir = os.getcwd() + "/tests/test_files/"
db = test_dir + "/exp-test/"
levelfile = db + "/levelfile"
levels = [15, 5, 25, 35]
exp_levels_pre = [5, 15, 25, 35]  # after ivdb init only
exp_levels_post = [5, 15, 25, 35, 45]  # after the arbmax addition
data = 'dname'


def test_init_none():
    r0 = r1 = r2 = False
    ig = isg_gen.IsoGeomGen()
    if ig.levels is None:
        r0 = True
    if ig.data is None:
        r1 = True
    if ig.db == os.getcwd() + "/tmp":
        r2 = True
    assert(all([r0, r1, r2]))


def test_init_input_list():
    r0 = r1 = r2 = False
    ig = isg_gen.IsoGeomGen(levels=levels, data=data, db=db)
    if ig.levels == exp_levels_pre:
        r0 = True
    if ig.data == data:
        r1 = True
    if ig.db == db:
        r2 = True
    assert(all([r0, r1, r2]))


def test_init_input_file():
    r0 = r1 = r2 = False
    ig = isg_gen.IsoGeomGen(levels=levelfile, data=data, db=db)
    if ig.levels == exp_levels_post:
        # level file has the arbmax in it already
        r0 = True
    if ig.data == data:
        r1 = True
    if ig.db == db:
        r2 = True
    assert(all([r0, r1, r2]))


def test_read_levels_list():
    ig = isg_gen.IsoGeomGen()
    ig.read_levels(levels)
    assert(ig.levels == exp_levels_pre)


def test_read_levels_file():
    ig = isg_gen.IsoGeomGen()
    ig.read_levels(levelfile)
    assert(ig.levels == exp_levels_post)


@pytest.mark.parametrize("input,error_str", [(1, "Type of levels"),
                                             ('fake', "does not exist")])
def test_read_levels_error(input, error_str):
    ig = isg_gen.IsoGeomGen()
    with pytest.raises(RuntimeError) as error_info:
        ig.read_levels(input)
    assert error_str in str(error_info)


def test_assign_levels():
    ig = isg_gen.IsoGeomGen()
    ig._IsoGeomGen__assign_levels(levels)
    assert(ig.levels == exp_levels_pre)
