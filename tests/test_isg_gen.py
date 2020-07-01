# test parent class

import os
import pytest

from IsogeomGenerator import isg_gen

levels = [1.0, 2.0, 5.0, 3.0]
levels_exp = [1.0, 2.0, 3.0, 5.0]
data = 'ww'
db = test_dir = os.getcwd() + "/tests/test_files/"


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
    if ig.levels == levels_exp:
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
