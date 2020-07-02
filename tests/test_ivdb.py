"""Tests for the ivdb module"""

import os
import pytest

from IsogeomGenerator import ivdb

test_dir = os.getcwd() + "/tests/test_files/"
db = test_dir + "/exp-test/"
levelfile = db + "/levelfile"
levels = [15, 5, 25, 35]
exp_levels = [5, 15, 25, 35]
data = 'dname'


def test_init_none():
    r0 = r1 = r2 = r3 = False
    iv = ivdb.IvDb()
    if iv.levels is None:
        r0 = True
    if iv.data is None:
        r1 = True
    if iv.db is None:
        r2 = True
    if iv.completed is False:
        r3 = True
    assert(all([r0, r1, r2, r3]))


def test_init_input_list():
    r0 = r1 = r2 = r3 = False
    iv = ivdb.IvDb(levels=levels, data=data, db=db)
    if iv.levels == exp_levels:
        r0 = True
    if iv.data == data:
        r1 = True
    if iv.db == db:
        r2 = True
    if iv.completed is False:
        r3 = True
    assert(all([r0, r1, r2, r3]))


def test_init_input_file():
    r0 = r1 = r2 = r3 = False
    iv = ivdb.IvDb(levels=levelfile, data=data, db=db)
    if iv.levels == exp_levels:
        r0 = True
    if iv.data == data:
        r1 = True
    if iv.db == db:
        r2 = True
    if iv.completed is False:
        r3 = True
    assert(all([r0, r1, r2, r3]))
