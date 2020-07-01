"""Tests for the ivdb module"""

import os
import pytest

from IsogeomGenerator import ivdb


def test_init():
    iv = ivdb.IvDb(levels=[1, 2], data='ww')
    assert(iv.levels == [1., 2.])
    assert(iv.data == 'ww')
    assert(iv.db is None)
