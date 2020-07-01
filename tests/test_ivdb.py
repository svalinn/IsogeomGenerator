"""Tests for the ivdb module"""

import os
import pytest

from IsogeomGenerator import ivdb


def test_init():
    iv = ivdb.IvDb()
    assert(iv.levels is None)
