#!/usr/bin/python

import os
import nose.tools
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal
import vol as v


def test_assign_levels():
    pass


def test_generate_levels_linear():
    pass


def test_generate_levels_log():
    pass


def test_generate_levels_ratio():
    pass


def test_generate_volumes():
    pass
    # assert number of files, assert names of files, assert not empty?
    # too difficult? assert file contents?
    # Test in different scenarios:
    #   1. both min and max are within range of data in the mesh
    #   2. min within range, max out of range
    #   3. min out of range, max within range
    #   4. single level


def test_create_geometry():
    # geom file creation starting from existing STL level files
    # no mat or viz tags
    pass


def test_create_geometry_full():
    # beginning to end geom creation (assigning levels)
    # no mat or viz tags
    pass


def test_separate_volumes():
    pass


def test_merge():
    pass


def test_add_mats():
    # full geom w/ mats
    pass


def test_viz_tag():
    # full geom w/ visualization tags
    pass
