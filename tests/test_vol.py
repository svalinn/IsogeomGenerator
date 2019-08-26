#!/usr/bin/python

import os
import nose.tools
from nose.plugins.skip import SkipTest
from nose.tools import *
import vol as v

ww_file = "cwwm.vtk"


def test_assign_levels():
    """Test if values are being assigned properly when passed an
    unsorted list of values."""

    levels = [0.1, 0.15, 0.05]

    g = v.IsoVolume()
    g.assign_levels(levels)

    # tests
    assert_equal([0.05, 0.1, 0.15], g.levels)
    assert_equal(0.05, g.minN)
    assert_equal(0.15, g.maxN)
    assert_equal(3, g.N)



def test_generate_levels_linear():
    """Test if levels are correctly generated in a linear fashion when
    specifying a min/max value and a number of levels.
    """

    minN = 0.1
    maxN = 1.0
    N = 5
    levels_expected = [0.1, 0.325, 0.55, 0.775, 1.0]

    g = v.IsoVolume()
    g.generate_levels(N, minN, maxN, log=False)

    # tests
    assert_equal(levels_expected, g.levels)
    assert_equal(minN, g.minN)
    assert_equal(maxN, g.maxN)
    assert_equal(N, g.N)


def test_generate_levels_log():
    """Test if levels are correctly generated in a log fashion when
    specifying a min/max value and a number of levels.
    """

    minN = 0.1
    maxN = 1.e3
    N = 5
    levels_expected = [0.1, 1., 10., 100., 1000.]

    g = v.IsoVolume()
    g.generate_levels(N, minN, maxN, log=True)

    # tests (account for rounding errors w/ np logspace)
    for i, val in enumerate(levels_expected):
        assert_almost_equal(val, g.levels[i])

    assert_almost_equal(minN, g.minN)
    assert_almost_equal(maxN, g.maxN)
    assert_equal(N, g.N)


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


def test_create_geometry_full_norm():
    # add normalization factor
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
