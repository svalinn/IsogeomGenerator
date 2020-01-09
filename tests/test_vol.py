#!/usr/bin/python

import os
from os import listdir
from os.path import isfile, join
import nose.tools
from nose.plugins.skip import SkipTest
from nose.tools import *
import vol as v
import filecmp
import shutil

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


def test_generate_levels_ratio_1():
    """Test if levels are generated properly when assigned a spacing
    ratio. First test: max value should be included in list of levels.
    """
    minN = 1.0
    maxN = 125.0
    N = 5.
    levels_expected = [1., 5., 25., 125.]

    g = v.IsoVolume()
    g.generate_levels(N, minN, maxN, ratio=True)
    assert_equal(levels_expected, g.levels)
    assert_equal(minN, g.minN)
    assert_equal(maxN, g.maxN)
    assert_equal(4, g.N)


def test_generate_levels_ratio_2():
    """Test if levels are generated properly when assigned a spacing
    ratio. Second test: max value should not be included in list of
    levels.
    """
    minN = 1.0
    maxN = 122.0
    N = 5.
    levels_expected = [1., 5., 25.]

    g = v.IsoVolume()
    g.generate_levels(N, minN, maxN, ratio=True)
    assert_equal(levels_expected, g.levels)
    assert_equal(minN, g.minN)
    assert_equal(25., g.maxN)
    assert_equal(3, g.N)


def test_generate_levels_ratio_3():
    """Test if levels are generated properly when assigned a spacing
    ratio. Third test: make sure that log=True is overrided.
    """
    minN = 1.0
    maxN = 125.0
    N = 5.
    levels_expected = [1., 5., 25., 125.]

    g = v.IsoVolume()
    g.generate_levels(N, minN, maxN, ratio=True, log=True)
    assert_equal(levels_expected, g.levels)
    assert_equal(minN, g.minN)
    assert_equal(maxN, g.maxN)
    assert_equal(4, g.N)


def test_generate_volumes_in_range():
    g = v.IsoVolume()
    g.generate_levels(4, 8.e-7, 1.7e-6, log=False)
    db = os.getcwd()+"/test-1"
    g.generate_volumes(ww_file, 'ww_n', dbname=db)

    comp_vols = os.getcwd() + "/vols-1/vols"
    comp_files = [f for f in listdir(comp_vols) if isfile(join(comp_vols, f))]
    gen_vols = db + "/vols"

    # check that files produced are the same
    res = filecmp.cmpfiles(comp_vols, gen_vols, comp_files)
    match_list = res[0]
    non_match = res[1]
    assert_equal(match_list, comp_files)
    assert_equal(non_match, [])

    shutil.rmtree(db)


def test_generate_volumes_min_in_range():
    g = v.IsoVolume()
    g.generate_levels(4, 8.e-7, 3.e-6, log=False)
    db = os.getcwd()+"/test-2"
    g.generate_volumes(ww_file, 'ww_n', dbname=db)

    comp_vols = os.getcwd() + "/vols-2/vols"
    comp_files = [f for f in listdir(comp_vols) if isfile(join(comp_vols, f))]
    gen_vols = db + "/vols"

    # check that files produced are the same
    res = filecmp.cmpfiles(comp_vols, gen_vols, comp_files)
    match_list = res[0]
    non_match = res[1]
    assert_equal(match_list, comp_files)
    assert_equal(non_match, [])

    shutil.rmtree(db)


def test_generate_volumes_max_in_range():
    g = v.IsoVolume()
    g.generate_levels(4, 5.e-7, 1.7e-6, log=False)
    db = os.getcwd()+"/test-3"
    g.generate_volumes(ww_file, 'ww_n', dbname=db)

    comp_vols = os.getcwd() + "/vols-3/vols"
    comp_files = [f for f in listdir(comp_vols) if isfile(join(comp_vols, f))]
    gen_vols = db + "/vols"

    # check that files produced are the same
    res = filecmp.cmpfiles(comp_vols, gen_vols, comp_files)
    match_list = res[0]
    non_match = res[1]
    assert_equal(match_list, comp_files)
    assert_equal(non_match, [])

    shutil.rmtree(db)


def test_generate_volumes_not_in_range():
    g = v.IsoVolume()
    g.generate_levels(4, 5.e-7, 3.e-6, log=False)
    db = os.getcwd()+"/test-4"
    g.generate_volumes(ww_file, 'ww_n', dbname=db)

    comp_vols = os.getcwd() + "/vols-4/vols"
    comp_files = [f for f in listdir(comp_vols) if isfile(join(comp_vols, f))]
    gen_vols = db + "/vols"

    # check that files produced are the same
    res = filecmp.cmpfiles(comp_vols, gen_vols, comp_files)
    match_list = res[0]
    non_match = res[1]
    assert_equal(match_list, comp_files)
    assert_equal(non_match, [])

    shutil.rmtree(db)


# def test_generate_volumes_single_level():
#     # test that two levels are generated if a single value is passed
#     # (code needs work to accomplish this)
#     pass


def test_create_geometry():
    """Test geom creation from existing STL files."""
    # all levels are within range of data

    g = v.IsoVolume()
    g.db = os.getcwd()+"/vols-1"


# def test_create_geometry_full():
#     # beginning to end geom creation (assigning levels)
#     # no mat or viz tags
#     pass
#
#
# def test_create_geometry_full_norm():
#     # add normalization factor
#     pass
#
#
# def test_separate_volumes():
#     pass
#
#
# def test_merge():
#     pass
#
#
# def test_add_mats():
#     # full geom w/ mats
#     pass
#
#
# def test_viz_tag():
#     # full geom w/ visualization tags
#     pass
