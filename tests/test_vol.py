import os
from os import listdir
from os.path import isfile, join
import filecmp
import shutil
import pytest

from IsogeomGenerator import vol


test_dir = os.getcwd() + "/tests/test_files/"
ww_file = test_dir + "cwwm.vtk"


def test_assign_levels():
    """assign predetermined levels (out of order)"""
    g = vol.IsoVolDatabase()
    levels = [0.1, 0.2, 0.05]
    exp = sorted(levels)
    g.assign_levels(levels)
    assert(g.levels == exp)


def test_read_levels():
    """Read levels from file"""
    levelfile = test_dir + "/vols-assign/levelfile"
    exp = [8e-07, 1.2e-06, 1.7e-06]
    g = vol.IsoVolDatabase()
    g.read_levels(levelfile)
    assert(g.levels == exp)


def test_read_levels_nofile():
    """Read levels from nonexisting file"""
    levelfile = ""
    exp = [8e-07, 1.2e-06, 1.7e-06]
    g = vol.IsoVolDatabase()
    with pytest.raises(RuntimeError) as error_info:
        g.read_levels(levelfile)


# Generate Levels parametrized tests:
# linear: (6, 5, 15, 'lin', [5., 7., 9., 11., 13., 15.]
# log: (6, 1, 1e5, 'log', [1, 10, 1e2, 1e3, 1e4, 1e5])
# ratio, max included: (5, 1, 625, 'ratio', [1., 5., 25., 125., 625.])
# ratio, max not included: (5, 1, 700, 'ratio', [1., 5., 25., 125., 625.])
@pytest.mark.parametrize("N,minN,maxN,mode,exp",
                         [(6, 5, 15, 'lin', [5., 7., 9., 11., 13., 15.]),
                          (6, 1, 1e5, 'log', [1, 10, 1e2, 1e3, 1e4, 1e5]),
                          (5, 1, 625, 'ratio', [1., 5., 25., 125., 625.]),
                          (5, 1, 700, 'ratio', [1., 5., 25., 125., 625.])])
def test_generate_levels(N, minN, maxN, mode, exp):
    """generate levels with different modes"""
    g = vol.IsoVolDatabase()
    g.generate_levels(N, minN, maxN, mode=mode)
    assert(g.levels == exp)


def test_generate_levels_error():
    """generate levels with invalid mode"""
    g = vol.IsoVolDatabase()
    with pytest.raises(RuntimeError) as error_info:
        g.generate_levels(6, 5, 1e5, mode='nonsense')


# Generate Volumes parametrized tests:
# Min and Max w/in data bounds: (4, 8.e-7, 1.7e-6, 1)
# Max out of data bounds: (4, 8.e-7, 3.e-6, 2)
# Min out of data bounds: (4, 5.e-7, 1.7e-6, 3)
# Min and Max out of data bounds: (4, 5.e-7, 3.e-6, 4)
@pytest.mark.parametrize("N,minN,maxN,id", [(4, 8.e-7, 1.7e-6, 1),
                                            (4, 8.e-7, 3.e-6, 2),
                                            (4, 5.e-7, 1.7e-6, 3),
                                            (4, 5.e-7, 3.e-6, 4)])
def test_generate_volumes(N, minN, maxN, id):
    """Generate all isovolume files with different bounds"""
    # Expected results:
    exp_vols_dir = test_dir + "/vols-{}/vols".format(id)
    common_files = [f for f in listdir(exp_vols_dir)
                    if isfile(join(exp_vols_dir, f))]
    exp_levelfile = test_dir + "/vols-{}/levelfile".format(id)
    # Generate the volumes
    g = vol.IsoVolDatabase()
    g.generate_levels(N, minN, maxN, mode='lin')
    db = test_dir + "/test-{}".format(id)
    if os.path.isdir(db):
        shutil.rmtree(db)
    # data out of range should produce warning
    with pytest.warns(None) as record:
        g.generate_volumes(ww_file, 'ww_n', dbname=db)
    gen_vols_dir = db + "/vols"
    levelfile = db + "/levelfile"
    # check that files produced are the same
    res = filecmp.cmpfiles(exp_vols_dir, gen_vols_dir, common_files)
    match_list = res[0]
    non_match = res[1]
    assert(match_list == common_files)
    assert(non_match == [])
    # check that the level files are the same
    res = filecmp.cmp(exp_levelfile, levelfile)
    assert(res)
    shutil.rmtree(db)


def test_generate_volumes_assign_levels():
    """Assign levels at time of generating levels"""
    # Expected results:
    exp_vols_dir = test_dir + "/vols-assign/vols"
    common_files = [f for f in listdir(exp_vols_dir)
                    if isfile(join(exp_vols_dir, f))]
    exp_levelfile = test_dir + "/vols-assign/levelfile"
    # Generate the volumes
    g = vol.IsoVolDatabase()
    levels = [8e-07, 1.2e-6, 1.7e-06]
    db = test_dir + "/test-assign"
    if os.path.isdir(db):
        shutil.rmtree(db)
    g.generate_volumes(ww_file, 'ww_n', dbname=db, levels=levels)
    gen_vols_dir = db + "/vols"
    levelfile = db + "/levelfile"
    # check that files produced are the same
    res = filecmp.cmpfiles(exp_vols_dir, gen_vols_dir, common_files)
    match_list = res[0]
    non_match = res[1]
    assert(match_list == common_files)
    assert(non_match == [])
    # check that the level files are the same
    res = filecmp.cmp(exp_levelfile, levelfile)
    assert(res)
    shutil.rmtree(db)


def test_generate_volumes_levelfile():
    """Read levels from file at time of generating levels"""
    # Expected results:
    exp_vols_dir = test_dir + "/vols-assign/vols"
    common_files = [f for f in listdir(exp_vols_dir)
                    if isfile(join(exp_vols_dir, f))]
    exp_levelfile = test_dir + "/vols-assign/levelfile"
    # Generate the volumes
    g = vol.IsoVolDatabase()
    levelfile_in = test_dir + "/vols-assign/levelfile"
    db = test_dir + "/test-read"
    if os.path.isdir(db):
        shutil.rmtree(db)
    g.generate_volumes(ww_file, 'ww_n', dbname=db, levelfile=levelfile_in)
    gen_vols_dir = db + "/vols"
    levelfile_out = db + "/levelfile"
    # check that files produced are the same
    res = filecmp.cmpfiles(exp_vols_dir, gen_vols_dir, common_files)
    match_list = res[0]
    non_match = res[1]
    assert(match_list == common_files)
    assert(non_match == [])
    # check that the level files are the same
    res = filecmp.cmp(exp_levelfile, levelfile_out)
    assert(res)
    shutil.rmtree(db)


def test_generate_volumes_genmode_pass():
    """Generate levels at time of generating volumes - all input included"""
    # Expected results:
    exp_vols_dir = test_dir + "/vols-1/vols"
    common_files = [f for f in listdir(exp_vols_dir)
                    if isfile(join(exp_vols_dir, f))]
    exp_levelfile = test_dir + "/vols-1/levelfile"
    # Generate the volumes
    g = vol.IsoVolDatabase()
    db = test_dir + "/test-gen"
    if os.path.isdir(db):
        shutil.rmtree(db)
    g.generate_volumes(ww_file, 'ww_n', dbname=db, N=4, minN=8.e-7,
                       maxN=1.7e-6, genmode='lin')
    gen_vols_dir = db + "/vols"
    levelfile = db + "/levelfile"
    # check that files produced are the same
    res = filecmp.cmpfiles(exp_vols_dir, gen_vols_dir, common_files)
    match_list = res[0]
    non_match = res[1]
    assert(match_list == common_files)
    assert(non_match == [])
    # check that the level files are the same
    res = filecmp.cmp(exp_levelfile, levelfile)
    assert(res)
    shutil.rmtree(db)


def test_generate_volumes_genmode_error():
    """Generate levels at time of generating volumes - missing input"""
    # Expected results:
    exp_vols_dir = test_dir + "/vols-1/vols"
    common_files = [f for f in listdir(exp_vols_dir)
                    if isfile(join(exp_vols_dir, f))]
    exp_levelfile = test_dir + "/vols-1/levelfile"
    # Generate the volumes
    g = vol.IsoVolDatabase()
    # check error
    with pytest.raises(RuntimeError) as error_info:
        g.generate_volumes(ww_file, 'ww_n', genmode='lin')


def test_generate_volumes_preset():
    """Set levels as member variable before generating vols"""
    # Expected results:
    exp_vols_dir = test_dir + "/vols-assign/vols"
    common_files = [f for f in listdir(exp_vols_dir)
                    if isfile(join(exp_vols_dir, f))]
    exp_levelfile = test_dir + "/vols-assign/levelfile"
    # Generate the volumes
    g = vol.IsoVolDatabase()
    levels = [8e-07, 1.2e-6, 1.7e-06]
    g.levels = levels
    db = test_dir + "/test-set"
    if os.path.isdir(db):
        shutil.rmtree(db)
    g.generate_volumes(ww_file, 'ww_n', dbname=db)
    gen_vols_dir = db + "/vols"
    levelfile = db + "/levelfile"
    # check that files produced are the same
    res = filecmp.cmpfiles(exp_vols_dir, gen_vols_dir, common_files)
    match_list = res[0]
    non_match = res[1]
    assert(match_list == common_files)
    assert(non_match == [])
    # check that the level files are the same
    res = filecmp.cmp(exp_levelfile, levelfile)
    assert(res)
    shutil.rmtree(db)


def test_generate_volumes_init():
    """Set levels in the init"""
    # Expected results:
    exp_vols_dir = test_dir + "/vols-assign/vols"
    common_files = [f for f in listdir(exp_vols_dir)
                    if isfile(join(exp_vols_dir, f))]
    exp_levelfile = test_dir + "/vols-assign/levelfile"
    # Generate the volumes
    levels = [8e-07, 1.2e-6, 1.7e-06]
    g = vol.IsoVolDatabase(levels)
    db = test_dir + "/test-init"
    if os.path.isdir(db):
        shutil.rmtree(db)
    g.generate_volumes(ww_file, 'ww_n', dbname=db)
    gen_vols_dir = db + "/vols"
    levelfile = db + "/levelfile"
    # check that files produced are the same
    res = filecmp.cmpfiles(exp_vols_dir, gen_vols_dir, common_files)
    match_list = res[0]
    non_match = res[1]
    assert(match_list == common_files)
    assert(non_match == [])
    # check that the level files are the same
    res = filecmp.cmp(exp_levelfile, levelfile)
    assert(res)
    shutil.rmtree(db)


def test_generate_volumes_no_levels():
    """Try to generate volumes without assigning levels first (error)"""
    g = vol.IsoVolDatabase()
    with pytest.raises(RuntimeError) as error_info:
        g.generate_volumes(ww_file, 'ww_n')


def test_generate_volumes_dir_exists():
    """Catch warning if database already exists"""
    # Generate the volumes
    g = vol.IsoVolDatabase()
    g.levels = [8e-07, 1.2e-6, 1.7e-06]
    # create exisit dir
    db = test_dir + "test-exist"  # already exists
    if os.path.isdir(db):
        shutil.rmtree(db)
    os.mkdir(db)
    with pytest.warns(None) as record:
        g.generate_volumes(ww_file, 'ww_n', dbname=db)
    new_db = test_dir + "/test-exist-1"
    assert(os.path.isdir(new_db))
    shutil.rmtree(db)
    shutil.rmtree(new_db)


def test__generate_vols():
    pass


def test__get_isovol():
    pass


def test__write_levels():
    pass


# Test Moab functions

# Create a useable IsoVolDatabase object for testing
isovol_obj = vol.IsoVolDatabase()
isovol_obj.completed = True
isovol_obj.data = 'wwn'
isovol_obj.db = ''  # decide which one to use
isovol_obj.levels = []  # depends on which test is decided


# tests for checking if supplied variables are properly handled
def test_create_geom_init_obj():
    """init with an object"""
    pass


def test_create_geom_pass_obj():
    """pass with an object"""
    pass


def test_create_geom_pass_vars():
    """no object, pass variables in function call"""
    pass


def test_create_geom_no_vars_no_obj():
    """no obj or vars provided - should raise error"""
    pass


def test_create_geom_init_vars_and_pass_obj():
    """set vars in the init and pass object in function call - raise warn"""
    pass


def test_create_geom_init_obj_and_pass_vars():
    """init with an object and pass vars - should raise warning"""
    pass


def test_create_geom_incomple_obj():
    """pass incomplete IsoVolDatabase obj - raise error"""
    pass


# test full geometry creation with different options
def test_create_geom_full():
    """test the complete geometry created"""
    pass


def test_create_geom_norm():
    """use a norm factor"""
    pass


def test_create_geom_tags():
    """provided tags are set correctly"""
    pass


def test_create_geom_viz():
    """test viz tags are set"""
    pass


def test__read_levels():
    pass


def test_read_database():
    pass


def test__read_database_error():
    pass


def test__separate_isvols():
    pass


def test_get_surf_tris():
    pass


def test__list_coords():
    pass


def test__list_coords_invert():
    pass


def test__get_matches_exact():
    pass


def test__get_matches_approx():
    pass


def test__compare_surfs():
    pass


def test__imprint_merge():
    pass


def test__make_family():
    pass


def test__tag_for_viz():
    pass


def test__set_tags():
    pass


def test__write_geom():
    pass


def test__write_geom_invalid_ext():
    pass
