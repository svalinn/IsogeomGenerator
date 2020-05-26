import os
from os import listdir, remove
from os.path import isfile, join
import filecmp
import shutil
import pytest

from IsogeomGenerator import vol


# Set up test files and expected results
test_dir = os.getcwd() + "/tests/test_files/"
test_mesh = test_dir + "test_mesh.vtk"
dataname = 'dname'
exp_db = test_dir + "/exp-test/"
exp_vols_dir = exp_db + "/vols"
common_files = [f for f in listdir(exp_vols_dir)
                if isfile(join(exp_vols_dir, f))]
exp_levelfile = exp_db + "/levelfile"
exp_levels = [5, 15, 25, 35]
exp_geom = test_dir + '/exp-isogeom.h5m'


def __compare_levelfiles(f1, f2):
    """compares the levelfile contents regardless of string format"""
    fs1 = open(f1, 'r')
    fs2 = open(f2, 'r')

    levels1 = []
    lines = fs1.readlines()
    for line in lines:
        levels1.append(float(line))

    levels2 = []
    lines = fs2.readlines()
    for line in lines:
        levels2.append(float(line))

    if levels1 == levels2:
        return True
    else:
        return False


def test_assign_levels():
    """assign predetermined levels (out of order)"""
    g = vol.IsoVolDatabase()
    levels = [0.1, 0.2, 0.05]
    exp = sorted(levels)
    g.assign_levels(levels)
    assert(g.levels == exp)


def test_read_levels():
    """Read levels from file"""
    g = vol.IsoVolDatabase()
    g.read_levels(exp_levelfile)
    assert(g.levels == exp_levels)


def test_read_levels_nofile():
    """Read levels from nonexisting file"""
    levelfile = ""
    g = vol.IsoVolDatabase()
    with pytest.raises(RuntimeError) as error_info:
        g.read_levels(levelfile)
    assert 'does not exist' in str(error_info)


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
    assert 'Level generation' in str(error_info)


# Generate Volumes parametrized tests:
#   (1) Min and Max w/in data bounds
#   (2) Mid level no data
#   (3) Min and max out of data bounds
@pytest.mark.parametrize("levels,id", [([5, 15, 25, 35], 1),
                                       ([5, 15, 25, 28, 35], 2),
                                       ([-5, 5, 15, 25, 35, 45], 3)])
def test_generate_volumes(levels, id):
    """Generate all isovolume files with different bounds.
    Some should raise warnings. Make sure generated files are correct."""
    # Generate the volumes
    g = vol.IsoVolDatabase()
    db = test_dir + "/test-gen-vols-{}".format(id)
    if os.path.isdir(db):
        shutil.rmtree(db)
    # data out of range should produce warning
    with pytest.warns(None) as warn_info:
        g.generate_volumes(test_mesh, dataname, dbname=db, levels=levels)
    # assert number of warnings raised
    if id == 2:
        # no data found warning
        assert(len(warn_info) == 1)
    elif id == 3:
        # levels are out of bounds
        assert(len(warn_info) == 2)
    else:
        assert(len(warn_info) == 0)
    # check that files produced are the same
    gen_vols_dir = db + "/vols"
    levelfile_out = db + "/levelfile"
    res = filecmp.cmpfiles(exp_vols_dir, gen_vols_dir, common_files)
    match_list = res[0]
    non_match = res[1]
    assert(match_list == common_files)
    assert(non_match == [])
    # check that the level files are the same
    res = __compare_levelfiles(exp_levelfile, levelfile_out)
    assert(res)
    shutil.rmtree(db)


@pytest.mark.filterwarnings("ignore:Level")
def test_generate_volumes_levels_outofbounds():
    """level values provided are not within data bounds, so raise error"""
    # Generate the volumes
    g = vol.IsoVolDatabase()
    # data out of range should produce warning
    with pytest.raises(RuntimeError) as error_info:
        g.generate_volumes(test_mesh, dataname, levels=[-5, 0, 45])
    assert "No data exists" in str(error_info)


def test_generate_volumes_levelfile():
    """Read levels from file at time of generating levels"""
    # Generate the volumes
    g = vol.IsoVolDatabase()
    levelfile_in = test_dir + "/exp-test/levelfile"
    db = test_dir + "/test-read-levels"
    if os.path.isdir(db):
        shutil.rmtree(db)
    g.generate_volumes(test_mesh, dataname, dbname=db, levelfile=levelfile_in)
    levelfile_out = db + "/levelfile"
    # check that levels are read and written correctly
    assert(g.levels == exp_levels)
    # check that the level files are the same
    res = __compare_levelfiles(exp_levelfile, levelfile_out)
    assert(res)
    shutil.rmtree(db)


def test_generate_volumes_genmode_pass():
    """Generate levels at time of generating volumes - all input included"""
    # Generate the volumes
    g = vol.IsoVolDatabase()
    db = test_dir + "/test-gen"
    if os.path.isdir(db):
        shutil.rmtree(db)
    g.generate_volumes(test_mesh, dataname, dbname=db, N=4, minN=5,
                       maxN=35, genmode='lin')
    levelfile_out = db + "/levelfile"
    # check that levels are read and written correctly
    assert(g.levels == exp_levels)
    # check that the level files are the same
    res = __compare_levelfiles(exp_levelfile, levelfile_out)
    assert(res)
    shutil.rmtree(db)


def test_generate_volumes_genmode_error():
    """Generate levels at time of generating volumes - missing input"""
    # Generate the volumes
    g = vol.IsoVolDatabase()
    # check error
    with pytest.raises(RuntimeError) as error_info:
        g.generate_volumes(test_mesh, dataname, genmode='lin')
    assert 'requires input values' in str(error_info)


def test_generate_volumes_init():
    """Set levels in the init"""
    # Generate the volumes
    levels = [5, 15, 25, 35]
    g = vol.IsoVolDatabase(levels)
    db = test_dir + "/test-init"
    if os.path.isdir(db):
        shutil.rmtree(db)
    g.generate_volumes(test_mesh, dataname, dbname=db)
    levelfile_out = db + "/levelfile"
    # check that levels are read and written correctly
    assert(g.levels == exp_levels)
    # check that the level files are the same
    res = __compare_levelfiles(exp_levelfile, levelfile_out)
    assert(res)
    shutil.rmtree(db)


def test_generate_volumes_no_levels():
    """Try to generate volumes without assigning levels first (error)"""
    g = vol.IsoVolDatabase()
    with pytest.raises(RuntimeError) as error_info:
        g.generate_volumes(test_mesh, dataname)
    assert 'levels must be provided' in str(error_info)


def test_generate_volumes_dir_exists():
    """Catch warning if database already exists"""
    # Generate the volumes
    g = vol.IsoVolDatabase()
    g.levels = [5, 15, 25, 35]
    # create exisit dir
    db = test_dir + "/test-exist"  # already exists
    if os.path.isdir(db):
        shutil.rmtree(db)
    os.mkdir(db)
    with pytest.warns(None) as warn_info:
        g.generate_volumes(test_mesh, dataname, dbname=db)
    assert(len(warn_info) == 1)
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

def __create_isvolobj(complete):
    """Create a useable IsoVolDatabase object for testing"""
    isovolobj = vol.IsoVolDatabase()
    isovolobj.completed = complete
    isovolobj.data = dataname
    isovolobj.db = exp_db
    isovolobj.levels = exp_levels
    return isovolobj


# tests for checking if supplied variables are properly handled
def test_init_obj():
    """init with an object"""
    ivo = __create_isvolobj(True)
    g = vol.IsoSurfGeom(isovoldbobj=ivo)
    assert(g.data == dataname)
    assert(g.db == exp_db)
    assert(g.levels == exp_levels)


def test_init_obj_multiple():
    """init with an object and redundant info. Test for correct variables and
    warnings are raised."""
    ivo = __create_isvolobj(True)
    with pytest.warns(None) as warn_info:
        g = vol.IsoSurfGeom(isovoldbobj=ivo, data='nonsense', dbname='fake_db')
    assert(len(warn_info) == 2)
    assert(g.data == dataname)
    assert(g.db == exp_db)
    assert(g.levels == exp_levels)


def test_init_obj_incomplete():
    """Try to init with an obj that was not completed. Should raise error."""
    ivo = __create_isvolobj(False)
    with pytest.raises(RuntimeError) as error_info:
        g = vol.IsoSurfGeom(isovoldbobj=ivo)
    assert "Incomplete IsoVolDatabase object" in str(error_info)


def test_create_geom():
    """Test geometry is created properly when using defaults"""
    ivo = __create_isvolobj(True)
    g = vol.IsoSurfGeom()
    g.create_geometry(isovoldbobj=ivo)
    geom_file = exp_db + '/isogeom.h5m'
    # compare h5m files, ignoring timestamps
    diffs = os.system('h5diff --exclude-path /tstt/history ' +
                      '{} {}'.format(exp_geom, geom_file))
    #assert(diffs == 0)
    remove(geom_file)


def test_create_geom_pass_obj():
    """test values are set properly when object is passed"""
    ivo = __create_isvolobj(True)
    g = vol.IsoSurfGeom()
    g.create_geometry(isovoldbobj=ivo)
    assert(g.data == dataname)
    assert(g.db == exp_db)
    assert(g.levels == exp_levels)
    remove(exp_db + '/isogeom.h5m')


def test_create_geom_init_and_pass_obj():
    """called with an object after already having set in the initself.
    Only test that variables are properly set and warnings are raised.
    """
    # create object that will be overwritten
    init_obj = vol.IsoVolDatabase()
    init_obj.completed = True
    init_obj.data = 'fake_name'
    init_obj.db = 'fake_db'
    init_obj.levels = [-1, 1]
    # correct object
    ivo = __create_isvolobj(True)
    g = vol.IsoSurfGeom(isovoldbobj=init_obj)
    with pytest.warns(None) as warn_info:
        g.create_geometry(isovoldbobj=ivo)
    assert(len(warn_info) == 3)  # initial warning, plus data and db variables
    assert(g.data == dataname)
    assert(g.db == exp_db)
    assert(g.levels == exp_levels)
    remove(exp_db + '/isogeom.h5m')


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
