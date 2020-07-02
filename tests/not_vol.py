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
    # assert(diffs == 0)
    remove(geom_file)


def test_create_geom_pass_obj():
    """test values are set properly when object is passed"""
    ivo = __create_isvolobj(True)
    g = vol.IsoSurfGeom()
    g.create_geometry(isovoldbobj=ivo)
    assert(g.data == dataname)
    assert(g.db == exp_db)
    assert(g.levels == exp_levels)
    assert(isfile(exp_db + '/isogeom.h5m'))
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
    assert(isfile(exp_db + '/isogeom.h5m'))
    remove(exp_db + '/isogeom.h5m')


def test_create_geom_pass_vars():
    """no object, pass variables in function call"""
    g = vol.IsoSurfGeom()
    g.create_geometry(data=dataname, dbname=exp_db)
    assert(g.data == dataname)
    assert(g.db == exp_db)
    assert(g.levels == exp_levels)
    assert(isfile(exp_db + '/isogeom.h5m'))
    remove(exp_db + '/isogeom.h5m')


def test_create_geom_no_data():
    g = vol.IsoSurfGeom()
    with pytest.raises(RuntimeError) as error_info:
        g.create_geometry()
    assert "Variable 'data' must be provided" in str(error_info)


def test_create_geom_db_nonexistent():
    g = vol.IsoSurfGeom()
    with pytest.raises(RuntimeError) as error_info:
        g.create_geometry(data=dataname)
    assert "does not contain an isovolume database" in str(error_info)
    assert(g.db == os.getcwd() + '/tmp')


def test_create_geom_init_vars_and_pass_obj():
    """set vars in the init and pass object in function call - raise warn"""
    # create object that will be overwritten
    g = vol.IsoSurfGeom(data='fake_data', dbname='fake_db')
    ivo = __create_isvolobj(True)
    with pytest.warns(None) as warn_info:
        g.create_geometry(isovoldbobj=ivo)
    assert(len(warn_info) == 2)
    assert(g.data == dataname)
    assert(g.db == exp_db)
    assert(g.levels == exp_levels)
    assert(isfile(exp_db + '/isogeom.h5m'))
    remove(exp_db + '/isogeom.h5m')


def test_create_geom_init_and_pass_vars():
    """init with an object and pass vars - should raise warning"""
    g = vol.IsoSurfGeom(data='fake_data', dbname='fake_db')
    with pytest.warns(None) as warn_info:
        g.create_geometry(data=dataname, dbname=exp_db)
    assert(len(warn_info) == 2)  # initial warning, plus data and db variables
    assert(g.data == dataname)
    assert(g.db == exp_db)
    assert(g.levels == exp_levels)
    assert(isfile(exp_db + '/isogeom.h5m'))
    remove(exp_db + '/isogeom.h5m')


def test_create_geom_init_obj_and_pass_vars():
    """init with an object and also pass variables - warnings"""
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
        g.create_geometry(data=dataname, dbname=exp_db)
    assert(len(warn_info) == 2)
    assert(g.data == dataname)
    assert(g.db == exp_db)
    assert(g.levels == exp_levels)
    assert(isfile(exp_db + '/isogeom.h5m'))
    remove(exp_db + '/isogeom.h5m')


def test__read_isovol():
    ivo = __create_isvolobj(True)
    g = vol.IsoSurfGeom()
    g.isovoldbobj = ivo
    g._IsoSurfGeom__read_isovol(None, None)
    assert(g.data == dataname)
    assert(g.db == exp_db)
    assert(g.levels == exp_levels)


def test__read_isovol_incomplete():
    ivo = __create_isvolobj(False)
    g = vol.IsoSurfGeom()
    g.isovoldbobj = ivo
    with pytest.raises(RuntimeError) as error_info:
        g._IsoSurfGeom__read_isovol(None, None)
    assert "Incomplete IsoVolDatabase object" in str(error_info)


def test__read_isovol_var_exists_init():
    ivo = __create_isvolobj(True)
    g = vol.IsoSurfGeom()
    g.isovoldbobj = ivo
    g.data = 'fake_data'
    g.db = 'fake_db'
    with pytest.warns(None) as warn_info:
        g._IsoSurfGeom__read_isovol(None, None)
    assert(len(warn_info) == 2)
    assert(g.data == dataname)
    assert(g.db == exp_db)
    assert(g.levels == exp_levels)


def test__read_isovol_var_exists_pass():
    ivo = __create_isvolobj(True)
    g = vol.IsoSurfGeom()
    g.isovoldbobj = ivo
    with pytest.warns(None) as warn_info:
        g._IsoSurfGeom__read_isovol('fake_data', 'fake_db')
    assert(len(warn_info) == 2)
    assert(g.data == dataname)
    assert(g.db == exp_db)
    assert(g.levels == exp_levels)


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
