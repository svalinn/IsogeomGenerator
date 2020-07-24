"""tests for the IsoSurfGeom module"""
from os import listdir, remove, getcwd, mkdir
from os.path import isfile, isdir, join
import pytest
from pymoab import core, types
import numpy as np

from IsogeomGenerator import isg, ivdb

# Set up test files and expected results
test_dir = getcwd() + "/tests/test_files/"
test_mesh = test_dir + "test_mesh.vtk"
data = 'dname'
levels = [15, 5, 25, 35]
exp_db = test_dir + "/exp-test/"
exp_vols_dir = exp_db + "/vols"
common_files = [f for f in listdir(exp_vols_dir)
                if isfile(join(exp_vols_dir, f))]
exp_levelfile = exp_db + "/levelfile"
exp_levels = [5, 15, 25, 35]
exp_geom = test_dir + '/exp-isogeom.h5m'


def __ivdb_obj(completed):
    iv = ivdb.IvDb(levels=exp_levels, data=data, db=exp_db)
    iv.completed = completed
    return iv


def test_init_none():
    r0 = r1 = r2 = r3 = r4 = False
    ig = isg.IsGm()
    if ig.levels is None:
        r0 = True
    if ig.data is None:
        r1 = True
    if ig.db == getcwd() + "/tmp":
        r2 = True
    if isinstance(ig.mb, type(core.Core())):
        r3 = True
    if ig.isovol_meshsets == {}:
        r4 = True
    assert(all([r0, r1, r2, r3, r4]))


def test_init_input():
    r0 = r1 = r2 = False
    ig = isg.IsGm(levels=levels, data=data, db=exp_db)
    if ig.levels == exp_levels:
        r0 = True
    if ig.data == data:
        r1 = True
    if ig.db == exp_db:
        r2 = True
    assert(all([r0, r1, r2]))


def test_init_ivdb():
    """test that info is taken from ivdb"""
    r0 = r1 = r2 = False
    iv = __ivdb_obj(True)
    ig = isg.IsGm(ivdb=iv)
    if ig.levels == exp_levels:
        r0 = True
    if ig.data == data:
        r1 = True
    if ig.db == exp_db:
        r2 = True
    assert(all([r0, r1, r2]))


def test_init_input_ivdb():
    """test that info from ivdb overwrites other input"""
    r0 = r1 = r2 = False
    iv = __ivdb_obj(True)
    ig = isg.IsGm(ivdb=iv, levels=[0, 2], data='nonsense', db='fake_db')
    if ig.levels == exp_levels:
        r0 = True
    if ig.data == data:
        r1 = True
    if ig.db == exp_db:
        r2 = True
    assert(all([r0, r1, r2]))


def test_init_input_file():
    ig = isg.IsGm(levels=exp_levelfile)
    assert(ig.levels == exp_levels)


def test_read_ivdb():
    """read info from ivdb obj"""
    iv = __ivdb_obj(True)
    ig = isg.IsGm()
    ig.read_ivdb(iv)
    r0 = r1 = r2 = False
    if ig.levels == exp_levels:
        r0 = True
    if ig.data == data:
        r1 = True
    if ig.db == exp_db:
        r2 = True
    assert(all([r0, r1, r2]))


def test_read_ivdb_incomplete():
    """raise error if incomplete ivdb obj"""
    iv = __ivdb_obj(False)
    ig = isg.IsGm()
    with pytest.raises(RuntimeError) as error_info:
        ig.read_ivdb(iv)
    assert "Incomplete IvDb object" in str(error_info)


def test_read_database():
    """check that meshsets are properly populated with read_database"""
    # create obj and read database
    ig = isg.IsGm(levels=levels, data=data, db=exp_db)
    ig.read_database()
    # expected meshset entity handles
    ehs = [12682136550675316737,
           12682136550675316738,
           12682136550675316739,
           12682136550675316740,
           12682136550675316741]
    # setup truth array
    res = np.full(len(ehs) + 1, False)
    # check that meshsets exist in the moab instance
    for r, eh in enumerate(ehs):
        try:
            # any moab call that will work if meshet exists, else
            # it will fail
            ig.mb.get_child_meshsets(eh)
        except RuntimeError:
            pass
        else:
            res[r] = True
    # check meshets and bound information are in dictionary
    exp_meshsets = {(0, ehs[0]): {'bounds': (None, 5.0)},
                    (1, ehs[1]): {'bounds': (5.0, 15.0)},
                    (2, ehs[2]): {'bounds': (15.0, 25.0)},
                    (3, ehs[3]): {'bounds': (25.0, 35.0)},
                    (4, ehs[4]): {'bounds': (35.0, None)}}
    if sorted(ig.isovol_meshsets) == sorted(exp_meshsets):
        res[-1] = True
    # assert all pass
    assert(all(res))


def test_read_database_numfiles_error():
    """read_database throws error if num levels and files mismatch"""
    # create obj and read database
    ig = isg.IsGm(levels=[300], data=data, db=exp_db)
    with pytest.raises(RuntimeError) as error_info:
        ig.read_database()
    assert "does not match number" in str(error_info)


def test_read_database_nolevels_error():
    """read_database throws error no levels are defined"""
    # create obj and read database
    ig = isg.IsGm()
    with pytest.raises(RuntimeError) as error_info:
        ig.read_database()
    assert "levels defined" in str(error_info)


def test_separate_isovols():
    # load mesh that needs separation
    ig = isg.IsGm()
    fs = ig.mb.create_meshset()
    ig.mb.load_file(test_dir + '/vol-files/separate-vols.stl', file_set=fs)
    # create useable meshset dict
    ig.isovol_meshsets[(0, fs)] = {}
    # separate the volumes
    ig.separate_isovols()
    # check there are two new surfaces
    r = np.full(3, False)
    num_surfs = len(ig.isovol_meshsets[(0, fs)]['surfs_EH'])
    if num_surfs == 2:
        r[0] = True
    # check that no triangles or vertices are shared between the surfs
    surf0 = ig.isovol_meshsets[(0, fs)]['surfs_EH'][0]
    verts0 = set(ig.mb.get_entities_by_type(surf0, types.MBVERTEX))
    tris0 = set(ig.mb.get_entities_by_type(surf0, types.MBTRI))
    surf1 = ig.isovol_meshsets[(0, fs)]['surfs_EH'][1]
    verts1 = set(ig.mb.get_entities_by_type(surf1, types.MBVERTEX))
    tris1 = set(ig.mb.get_entities_by_type(surf1, types.MBTRI))
    common_verts = verts0 & verts1
    common_tris = tris0 & tris1
    if len(common_verts) == 0:
        r[1] = True
    if len(common_tris) == 0:
        r[2] = True
    assert(all(r))


def test_separate_isovols_single():
    # load mesh that does not need separation
    ig = isg.IsGm()
    fs = ig.mb.create_meshset()
    ig.mb.load_file(test_dir + '/vol-files/single-box-1.stl', file_set=fs)
    # create useable meshset dict
    ig.isovol_meshsets[(0, fs)] = {}
    # separate the volumes
    ig.separate_isovols()
    # check there are two new surfaces
    r = np.full(3, False)
    num_surfs = len(ig.isovol_meshsets[(0, fs)]['surfs_EH'])
    if num_surfs == 1:
        r[0] = True
    # check number of triangles and vertices in surfaces (8 verts, 12 tris)
    surf = ig.isovol_meshsets[(0, fs)]['surfs_EH'][0]
    verts = ig.mb.get_entities_by_type(surf, types.MBVERTEX)
    tris = ig.mb.get_entities_by_type(surf, types.MBTRI)
    if len(verts) == 8:
        r[1] = True
    if len(tris) == 12:
        r[2] = True
    assert(all(r))


def test_imprint_merge():
    pass


def test_make_family():
    pass


def test_tag_for_viz():
    pass


def test_set_tags():
    pass


def test_write_geometry():
    pass


def test_get_surf_triangles():
    """get triangles when one coord is not good"""
    # setup IsGm instance
    ig = isg.IsGm()
    ig.mb = core.Core()
    # make good and bad coords/verts
    all_coords = np.array((0., 0., 0.,
                           1., 0., 0.,
                           .5, 1., 0.,
                           2., 1., 0.), dtype='float64')
    ig.mb.create_vertices(all_coords)
    rs = ig.mb.get_root_set()
    verts_all = ig.mb.get_entities_by_type(rs, types.MBVERTEX)
    verts_good = verts_all[0:3]
    verts_bad = verts_all[1:]
    # make good triangle with coords:
    tri_good = ig.mb.create_element(types.MBTRI, verts_good)
    # make bad triangle:
    tri_bad = ig.mb.create_element(types.MBTRI, verts_bad)
    # test get surf tris
    tris_out = ig._IsGm__get_surf_triangles(verts_good)
    assert(tris_out == tri_good)


def test_get_surf_triangles_all():
    """get triangles when all coords are good"""
    # setup IsGm instance
    ig = isg.IsGm()
    ig.mb = core.Core()
    # make two sets of verts
    all_coords = np.array((0., 0., 0.,
                           1., 0., 0.,
                           .5, 1., 0.,
                           2., 1., 0.), dtype='float64')
    ig.mb.create_vertices(all_coords)
    rs = ig.mb.get_root_set()
    verts_all = ig.mb.get_entities_by_type(rs, types.MBVERTEX)
    verts1 = verts_all[0:3]
    verts2 = verts_all[1:]
    # make two triangles
    tri1 = ig.mb.create_element(types.MBTRI, verts1)
    tri2 = ig.mb.create_element(types.MBTRI, verts2)
    # test get surf tris
    tris_out = ig._IsGm__get_surf_triangles(verts_all)
    assert(sorted(tris_out) == sorted([tri1, tri2]))


def test_list_coords():
    """test list coords returns correct dictionary"""
    # setup IsGm instance
    ig = isg.IsGm()
    ig.mb = core.Core()
    # create some vertices
    coords = np.array((1., 2., 3.), dtype='float64')
    ig.mb.create_vertices(coords)
    rs = ig.mb.get_root_set()
    eh = list(ig.mb.get_entities_by_type(rs, types.MBVERTEX))[0]
    exp_coords_dict = {eh: (1., 2., 3.)}
    # list coords
    coords_out = ig._IsGm__list_coords(rs)
    assert(exp_coords_dict == coords_out)


def test_list_coords_invert():
    """test list coords correctly returns inverted dictionary"""
    # setup IsGm instance
    ig = isg.IsGm()
    ig.mb = core.Core()
    # create some vertices
    coords = np.array((1., 2., 3.), dtype='float64')
    ig.mb.create_vertices(coords)
    rs = ig.mb.get_root_set()
    eh = list(ig.mb.get_entities_by_type(rs, types.MBVERTEX))[0]
    exp_coords_dict = {(1., 2., 3.): eh}
    # list coords
    coords_out = ig._IsGm__list_coords(rs, invert=True)
    assert(exp_coords_dict == coords_out)


def test_get_matches():
    # create instance
    ig = isg.IsGm()
    # set up verts to test:
    vertsA = {10: (1., 2., 3.), 20: (4., 5., 6.), 30: (7., 8., 9.)}
    vertsB = {(7., 8., 9.): 40, (2., 1., 3.): 50}
    merge_tol = 1.e-5

    a_eh, a_coords, b_eh, b_coords, match_dict = \
        ig._IsGm__get_matches(vertsA, vertsB, merge_tol)

    # expected values
    a_eh_exp = [30]
    a_coords_exp = [(7., 8., 9.)]
    b_eh_exp = [40]
    b_coords_exp = [(7., 8., 9.)]
    match_dict_exp = {40: 30}

    res = np.full(5, False)
    if a_eh == a_eh_exp:
        res[0] = True
    if a_coords == a_coords_exp:
        res[1] = True
    if b_eh == b_eh_exp:
        res[2] = True
    if b_coords == b_coords_exp:
        res[3] = True
    if match_dict == match_dict_exp:
        res[4] = True

    assert(all(res))


def test_get_matches_approx():
    # create instance
    ig = isg.IsGm()
    # set up verts to test:
    vertsA = {10: (1., 2., 3.), 20: (4., 5., 6.), 30: (7., 8., 9.)}
    vertsB = {(7., 8., 9.01): 40, (2., 1., 3.): 50}
    merge_tol = 1.e-1

    a_eh, a_coords, b_eh, b_coords, match_dict = \
        ig._IsGm__get_matches(vertsA, vertsB, merge_tol)

    # expected values
    a_eh_exp = [30]
    a_coords_exp = [(7., 8., 9.)]
    b_eh_exp = [40]
    b_coords_exp = [(7., 8., 9.01)]
    match_dict_exp = {40: 30}

    res = np.full(5, False)
    if a_eh == a_eh_exp:
        res[0] = True
    if a_coords == a_coords_exp:
        res[1] = True
    if b_eh == b_eh_exp:
        res[2] = True
    if b_coords == b_coords_exp:
        res[3] = True
    if match_dict == match_dict_exp:
        res[4] = True

    assert(all(res))


def test_compare_surfs():
    pass
