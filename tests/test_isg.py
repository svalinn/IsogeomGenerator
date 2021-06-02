"""tests for the IsoSurfGeom module"""
from os import listdir, remove, getcwd, mkdir
from os.path import isfile, isdir, join
import pytest
from pymoab import core, types
import numpy as np
import itertools
import warnings

from IsogeomGenerator import isg, ivdb

# Set up test files and expected results
test_dir = getcwd() + "/tests/test_files/"
test_mesh = test_dir + "test_mesh.vtk"
data = 'dname'
levels = [15, 5, 25, 35, 45]
exp_db = test_dir + "/exp-test/"
exp_vols_dir = exp_db + "/vols"
common_files = [f for f in listdir(exp_vols_dir)
                if isfile(join(exp_vols_dir, f))]
exp_levelfile = exp_db + "/levelfile"
exp_levels = [5, 15, 25, 35, 45]
exp_geom = test_dir + '/exp-isogeom.h5m'
# geometric extent info
exp_ext_min = -10.
exp_ext_max = 10.
exts = [np.full(3, exp_ext_min), np.full(3, exp_ext_max)]


def __ivdb_obj(completed):
    # manually generated a usuable ivdb object
    iv = ivdb.IvDb(levels=levels, data=data, db=exp_db)
    iv.xmin = iv.ymin = iv.zmin = exp_ext_min
    iv.xmax = iv.ymax = iv.zmax = exp_ext_max
    iv.completed = completed
    return iv


def test_init_none():
    r = np.full(6, False)
    ig = isg.IsGm()
    if ig.levels is None:
        r[0] = True
    if ig.data is None:
        r[1] = True
    if ig.db == getcwd() + "/tmp":
        r[2] = True
    if isinstance(ig.mb, type(core.Core())):
        r[3] = True
    if ig.isovol_meshsets == {}:
        r[4] = True
    if ig.xmin == ig.xmax == ig.ymin == ig.ymax == ig.zmin == ig.zmax == None:
        r[5] = True
    assert(all(r))


def test_init_input():
    r = np.full(5, False)
    ig = isg.IsGm(levels=levels, data=data, db=exp_db, extents=exts)
    if ig.levels == exp_levels:
        r[0] = True
    if ig.data == data:
        r[1] = True
    if ig.db == exp_db:
        r[2] = True
    if ig.xmin == ig.ymin == ig.zmin == exp_ext_min:
        r[3] = True
    if ig.xmax == ig.ymax == ig.zmax == exp_ext_max:
        r[4] = True
    assert(all(r))


def test_init_ivdb():
    """test that info is taken from ivdb"""
    r = np.full(5, False)
    iv = __ivdb_obj(True)
    ig = isg.IsGm(ivdb=iv)
    if ig.levels == exp_levels:
        r[0] = True
    if ig.data == data:
        r[1] = True
    if ig.db == exp_db:
        r[2] = True
    if ig.xmin == ig.ymin == ig.zmin == exp_ext_min:
        r[3] = True
    if ig.xmax == ig.ymax == ig.zmax == exp_ext_max:
        r[4] = True
    assert(all(r))


def test_init_input_ivdb():
    """test that info from ivdb overwrites other input"""
    r = np.full(5, False)
    iv = __ivdb_obj(True)
    ig = isg.IsGm(ivdb=iv, levels=[0, 2], data='nonsense', db='fake_db')
    if ig.levels == exp_levels:
        r[0] = True
    if ig.data == data:
        r[1] = True
    if ig.db == exp_db:
        r[2] = True
    if ig.xmin == ig.ymin == ig.zmin == exp_ext_min:
        r[3] = True
    if ig.xmax == ig.ymax == ig.zmax == exp_ext_max:
        r[4] = True
    assert(all(r))


def test_init_input_file():
    ig = isg.IsGm(levels=exp_levelfile)
    assert(ig.levels == exp_levels)


def test_read_ivdb():
    """read info from ivdb obj"""
    iv = __ivdb_obj(True)
    ig = isg.IsGm()
    ig.read_ivdb(iv)
    r = np.full(5, False)
    if ig.levels == exp_levels:
        r[0] = True
    if ig.data == data:
        r[1] = True
    if ig.db == exp_db:
        r[2] = True
    if ig.xmin == ig.ymin == ig.zmin == exp_ext_min:
        r[3] = True
    if ig.xmax == ig.ymax == ig.zmax == exp_ext_max:
        r[4] = True
    assert(all(r))


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


def test_separate_isovols_exterior():
    """test that disjoint volumes are properly separated"""
    # load mesh that needs separation
    ig = isg.IsGm()
    fs = ig.mb.create_meshset()
    ig.mb.load_file(test_dir + '/vol-files/separate-vols.stl', file_set=fs)
    # create useable meshset dict
    ig.isovol_meshsets[(0, fs)] = {}
    # manually set the geometric extents
    # these are chosen such that the volume file aligns on the x plane
    # geometric extents (-10, 15). The volume file y and z are -5 to 5,
    # so if this is considered to be one volume in a larger geometry,
    # only one the surfaces on the x planes are considered exterior.
    ig.xmin = -10.
    ig.xmax = 15.
    ig.ymin = ig.zmin = -15.
    ig.ymax = ig.zmax = 15.
    # separate the volumes
    ig.separate_isovols()
    # check there are four new surfaces
    r = np.full(2, False)
    num_surfs = len(ig.isovol_meshsets[(0, fs)]['surfs_EH'])
    if num_surfs == 4:
        r[0] = True
    # check that no triangles are shared between the each of the surfaces
    surf0 = ig.isovol_meshsets[(0, fs)]['surfs_EH'][0]
    tris0 = set(ig.mb.get_entities_by_type(surf0, types.MBTRI))
    surf1 = ig.isovol_meshsets[(0, fs)]['surfs_EH'][1]
    tris1 = set(ig.mb.get_entities_by_type(surf1, types.MBTRI))
    surf2 = ig.isovol_meshsets[(0, fs)]['surfs_EH'][2]
    tris2 = set(ig.mb.get_entities_by_type(surf2, types.MBTRI))
    surf3 = ig.isovol_meshsets[(0, fs)]['surfs_EH'][3]
    tris3 = set(ig.mb.get_entities_by_type(surf3, types.MBTRI))
    common_tris = [list(tris0 & tris1), list(tris0 & tris2),
                   list(tris0 & tris3), list(tris1 & tris2),
                   list(tris1 & tris3), list(tris2 & tris3)]
    if not common_tris == 0:
        r[1] = True
    assert(all(r))


def test_separate_isovols_single_exterior():
    """test a single vol with an exterior surface is split in separation"""
    # load mesh that does not need separation
    ig = isg.IsGm()
    fs = ig.mb.create_meshset()
    ig.mb.load_file(test_dir + '/vol-files/single-box-1.stl', file_set=fs)
    # create useable meshset dict
    ig.isovol_meshsets[(0, fs)] = {}
    # manually set the geometric extents
    # these are chosen such that the volume file aligns on the -x plane
    # geometric extents (-5). The volume file x, y, and z are -5 to 5,
    # so if this is considered to be one volume in a larger geometry,
    # only one the surfaces on the -x plane is considered exterior.
    ig.xmin = -5.
    ig.xmax = 15.
    ig.ymin = ig.zmin = -15.
    ig.ymax = ig.zmax = 15.
    # separate the volumes
    ig.separate_isovols()
    # check there are two new surfaces
    r = np.full(2, False)
    num_surfs = len(ig.isovol_meshsets[(0, fs)]['surfs_EH'])
    if num_surfs == 2:
        r[0] = True
    # check number of triangles and vertices in surfaces (8 verts, 12 tris)
    # check that no triangles are shared between the each of the surfaces
    surf0 = ig.isovol_meshsets[(0, fs)]['surfs_EH'][0]
    tris0 = set(ig.mb.get_entities_by_type(surf0, types.MBTRI))
    surf1 = ig.isovol_meshsets[(0, fs)]['surfs_EH'][1]
    tris1 = set(ig.mb.get_entities_by_type(surf1, types.MBTRI))
    common_tris = tris0 & tris1
    if len(common_tris) == 0:
        r[1] = True
    assert(all(r))


def test_separate_isovols_single_interior():
    """test a single interior vol is unchanged when it goes through separation"""
    # load mesh that does not need separation
    ig = isg.IsGm()
    fs = ig.mb.create_meshset()
    ig.mb.load_file(test_dir + '/vol-files/single-box-1.stl', file_set=fs)
    # create useable meshset dict
    ig.isovol_meshsets[(0, fs)] = {}
    # manually set the geometric extents so that no surface is on the
    # exterior
    ig.xmin = -15.
    ig.xmax = 15.
    ig.ymin = ig.zmin = -15.
    ig.ymax = ig.zmax = 15.
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


def __setup_geom():
    """function for other tests to create a useable isogeom object"""
    # load two coincident volumes that need merging
    ig = isg.IsGm()
    fs1 = ig.mb.create_meshset()
    ig.mb.load_file(test_dir + '/vol-files/single-box-1.stl', file_set=fs1)
    fs2 = ig.mb.create_meshset()
    ig.mb.load_file(test_dir + '/vol-files/single-box-2.stl', file_set=fs2)
    # create useable meshset dict
    # assign arbitrary values for levels 0, 5, 10
    iv1 = (0, fs1)
    iv2 = (1, fs2)
    ig.isovol_meshsets[iv1] = {}
    ig.isovol_meshsets[iv1]['bounds'] = (0., 5.)
    ig.isovol_meshsets[iv2] = {}
    ig.isovol_meshsets[iv2]['bounds'] = (5., 10.)
    ig.levels = [5., 10.]
    ig.data = 'dataname'
    ig.xmin = -5.
    ig.xmax = 15.
    ig.ymin = ig.zmin = -5.
    ig.ymax = ig.zmax = 5.
    ig.separate_isovols()  # use this to get surface EHs
    return ig


def test_imprint_merge():
    """test mesh imprint and merge capability"""
    # get setup
    ig = __setup_geom()
    ivs = sorted(ig.isovol_meshsets.keys())  # isovol info: (vol id, EH)
    iv1 = ivs[0]
    iv2 = ivs[1]
    fs1 = list(iv1)[1]
    fs2 = list(iv2)[1]
    # imprint and merge
    norm = 1.5
    ig.imprint_merge(norm)
    # checks
    r = np.full(2, False)  # truth array for checks
    surfs_1 = ig.isovol_meshsets[iv1]['surfs_EH']
    surfs_2 = ig.isovol_meshsets[iv2]['surfs_EH']
    all_surfs = list(set(surfs_1).union(set(surfs_2)))
    common_surf = list(set(surfs_1) & set(surfs_2))[0]
    # check the tags (all surfaces should have tags now)
    # sense tag:
    #   common surf should have both volumes [vol1, vol2]
    #   other surfaces should have only one vol (bwd sense is 0): [vol, 0]
    # value tag:
    #   common surface should be 7.5 (shared bound = 5, norm=1.5 -> 5*1.5=7.5)
    #   other surfaces should be 0
    val_exp = 7.5
    sense_exp = [fs1, fs2]
    tmp_val = [False, False, False]
    tmp_sense = [False, False, False]
    ivs = sorted(ig.isovol_meshsets.keys())
    vols = [fs1, fs2]  # get EHs
    for i, surf in enumerate(all_surfs):
        val_out = ig.mb.tag_get_data(ig.val_tag, surf)[0][0]
        sense_out = list(ig.mb.tag_get_data(ig.sense_tag, surf)[0])
        if surf == common_surf:
            if val_out == val_exp:
                tmp_val[i] = True
            if sense_out == sense_exp:
                tmp_sense[i] = True
        else:
            if val_out == 0.0:
                tmp_val[i] = True
            if (sense_out[0] in vols) and (sense_out[1] == np.uint64(0)):
                tmp_sense[i] = True
                vols.remove(sense_out[0])
    if all(tmp_val):
        r[0] = True
    if all(tmp_sense):
        r[1] = True
    assert(all(r))


def test_make_family():
    """test tags are added properly"""
    # get setup
    ig = __setup_geom()
    ivs = sorted(ig.isovol_meshsets.keys())  # isovol info: (vol id, EH)
    iv1 = ivs[0]
    iv2 = ivs[1]
    fs1 = list(iv1)[1]
    fs2 = list(iv2)[1]
    norm = 1.5
    ig.imprint_merge(norm)  # do this to get shared surfaces
    # run make_family
    ig.make_family()
    # checks:
    r = np.full(9, False)
    # Check Tags:
    #   geom dim: vols=3, surfs=2, curves=1
    #   category: 'Volume', 'Surface', 'Curve'
    #   global id: 1..N for each vol, surf, curve
    all_vols = [fs1, fs2]
    surfs_1 = ig.isovol_meshsets[iv1]['surfs_EH']
    surfs_2 = ig.isovol_meshsets[iv2]['surfs_EH']
    all_surfs = list(set(surfs_1 + surfs_2))
    common_surf = list(set(surfs_1) & set(surfs_2))[0]
    # tag EHs
    geom_dim = ig.mb.tag_get_handle('GEOM_DIMENSION')
    category = ig.mb.tag_get_handle('CATEGORY')
    global_id = ig.mb.tag_get_handle('GLOBAL_ID')
    # volume tags
    vol_dim = ig.mb.tag_get_data(geom_dim, all_vols)
    vol_category = ig.mb.tag_get_data(category, all_vols)
    vol_id = ig.mb.tag_get_data(global_id, all_vols)
    if list(vol_dim) == [3, 3]:
        r[0] = True
    if list(vol_category) == ['Volume', 'Volume']:
        r[1] = True
    if sorted(list(vol_id)) == [1, 2]:
        r[2] = True
    # surface tags
    surf_dim = ig.mb.tag_get_data(geom_dim, all_surfs)
    surf_category = ig.mb.tag_get_data(category, all_surfs)
    surf_id = ig.mb.tag_get_data(global_id, all_surfs)
    if list(surf_dim) == [2, 2, 2]:
        r[3] = True
    if list(surf_category) == ['Surface', 'Surface', 'Surface']:
        r[4] = True
    if sorted(list(surf_id)) == [1, 2, 3]:
        r[5] = True
    # Check Parent/Child Relationships
    #   1) each vol is parent to two surfs
    #   2) one surf (shared surface) has 2 vol parents
    #   3) two surfs each have one vol parent
    # 1) volumes have two child surfaces
    tmp1 = np.full(len(all_vols), False)
    tmp2 = np.full(len(all_vols), False)
    for i, vol in enumerate(all_vols):
        child_surfs = list(ig.mb.get_child_meshsets(vol))
        if all(child in all_surfs for child in child_surfs):
            tmp1[i] = True
        if len(child_surfs) == 2:
            tmp2[i] = True
    if all(tmp1) and all(tmp2):
        r[6] = True
    # 2) common surf has 2 vol parents
    parent_vols = list(ig.mb.get_parent_meshsets(common_surf))
    if sorted(parent_vols) == sorted(all_vols):
        r[7] = True
    # 3) other two surfaces have one vol parent
    all_surfs.remove(common_surf)
    tmp1 = np.full(len(all_surfs), False)
    tmp2 = np.full(len(all_surfs), False)
    for i, surf in enumerate(all_surfs):
        parent_vols = list(ig.mb.get_parent_meshsets(surf))
        if all(parent in all_vols for parent in parent_vols):
            tmp1[i] = True
            all_vols.remove(parent_vols)  # no shared parent vols
        if len(parent_vols) == 1:
            tmp2[i] = True
    if all(tmp1) and all(tmp2):
        r[8] = True
    assert(all(r))


def test_tag_for_viz():
    """test visualization tags are added to triangles"""
    # load volume
    ig = isg.IsGm()
    fs = ig.mb.create_meshset()
    ig.mb.load_file(test_dir + '/vol-files/single-box-1.stl', file_set=fs)
    iv = (0, fs)
    ig.isovol_meshsets[iv] = {}
    ig.separate_isovols()  # use this to get surface EHs
    # set val_tag (used for viz)
    val = 5.0
    surf = ig.isovol_meshsets[iv]['surfs_EH'][0]
    ig.mb.tag_set_data(ig.val_tag, surf, val)
    # tag for viz
    ig.tag_for_viz()
    # test viz tags
    all_tris = ig.mb.get_entities_by_type(surf, types.MBTRI)
    vals_out = list(ig.mb.tag_get_data(ig.val_tag, all_tris))
    vals_exp = list(np.full(12, 5.0))
    assert(vals_exp == vals_out)


@pytest.mark.parametrize("tagname, tagval, expname, expval, exptype",
                         [('int_tag', 1, 'int_tag', [1], np.int32),
                          ('float_tag', 2.0, 'float_tag', [2.0], np.float64),
                          ('list_tag', [3., 4.], 'list_tag',
                           [3., 4.], np.float64),
                          ('str_tag', 'val', 'str_tag', ['val'], np.string_),
                          (1.0, 'convert', '1.0', ['convert'], np.string_)])
def test_set_tags(tagname, tagval, expname, expval, exptype):
    """test tags are set on root set for various lengths and data types."""
    ig = isg.IsGm()
    tags = {tagname: tagval}
    # set tags
    ig.set_tags(tags)
    # check tag names, data, and type
    rs = ig.mb.get_root_set()
    r = np.full(3, False)
    try:
        # try to get tag by expected name
        tag_out = ig.mb.tag_get_handle(expname)
        r[0] = True
    except:
        # fails if tag does not exist
        r[0] = False
    # only test the rest if tag exists
    if r[0]:
        data_out = ig.mb.tag_get_data(tag_out, rs)
        # check data value
        if list(data_out[0]) == expval:
            r[1] = True
        # check data type
        if type(data_out[0][0]) is exptype:
            r[2] = True
    assert(all(r))


def test_write_geometry():
    """test that a file of the correct name is created"""
    # write a file
    ig = isg.IsGm()
    sname = 'write-test.h5m'
    ig.write_geometry(sname, test_dir)
    # check that file exists
    exp_file = test_dir + '/' + sname
    r = False
    if isfile(exp_file):
        r = True
        remove(exp_file)
    assert(r)


def test_write_geometry_ext():
    """generated file should have a different name than supplied + warnings"""
    r = np.full(4, False)
    # write file with incorrect extension
    ig = isg.IsGm()
    sname = 'write-test.bad'
    # check for a warning
    with warnings.catch_warnings(record=True) as w:
        ig.write_geometry(sname, test_dir)
        warnings.simplefilter("always")
        if len(w) == 1:
            r[0] = True
        if "File will be saved as type .h5m" in str(w[-1].message):
            r[1] = True
    # check that file exists
    good_file = test_dir + '/write-test.h5m'
    bad_file = test_dir + '/' + sname
    if isfile(good_file):
        # check that name was changed
        r[2] = True
        remove(good_file)
    if not isfile(bad_file):
        # check that bad name was not used
        r[3] = True
    else:
        # file exists, needs removed
        remove(bad_file)
    assert(all(r))


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


def test_compare_surfs():
    """test that new surf is correctly generated when comparing two"""
    # get setup
    ig = __setup_geom()
    ivs = sorted(ig.isovol_meshsets.keys())  # isovol info: (vol id, EH)
    iv1 = ivs[0]
    iv2 = ivs[1]
    fs1 = list(iv1)[1]
    fs2 = list(iv2)[1]
    # compare surfs
    norm = 1.5
    ig._IsGm__compare_surfs(iv1, iv2, norm)
    # checks
    r = np.full(5, False)  # truth array for checks
    # check number of surfaces in each volume (should be two each)
    surfs_1 = ig.isovol_meshsets[iv1]['surfs_EH']
    surfs_2 = ig.isovol_meshsets[iv2]['surfs_EH']
    if len(surfs_1) == 2:
        r[0] = True
    if len(surfs_2) == 2:
        r[1] = True
    # check there is a common surface
    common_surf = set(surfs_1) & set(surfs_2)
    if len(common_surf) == 1:
        r[2] = True
    # check the tags on the merged surface
    # sense tag:
    sense_out = list(ig.mb.tag_get_data(ig.sense_tag, common_surf)[0])
    sense_exp = [fs1, fs2]
    if sense_out == sense_exp:
        r[3] = True
    # value tag:
    # shared surface should be 7.5 (shared bound = 5, norm=1.5 -> 5*1.5=7.5)
    val_out = ig.mb.tag_get_data(ig.val_tag, common_surf)[0][0]
    val_exp = 7.5
    if val_out == val_exp:
        r[4] = True
    assert(all(r))


def test_compare_surfs_no_val():
    """no matching value - should throw warning"""
    # get setup
    ig = __setup_geom()
    iv = sorted(ig.isovol_meshsets.keys())  # isovol info: (vol id, EH)
    fs1 = list(iv[0])[1]
    fs2 = list(iv[1])[1]
    # change level info on one so there is no common value
    ig.isovol_meshsets[iv[1]]['bounds'] = (6., 10.)
    # compare surfs
    norm = 1.5
    r = np.full(3, False)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ig._IsGm__compare_surfs(iv[0], iv[1], norm)
        # check warning
        if len(w) == 1:
            r[0] = True
        if "No matching value" in str(w[-1].message):
            r[1] = True
    # common surface should be assigned a value of 0.0
    surfs_1 = ig.isovol_meshsets[iv[0]]['surfs_EH']
    surfs_2 = ig.isovol_meshsets[iv[1]]['surfs_EH']
    common_surf = set(surfs_1) & set(surfs_2)
    val_out = ig.mb.tag_get_data(ig.val_tag, common_surf)[0][0]
    val_exp = 0.0
    if val_out == val_exp:
        r[2] = True
    assert(all(r))


def test_calc_centroid():
    """test centroid calculation"""
    ig = isg.IsGm()
    # triangle defined by coordinates (3, 0, 0), (0, 3, 0), (0, 0, 3)
    coords = [3., 0., 0., 0., 3., 0., 0., 0., 3.]
    # expected centroid:
    # x = (3 + 0 + 0)/3 = 1
    # y = (0 + 3 + 0)/3 = 1
    # z = (0 + 0 + 3)/3 = 1
    exp_centroid = [1., 1., 1.]
    centroid = ig._IsGm__calc_centroid(coords)
    assert(exp_centroid == list(centroid))


def test_calc_centroid_error():
    """test that error is raised if wrong number of values in centroid calc"""
    ig = isg.IsGm()
    # incorrect list of only 2 coords (6 values)
    coords = [3., 0., 0., 0., 3., 0.]
    with pytest.raises(RuntimeError) as error_info:
        ig._IsGm__calc_centroid(coords)
    assert "Cannot calculate centroid" in str(error_info)


@pytest.mark.parametrize("point", [([-5, 0, 0]), ([5, 0, 0]),
                                   ([0, -5, 0]), ([0, 5, 0]),
                                   ([0, 0, -5]), ([0, 0, 5])])
def test_check_exterior_true(point):
    """test that point is correctly identified as being on the surface"""
    ig = isg.IsGm()
    # define the geometric extents to be -5 to 5 in all directions
    ig.xmin = ig.ymin = ig.zmin = -5.
    ig.xmax = ig.ymax = ig.zmax = 5.
    # test points on each surface (every point should be on a surface)
    result = ig._IsGm__check_exterior(point)
    assert(result)


def test_check_exterior_false():
    """test that point is correctly identified as being on the interior"""
    ig = isg.IsGm()
    # define the geometric extents to be -5 to 5 in all directions
    ig.xmin = ig.ymin = ig.zmin = -5.
    ig.xmax = ig.ymax = ig.zmax = 5.
    point = [0, 0, 0]  # interior
    result = ig._IsGm__check_exterior(point)
    assert(not result)
