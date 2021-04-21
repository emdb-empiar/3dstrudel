import pytest
import numpy as np
from threed_strudel.utils import bio_utils
from threed_strudel.parse.map_parser import MapParser
from threed_strudel.chop.chop_map import ChopMap


@pytest.fixture
def map_obj():
    map_path = "data/in/6272.mrc"
    map_obj_ = MapParser(map_path)
    return map_obj_


@pytest.fixture
def model_obj():
    model_path = "data/in/6272.pdb"
    model_obj_ = bio_utils.load_structure(model_path)
    return model_obj_


def test_chop_cube(map_obj, model_obj):
    model = model_obj[0]['A'][158]
    structure = bio_utils.residues2structure(model)
    bio_utils.save_model(structure, 'data/out/chop/res_158.cif')

    zero_out_map, shifts = ChopMap.chop_cube(model, map_obj, 4, True)
    zero_out_map.write_map('data/out/chop/res_158_zero_orig.mrc')

    out_map, shifts = ChopMap.chop_cube(model, map_obj, 4, False)
    out_map.write_map('data/out/chop/res_158.mrc')

    new = MapParser('data/out/chop/res_158.mrc')
    new_zero = MapParser('data/out/chop/res_158_zero_orig.mrc')
    ref_zero = MapParser('data/in/chop/ref_res_158_zero_orig.mrc')
    ref = MapParser('data/in/chop/ref_res_158.mrc')

    assert (np.array_equal(new_zero.data, ref_zero.data))
    assert (np.array_equal(new_zero.origin, ref_zero.origin))
    assert (np.array_equal(new_zero.n_start, ref_zero.n_start))

    assert (np.array_equal(new.data, ref.data))
    assert (np.array_equal(new.origin, ref.origin))
    assert (np.array_equal(new.n_start, ref.n_start))


def test_chop_map_list(map_obj, model_obj):
    models = [model_obj[0]['A'][160], model_obj[0]['A'][161]]
    chop = ChopMap()
    maps, shifts = chop.chop_cube_list(models, map_obj, 4)
    assert len(maps) == len(models) and len(shifts) == len(models)
    assert all([isinstance(el, MapParser) for el in maps])


def test_chop_soft_radius():
    model = bio_utils.load_structure('data/in/chop/ref_res_158.cif')
    map_obj = MapParser('data/in/chop/ref_res_158.mrc')
    out_map = ChopMap.chop_soft_radius(model, map_obj, hard_radius=2, soft_radius=1)
    out_map.write_map('data/out/chop/soft_res_158.mrc')

    ref = MapParser('data/in/chop/soft_res_158.mrc')
    assert (np.array_equal(out_map.data, ref.data))


def test_find_near_atoms(model_obj):

    atom = model_obj[0]['A'][158]['CZ']
    print(atom)
    near_atoms = ChopMap.find_near_atoms(atom, model_obj, distance=4)
    result = [(a.parent.id[1], a.id) for a in near_atoms]
    result.sort()
    assert result == [(188, 'CD1'), (188, 'CG1'), (326, 'CD1')]


def test_chop_soft_radius_watershed():
    mod_p = '/Volumes/data/Work/covid/new_13_may/11007/bundle_30178_new_chop/input/7btf.cif'
    map_p = '/Volumes/data/Work/covid/new_13_may/11007/bundle_30178_new_chop/input/emd_30178.map'
    model_obj = bio_utils.load_structure(mod_p)
    map_obj = MapParser(map_p)
    local_model = model_obj[0]['A'][764]
    chop = ChopMap()

    cube_map_obj, shifts = chop.chop_cube(local_model, map_obj, 4, zero_origin=False)
    cube_map_obj.grid_resample_emda(0.25)
    side_chain = bio_utils.del_main_chain(local_model)
    import time

    t1 = time.time()
    # map_obj_w, mask, outer_mask = chop.chop_soft_radius_watershed(side_chain, cube_map_obj, model_obj, radius=2, soft_radius=1, asymmetric_delta=0.5)
    map_obj_w = chop.chop_soft_radius_watershed(side_chain, cube_map_obj, model_obj, radius=2,
                                                                  soft_radius=1, asymmetric_delta=0.5)
    print(f"New time: {time.time()-t1}")

    map_obj_w.write_map('/Volumes/data/Work/covid/new_13_may/11007/bundle_30178_new_chop/input/784_5_ccano_asim-0.5_tmp.mrc')
    # outer_mask.write_map('/Volumes/data/Work/covid/new_13_may/11007/bundle_30178_new_chop/input/784_5_ccano_asim-0.5_out_mask_nozero5.mrc')
    # mask.write_map(
    #     '/Volumes/data/Work/covid/new_13_may/11007/bundle_30178_new_chop/input/784_3_ccano_asim-1_fin_mask.mrc')

    # t1 = time.time()
    # map_obj_w_old = chop.chop_soft_radius_watershed_old(side_chain, cube_map_obj, model_obj, radius=2,
    #                                                             soft_radius=1)
    # print(f"Old time: {time.time()-t1}")
    # t1 = time.time()
    # map_obj_w_old.write_map('/Volumes/data/Work/covid/new_13_may/11007/bundle_30178_new_chop/input/784_2_old.mrc')
    #
    # map_obj_w_old = chop.chop_soft_radius(side_chain, cube_map_obj, hard_radius=2,
    #                                                     soft_radius=1)
    # print(f"no env time: {time.time() - t1}")
    # map_obj_w_old.write_map('/Volumes/data/Work/covid/new_13_may/11007/bundle_30178_new_chop/input/784_2_soft.mrc')
#
#
# def test_chop_soft_radius_watershed_var2(model_obj):
#
#     map_obj = MapParser('data/in/chop/ref_res_158.mrc')
#     local_model = model_obj[0]['A'][158]
#     chop = ChopMap()
#     map_obj, mask, outer_mask = chop.chop_soft_radius_watershed_var2(local_model, map_obj, model_obj, radius=2, soft_radius=1)
#     map_obj.write_map('data/out/chop/watershed/res9_158.mrc')
#     mask.write_map('data/out/chop/watershed/mask9_158.mrc')
#     outer_mask.write_map('data/out/chop/watershed/outer_mask9_158.mrc')