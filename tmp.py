from threed_strudel.chop.chop_map import ChopMap

chop = ChopMap()

# m_path = '/Volumes/data/debug_strudel/debug_chop/emd_4225.map'
# cube_path = '/Volumes/data/debug_strudel/debug_chop/emd_4225_cube1.map'
# mod_path = '/Volumes/data/debug_strudel/debug_chop/6fbs.cif'
# chop.chop_cube(mod_path, m_path, out_map_path=cube_path, zero_origin=False)
#
# ala = '/Volumes/data/debug_strudel/debug_chop/ala-93-A.cif'
# ala_path = '/Volumes/data/debug_strudel/debug_chop/ala_soft1.mrc'
# chop.chop_soft_radius_watershed(ala, m_path, out_map=ala_path)
# import json
# import os
# def read_json(json_path):
#     with open(json_path) as j:
#         data = json.load(j)
#     return data
#
# data = read_json('/Volumes/data/debug_strudel/debug_chop/valid/segments_list.json')
#
# for pair in data:
#     res_name = os.path.basename(pair[1]).split('.')[0]
#     name_lst = res_name.split('-')
#     res_type = name_lst[0]
#     res_nr = name_lst[1]
#     chain = name_lst[2]
#     try:
#         int(res_nr)
#     except ValueError:
#         print(pair)
ll = 'asn--1-A_soft.mrc'.split('-')
ll1 = []

if ll[1] == '':
    ll[2] = '-' + ll[2]
    del ll[1]
print(ll)