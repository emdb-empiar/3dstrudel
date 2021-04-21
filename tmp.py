from threed_strudel.chop.chop_map import ChopMap
model = '/Users/andrei/PycharmProjects/strudel/threed_strudel/unitest/data/in/validate/6272.pdb'
map_ = '/Users/andrei/PycharmProjects/strudel/threed_strudel/unitest/data/in/validate/6272.mrc'

out_map = '/Users/andrei/PycharmProjects/strudel/threed_strudel/unitest/data/in/validate/6272_cut.mrc'

ChopMap.chop_cube(model, map_, zero_origin=False, out_map_path=out_map)