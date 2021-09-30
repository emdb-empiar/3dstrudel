from threed_strudel.classify.rotamer_utils import read_rotamer_names


inp_bounds = '/Applications/phenix-1.13-2998/modules/cctbx_project/mmtbx/rotamer/rotamer_names.props'
out_j = '/Users/andrei/PycharmProjects/strudel/threed_strudel/lib/strudel/rotamer_bounds_lib.json'

read_rotamer_names(inp_bounds, out_j)