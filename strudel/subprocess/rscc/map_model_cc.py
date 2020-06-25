from __future__ import division
import sys
import iotbx.pdb
from libtbx import group_args, easy_pickle
import mmtbx.utils
import mmtbx.maps.map_model_cc
import json
import os


# Usage
# python map_model_cc.py map.mrc residue.pdb resolution=*.* atom_radius=3.0 scattering_table=electron out_file=file.json


master_params_str = """\
    map_file_name = None
        .type = str
        .help = Map file name
    model_file_name = None
        .type = str
        .help = Model file name
    pkl_file_name = None
        .type = in_dir
        .help = File name for pickle file with results
        .expert_level = 3
    include scope mmtbx.maps.map_model_cc.master_params
    """


def master_params():
    return iotbx.phil.parse(master_params_str, process_includes=True)


def get_inputs(args, master_params):

    inputs = mmtbx.utils.process_command_line_args(args=args,
                                                   master_params=master_params)
    # Model
    pdb_file_name = inputs.pdb_file_names[0]
    pdb_hierarchy = iotbx.pdb.input(file_name=pdb_file_name).construct_hierarchy()
    # Map
    ccp4_map_object = inputs.ccp4_map
    crystal_symmetry = inputs.crystal_symmetry
    #
    return group_args(params=inputs.params.extract(),
                      pdb_file_name=pdb_file_name,
                      pdb_hierarchy=pdb_hierarchy,
                      ccp4_map_object=ccp4_map_object,
                      crystal_symmetry=crystal_symmetry)


def run(args):
    # Get inputs
    out_basename = 'residue_cc'
    for index, arg in enumerate(args):
        if arg.split('=')[0] == 'out_file':
            out_basename = arg.split('=')[1]
            del args[index]

    inputs = get_inputs(args=args,
                        master_params=master_params())

    # Run task in 4 separate steps
    task_obj = mmtbx.maps.map_model_cc.map_model_cc(map_data=inputs.ccp4_map_object.data.as_double(),
                                                    pdb_hierarchy=inputs.pdb_hierarchy,
                                                    crystal_symmetry=inputs.crystal_symmetry,
                                                    params=inputs.params.map_model_cc)
    task_obj.validate()
    task_obj.run()
    results = task_obj.get_results()

    head, tail = os.path.split(out_basename)
    out_basename = os.path.join(head, tail.split('.')[0])
    out_json = out_basename + '.json'
    out_txt = out_basename + '.txt'

    with open(out_txt, 'w') as txt_file:
        txt_file.write("chain ID res    Nr    CC      <B>   <occ>\n")
        data = []
        for r in results.cc_per_residue:
            fmt = '{}        {}  {}  {:2.5f} {:8.3f} {:4.2f}\n'
            txt_file.write(fmt.format(r.chain_id, r.resname, r.resseq, r.cc, r.b_iso_mean, r.occ_mean))
            data.append([r.chain_id, r.resname, r.resseq, r.cc])

    with open(out_json, 'w') as j:
        json.dump(data, j, indent=4)
    return None


if __name__ == "__main__":
    run(args=sys.argv[1:])
