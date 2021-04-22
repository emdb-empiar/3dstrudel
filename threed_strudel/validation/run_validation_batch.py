"""
mapMotifValidation.py

This module runs map versus amino acid residues map motif library validation.

Copyright [2013] EMBL - European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the
"License"); you may not use this file except in
compliance with the License. You may obtain a copy of
the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on
an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
"""

__author__ = 'Andrei Istrate'
__email__ = 'andrei@ebi.ac.uk'
__date__ = '2018-05-29'

import os
import sys
import argparse
import subprocess
from shutil import copy
from distutils.spawn import find_executable
from threed_strudel import configure as config
from threed_strudel.utils import functions as func



def get_motif_resolution_ranges(library_path):
    lib_list = []
    folders = os.listdir(library_path)
    for folder in folders:
        if folder.startswith('motifs_'):
            res_range = folder.split('_')[-1].split('-')
            res_range = [float(i) for i in res_range]
            lib_list.append([folder, res_range])
    return lib_list


def main():
    parser = argparse.ArgumentParser(description='Map Validation')
    parser.add_argument("-i", "--json_input", dest="json_inp", required=True,
                        help="List of models maps resolution and out folder in json format. "
                             "e.g.: [[run_flag(0 or 1), map, model, resolution, out_folder]]")
    parser.add_argument("-o", "--out", dest="out", required=True, help="Output directory")
    parser.add_argument("-l", "--lib", dest="lib", required=True, help="Strudel motif libraries path")
    parser.add_argument("-M", "--memory", dest="mem", required=True, help="Maximum memory")
    parser.add_argument("-np", "--num_proc", dest="np", required=True, help="Number of processors")
    # parser.add_argument("-process_log", "--process_log", dest="process_log", default=None, required=False, help="Log file")
    # parser.add_argument("-v", "--voxel", dest="voxel", required=False, default=0.25, type=float, help="Segments voxel size")
    # parser.add_argument("-r", "--recompute", dest="recompute_scores", action='store_true',
    #                     help="Recalculate correlations")
    # parser.add_argument("-rs", "--recompute_segments", dest="recompute_segments", action='store_true',
    #                     help="Repeat segmenting of the input model and map if found in the output")
    # parser.add_argument("-vc", "--verbose_chimerax", dest="v_c", action='store_true',
    #                     help="Run ChimeraX in verbose mode")
    # parser.add_argument("-wl", "--warning_level", dest="warning_level", default='info',
    #                     help="Log file warning level [info[debug]")

    if not find_executable(config.CHIMERA_PATH):
        print('ChimeraX path not set!\nPlease run strudel_setChimeraX.py to set the ChimeraX path')
        sys.exit()


    print(' '.join(sys.argv))
    args = parser.parse_args()
    # if args.log:
    #     log_file = args.log
    # else:
    #     log_file = os.path.join(args.out, 'map_motif_validation.process_log')

    lib_list = get_motif_resolution_ranges(args.lib)

    entry_list = func.read_json(args.json_inp)
    if not os.path.exists(args.out):
        os.makedirs(args.out)
    tmp = os.path.join(args.out, 'tmp')
    if not os.path.exists(tmp):
        os.makedirs(tmp)
    for entry in entry_list:
        if entry[0]:
            resolution = entry[3]
            for lib in lib_list:
                if lib[1][0] < resolution <= lib[1][1]:
                    lib_path = os.path.join(args.lib, lib[0])
                    out_path = os.path.join(args.out, entry[4])
                    if not os.path.exists(out_path):
                        os.makedirs(out_path)
                    command = f'bsub -o o_{lib[0]} -e e_{lib[0]} -M {args.mem} -n {args.np} strudel_mapMotifValidation.py ' \
                              f'-p {entry[2]} -m {entry[1]} -l {lib_path} -o out -log vs_{lib[0]}.log -np {args.np}'
                    print(f'running {command}')

                    # subprocess.call(command, cwd=out_path, shell=True)


if __name__ == '__main__':
    main()