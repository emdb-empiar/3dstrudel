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
import csv
import json
import subprocess
import numpy as np
import psutil
import operator
import shutil
import logging
from datetime import datetime
import time
import argparse
from multiprocessing import Process, Manager
from distutils.spawn import find_executable

from threed_strudel.chop.chop_map import ChopMap
from threed_strudel.utils import functions as func
import threed_strudel.configure as config
from threed_strudel.utils import bio_utils
from threed_strudel import nomenclature

log = logging.getLogger(__name__)


class DictKeys:
    """
    Class for storing the fields names in the output csv file
    """
    ID = 'id'
    OUTLIER = 'outlier'
    DELTA = 'relative_difference'
    RES_TYPE = 'residue_type'
    RES_NR = 'residue_nr'
    CHAIN = 'chain'
    RES_MAP = 'residue_map'
    RES_MODEL = 'residue_model'
    M_TOP_TYPE = 'top_motif_type'
    M_TOP_NAME = 'top_motif_name'
    M_TOP_MATRIX = 'top_motif_matrix'
    TOP_CC = 'top_rscc'
    SAME_TYPE_NAME = 'same_type_motif_name'
    SAME_TYPE_CC = 'same_type_motif_rscc'
    SAME_TYPE_MATRIX = 'same_type_motif_matrix'
    # FULL NAME EXAMPLE "ala_cc" suffix for top correlation of specific residue type library motifs
    TYPE_TOP_CC = '_top_motif_rscc'
    # FULL NAME EXAMPLE "ala_motif_name"
    TYPE_TOP_NAME = '_top_motif_name'
    TYPE_TOP_MATRIX = '_top_motif_matrix'
    SCORES = 'scores'
    ALL_TOP_MOTIFS = 'all_top_motifs'
    TOP_SAME_TYPE_DIFF = 'top_same_dif'
    COMPLETE_DATA = 'complete_data'



def dict_list_to_csv(dict_list, csv_path):
    """
    Dumps a list of dictionaries with identical keys to a csv file
    :param dict_list: list of dictionaries
    :param csv_path: path
    """
    with open(csv_path, mode='w') as csv_file:
        fieldnames = [key for key in dict_list[0].keys()]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for item in dict_list:
            writer.writerow(item)


def json_to_full_csv(json_path, csv_path):
    """
    Converts validation json output file to csv format
    :param json_path: input json file path
    :param csv_path: output csv path
    """
    k = DictKeys()
    with open(json_path) as j:
        data = json.load(j)

    out_data = []
    for item in data:
        tmp = {}
        for key, value in item.items():
            if key != k.SCORES:
                tmp[key] = value
            else:
                for el in value:
                    tmp[el[0]] = el[1]
        out_data.append(tmp)

    dict_list_to_csv(out_data, csv_path)


def csv_to_top_csv(csv_path, out_csv_path, outlier_diff=0.05):
    """
    Prepares the input file for Strudel Score tool.
    Keeps only the top scoring motif data for each residue type
    :param csv_path: input scv file path
    :param out_csv_path: output scv file path
    :param outlier_diff: difference for outliers, diff = (top_score-same_type_score) / top_score
    """
    data_frame = []
    k = DictKeys()
    dd = DictKeys.__dict__
    k_values = [dd[k] for k in dd.keys() if not k.startswith('__')]

    with open(csv_path, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)

        for row in csv_reader:
            codes = nomenclature.AA_RESIDUES_LIST
            row_d = {}
            tmp = {}
            if any(['None' in v for v in row.values()]):
                row_d[k.COMPLETE_DATA] = 0
            else:
                row_d[k.COMPLETE_DATA] = 1
            for key, value in row.items():
                code = key[:3]
                if code.upper() in codes:
                    corr, _, matrix = value.partition('_m_')
                    if corr == 'None':
                        corr = None
                    if matrix == 'None':
                        matrix = None
                    if code not in tmp.keys():
                        if 'None' not in value:
                            tmp[code] = [float(corr), key, matrix]
                        else:
                            tmp[code] = [None, None, None]
                    elif tmp[code][0] is not None and 'None' not in value:
                        if tmp[code][0] < float(corr):
                            tmp[code][0] = float(corr)
                            tmp[code][1] = key
                            tmp[code][2] = matrix

                elif key in k_values:
                    try:
                        row_d[key] = int(value)
                    except ValueError:
                        row_d[key] = value
                else:
                    raise Exception(f'Key {key} is unknown to the data structure cannot continue')

            max_correlation = -1
            row_d[k.M_TOP_TYPE] = None
            row_d[k.M_TOP_NAME] = None
            row_d[k.TOP_CC] = None
            row_d[k.SAME_TYPE_NAME] = None
            row_d[k.SAME_TYPE_CC] = None
            row_d[k.SAME_TYPE_MATRIX] = None
            for key, value in tmp.items():
                if key == row_d[k.RES_TYPE]:
                    row_d[k.SAME_TYPE_NAME] = value[1]
                    row_d[k.SAME_TYPE_CC] = value[0]
                    row_d[k.SAME_TYPE_MATRIX] = value[2]
                if value[0] is not None:
                    if value[0] > max_correlation:
                        max_correlation = value[0]
                        row_d[k.M_TOP_TYPE] = key
                        row_d[k.M_TOP_NAME] = value[1]
                        row_d[k.TOP_CC] = value[0]
                        row_d[k.M_TOP_MATRIX] = value[2]

            diff = (row_d[k.TOP_CC] - row_d[k.SAME_TYPE_CC]) / row_d[k.TOP_CC]
            row_d[k.DELTA] = round(diff, 4)
            if diff > outlier_diff:
                row_d[k.OUTLIER] = 1
            else:
                row_d[k.OUTLIER] = 0

            for key, value in tmp.items():
                row_d[key + k.TYPE_TOP_CC] = value[0]
                row_d[key + k.TYPE_TOP_NAME] = value[1]
                row_d[key + k.TYPE_TOP_MATRIX] = value[2]
            data_frame.append(row_d)

    data_frame.sort(key=operator.itemgetter(k.CHAIN, k.RES_NR))
    dict_list_to_csv(data_frame, out_csv_path)


def csv_to_top_scores_only_csv(csv_path, out_csv_path, outlier_diff=0.05):
    """
    Keeps only minimum data for the web version of the Strudel Score.
    :param csv_path: input scv file path
    :param out_csv_path: output scv file path
    :param outlier_diff: difference for outliers, diff = (top_score-same_type_score) / top_score
    """
    data_frame = []
    k = DictKeys()
    dd = DictKeys.__dict__
    k_values = [dd[k] for k in dd.keys() if not k.startswith('__')]

    with open(csv_path, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)

        for i, row in enumerate(csv_reader):
            codes = nomenclature.AA_RESIDUES_LIST
            row_d = {}
            tmp = {}
            row_d[k.ID] = i
            if any(['None' in v for v in row.values()]):
                row_d[k.COMPLETE_DATA] = 0
            else:
                row_d[k.COMPLETE_DATA] = 1
            for key, value in row.items():
                code = key[:3]
                if code.upper() in codes:
                    corr, _, matrix = value.partition('_m_')
                    if code not in tmp.keys():
                        if 'None' not in value:
                            tmp[code] = [float(corr), key]
                        else:
                            tmp[code] = [None, None]

                    elif tmp[code][0] is not None and 'None' not in value:
                        if tmp[code][0] < float(corr):
                            tmp[code][0] = float(corr)
                            tmp[code][1] = key

                elif key in k_values:
                    try:
                        row_d[key] = int(value)
                    except ValueError:
                        row_d[key] = value
                else:
                    raise Exception(f'Key {key} is unknown to the data structure cannot continue')

            max_correlation = -1
            row_d[k.M_TOP_TYPE] = None
            row_d[k.TOP_CC] = None
            for key, value in tmp.items():
                if key == row_d[k.RES_TYPE]:
                    row_d[k.SAME_TYPE_CC] = value[0]
                if value[0] is not None:
                    if value[0] > max_correlation:
                        max_correlation = value[0]
                        row_d[k.M_TOP_TYPE] = key
                        row_d[k.TOP_CC] = value[0]

            diff = (row_d[k.TOP_CC] - row_d[k.SAME_TYPE_CC]) / row_d[k.TOP_CC]
            row_d[k.DELTA] = round(diff, 4)
            if diff > outlier_diff:
                row_d[k.OUTLIER] = 1
            else:
                row_d[k.OUTLIER] = 0
            for key, value in tmp.items():
                row_d[key] = value[0]
            data_frame.append(row_d)
    data_frame.sort(key=operator.itemgetter(k.CHAIN, k.RES_NR))
    dict_list_to_csv(data_frame, out_csv_path)


class ComputeScores:
    """
    Class for correlations calculation between each amino acid residue map and each map motif
    """
    def __init__(self):
        self.residues = nomenclature.AA_RESIDUES_LIST
        self.chimera_path = self.find_chimerax()
        basedir = os.path.dirname(os.path.abspath(__file__))
        self.score_chimera_script = os.path.join(basedir, 'score_chimerax.py')
        self.segmented_pairs = 'segments_list.json'
        self.scores_file = 'scores.json'
        self.score_log = 'scores.log'
        self.work_grid_sampling = 0.25
        self.chopped_pairs = Manager().list()
        self.chopped_res_names = Manager().list()
        self.chop = ChopMap()
        self.segments = None
        self.work_dir = None
        self.lib = None
        self.in_map = None
        self.in_model = None
        self.nr_residues = None
        self.input = None
        self.out_dir = None

    def find_chimerax(self):
        path = config.CHIMERA_PATH
        out = find_executable(path)
        if out is None:
            raise Exception(f'Could not find ChimeraX\n please edit {os.path.abspath(config.__file__)}')
        else:
            return path

    def set_paths(self, work_dir, lib, in_map, in_model):
        """
        Creates the directory tree for the output and sets the input files paths
        :param work_dir: output directory
        :param lib: map motif library path
        :param in_map: input map path
        :param in_model: input model path
        """
        self.work_dir = os.path.abspath(work_dir)
        if os.path.exists(lib):
            self.lib = os.path.abspath(lib).rstrip('/')
            self.out_dir = os.path.join(self.work_dir, 'vs_' + os.path.basename(self.lib))
            self.lib = os.path.join(self.lib, 'motifs')
            self.check_motif_lib(self.lib)
        else:
            log.error('Motif library %s not found', lib)
            sys.exit()
        self.in_map = os.path.abspath(in_map)
        self.in_model = os.path.abspath(in_model)

        self.segments = os.path.join(self.work_dir, 'segments')
        self.input = os.path.join(self.work_dir, 'input')

        for path in [self.segments, self.input, self.out_dir]:
            try:
                os.makedirs(path)
            except FileExistsError:
                pass

        shutil.copy(self.in_map, self.input)
        shutil.copy(self.in_model, self.input)

    @staticmethod
    def check_motif_lib(motif_lib_dir):
        """
        Checks the integrity of the motif library
        :param motif_lib_dir: motif library directory
        """
        log.info(f'Checking @{os.path.basename(motif_lib_dir)} motif library')

        files = os.listdir(motif_lib_dir)
        if len(files) == 0:
            log.error('Empty motif library')
        pdb_names = [i.split('.')[0] for i in files if i.endswith('.cif') and not i.startswith('.')]
        map_names = [i.split('.')[0] for i in files if (i.endswith('.mrc') or
                                                        i.endswith('.map')) and not i.startswith('.')]
        pdb_names.sort()
        map_names.sort()
        if len(pdb_names) != len(map_names):
            log.warning('There are missing files in the map motif library: %s', motif_lib_dir)
        else:
            log.info(f'OK')
        for name in pdb_names:
            if name not in map_names:
                log.warning('Missing map file for %s motif', name)
        for name in map_names:
            if name not in pdb_names:
                log.warning('Missing model file for %s motif', name)

    def chop_structure(self, n_cores=4, voxel=0.25, replace=True):
        """
        Chops the input map and model into amino-acid fragments
        :param n_cores: number of computing cores
        :param voxel: library voxel size
        :param replace: replace or not fragments chopped in during previous runs
        :return: segments list file path (json format)
        """
        try:
            in_model = bio_utils.load_structure(self.in_model)[0]
            log.info('Failed to load structure using "auth" cif records. Attempting to use "label" records')
        except:
            in_model = bio_utils.load_structure_label_id(self.in_model)[0]

        del_atom_list = []
        for chain in in_model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name().upper().startswith('H'):
                        del_atom_list.append((chain.get_id(), residue.get_id(), atom.get_id()))

        for chain, res, atom in del_atom_list:
            del in_model[chain][res][atom]

        residues = [r for r in in_model.get_residues() if r.get_resname().upper() in self.residues]
        self.nr_residues = len(residues)
        # residues = residues[:]
        log.info('Found %s residues in the input model', len(residues))
        log.info('Running map chopping using %s cores', n_cores)
        map_obj_list, shifts_list = self.chop.chop_cube_list(residues, self.in_map, 4, zero_origin=False)
        split_residues = func.split_list(residues, n_cores)
        split_map_obj_list = func.split_list(map_obj_list, n_cores)
        split_shifts_list = func.split_list(shifts_list, n_cores)
        thread_list = []
        for i in range(len(split_residues)):
            t = Process(target=self._slave,
                        args=(split_residues[i], split_map_obj_list[i],
                              split_shifts_list[i], in_model, self.segments, voxel, replace))
            thread_list.append(t)
        for thread in thread_list:
            thread.start()
        for thread in thread_list:
            thread.join()

        out_pairs = []
        for i in self.chopped_pairs:
            out_pairs.append(i)

        segments_list_path = os.path.join(self.segments, self.segmented_pairs)

        if not os.path.exists(segments_list_path):
            func.save_json(segments_list_path, out_pairs)
        else:
            lst = func.read_json(segments_list_path)
            if len(out_pairs) > len(lst):
                p = psutil.Process()
                new_lst = os.path.join(self.segments, 'segments_' + str(p.pid) + '.json')
                func.save_json(new_lst, out_pairs)
                return new_lst

        return segments_list_path

    def _slave(self, residues, map_obj_list, shifts_list, whole_model, work_dir, voxel, replace):
        """
        Chopping worker
        :param residues: list of biopython residues objects
        :param map_obj_list: strudel map objects list
        :param shifts_list: list of translation matrices
        :param whole_model: input model as biopython object
        :param work_dir: work directory
        :param voxel: voxel size
        :param replace: replace or not fragments chopped in during previous runs
        """

        for i in range(len(residues)):
            residue = residues[i]
            matrix = shifts_list[i]

            res_type = residue.get_resname().lower()
            chain_id = residue.parent.id
            res_nr = residue.id[1]
            name_prefix = '{}-{}-{}'.format(res_type, res_nr, chain_id)

            # cube_map_path = os.path.join(work_dir, name_prefix + '_cube.mrc')
            res_cube_map_path = os.path.join(work_dir, name_prefix + '_cube_res.mrc')
            fin_map_path = os.path.join(work_dir, name_prefix + '_soft.mrc')
            residue_path = os.path.join(work_dir, name_prefix + '.cif')


            if not os.path.exists(fin_map_path) or not os.path.exists(residue_path) or replace:
                # map_obj_list[i].write_map(cube_map_path)
                cube_map_obj = map_obj_list[i]
                # cube_map_obj.write_map(cube_map_path)

                cube_map_obj.grid_resample_emda(voxel)

                # cube_map_obj.write_map(res_cube_map_path)

                bio_utils.shift_coord(matrix, residue)
                struct = bio_utils.residues2structure(residue)
                bio_utils.save_model(struct, residue_path)
                side_chain = bio_utils.del_main_chain(residue)
                # fin_map = self.chop.chop_soft_radius(side_chain, res_cube_map_path, hard_radius=2, soft_radius=1,)
                fin_map = self.chop.chop_soft_radius_watershed(side_chain, cube_map_obj, whole_model,
                                                               radius=2, soft_radius=1, )
                if np.isnan(np.sum(fin_map.data)):
                    fin_map.data[np.isnan(fin_map.data)] = 0
                    log.warning("NaN values in {}\n"
                                "NaN values were replaced by 0 in further calculations.".format(fin_map_path))
                fin_map.write_map(fin_map_path)

            self.chopped_pairs.append((fin_map_path, residue_path))
            self.chopped_res_names.append(name_prefix)

    def compute_correlations(self, pairs_json, n_cores=4, recompute=False, verbose=False):
        """
        Calculates input map fragments versus library map motifs correlations
        Runs Chimera as subprocess
        :param pairs_json: json file with chopped maps amd models paths
        :param n_cores: number of cores
        :param recompute: recompute or not known correlations from previous runs
        :param verbose: to run ChimeraX in verbose mode or not
        :return: output json file path
        """
        log.info('Running correlations calculations with ChimeraX')
        # Update keys in the chimeraX script
        updated_script = self.update_chimerax_script(self.score_chimera_script)
        if recompute:
            recompute = '-r'
        else:
            recompute = ''

        if verbose:
            v_flag = ''
        else:
            v_flag = '--silent'

        pairs_json = os.path.abspath(pairs_json)
        json_out = os.path.join(self.out_dir, self.scores_file)
        json_out = os.path.abspath(json_out)
        log_path = os.path.join(self.out_dir, self.score_log)
        if os.path.exists(log_path):
            shutil.copy(log_path, log_path + '_back')
            os.remove(log_path)

        command = f'{self.chimera_path} --nogui {v_flag} {updated_script}' \
                  f' -p {pairs_json} -sl {self.lib} -np {n_cores} -o {json_out} {recompute} ' \
                  f'-l {self.score_log} > chimera.log'

        subprocess.call(command, cwd=self.out_dir, shell=True)
        try:
            with open(log_path, 'r') as f:
                lines = f.readlines()
                try:
                    ln = len(lines[0].split()[0])
                    lines[0] = lines[0][ln+2:]
                except IndexError:
                    pass
                try:
                    lines[-1] = lines[-1].rstrip()
                except IndexError:
                    pass
                log.info(''.join(lines))
        except FileNotFoundError:
            pass
        return json_out

    @staticmethod
    def update_chimerax_script(script_path):
        """
        Updates the dictionary keys in the Chimerax script
        :param script_path: Chimerax script path
        :return:
        """
        dd = DictKeys.__dict__
        with open(script_path, 'r') as f:
            lines = f.readlines()
        for key, value in dd.items():
            if not key.startswith('__'):
                for i, line in enumerate(lines):
                    if line.startswith('    ' + key):
                        lines[i] = f"    {key} = '{value}'\n"
        head, tail = os.path.split(script_path)
        updated_script_path = os.path.join(head, 'updated_' + tail)
        with open(updated_script_path, 'w') as of:
            for line in lines:
                of.write(line)
        return updated_script_path

    @staticmethod
    def check_scores_completeness(segmented_residues, top_csv_file, scoring_log):
        k = DictKeys()
        scored = []
        complete = True
        with open(top_csv_file, mode='r') as csv_file:
            csv_reader = csv.DictReader(csv_file)

            for row in csv_reader:
                name = f'{row[k.RES_TYPE]}-{row[k.RES_NR]}-{row[k.CHAIN]}'
                scored.append(name)
                if not int(row[k.COMPLETE_DATA]):
                    complete = False
                    log.warning(f'The data for residue {name} is not complete.'
                                f' Check {scoring_log} for details')
            not_scored = []
            for name in segmented_residues:
                if name not in scored:
                    not_scored.append(name)
                    complete = False

            if len(not_scored) == 1:
                log.warning(f'Residue {not_scored[0]} was not scored. '
                            f'Please restart the program with the same input parameters.')
            elif len(not_scored) > 1:
                log.warning(f'Residues {", ".join(not_scored)} were not scored.\n'
                            f'Please restart the program with the same input parameters.')
        if complete:
            log.info("Validation completed successfully")
        return complete
                    

def main():
    parser = argparse.ArgumentParser(description='Map Validation')
    parser.add_argument("-p", "--model", dest="in_model", required=True, help="Atomic model path")
    parser.add_argument("-m", "--map", dest="in_map", required=True, help="Map file path")
    parser.add_argument("-l", "--lib", dest="lib", required=True, help="Strudel motif library path")
    parser.add_argument("-np", "--n_processors", dest="np", required=False, default=2, help="Number of processors")
    parser.add_argument("-o", "--out", dest="out", required=True, help="Output directory")
    parser.add_argument("-log", "--log", dest="log", default=None, required=False, help="Log file")
    parser.add_argument("-v", "--voxel", dest="voxel", required=False, default=0.25, type=float, help="Segments voxel size")
    parser.add_argument("-r", "--recompute", dest="recompute_scores", action='store_true',
                        help="Recalculate correlations")
    parser.add_argument("-rs", "--recompute_segments", dest="recompute_segments", action='store_true',
                        help="Repeat segmenting of the input model and map if found in the output")
    parser.add_argument("-vc", "--verbose_chimerax", dest="v_c", action='store_true',
                        help="Run ChimeraX in verbose mode")
    parser.add_argument("-wl", "--warning_level", dest="warning_level", default='info',
                        help="Log file warning level [info[debug]")
    parser.add_argument("-diff", "--outlier_difference", dest="diff", required=False, default=0.05, type=float,
                        help="Outliers threshold, (calculated as: (top_score-same_type_score) / top_score")
    parser.add_argument("-s", "--keep_segments", dest="segm", action='store_true',
                        help="Keep residues segments after finishing")

    args = parser.parse_args()

    if not os.path.exists(args.out):
        os.makedirs(args.out)

    logging.basicConfig(filename=args.log, level=logging.INFO, format='%(levelname)s:  %(message)s')
    date_time = datetime.now().strftime("%H:%M:%S %Y-%m-%d")
    logging.info('\n{:_^100}'.format('MapMotifValidation') + '\n\nRun command:\n' +
                 ' '.join(sys.argv) + '\n\nStarted: {}\n'.format(date_time))

    args.np = int(args.np)
    start = time.time()

    compute = ComputeScores()
    compute.set_paths(work_dir=args.out, lib=args.lib, in_map=args.in_map, in_model=args.in_model)
    chopped_pairs = compute.chop_structure(n_cores=args.np, voxel=args.voxel, replace=args.recompute_segments)

    # correlations_log = f'score_{datetime.now().strftime("%Y-%m-%d")}.log'
    json_file_path = compute.compute_correlations(chopped_pairs, args.np, args.recompute_scores, verbose=args.v_c)

    prefix = os.path.splitext(json_file_path)[0]
    if os.path.exists(prefix + '.csv'):
        csv_to_top_csv(prefix+'.csv', prefix+'_top.csv', args.diff)
        compute.check_scores_completeness(compute.chopped_res_names, prefix + '_top.csv', compute.score_log)
        csv_to_top_scores_only_csv(prefix + '.csv', prefix + '_top_for_web.csv')

    if not args.segm:
        logging.info('Cleaning temporary files...')
        shutil.rmtree(compute.segments)
        out_files = os.listdir(compute.out_dir)
        tmp_files = [f for f in out_files if f.endswith('_bak')]
        for tmp in tmp_files:
            os.remove(os.path.join(compute.out_dir, tmp))

    logging.info('Elapsed: {}\n{:_^100}'.format(func.report_elapsed(start), ''))


if __name__ == '__main__':
    main()