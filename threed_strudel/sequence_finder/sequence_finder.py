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
from Bio.PDB.Polypeptide import PPBuilder
# from strudel.chop.chopMap import ChopMap
from threed_strudel.utils import functions as func
import threed_strudel.configure as config
from threed_strudel.utils import bio_utils
from threed_strudel import nomenclature
from threed_strudel.chop.segment_residue_map import ExtractMap

log = logging.getLogger(__name__)


class DictKeys:
    """
    Class for storing the fields names in the output csv file
    """
    RES_TYPE = 'residue_type'
    RES_NR = 'residue_nr'
    CHAIN = 'chain'
    RES_MAP = 'residue_map'
    SCORES = 'scores'
    #scores dict
    TYPE = 'type'
    CORR = 'corr'
    MAP = 'map'
    MATRIX = 'matrix'


class ComputeScores:
    """
    Class for correlations calculation between each amino acid residue map and each map motif
    """
    def __init__(self):
        self.residues = nomenclature.AA_RESIDUES_LIST
        self.aa_3to1 = nomenclature.AA_3LET_TO_1LET
        self.aa_1to3 = {v: k for k, v in self.aa_3to1.items()}
        self.chimera_path = self.find_chimerax()
        basedir = os.path.dirname(os.path.abspath(__file__))
        self.score_chimera_script = os.path.join(basedir, 'score_chimerax.py')
        self.segmented_pairs = 'segments_list.json'
        self.scores_file = 'scores.json'
        self.work_grid_sampling = 0.25
        self.chopped_pairs = Manager().list()
        # self.chop = ChopMap()
        self.extract = None
        self.segments = None
        self.work_dir = None
        self.lib = None
        self.in_map = None
        self.in_model = None
        self.input = None
        self.out_dir = None
        self.all_libs = '/Volumes/data/libs/libs'
        self.scores_paths = {}
        self.sequence_dir = None
        self.sequence_dict = {}
        self.keys = DictKeys()

    def find_chimerax(self):
        path = config.CHIMERA_PATH
        out = find_executable(path)
        if out is None:
            raise Exception(f'Could not find ChimeraX\n please edit {os.path.abspath(config.__file__)}')
        else:
            return path

    def set_dir_tree(self, work_dir):
        self.work_dir = os.path.abspath(work_dir)
        self.segments = os.path.join(self.work_dir, 'segments')
        self.input = os.path.join(self.work_dir, 'input')
        self.sequence_dir = os.path.join(self.work_dir, 'sequence')

        for path in [self.segments, self.input, self.out_dir, self.sequence_dir]:
            try:
                os.makedirs(path)
            except FileExistsError:
                pass
            except TypeError:
                pass

    def set_paths(self, work_dir, in_map, in_model, lib=None):
        """
        Creates the directory tree for the output and sets the input files paths
        :param work_dir: output directory
        :param lib: map motif library path
        :param in_map: input map path
        :param in_model: input model path
        """

        # print("WDD",self.work_dir)

        if os.path.exists(lib):
            self.lib = os.path.abspath(lib).rstrip('/')
            self.out_dir = os.path.join(self.work_dir, 'vs_' + os.path.basename(self.lib))
            self.lib = os.path.join(self.lib, 'motifs')
            self.check_motif_lib(self.lib)
        else:
            log.info('Motif library %s not found', lib)
        self.in_map = os.path.abspath(in_map)
        self.in_model = os.path.abspath(in_model)
        try:
            os.makedirs(self.out_dir)
        except FileExistsError:
            pass

        shutil.copy(self.in_map, self.input)
        shutil.copy(self.in_model, self.input)

    def set_results_paths(self, project_dir):
        folders = os.listdir(project_dir)

        tmp = os.path.join(project_dir, 'segments')
        if os.path.exists(tmp):
            self.segments = tmp
        for folder in folders:
            if folder.startswith('vs_'):
                path = os.path.join(project_dir, folder, 'scores.json')
                if os.path.exists(path):
                    lib = folder.split('_')[-1]
                    self.scores_paths[lib] = path
        self.all_libs = list(sorted(self.scores_paths.keys()))
        inp = os.path.join(project_dir, 'input')
        if os.path.exists(inp):
            in_files = os.listdir(inp)
            for ext in ['.cif', '.mmcif', '.pdb', '.ent']:
                for f in in_files:
                    if f.endswith(ext):
                        self.in_model = os.path.join(inp, f)
                        break
            for ext in ['.map', '.mrc']:
                for f in in_files:
                    if f.endswith(ext):
                        self.in_map = os.path.join(inp, f)
                        break
        seq_dir = os.path.join(project_dir, 'sequence')
        if os.path.exists(seq_dir):
            in_files = os.listdir(seq_dir)
            for f in in_files:
                if not f.startswith('.'):
                    seq_tmp = self.load_sequences(os.path.join(seq_dir, f))
                    if len(seq_tmp.keys()) > 0:
                        for k, v in seq_tmp.items():
                            self.sequence_dict[k] = v
        else:
            os.makedirs(seq_dir)
            self.sequence_dir = seq_dir

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

    def segment_structure(self, n_cores=4, voxel=0.2, replace=True):
        """
        Chops the input map and model into amino-acid fragments
        :param n_cores: number of computing cores
        :param voxel: library voxel size
        :param replace: replace or not fragments chopped in during previous runs
        :return: segments list file path (json format)
        """
        # try:
        #     in_model = bioUtils.load_structure(self.in_model)[0]
        #     process_log.info('Failed to load structure using "auth" cif records. Attempting to use "label" records')
        # except:
        #     in_model = bioUtils.load_structure_label_id(self.in_model)[0]

        # del_atom_list = []
        # for chain in in_model:
        #     for residue in chain:
        #         for atom in residue:
        #             if atom.get_name().upper().startswith('H'):
        #                 del_atom_list.append((chain.get_id(), residue.get_id(), atom.get_id()))
        #
        # for chain, res, atom in del_atom_list:
        #     del in_model[chain][res][atom]

        segments_list_path = os.path.join(self.segments, self.segmented_pairs)
        if os.path.exists(segments_list_path):
            return segments_list_path

        self.extract = ExtractMap(map_path=self.in_map, model_path=self.in_model, work_voxel_size=voxel)
        residues = [r for r in self.extract.in_model.get_residues() if r.get_resname().upper() in self.residues]
        # residues = residues[:1]
        log.info('Found %s residues in the input model', len(residues))
        log.info('Running map chopping using %s cores', n_cores)
        # map_obj_list, shifts_list = self.chop.chop_cube_list(residues, self.in_map, 4, zero_origin=False)

        split_residues = func.split_list(residues, n_cores)
        # split_map_obj_list = func.split_list(map_obj_list, n_cores)
        # split_shifts_list = func.split_list(shifts_list, n_cores)
        thread_list = []
        for i in range(len(split_residues)):
            t = Process(target=self._slave,
                        args=(split_residues[i], self.segments, replace))
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

    def _slave(self, residues, work_dir, replace):
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

        for residue in residues:
            res_type = residue.get_resname().lower()
            chain_id = residue.parent.id
            res_nr = residue.id[1]
            name_prefix = '{}-{}-{}'.format(res_type, res_nr, chain_id)

            res_cube_map_path = os.path.join(work_dir, name_prefix + '_cube_res.mrc')
            res_trace_map_path = os.path.join(work_dir, name_prefix + '_trace_res.mrc')
            res_walk_map_path = os.path.join(work_dir, name_prefix + '_walk_res.mrc')
            fin_map_path = os.path.join(work_dir, name_prefix + '_soft.mrc')

            if not os.path.exists(fin_map_path) or replace:
                ca_map, trace_map, walk_map, residue_map = self.extract.segment_residue_density(residue)
                if np.isnan(np.sum(residue_map.data)):
                    log.error("NaN values in {}".format(fin_map_path))
                residue_map.write_map(fin_map_path)
                ca_map.write_map(res_cube_map_path)
                trace_map.write_map(res_trace_map_path)
                walk_map.write_map(res_walk_map_path)
            self.chopped_pairs.append(fin_map_path)

    def compute_correlations(self, pairs_json, sd_level, n_cores=4, recompute=False, verbose=False):
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
        updated_script = self.score_chimera_script
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

        if sd_level:
            sd = f'-sd {sd_level}'
        else:
            sd = ''

        command = f'{self.chimera_path} --nogui {v_flag} {updated_script}' \
                  f' -s {pairs_json} -m {self.in_model} -l {self.lib} -np {n_cores} {sd} -o {json_out} {recompute} > chimera.process_log'

        subprocess.call(command, cwd=self.out_dir, shell=True)
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

    def score_sequence(self, seq, score_list, cut_off):
        if len(seq) != len(score_list):
            raise Exception("The length of the sequence_dir should be equal to the number of residues in the scores list")
        seq_score = 0
        for nr, char in enumerate(seq):
            try:
                code = self.aa_1to3[char.upper()].lower()
            except KeyError:
                continue
            try:
                corr = score_list[nr][self.keys.SCORES][code][1]
            except KeyError:
                print(f'Code {code} not valid')
                corr = 0
            # all_c = [i[1] for i in score_list[nr][self.keys.SCORES].values()]
            # quart = np.percentile(all_c, 80)
            # aver = sum(all_c)/len(all_c)
            if corr < cut_off:
            # if corr < quart:
                continue
            else:
                seq_score += corr
        # print(f'SEQ: {seq} SCORE: {seq_score}')
        return seq_score

    def score_sequence_(self, seq, score_list, cut_off=0.9):
        hgh = 'FTYWPHRL'
        if len(seq) != len(score_list):
            raise Exception("The length of the sequence_dir should be equal to the number of residues in the scores list")
        seq_score = 0
        for nr, char in enumerate(seq):
            try:
                code = self.aa_1to3[char.upper()].lower()
            except KeyError:
                continue
            try:
                corr = score_list[nr][self.keys.SCORES][code+'-class'][1]
            except KeyError:
                corr = 0
            all_c = [i[1] for i in score_list[nr][self.keys.SCORES].values()]
            aver = sum(all_c)/len(all_c)
            if not char.upper() in hgh:
                corr = corr * 0.9

            seq_score += corr

        return seq_score

    @staticmethod
    def load_sequences(file_path):
        with open(file_path, 'r') as fileobject:
            lines = fileobject.readlines()
        sequences = {}
        sequence = ''
        name = ''
        for line in lines:
            if not name and line[0] == '>':
                name = line[1:].rstrip()
            elif name and line[0] == '>':
                sequences[name] = sequence
                name = line[1:].rstrip()
                sequence = ''
            else:
                for char in line:
                    if char != '\n':
                        sequence += char
        sequences[name] = sequence
        return sequences

    @staticmethod
    def fragment_sequence(sequence, length):
        seq_lst = []
        for i in range(len(sequence) - length + 1):
            seq_lst.append(sequence[i:length + i])
        return seq_lst

    @staticmethod
    def get_sequence(model_path):
        model = bio_utils.load_structure(model_path)[0]
        ppb = PPBuilder()
        out_pp = []
        for pp in ppb.build_peptides(model):
            out_pp.append(pp.get_sequence())
        return out_pp

    def get_cutt_off(self, score_lst, percentile):
        all_c = []
        for nr in range(len(score_lst)):
            all_c += [i[1] for i in score_lst[nr][self.keys.SCORES].values()]
        print('ALL corr', len(all_c))
        return np.percentile(all_c, percentile)


    def test(self, score_lst, cs, cut_off):

        print("Sequence length: {}".format(len(score_lst)))
        rev_score_lst = score_lst[-1::-1]

        seq = self.sequence_dict
        # print(seq)
        keys = seq.keys()
        # print(keys)
        # print('Correct sequence: {}'.format(seq[keys[0]][start:end]))
        print('\nTop score sequences:')
        seq_scores = []
        length = len(score_lst)
        nr_s = 0
        tot_res = 0
        for full_seq in seq.values():
            # print(full_seq)
            nr_s += 1
            tot_res += len(full_seq)
            fragments = self.fragment_sequence(full_seq, length)
            # fragments.append('YAGDNFVRYTGDTISMAQTQLFAWEA')
            for fragment in fragments:
                f_score = self.score_sequence(fragment, score_lst, cut_off=cut_off)
                seq_scores.append((f_score, fragment, 'Forward'))
                f_score = self.score_sequence(fragment, rev_score_lst, cut_off=cut_off)
                seq_scores.append((f_score, fragment, 'Reversed'))

        seq_scores.sort(key=lambda z: z[0], reverse=True)
        # print("Total Sequences: {}".format(len(seq_scores)))
        print(f'Total sequences: {len(seq.values())}. Total fragments: {len(seq_scores)}')
        z = []
        for i, pair in enumerate(seq_scores[:]):

            scores = np.array([x[0] for x in seq_scores[:i] + seq_scores[i + 1:]])
            aver = scores.sum() / len(scores)

            Z = (pair[0] - aver) / scores.std()
            z.append(Z)
            # cs = 'YAGDNFVRYTGDTISMAQTQLFAWEA'
            # cs = 'LSRPERPDLVFEEEDLPYEEEIMRNQ'
            if pair[1] == cs:
                print('Correct sequencs:', cs)
                print('CORRECT Z= {}'.format(Z))

            if i < 10:
                seq_code = ''
                for key, value in seq.items():
                    if pair[1] in value:
                        seq_code = key
                print('Score= {:.4f}, Z= {:.4f}, Fragment: {} Parent seq name: {}, Direction: {}'.format(pair[0], Z,
                                                                                                         pair[1],
                                                                                                         seq_code,
                                                                                                         pair[2]))
            elif i > 1000:
                break
        all_scores = np.array([x[0] for x in seq_scores])
        aver = all_scores.sum() / len(all_scores)

        z = np.array(z)
        aver_z = z.sum() / len(z)

        print('Average score: {}'.format(aver))
        print('Average Z score: {}'.format(aver_z))
        print(nr_s)
        print(tot_res)

    def find_fragment(self, fragment_scores_lst, cut_off, correct_seq=None, verbose=True):
        # cut_off = self.get_cutt_off(fragment_scores_lst, 75)
        length = len(fragment_scores_lst)
        if verbose:
            print(f'Correct sequence: {correct_seq}')
            print("Sequence length: {}".format(len(fragment_scores_lst)))

        # rev_score_lst = fragment_scores_lst[-1::-1]
        seq = self.sequence_dict

        seq_scores = []
        nr_s = 0
        tot_res = 0
        for full_seq in seq.values():
            nr_s += 1
            tot_res += len(full_seq)
            fragments = self.fragment_sequence(full_seq, length)
            for fragment in fragments:
                f_score = self.score_sequence(fragment, fragment_scores_lst, cut_off=cut_off)
                seq_scores.append((f_score, fragment, 'Forward'))
                # f_score = self.score_sequence(fragment, rev_score_lst, cut_off=cut_off)
                # seq_scores.append((f_score, fragment, 'Reversed'))

        seq_scores.sort(key=lambda z: z[0], reverse=True)
        if verbose:
            print(f'Total sequences: {len(seq.values())}. Total fragments: {len(seq_scores)}')

        place = -1
        cor_Z = 0
        sc = 0
        for i, pair in enumerate(seq_scores[:]):

            scores = np.array([x[0] for x in seq_scores[:i] + seq_scores[i + 1:]])
            aver = scores.sum() / len(scores)

            Z = (pair[0] - aver) / scores.std()

            if i < 1000:
                seq_code = ''
                for key, value in seq.items():
                    if pair[1] in value:
                        seq_code = key
                if verbose:
                    print('Score= {:.4f}, Z= {:.4f}, Fragment: {} Parent seq name: {}, Direction: {}'.format(pair[0], Z,
                                                                                                         pair[1],
                                                                                                         seq_code,
                                                                                                         pair[2]))
                if pair[1] == correct_seq and place == -1:
                    place = i
                    cor_Z = Z
                    sc = pair[0]
            else:
                break

        return place, cor_Z, sc

    def find_fragment_per_chain_not_finished(self, fragment_scores_lst, cut_off, correct_seq=None, verbose=True):
        # cut_off = self.get_cutt_off(fragment_scores_lst, 75)
        length = len(fragment_scores_lst)
        if verbose:
            print(f'Correct sequence: {correct_seq}')
            print("Sequence length: {}".format(len(fragment_scores_lst)))

        rev_score_lst = fragment_scores_lst[-1::-1]
        seq = self.sequence_dict

        place = -1
        cor_Z = 0
        sc = 0
        for full_seq in seq.values():
            seq_scores = []

            fragments = self.fragment_sequence(full_seq, length)
            for fragment in fragments:
                f_score = self.score_sequence(fragment, fragment_scores_lst, cut_off=cut_off)
                seq_scores.append((f_score, fragment, 'Forward'))
                f_score = self.score_sequence(fragment, rev_score_lst, cut_off=cut_off)
                seq_scores.append((f_score, fragment, 'Reversed'))

            seq_scores.sort(key=lambda z: z[0], reverse=True)
            if verbose:
                print(f'Total sequences: {len(seq.values())}. Total fragments: {len(seq_scores)}')


            for i, pair in enumerate(seq_scores[:]):

                scores = np.array([x[0] for x in seq_scores[:i] + seq_scores[i + 1:]])
                aver = scores.sum() / len(scores)

                Z = (pair[0] - aver) / scores.std()

                if i < 1000:
                    seq_code = ''
                    for key, value in seq.items():
                        if pair[1] in value:
                            seq_code = key
                    if verbose:
                        print('Score= {:.4f}, Z= {:.4f}, Fragment: {} Parent seq name: {}, Direction: {}'.format(pair[0], Z,
                                                                                                             pair[1],
                                                                                                            seq_code,
                                                                                                           pair[2]))
                    if pair[1] == correct_seq and place == -1:
                        place = i
                        cor_Z = Z
                        sc = pair[0]
                else:
                    break

        return place, cor_Z, sc

    def score_fragments_old(self, scores_lst, length, cut_off, verbose=True):
        all_frag = []
        # c_off = self.get_cutt_off(scores_lst, 75)
        transl = nomenclature.AA_3LET_TO_1LET
        l = len(scores_lst)

        for i in range(l - length):
            seq = ''
            tmp = scores_lst[i:i+length]
            for e in tmp:
                code = transl[e["residue_type"].upper()]
                seq += code
            place, Z, sc = self.find_fragment(tmp, cut_off, seq, verbose)
            all_frag.append({'place': place, 'z': Z, 'seq': seq})
            # print(place, Z, seq, sc)
        return all_frag

    def find_record(self, score_list, chain, r_nr):
        for record in score_list:
            if record['chain'] == chain and record['residue_nr'] == r_nr:
                return record


    def score_fragments(self, scores_lst, length, cut_off, model_path, verbose=True):
        model = bio_utils.load_structure(model_path)[0]
        transl = nomenclature.AA_3LET_TO_1LET
        all_frag = []
        not_found = ''
        for ch in model.get_chains():
            if len(ch) < length:
                continue
            else:
                l = len(ch)
                all_res = [r for r in ch.get_residues()]
                for i in range(l - length):
                    fragm_scores = []
                    for res in all_res[i:i + length]:
                        c = res.parent.id
                        r_nr = res.id[1]
                        rec = self.find_record(scores_lst, c, r_nr)
                        if rec is not None:
                            fragm_scores.append(rec)
                        else:
                            not_found = f'{c}-{r_nr}'
                            break
                    seq = ''
                    for e in fragm_scores:
                        code = transl[e["residue_type"].upper()]
                        seq += code
                    if len(fragm_scores) == length:
                        place, Z, sc = self.find_fragment(fragm_scores, cut_off, seq, verbose)
                        all_frag.append({'place': place, 'z': Z, 'seq': seq, 'not_found': None})
                    else:
                        all_frag.append({'place': -1, 'z': 0, 'seq': '', 'not_found': not_found})

        return all_frag


def main():
    parser = argparse.ArgumentParser(description='Map Validation')
    parser.add_argument("-M", "--mode", dest="mode", required=True, help="Running mode")
    parser.add_argument("-s", "--sequence_dir", dest="seq", required=False, help="Atomic model path")
    parser.add_argument("-p", "--model", dest="in_model", required=True, help="Atomic model path")
    parser.add_argument("-m", "--map", dest="in_map", required=False, help="Map file path")
    parser.add_argument("-l", "--lib", dest="lib", required=False, help="Strudel motif library path")
    parser.add_argument("-np", "--n_processors", dest="np", required=False, default=2, help="Number of processors")
    parser.add_argument("-o", "--out", dest="out", required=True, help="Output directory")
    parser.add_argument("-process_log", "--process_log", dest="process_log", default=None, required=False, help="Log file")
    parser.add_argument("-v", "--voxel", dest="voxel", required=False, default=0.25, type=float, help="Segments voxel size")
    parser.add_argument("-r", "--recompute", dest="recompute_scores", action='store_true',
                        help="Recalculate correlations")
    parser.add_argument("-rs", "--recompute_segments", dest="recompute_segments", action='store_true',
                        help="Repeat segmenting of the input model and map if found in the output")
    parser.add_argument("-vc", "--verbose_chimerax", dest="v_c", action='store_true',
                        help="Run ChimeraX in verbose mode")
    parser.add_argument("-wl", "--warning_level", dest="warning_level", default='info',
                        help="Log file warning level [info[debug]")
    parser.add_argument("-sd", "--sd_level", dest="sd_level", required=False,
                        help="SD level for correlations calculations")
    parser.add_argument("-co", "--cut_off", dest="cut_off", type=int, required=False,
                        help="Cut-off quantile")


    args = parser.parse_args()
    if args.log:
        log_file = args.log
    else:
        log_file = os.path.join(args.out, 'sequence_finder.process_log')

    if not os.path.exists(args.out):
        os.makedirs(args.out)
    # print(args.out)

    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(levelname)s:  %(message)s')
    date_time = datetime.now().strftime("%H:%M:%S %Y-%m-%d")
    logging.info('\n{:_^100}'.format('SequenceFinder') + '\n\nStarted: {}\n'.format(date_time))

    args.np = int(args.np)
    start = time.time()

    compute = ComputeScores()
    compute.set_dir_tree(work_dir=args.out)
    if args.lib:
        lib = args.lib
    else:
        lib = None

    if args.mode == 'C':
        compute.set_paths(work_dir=args.out, in_map=args.in_map, in_model=args.in_model, lib=lib)
        chopped_pairs = compute.segment_structure(n_cores=args.np, voxel=args.voxel, replace=args.recompute_segments)
        # chopped_pairs = os.path.join(args.out, 'segments', 'segments_list.json')
        json_file_path = compute.compute_correlations(chopped_pairs, args.sd_level, args.np, args.recompute_scores, verbose=args.v_c)
        #
        # prefix = os.path.splitext(json_file_path)[0]
    elif args.mode == 'F':
        # print(args.seq)
        if args.seq and os.path.exists(args.seq):
            shutil.copy(args.seq, compute.sequence_dir)
        compute.set_results_paths(args.out)
        # print("CORRECT SEQUENCE\n")
        # seq = compute.get_sequence(compute.in_model)
        # print(seq)
        for l, p in compute.scores_paths.items():
            scores_lst = func.read_json(p)
            tmp = []
            count = []
            for entry in scores_lst:
                if entry['residue_nr'] not in count:
                    tmp .append(entry)
                    count.append(entry['residue_nr'])
            # print("COUNTER", count)
            tmp = tmp[:]

            seq = ''
            transl = nomenclature.AA_3LET_TO_1LET
            for e in tmp:
                code = transl[e["residue_type"].upper()]
                seq += code
            print('Correct sequencs:', seq)
            c_off = compute.get_cutt_off(tmp, 75)
            print(c_off)
            # compute.test(tmp1, seq, cut_off=c_off)

            # result = compute.find_fragment(tmp[10:25], c_off, seq, True)
            # print(result)
            # o_j_p = '/Volumes/data/Work/test_segm/out_6/vs_2.5-3_formated/res_jsons_1000_60'
            # o_j_p = '/Volumes/data/find_seq_20584/out/vs_motifs_3.2-3.5/jsons_75'
            # o_j_p = '/Volumes/data/find_seq/untitled_folder/0711_out_sd1/jsons_75'
            compute.sequence_dict['scored'] = seq
            o_j_p = os.path.join(args.out, f'jsons_{args.cut_off}')
            if not os.path.exists(o_j_p):
                os.makedirs(o_j_p)
            for i in range(5, 40):
                o_j = os.path.join(o_j_p, f'{i}.json')
                res = compute.score_fragments(tmp, i, c_off, args.in_model, False)
                func.save_json(o_j, res)

                r_sc = [1 for i in res if i['place'] == 0]
                print(len(r_sc)/len(res) * 100)




    logging.info('\nElapsed: {}\n{:_^100}'.format(func.report_elapsed(start), ''))


if __name__ == '__main__':
    main()