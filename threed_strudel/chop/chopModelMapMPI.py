"""
chopModelMap.py

Class for Chopping an atomic residue (and EM map) into residues and residue side chains which are saved
as separate files.
The chopped side chains are used to chop the corresponding map.
This code uses the chopMap module for map chopping.

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


import copy
import os
import numpy as np
from datetime import datetime
import time
import json
from multiprocessing import Process, Array
try:
    from mpi4py import MPI
    mpi_avail = True
except ImportError:
    mpi_avail = False

from random import shuffle

from threed_strudel.chop.chop_map import ChopMap, MapParser
import threed_strudel.utils.functions as func
from threed_strudel.utils import bio_utils
import threed_strudel.nomenclature as nomenclature


# noinspection PyTypeChecker
class ChopModelMap:
    """
    Class for chopping residues and residue side chains out of proteins. The chopped
    side chains are used as guide for density map chopping
    """

    def __init__(self, rank=0, parallelism="shared"):
        self.rank = rank
        self.parallelism = parallelism
        if self.parallelism == "mpi" and mpi_avail:
            self.comm = MPI.COMM_WORLD
        else:
            self.comm = None
            self.parallelism = "shared"
        self.start_time = None
        self.min_rscc = None
        self.no_charged_check = None
        self.chop_log = None
        self.allowed_b = 0.0
        self.entry_code = ''
        self.work_dir = ''
        self.tmp_dir = None
        self.tmp_log_dir = None
        self.master_log_path = None
        self.input_model = None
        self.input_model_path = ''
        self.map_file_path = ''
        self.map_object = None
        self.rscc_list = []
        self.check_rscc = None
        self.side_chain_end = ''
        self.residue_end = ''
        self.shift_side_end = ''
        self.shift_residue_end = ''
        self.cube_res_end = ''
        self.cube_end = ''
        self.soft_map_end = ''
        self.hard_map_end = ''
        self.cube_radius = None
        self.final_voxel = None
        self.chop_radius = None
        self.chop_soft_radius = None
        # self.script_path = os.path.dirname(os.path.abspath(__file__))
        self.charged_res = nomenclature.CHARGED_RESIDUES
        # self.inclusion_threshold = None
        self.check_b_values = None
        # self.known_inclusion_levels = os.path.join(self.script_path, 'known_inclusion_levels.txt')
        self.inclusion_charged_tolerance = 3
        self.heavy_atoms = nomenclature.HEAVY_ATOMS

    def set_env(self, work_dir,
                entry_code,
                model_path,
                map_path,
                log_path,
                rscc_list=None,
                allowed_b=None,
                filter_identical_chains=True,
                warning_level='info'):
        """
        Set the environment for residue (and density map) chopping
        :param warning_level:
        :param inclusion_fraction: fraction of atoms inside the threshold
        :param allowed_b:
        :param rscc_list:
        :param log_path:
        :param work_dir: The directory for output files
        :param entry_code: EMD or pdb code of the entry
        :param model_path: Path to the residue pdb or cif file
        :param map_path: Path to the residue map file
        :param warning_level: warning level
        """
        self.start_time = time.time()
        s_env_t = time.time()
        self.work_dir = work_dir
        self.create_dir(work_dir)
        self.entry_code = entry_code
        self.input_model_path = model_path
        self.map_file_path = map_path
        self.master_log_path = log_path
        self.rscc_list = rscc_list

        self.tmp_log_dir = os.path.join(work_dir, 'tmp_log')
        self.create_dir(self.tmp_log_dir)
        tmp_log_path = os.path.join(self.tmp_log_dir, '{}{}{}'.format('rank_', self.rank, '.log'))

        if self.rank == 0:
            self.chop_log = func.setup_logger('chop_log', log_path, warning_level=warning_level)
            date_time = datetime.now().strftime("%H:%M:%S %Y-%m-%d")
            text = '\n{:_^100}'.format('ChopModelMap') + '\n\nStarted: {}'.format(date_time)
            self.chop_log.info('%s\n', text)
            self.chop_log.info('Reading the input files and preparing the environment..\n')
        else:
            if os.path.exists(tmp_log_path):
                os.remove(tmp_log_path)
            self.chop_log = func.setup_logger('chop_log', tmp_log_path, warning_level=warning_level)

        if self.parallelism == "mpi":
            if self.rank == 0:
                start = time.time()
                self.map_object = MapParser(map_path)
                self.chop_log.info('Map loaded in: %s', func.report_elapsed(start))
                self.input_model = bio_utils.load_structure(model_path)
                if filter_identical_chains:
                    self.input_model = self.delete_identical_chains(self.input_model)
            self.input_model = self.comm.bcast(self.input_model, root=0)
        else:
            self.map_object = MapParser(map_path)
            self.input_model = bio_utils.load_structure(model_path)
            if filter_identical_chains:
                self.input_model = self.delete_identical_chains(self.input_model)

        if allowed_b is None:
            self.allowed_b = self.calc_b_value_cut_off()
            if self.rank == 0:
                self.chop_log.info('\nNo B-factor cut-off was specified\n'
                                   'The cut-off value vas set to median + 2 * 1Qu: %s',
                                   self.allowed_b)

        self.set_out_file_suffixes()
        if self.rank == 0:
            self.chop_log.info("Spent on set_paths %s", func.report_elapsed(s_env_t))

    def set_map_chop_parameters(self, cube_radius=5.0, final_voxel=0.5, chop_radius=3.0, chop_soft_radius=2.0,
                                check_rscc=True, min_rscc=0.7, check_b_values=True, no_charged_check=True):
        """
        Set custom parameters for map chopping
        :param check_b_values:
        :param check_rscc:
        :param check_inclusion:
        :param cube_radius: The distance between the cube edge the residue
        :param final_voxel: Final voxel size of the chopped maps
        :param chop_radius: Radius around the molecule for map chopping
        :param chop_soft_radius: Soft hard_radius around the molecule for map chopping
        :param resolution:
        """
        self.cube_radius = cube_radius
        self.final_voxel = final_voxel
        self.chop_radius = chop_radius
        self.chop_soft_radius = chop_soft_radius
        self.check_rscc = check_rscc
        self.min_rscc = min_rscc
        self.check_b_values = check_b_values
        self.no_charged_check = no_charged_check

        # Set RSCC data
        if check_rscc and self.rscc_list is None:
            if self.rank == 0:
                raise Exception('Please provide the map resolution or a list of per residue RSCC')

    def set_out_file_suffixes(self):
        """
        Sets suffixes for output files
        """
        self.side_chain_end = '-{}_side.cif'.format(self.entry_code)
        self.shift_side_end = '-{}_side_shift.cif'.format(self.entry_code)
        self.residue_end = '-{}_residue.cif'.format(self.entry_code)
        self.shift_residue_end = '-{}_residue_shift.cif'.format(self.entry_code)
        self.cube_end = '-{}_cube.mrc'.format(self.entry_code)
        self.cube_res_end = '-{}_cube_resampled.mrc'.format(self.entry_code)
        self.soft_map_end = '-{}_soft.mrc'.format(self.entry_code)
        self.hard_map_end = '-{}_hard.mrc'.format(self.entry_code)

    def delete_identical_chains(self, structure, delta=0.15):
        """
        Deletes duplicating chains in a structure residue based on pairwise RMS (delta)
        :param delta: maximum pairwise RMS to consider chains identical
        :param structure: Biopython structure object
        :return: Biopython structure object
        """
        # if structure.level == 'S':
        #     structure = structure[0]
        chain_ids = [chain.get_id() for chain in structure.get_chains()]
        if self.rank == 0:
            start = time.time()
            self.chop_log.info('\nSearching identical chains..')
            self.chop_log.info('All chains:')
            self.chop_log.info(' '.join(chain_ids) + '\n')

        classes, max_distances = bio_utils.classify_chains(structure, delta=delta)
        duplicates = [cl for cl in classes if len(cl) > 1]
        # if self.rank == 0:
        #     self.chop_log.info(' '.join(classes))
        if len(duplicates) > 0:
            if self.rank == 0:
                self.chop_log.info('Identical chains:')
                for cl in duplicates:
                    self.chop_log.info(' '.join(cl))

            unique_best_chains = bio_utils.select_best_chains(structure, self.map_object, classes)
            if self.rank == 0:
                self.chop_log.info('Selected as best fitting in the map unique chains:')
                self.chop_log.info(' '.join(unique_best_chains))
            for chain_id in chain_ids:
                if chain_id not in unique_best_chains:
                    del structure[0][chain_id]
        elif self.rank == 0:
            self.chop_log.info("No identical chains were found")
            self.chop_log.info('Elapsed: %s', func.report_elapsed(start))
        return structure

    @staticmethod
    def get_residue_list(structure, res):
        """
        Create a list of BIO.PDB residue objects containing all the residues
        of the given type in the input residue
        :param structure: BIO.PDB structure object
        :param res: Residue code
        :return: A list of BIO.PDB residue objects
        """
        del_atom_list = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.get_name().upper().startswith('H'):
                            del_atom_list.append((chain.get_id(), residue.get_id(), atom.get_id()))
        for chain, resi, atom in del_atom_list:
            del structure[0][chain][resi][atom]
        residues = []
        for residue in structure.get_residues():
            if residue.get_resname() == res.upper():
                residues.append(residue)

        return residues

    def calc_b_value_cut_off(self):
        """
        Computes B-factor statistics
        :return: median + 2 * 1st Quartile
        """
        b = []
        for atom in [a for a in self.input_model.get_atoms()]:
            if atom.parent.resname.rstrip() != 'HOH':
                b.append(atom.bfactor)
        np_b = np.array(b)
        median = np.median(np_b)
        qu_1 = np.percentile(np_b, 25, interpolation='lower')
        b_cut = median + 2 * qu_1

        return int(b_cut)

    def check_b_factors(self, residue):
        """
        Analyse the B-factors of a BIO.PDB residue object in order to asses the local quality
        of the density map. Returns true if the B-factors are consistent with high quality map.
        :param residue: BIO.PDB residue object
        :return: True or False
        """
        valid = True
        for atom in residue:
            if atom.bfactor > self.allowed_b:
                valid = False
                break
        return valid

    def run_model_quality_checks(self, residue, chop_map):
        """
        Runs the residue residue quality checks
        :param residue: residue object
        :param chop_map: map dir
        :return: False or True
        """
        good_b = True
        good_rscc = True
        complete = True
        chain_id = residue.parent.id
        res_nr = residue.id[1]
        res_name = residue.get_resname()

        nr_atoms = 0
        single_conf = True
        for atom in residue:
            if atom.is_disordered():
                single_conf = False
            #     break
            # try:
            #     occupancy = atom.last_occupancy
            # except AttributeError:
            #     occupancy = 1
            # if occupancy < 1:
            #     single_conf = False
                self.chop_log.info('The %s %s %s residue has multiple conformations',
                                   residue.get_resname().upper(), residue.id[1], residue.parent.id)
                break
            nr_atoms += 1

        # Check heavy atom completeness
        if nr_atoms < self.heavy_atoms[res_name]:
            self.chop_log.info('Missing heavy atoms in %s %s %s (%s out of %s)', res_name,
                               res_nr, chain_id, nr_atoms, self.heavy_atoms[res_name])
            complete = False
        elif nr_atoms > self.heavy_atoms[res_name]:
            self.chop_log.info('Terminal residue will be skipped %s %s %s (%s out of %s)', res_name,
                               res_nr, chain_id, nr_atoms, self.heavy_atoms[res_name])
            complete = False

        # Check local map quality based on B-factors and atom inclusion
        if chop_map:
            if res_name.upper() not in self.charged_res or not self.no_charged_check:
                # Check B-factors
                if self.check_b_values:
                    good_b = self.check_b_factors(residue)
                    if not good_b:
                        self.chop_log.info('Residue %s %s %s has high B-factors', res_name, res_nr, chain_id)
                # Check map residue real space correlation coefficients
                if self.check_rscc:
                    rscc = self.find_residue_rscc(self.rscc_list, chain_id, res_name, res_nr)
                    if rscc < self.min_rscc:
                        self.chop_log.info('Residue %s %s %s has RSCC lover than %s',
                                           res_name, res_nr, chain_id, self.min_rscc)
                        good_rscc = False

        if good_rscc and good_b and complete and single_conf:
            return True
        else:
            return False

    @staticmethod
    def find_residue_rscc(rscc_list, chain, res_type, res_nr):
        """
        Find the residue RSCC in the list
        :param rscc_list:
        :param chain:
        :param res_type:
        :param res_nr:
        :return:
        """
        result = None
        for residue in rscc_list:
            if residue[0] == chain and residue[1] == res_type and int(residue[2]) == res_nr:
                result = residue[3]
                break
        return result

    def build_all_residue_object_list(self, residue_list):
        """
         Create a list of BIO.PDB residue objects containing all the residues
         to be chopped in the input residue
         :param residue_list: list of residues names
         :return: A list of BIO.PDB residue objects
         """
        all_res_obj_list = []
        for res in residue_list:
            # Create a list of all current type residue objects
            residues = self.get_residue_list(self.input_model, res)
            all_res_obj_list = all_res_obj_list + residues
        all_res_obj_list.sort(key=lambda z: z.get_resname())

        return all_res_obj_list

    def create_dir_structure(self, residue_list):
        """
        Creates the chopping directory structure
        :param residue_list: list of residue names
        """
        for res in residue_list:
            res_dir = self.work_dir + '/' + res.lower()
            self.create_dir(res_dir)
            side_chain_dir = res_dir + '/' + 'side_chain'
            self.create_dir(side_chain_dir)
            model_map_dir = res_dir + '/' + 'models_and_maps'
            self.create_dir(model_map_dir)
            lq_model_map_dir = res_dir + '/' + 'lq_models_and_maps'
            self.create_dir(lq_model_map_dir)

    def chop_worker_task(self, residue_obj_list, in_map, stat_storage,
                         chop_map=True, chopping_mode='soft', save_tmp_files=False):
        """
        This function is a wrapper for the code which will be executed in parallel
        :param stat_storage: list of containers for chopping statistics
        :param in_map: input map as MapObject or map file in_dir
        :param self: The outer scape self object
        :param residue_obj_list: list of residues objects
        :param chop_map: chop both residue and map
        :param chopping_mode: 'soft' or 'hard'
        """
        chop = ChopMap()
        resampling_time = 0
        res_order_dict = stat_storage[0]
        hq = stat_storage[1]
        lq = stat_storage[2]
        tmp_dir = '/tmp/chop'
        self.tmp_dir = tmp_dir
        try:
            os.makedirs(tmp_dir)
        except FileExistsError:
            pass
        map_obj_list, shifts_list = chop.chop_cube_list(residue_obj_list,
                                                        in_map,
                                                        self.cube_radius)

        for i in range(len(residue_obj_list)):
            c_map_obj = map_obj_list[i]
            residue = residue_obj_list[i]
            matrix = shifts_list[i]

            res_type = residue.get_resname().lower()
            chain_id = residue.parent.id
            res_nr = residue.id[1]
            res_dir = self.work_dir + '/' + res_type

            # Run local residue quality checks
            good_model = self.run_model_quality_checks(residue, chop_map)
            for atom in residue:
                atom.disordered_flag = 0
            # Increment ctype objects
            if good_model:
                hq[res_order_dict[res_type]] = hq[res_order_dict[res_type]] + 1
                model_map_dir = res_dir + '/' + 'models_and_maps'
            else:
                lq[res_order_dict[res_type]] = lq[res_order_dict[res_type]] + 1
                model_map_dir = res_dir + '/' + 'lq_models_and_maps'
            side_pdb_dir = res_dir + '/' + 'side_chain'
            # Save the residue as pdb file
            name_prefix = '{}-{}-{}'.format(res_type, res_nr, chain_id)
            # Delete the main chain atoms
            side_chain = residue
            tmp_res_path = os.path.join(tmp_dir, name_prefix + '.cif')
            struct = bio_utils.residues2structure(residue)
            bio_utils.save_model(struct, tmp_res_path)
            residue = bio_utils.load_structure(tmp_res_path)
            os.remove(tmp_res_path)
            side_chain = self.del_main_chain(side_chain)

            # Save the side chain
            cube_map_path = os.path.join(tmp_dir, name_prefix + self.cube_end)
            if save_tmp_files:
                c_map_obj.write_map(cube_map_path)
            # Resampled map path
            cube_new_grid_path = os.path.join(tmp_dir, name_prefix + self.cube_res_end)
            start = time.time()

            # chop.grid_resample(cube_map_path, cube_new_grid_path, self.final_voxel)

            c_map_obj.grid_resample_emda(self.final_voxel)
            if save_tmp_files:
                c_map_obj.write_map(cube_new_grid_path)

            end = time.time()
            resampling_time += (end - start)
            bio_utils.shift_coord(matrix, side_chain)
            if chopping_mode.lower() == 'hard':
                fin_map = os.path.join(model_map_dir, name_prefix + self.hard_map_end)
                # chop.chop_hard_radius(side_chain, cube_new_grid_path, fin_map, self.chop_radius)
                chop.chop_hard_radius(side_chain, c_map_obj, fin_map, self.chop_radius)
            elif chopping_mode.lower() == 'soft':
                fin_map = os.path.join(model_map_dir, name_prefix + self.soft_map_end)
                # chop.chop_soft_radius(side_chain, cube_new_grid_path, fin_map, self.chop_radius,
                #                       self.chop_soft_radius)
                chop.chop_soft_radius(side_chain, c_map_obj, fin_map, self.chop_radius,
                                      self.chop_soft_radius)
            else:
                raise Exception('The chopping_mode parameter can be hard or soft')

            # try:
            #     os.remove(cube_map_path)
            #     os.remove(cube_new_grid_path)
            # except FileNotFoundError:
            #     pass

            fin_res = os.path.join(model_map_dir, name_prefix + self.shift_residue_end)
            side_path = os.path.join(side_pdb_dir, name_prefix + self.side_chain_end)
            side_struct = bio_utils.residues2structure(side_chain)
            bio_utils.save_model(side_struct, side_path)
            bio_utils.shift_coord(matrix, residue)
            # res_struct = bioUtils.residues2structure(residue)
            bio_utils.save_model(residue, fin_res)
        self.chop_log.debug("Time spent on resampling: %s s.", resampling_time)
        if self.rank != 0:
            self.comm.send([res_order_dict, hq, lq], dest=0, tag=1)
        return [res_order_dict, hq, lq]

    def chop_model_map_parallel(self, residue_list, chop_map=True, chopping_mode='soft',
                                nr_cores=2):
        """
        Similar to chop_model_map but runs in parallel
        Chop the specified residues and their side chains out of the input residue. Each residue and
        the corresponding side chain is saved in separate pdb files. The residues are classified as
        having high quality (highq_residue) or low quality (lowq_residue) local density map,
        superimposed and saved in multimodel pdb files.
        :param nr_cores: Number of cores
        :param residue_list: The list of residue type to be chopped (ex. [ASP, TYR, ARG]).
        :param chop_map: Boolean. To perform map chopping or not.
        :param chopping_mode: 'soft' or 'hard'
        """
        start = time.time()
        chop = ChopMap()
        residue_list.sort()
        all_res_obj_list = []

        if self.rank == 0:
            # Cut the map around initial residue to save memory
            map_obj_list, shifts_l = chop.chop_cube_list([self.input_model], self.map_object, 4)
            if not map_obj_list:
                return None

            bio_utils.shift_coord(shifts_l[0], self.input_model)

            # Update the main map object
            self.map_object = map_obj_list[0]

            # Save the residue and map which will be actually used for chopping
            bio_utils.save_model(self.input_model, os.path.join(self.work_dir, self.entry_code + '_shifted.cif'))
            cut_map_path = os.path.join(self.work_dir, self.entry_code + '_model_part.mrc')
            self.map_object.write_map(cut_map_path)
            self.map_file_path = cut_map_path

            # Get the lists of residues for chopping
            all_res_obj_list = self.build_all_residue_object_list(residue_list)
            split_all_res_obj_list = func.split_sorted_list(all_res_obj_list, nr_cores)
            for sublist in split_all_res_obj_list:
                shuffle(sublist)
            # Create the chopping directory structure
            self.create_dir_structure(residue_list)

        else:
            split_all_res_obj_list = None

        # Dictionary to track the residue position in the list ("ASP": index)
        res_order_dict = {}
        for index, r in enumerate(residue_list):
            res_order_dict[r.lower()] = index
        if self.parallelism == "mpi":
            if self.rank == 0:
                hq = [0 for _ in range(len(residue_list))]
                lq = [0 for _ in range(len(residue_list))]
                stat_storage = []
                for i in range(nr_cores):
                    stat_storage.append(copy.deepcopy([res_order_dict, hq, lq]))
            else:
                stat_storage = None
        else:
            # Create ctype objects which are available to all Processes
            hq = Array('i', [0 for _ in range(len(residue_list))])
            lq = Array('i', [0 for _ in range(len(residue_list))])
            stat_storage = [res_order_dict, hq, lq]

        if self.rank == 0:
            self.chop_log.info('Spent on prep in chop_model_map: %s', func.report_elapsed(start))
            self.chop_log.info('All ready for chopping\nElapsed: %s', func.report_elapsed(self.start_time))
            if self.parallelism == 'mpi':
                mem = 'distributed'
            else:
                mem = 'shared'
            self.chop_log.info('\nStarting chopping on %s cores with %s memory parallelism..', nr_cores, mem)

            txt = ''
            for i, r in enumerate(residue_list):
                if i % 10 == 0 and i > 0:
                    txt = txt + '\n          '
                else:
                    txt = txt + r + ' '

            self.chop_log.info('Chopping: %s', txt)
            self.chop_log.info('%s residues to chop', len(all_res_obj_list))

        # Setup and run parallel processes
        start = time.time()

        if self.parallelism == 'mpi':
            if self.rank == 0:
                self.chop_log.info('\nBroadcasting data..')
            # distributed memory
            split_all_res_obj_list = self.comm.bcast(split_all_res_obj_list, root=0)
            self.map_object = self.comm.bcast(self.map_object, root=0)
            self.input_model = self.comm.bcast(self.input_model, root=0)
            self.map_file_path = self.comm.bcast(self.map_file_path, root=0)
            # self.inclusion_threshold = self.comm.bcast(self.inclusion_threshold, root=0)
            stat_storage = self.comm.scatter(stat_storage, root=0)
            self.chop_log.info("\nRank: %s spent %s on broadcasting", self.rank, func.report_elapsed(start))
            start = time.time()
            # Run
            stat_storage = self.chop_worker_task(split_all_res_obj_list[self.rank],
                                                 self.map_object,
                                                 stat_storage,
                                                 chop_map=chop_map,
                                                 chopping_mode=chopping_mode)
            self.chop_log.info('Rank %s finished chopping of %s residues in %s',
                               self.rank, len(split_all_res_obj_list[self.rank]), func.report_elapsed(start))
            self.chop_log.info('%s', '{:*^80}'.format(''))
            stat_storage = self.comm.gather(stat_storage, root=0)
        else:
            # shared memory
            if __name__ == '__main__':
                thread_list = []
                for residue_obj_list in split_all_res_obj_list:
                    t = Process(target=self.chop_worker_task,
                                args=(residue_obj_list, self.map_object, stat_storage, chop_map, chopping_mode))
                    thread_list.append(t)
                for thread in thread_list:
                    thread.start()
                for thread in thread_list:
                    thread.join()
                self.chop_log.info('Finished chopping\nElapsed: %s', func.report_elapsed(start))
                self.chop_log.info('%s', '\n{:*^80}\n'.format(''))

        # Build a dictionary to store chopping statistics ( {res_name: [hq_nr, lq_nr]} )
        if self.rank == 0:
            # Collect the statistics from all cores
            if self.parallelism == "mpi":
                tmp = stat_storage[0]
                for storage in stat_storage[1:]:
                    for i in range(len(tmp[1])):
                        tmp[1][i] += storage[1][i]
                        tmp[2][i] += storage[2][i]
                stat_storage = tmp
            hq = stat_storage[1]
            lq = stat_storage[2]

            # Add workers logs to the master log
            with open(self.master_log_path, 'a') as master_log:
                tmp_files = os.listdir(self.tmp_log_dir)
                for file in tmp_files:
                    tmp_file_path = os.path.join(self.tmp_log_dir, file)
                    with open(tmp_file_path, 'r') as tmp_log:
                        tmp_content = tmp_log.read()
                    master_log.write(tmp_content)

            statistics_dict = {}
            for i in range(len(residue_list)):
                statistics_dict[residue_list[i]] = [hq[i], lq[i]]

            text = self.format_statistics(statistics_dict)
            self.chop_log.info(text)

            date_time = datetime.now().strftime("%H:%M:%S %Y-%m-%d")
            self.chop_log.info("\nChopModelMap process successfully finished  %s", date_time)
            self.chop_log.info("Elapsed: %s", func.report_elapsed(self.start_time))
            self.chop_log.info('%s', '{:_^100}'.format(''))
            return statistics_dict
        else:
            # try:
            #     rmtree(self.tmp_dir)
            # except FileNotFoundError:
            #     pass
            return None

    @staticmethod
    def format_statistics(statistics_dict):
        """
        Formats the chopping statistics for the log output
        :param statistics_dict: Choping statistics dictionary ex: {res_name: [hq_nr, lq_nr]}
        :return: Formatted text
        """
        line1 = '{:<10}'.format('Residue')
        line2 = '{:<10}'.format(' High Q')
        line3 = '{:<10}'.format(' Low  Q')
        total_hq = 0
        total_lq = 0
        for key in sorted(statistics_dict.keys()):
            line1 = line1 + '{:>6}'.format(key)
            line2 = line2 + '{:>6}'.format(statistics_dict[key][0])
            line3 = line3 + '{:>6}'.format(statistics_dict[key][1])
            total_hq = total_hq + statistics_dict[key][0]
            total_lq = total_lq + statistics_dict[key][1]

        text = '\nResidue quality statistics:\n{}\n{}\n{}'.format(line1, line2, line3)
        text2 = '\n\n{:<18}{}\n{:<18}{}\n{:<18}{}'.format('Total:', total_hq + total_lq, 'High quality:', total_hq,
                                                          'Low quality:', total_lq)
        text = text + text2
        return text

    @staticmethod
    def del_main_chain(residue):
        """
        Delete main chain atoms
        :param residue: BIO.PDB object
        :return: BIO.PDB object
        """
        for atom in residue:
            if atom.get_name() == "O":
                del residue[atom.get_name()]
        for atom in residue:
            if atom.get_name() == "C":
                del residue[atom.get_name()]
        # for atom in residue:
        #    if atom.get_name() == "CA":
        #        del residue[atom.get_name()]
        for atom in residue:
            if atom.get_name() == "N":
                del residue[atom.get_name()]
        return residue

    @staticmethod
    def create_dir(path):
        """
        Create new directory
        :param path: new directory in_dir
        """
        try:
            os.makedirs(path)
        except FileExistsError:
            pass

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description='Chops a map and an atomic model into residues and corresponding maps')
    parser.add_argument("-p", "--pdb_cif", dest="in_model", required=True, help="Input model in pdb or cif format")
    parser.add_argument("-m", "--map", dest="in_map", required=True, help="Input map")

    parser.add_argument("-o", "--chop_dir", dest="chop_dir", required=True, help="Where to output the generated files")
    parser.add_argument("-c", "--chop_par", dest="chop_par", nargs='*', required=False,
                        help="Chopping parameters. Example: cube_radius=4 final_voxel=0.5 chop_radius=2.0 "
                             "chop_soft_radius=1.0 inclusion_fraction=90 chop_map=True chopping_mode=soft "
                             "check_rscc=True")
    parser.add_argument("-r", "--rscc_file", dest="rscc_file", required=False,
                        help='Per residue RSCC values in json format EX: [ ["A", "MET", "1", 0.827],..]')
    parser.add_argument("-re", "--resolution", dest="res", required=False,
                        help="Map resolution required if RSCC were not provided")
    parser.add_argument("-rl", "--residue_list", dest="res_list", required=False, nargs='*',
                        help="Residue names to be chopped Ex: ASP ARG ..")
    parser.add_argument("-np", "--n_processors", dest="np", required=False,
                        help="Number of processors in chopping mode")
    parser.add_argument("-mpi", "--mpi", dest="mpi", action='store_true',
                        help="Use distributed memory parallelism")
    parser.add_argument("-filter_chain", "--filter_chain", dest="filter", action='store_true', default=False,
                        help="Do not chop identical chains")

    default_parameters = {"cube_radius": 4, "final_voxel": 0.25, "chop_radius": 2.0, "chop_soft_radius": 1.0,
                          "chop_map": True, "chopping_mode": "soft", "check_rscc": True}
    args = parser.parse_args()
    if args.mpi:
        sys.excepthook = func.global_except_hook
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
    else:
        rank = 0
        size = 1

    #config = Configure()
    chop = ChopModelMap(rank, parallelism='mpi')
    chop_dir = args.chop_dir

    if args.rscc_file:
        with open(args.rscc_file, 'r') as j:
            rscc_list = json.load(j)
    else:
        rscc_list = None

    # Set chopping parameters
    parameters = {}
    if args.chop_par:
        for parameter in args.chop_par:
            pair = parameter.split('=')
            if pair[1].lower == 'true':
                pair[1] = True
            if pair[1].lower() == 'false':
                pair[1] = False
            if pair[0] not in ['chop_map', 'chopping_mod', 'check_rscc', 'chop_map']:
                pair[1] = float(pair[1])
            parameters[pair[0]] = pair[1]

        for key, value in default_parameters.items():
            if key not in parameters.keys():
                parameters[key] = value
    else:
        parameters = default_parameters

    if args.res_list:
        res_list = args.res_list
    else:
        res_list = nomenclature.AA_RESIDUES_LIST

    if args.np:
        np = args.np
    else:
        np = 2
    if args.res:
        args.res = float(args.res)

    model_name = os.path.basename(args.in_model).split('.')[0]
    map_name = os.path.basename(args.in_map).split('.')[0]
    name = str(model_name) + '-' + str(map_name)
    log_path = os.path.join(chop_dir, 'chop_' + name + '.log')

    chop.set_env(args.chop_dir, name, args.in_model, args.in_map, log_path,
                 rscc_list=rscc_list,
                 allowed_b=None,
                 warning_level='debug',
                 filter_identical_chains=args.filter)
    chop.set_map_chop_parameters(cube_radius=parameters["cube_radius"],
                                 final_voxel=parameters["final_voxel"],
                                 chop_radius=parameters["chop_radius"],
                                 chop_soft_radius=parameters["chop_soft_radius"],
                                 check_rscc=parameters["check_rscc"],
                                 )

    statistics = chop.chop_model_map_parallel(res_list, chop_map=parameters["chop_map"],
                                              chopping_mode='soft',
                                              nr_cores=size)

if __name__ == '__main__':
    main()