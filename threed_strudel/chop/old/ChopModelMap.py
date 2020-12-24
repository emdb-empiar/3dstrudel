"""
ChopModelMap.py

Chop an atomic residue into residues and residue side chains which are saved
as separate pdb files.
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

from .ChopMap import ChopMap
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBIO import PDBIO
from shutil import copy2
import copy
import os
import numpy as np
import mrcfile
from scipy.interpolate import RegularGridInterpolator
import logging
from datetime import datetime, timedelta
import time
import argparse
import json
from threed_strudel.core.Configure import Configure
import subprocess


class ChopModelMap:
    """
    Class for chopping residues and residue side chains out of proteins. The chopped
    side chains are used as guide for density map chopping
    """

    def __init__(self):
        self.config = Configure()
        self.cctbx_python = self.config.cctbx_python
        self.rscc_script_path = self.config.rscc_script_path
        self.min_rscc = 0.75
        self.chop_log = None
        self.allowed_b = 0.0
        self.entry_code = ''
        self.work_dir = ''
        self.input_model = None
        self.input_model_path = ''
        self.map_file_path = ''
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
        self.final_grid = None
        self.chop_radius = None
        self.chop_soft_radius = None
        self.script_path = os.path.dirname(os.path.abspath(__file__))
        self.charged_res = ['ARG', 'ASP', 'GLU', 'LYS']
        self.inclusion_level = None
        self.inclusion_fraction = None
        self.check_inclusion = None
        self.known_inclusion_levels = os.path.join(self.script_path, 'known_inclusion_levels.txt')
        self.inclusion_charged_tolerance = 3
        self.heavy_atoms = {'ARG': 11, 'LYS': 9, 'MET': 8, 'GLU': 9, 'GLN': 9, 'ASP': 8, 'ASN': 8,
                            'ILE': 8, 'LEU': 8, 'HIS': 10, 'TRP': 14, 'TYR': 12, 'PHE': 11, 'PRO': 7,
                            'THR': 7, 'VAL': 7, 'SER': 6, 'CYS': 6, 'ALA': 5, 'GLY': 4}

    def set_env(self, work_dir, entry_code, input_model, map_file_dir,
                log_path, rscc_list=None, resolution=None, inclusion_lvl=None,
                allowed_b=None, inclusion_fraction=None, warning_level='info'):
        """
        Set the environment for residue (and density map) chopping
        :param warning_level:
        :param inclusion_fraction:
        :param allowed_b:
        :param inclusion_lvl:
        :param resolution:
        :param rscc_list:
        :param log_path:
        :param work_dir: The directory for output files
        :param entry_code: EMD or pdb code of the entry
        :param input_model: Path to the residue pdb or cif file
        :param map_file_dir: Path to the residue map file
        :param warning_level: warning level
        """

        self.inclusion_level = inclusion_lvl
        self.inclusion_fraction = inclusion_fraction
        self.work_dir = work_dir
        self.create_dir(work_dir)
        self.entry_code = entry_code
        if input_model.split('.')[-1] == 'pdb':
            parser = PDBParser(PERMISSIVE=1)
            self.input_model = parser.get_structure(self.entry_code, input_model)
        elif input_model.split('.')[-1] == 'cif':
            parser = MMCIFParser()
            self.input_model = parser.get_structure(self.entry_code, input_model)
        else:
            raise Exception('Please provide a the input residue in PDB or CIF format')
        self.input_model_path = input_model
        self.map_file_path = map_file_dir
        # Setup log
        try:
            self.chop_log.removeHandler()
        except:
            pass
        self.chop_log = self.config.setup_logger('chop_log', log_path, warning_level=warning_level)

        if allowed_b is None:
            self.allowed_b = self.calc_b_cut_off()
            self.chop_log.info('No B-factor cut-off was specified\nThe cut-off value vas set to median + 2 * 1Qu: %s\n',
                               self.allowed_b)

        self.rscc_list = rscc_list
        if rscc_list is None:
            if resolution:
                self.rscc_list = self.calculate_rscc(resolution)
            else:
                raise Exception('Please provide the map resolution or a list of per residue RSCC')

        self.side_chain_end = '-{}_side.pdb'.format(self.entry_code)
        self.shift_side_end = '-{}_side_shift.pdb'.format(self.entry_code)

        self.residue_end = '-{}_residue.pdb'.format(self.entry_code)
        self.shift_residue_end = '-{}_residue_shift.pdb'.format(self.entry_code)
        self.cube_end = '-{}_cube.mrc'.format(self.entry_code)
        self.cube_res_end = '-{}_cube_resampled.mrc'.format(self.entry_code)
        self.soft_map_end = '-{}_soft.mrc'.format(self.entry_code)
        self.hard_map_end = '-{}_hard.mrc'.format(self.entry_code)

    def calculate_rscc(self, resolution):
        if resolution < 4:
            atom_radius = 2.0
        else:
            atom_radius = 2.5
        out_file = os.path.join(self.work_dir, self.entry_code + '_RSCC.json')
        rscc_parameters = ' resolution={} atom_radius={} scattering_table=electron'.format(resolution, atom_radius)
        command = '{} {} {} {} {} out_file={}'.format(self.cctbx_python, self.rscc_script_path,
                                                      self.input_model_path, self.map_file_path,
                                                      rscc_parameters, out_file)
        subprocess.call(command, shell=True, cwd=self.work_dir)
        try:
            with open(out_file, 'r') as j_file:
                rscc = json.load(j_file)
                return rscc
        except FileNotFoundError:
            self.chop_log.info("RSCC calculations failed")
            return []

    @staticmethod
    def find_residue_rscc(rscc_list, chain, res_type, res_nr):
        result = None
        for residue in rscc_list:
            if residue[0] == chain and residue[1] == res_type and int(residue[2]) == res_nr:
                result = residue[3]
                break
        return result

    def set_map_chop_parameters(self, cube_radius=5.0, final_grid=0.5, chop_radius=3.0, chop_soft_radius=2.0,
                                check_inclusion=False, check_rscc=True):
        """
        Set custom parameters for map chopping
        :param inclusion_fraction: fraction of atoms inside the threshold
        :param cube_radius: The distance between the cube edge the residue
        :param final_grid: Final voxel size of the chopped maps
        :param chop_radius: Radius around the molecule for map chopping
        :param chop_soft_radius: Soft hard_radius around the molecule for map chopping
        """
        self.cube_radius = cube_radius
        self.final_grid = final_grid
        self.chop_radius = chop_radius
        self.chop_soft_radius = chop_soft_radius
        self.check_inclusion = check_inclusion
        self.check_rscc = check_rscc
        if check_inclusion:
            if self.inclusion_level is None:
                lvl = self.check_known_inclusion_levels(self.inclusion_fraction)
                if lvl is not None:
                    self.inclusion_level = lvl
                else:
                    self.inclusion_level = self.find_threshold(self.inclusion_fraction)
                    self.save_inclusion_level(self.inclusion_fraction)

        date_time = datetime.now().strftime("%H:%M:%S %Y-%m-%d")
        text = '{:_^100}'.format('ChopModelMap') + '\n\nStarted: {}'.format(date_time)
        self.chop_log.info('%s\n\n', text)
        self.chop_log.info('Atom inclusion will be used for local map quality checks.\n'
                           'The threshold level is chosen such that %s%s of the atoms are inside the map\n'
                           'Threshold level = %s\n\n', self.inclusion_fraction, '%', self.inclusion_level)

    def check_known_inclusion_levels(self, inclusion_fraction):
        """
        Checks for known threshold levels
        :param inclusion_fraction: atom inclusion fraction
        :return: threshold level for the specified atom inclusion fraction
        """
        lvl = None
        with open(self.known_inclusion_levels, 'r') as f:
            lines = f.readlines()
            for line in lines:
                fields = line.split()
                if os.path.basename(self.map_file_path) in line and os.path.basename(self.input_model_path) in line \
                        and str(inclusion_fraction) in line:
                    lvl = float(fields[-1])
        return lvl

    def save_inclusion_level(self, inclusion_fraction):
        """
        Saves the calculated threshold level
        :param inclusion_fraction:
        """
        with open(self.known_inclusion_levels, 'a') as f:
            line = '{}  {}  {}  {}\n'.format(os.path.basename(self.map_file_path),
                                             os.path.basename(self.input_model_path),
                                             inclusion_fraction, self.inclusion_level)
            f.write(line)

    def interpolator(self):
        """
        Setup scipy regulargrid interpolation function
        :return: interpolation function
        """
        with mrcfile.open(self.map_file_path, mode='r+', permissive=True) as mrc:
            nx, ny, nz = mrc.data.shape
            x = range(nx)
            y = range(ny)
            z = range(nz)
            interpolator = RegularGridInterpolator((x, y, z), mrc.data)
        return interpolator

    def find_map_parameters(self):
        """
        Calculates map parameters which are used for atom inclusion
        :return: map parameters
        """
        with mrcfile.open(self.map_file_path, mode='r+', permissive=True) as mrc:
            x_voxel_size = mrc.header.cella.x / mrc.header.nx
            y_voxel_size = mrc.header.cella.y / mrc.header.ny
            z_voxel_size = mrc.header.cella.z / mrc.header.nz
            nxstart = mrc.header.nxstart
            nystart = mrc.header.nystart
            nzstart = mrc.header.nzstart
            map_max = mrc.data.max()
            map_min = mrc.data.min()
        return (x_voxel_size, y_voxel_size, z_voxel_size), (nxstart, nystart, nzstart), (map_min, map_max)

    def find_threshold(self, inclusion_fr, delta=2):
        """
        Calculates the threshold level for which the specified fraction of atoms is inside the map
        :param inclusion_fr: atom inclusion fraction (%)
        :param delta: atom inclusion fraction precision
        :return: threshold level
        """
        interpolator = self.interpolator()
        voxel_size, nstart, lvl_range = self.find_map_parameters()
        upper = lvl_range[1]
        lower = 0
        current_lvl = (upper - lower) / 2
        while True:
            included = 0
            not_included = 0
            for atom in self.input_model.get_atoms():
                atom_coord = atom.coord
                x_index = atom_coord[2] / voxel_size[0] - nstart[0]
                y_index = atom_coord[1] / voxel_size[1] - nstart[1]
                z_index = atom_coord[0] / voxel_size[2] - nstart[2]
                if interpolator([x_index, y_index, z_index]) > current_lvl:
                    included += 1
                else:
                    not_included += 1
            current_incl = included / (included + not_included) * 100
            if current_incl < inclusion_fr:
                upper = current_lvl
                current_lvl = current_lvl - (upper - lower) / 2
            elif current_incl > inclusion_fr + delta:
                lower = current_lvl
                current_lvl = current_lvl + (upper - lower) / 2
            else:
                final_lvl = round(current_lvl, 8)
                break
        return final_lvl

    def atom_inclusion(self, model):
        """
        Counts the number of atoms within and outside map at a given level
        :param model: Biopython structure object
        :param in_map: map file
        :param level: map level
        :return: [included(nr), not_included(nr)]
        """
        with mrcfile.open(self.map_file_path, mode='r+', permissive=True) as mrc:
            inclusion = []
            x_voxel_size = mrc.header.cella.x / mrc.header.nx
            y_voxel_size = mrc.header.cella.y / mrc.header.ny
            z_voxel_size = mrc.header.cella.z / mrc.header.nz
            x = range(mrc.header.nx)
            y = range(mrc.header.ny)
            z = range(mrc.header.nz)
            a = RegularGridInterpolator((x, y, z), mrc.data)

            for atom in model.get_atoms():
                atom_coord = atom.coord
                x_index = atom_coord[2] / x_voxel_size - mrc.header.nxstart
                y_index = atom_coord[1] / y_voxel_size - mrc.header.nystart
                z_index = atom_coord[0] / z_voxel_size - mrc.header.nzstart
                if a([x_index, y_index, z_index]) > self.inclusion_level:
                    inclusion.append(1)
                else:
                    inclusion.append(0)

        included = inclusion.count(1)
        not_included = inclusion.count(0)
        return included, not_included

    def get_residue_list(self, structure, res):
        """
        Create a list of BIO.PDB residue objects containing all the residues
        of the given type in the input residue
        :param structure: BIO.PDB structure object
        :param res: Residue code
        :return: A list of BIO.PDB residue objects
        """
        residues = []
        for residue in structure.get_residues():
            if residue.get_resname() == res.upper():
                single_conf = True
                for atom in residue:
                    try:
                        if atom.last_occupancy < 1:
                            single_conf = False
                    except AttributeError:
                        pass
                if single_conf:
                    residues.append(residue)
                else:
                    self.chop_log.info('The %s %s %s residue has multiple conformations. Will not be chopped.',
                                       res.upper(), residue.id[1], residue.parent.id)
        return residues

    def calc_b_cut_off(self):
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
        # min_b = 1000
        # max_b = 0
        for atom in residue:
            if atom.bfactor > self.allowed_b:
                valid = False
                break
        #    if atom.bfactor > max_b:
        #        max_b = atom.bfactor
        #    if atom.bfactor < min_b:
        #        min_b = atom.bfactor
        # if max_b * 0.33 > min_b:
        #    valid = False
        return valid

    def run_model_quality_checks(self, residue, chop_map):

        included = True
        good_b = True
        good_rscc = True
        complete = True
        chain_id = residue.parent.id
        res_nr = residue.id[1]
        res_name = residue.get_resname()

        # Check heavy atom completeness
        nr_atoms = 0
        for atom in residue:
            nr_atoms += 1
        if nr_atoms < self.heavy_atoms[res_name]:
            self.chop_log.info('Missing heavy atoms in %s %s %s (%s out of %s)', res_name,
                               res_nr, chain_id, nr_atoms, self.heavy_atoms[res_name])
            complete = False

        # Check local map quality based on B-factors and atom inclusion
        if chop_map:
            # Check atom inclusion
            if self.check_inclusion:
                not_included = self.atom_inclusion(residue)[1]
                if not_included > 0:
                    self.chop_log.info('Residue %s %s %s has %s atoms outside map at %s level', res_name,
                                       res_nr, chain_id, not_included, self.inclusion_level)
                else:
                    self.chop_log.debug('Residue %s %s %s has %s atoms outside map at %s level', res_name,
                                        res_nr, chain_id, not_included, self.inclusion_level)
                if res_name not in self.charged_res and not_included > 0:
                    included = False
                if res_name in self.charged_res and not_included > self.inclusion_charged_tolerance:
                    included = False
            # Check B-factors
            good_b = self.check_b_factors(residue)
            if not good_b:
                self.chop_log.info('Residue %s %s %s has high B-factors', res_name, res_nr, chain_id)
            # Check map residue real space correlation coefficients
            if self.check_rscc:
                rscc = self.find_residue_rscc(self.rscc_list, chain_id, res_name, res_nr)
                if rscc < self.min_rscc:
                    self.chop_log.info('Residue %s %s %s has RCSS lover than %s',
                                       res_name, res_nr, chain_id, self.min_rscc)
                    good_rscc = False

        if good_rscc and good_b and included and complete:
            return True
        else:
            return False

    def chop_model_map(self, residue_list, chop_map=True, chopping_mode='soft'):
        """
        Chop the specified residues and their side chains out of the input residue. Each residue and
        the corresponding side chain is saved in separate pdb files. The residues are classified as
        having high quality (highq_residue) or low quality (lowq_residue) local density map,
        superimposed and saved in multimodel pdb files.
        :param residue_list: The list of residue type to be chopped (ex. [ASP, TYR, ARG]).
        :param chop_map: Boolean. To perform map chopping or not.
        :param chopping_mode: 'soft' or 'hard'
        :param check_inclusion: Boolean. To check atom inclusion or not
        """
        for res in residue_list:
            self.chop_log.info('\nChopping: %s', res)
            res_dir = self.work_dir + '/' + res.lower()
            self.create_dir(res_dir)
            hq_residue_dir = res_dir + '/' + 'highq_residue'
            self.create_dir(hq_residue_dir)
            lq_residue_dir = res_dir + '/' + 'lowq_residue'
            self.create_dir(lq_residue_dir)
            side_chain_dir = res_dir + '/' + 'side_chain'
            self.create_dir(side_chain_dir)
            model_map_dir = res_dir + '/' + 'models_and_maps'
            self.create_dir(model_map_dir)
            lq_model_map_dir = res_dir + '/' + 'lq_models_and_maps'
            self.create_dir(lq_model_map_dir)

            # Create a list of all current type residue objects
            residues = self.get_residue_list(self.input_model, res)
            total_residues = len(residues)
            chopped_res = 0
            lq = 0
            hq = 0
            for residue in residues:
                chopped_res += 1
                chain_id = residue.parent.id
                res_nr = residue.id[1]
                # Run local residue quality checks
                good_model = self.run_model_quality_checks(residue, chop_map)
                if good_model:
                    highq = True
                    hq += 1
                else:
                    highq = False
                    lq += 1

                out_name = res.lower() + '-' + str(res_nr) + '-' + chain_id + self.residue_end
                if highq:
                    out_path = hq_residue_dir + '/' + out_name
                    self.save_pdb(residue, out_path)
                else:
                    out_path = lq_residue_dir + '/' + out_name
                    self.save_pdb(residue, out_path)

                # Delete the main chain atoms
                side_chain = copy.deepcopy(residue)
                side_chain = self.del_main_chain(side_chain)
                # Save the side chain
                prefix = res.lower() + '-' + str(res_nr) + '-' + chain_id
                out_name = prefix + self.side_chain_end
                out_path = side_chain_dir + '/' + out_name
                self.save_pdb(side_chain, out_path)

                if chop_map:
                    chop = ChopMap()
                    self.chop_log.info('Chopping %s out of %s residues', chopped_res, total_residues)
                    prefix = side_chain_dir + '/' + res.lower() + '-' + str(res_nr) + '-' + chain_id
                    cube_map_name = prefix + self.cube_end
                    matrix = chop.chop_cube(side_chain, self.map_file_path, cube_map_name, self.cube_radius)

                    shift_side = prefix + self.shift_side_end
                    chop.shift_coord(matrix, side_chain)
                    self.save_pdb(side_chain, shift_side)

                    cube_newgrid_name = prefix + self.cube_res_end
                    chop.grid_resample(cube_map_name, cube_newgrid_name, self.final_grid)
                    if chopping_mode.lower() == 'hard':
                        fin_map = prefix + self.hard_map_end
                        chop.chop_hard_radius(side_chain, cube_newgrid_name, fin_map, self.chop_radius)
                    elif chopping_mode.lower() == 'soft':
                        fin_map = prefix + self.soft_map_end
                        chop.chop_soft_radius(side_chain, cube_newgrid_name, fin_map, self.chop_radius,
                                              self.chop_soft_radius)
                    else:
                        raise Exception('The chopping_mode parameter can be hard or soft')

                    if highq:
                        copy2(fin_map, model_map_dir)
                    else:
                        copy2(fin_map, lq_model_map_dir)
                    if highq:
                        prefix = '{}/{}-{}-{}'.format(hq_residue_dir, res.lower(), res_nr, chain_id)
                        shift_residue = prefix + self.shift_residue_end
                        chop.shift_coord(matrix, residue)
                        self.save_pdb(residue, shift_residue)
                        copy2(shift_residue, model_map_dir)
                    else:
                        prefix = '{}/{}-{}-{}'.format(lq_residue_dir, res.lower(), res_nr, chain_id)
                        shift_residue = prefix + self.shift_residue_end
                        chop.shift_coord(matrix, residue)
                        self.save_pdb(residue, shift_residue)
                        copy2(shift_residue, lq_model_map_dir)

            self.chop_log.info('\n%s residues have high quality maps\n%s residues have low quality maps', hq, lq)
            # Superimpose the chopped residues
            self.superimpose(hq_residue_dir, lq_residue_dir)

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
        # for atom in residue:
        #    if atom.get_name() == "C":
        #        del residue[atom.get_name()]
        # for atom in residue:
        #    if atom.get_name() == "CA":
        #        del residue[atom.get_name()]
        for atom in residue:
            if atom.get_name() == "N":
                del residue[atom.get_name()]
        return residue

    def superimpose(self, residue_dir, lq_residue_dir):
        """
        Superimpose the residues with high and low quality local map.
        :param residue_dir: in_dir to high quality residues
        :param lq_residue_dir: in_dir to low quality residues
        """
        files = os.listdir(residue_dir)
        hq_pdb = [i for i in files if i.endswith('.pdb')]
        if len(hq_pdb) != 0:
            hq_pdb.sort(key=lambda x: (int(x.split('-')[1]), x.split('-')[2]))
            sup = residue_dir + '/' + 'super'
            self.create_dir(sup)
            for model in hq_pdb:
                aligned = self.align_n_ca_c(residue_dir + '/' + hq_pdb[0], residue_dir + '/' + model)
                self.save_pdb(aligned[0], sup + '/' + model.split('_')[0] + '.pdb')

            name = '{}_{}_highq.pdb'.format(hq_pdb[0].split('-')[0], self.entry_code)
            self.concatenate_pdb(sup, name)

        files = os.listdir(lq_residue_dir)
        lq_pdb = [i for i in files if i.endswith('.pdb')]
        if len(lq_pdb) != 0:
            lq_pdb.sort(key=lambda x: (int(x.split('-')[1]), x.split('-')[2]))
            sup = lq_residue_dir + '/' + 'super'
            self.create_dir(sup)
            if len(hq_pdb) != 0:
                for model in lq_pdb:
                    aligned = self.align_n_ca_c(residue_dir + '/' + hq_pdb[0], lq_residue_dir + '/' + model)
                    self.save_pdb(aligned[0], sup + '/' + model.split('_')[0] + '.pdb')
            else:
                for model in lq_pdb:
                    aligned = self.align_n_ca_c(lq_residue_dir + '/' + lq_pdb[0], lq_residue_dir + '/' + model)
                    self.save_pdb(aligned[0], sup + '/' + model.split('_')[0] + '.pdb')
            name = '{}_{}_highq.pdb'.format(lq_pdb[0].split('-')[0], self.entry_code)
            self.concatenate_pdb(sup, name)

    @staticmethod
    def concatenate_pdb(model_dir, out_file):
        """
        Merge all pdb files in a directory
        :param model_dir: A directory containing pdb files
        :param out_file: The name for the output multimodel pdb file
        """
        files = os.listdir(model_dir)
        pdb_files = [i for i in files if i.endswith('.pdb')]
        nr = 0
        lines = []
        with open(model_dir + '/' + out_file, 'w') as reference:
            for model in pdb_files:
                nr += 1
                with open(model_dir + '/' + model, 'r') as model_text:
                    lines.append('MODEL     ' + str(nr) + '\n')
                    for line in model_text:
                        if line != 'END\n':
                            lines.append(line)
                        else:
                            lines.append('ENDMDL\n')
            for line in lines:
                reference.write(line)
            reference.write('END\n')

    @staticmethod
    def align_n_ca_c(model_fixed, model_moving):
        """
        Superimpose 2 models containing single residues using N, CA, C atom selection
        :param model_fixed: pdb file
        :param model_moving: pdb file
        :return: Bio.PDB residue object, rotation-translation matrix
        """
        parser = PDBParser(PERMISSIVE=1)
        name = model_fixed.split('/')[-1].split('.')[0]
        str_fixed = parser.get_structure(name, model_fixed)
        name = model_moving.split('/')[-1].split('.')[0]
        str_moving = parser.get_structure(name, model_moving)
        fixed = []
        moving = []
        for atom in [a for a in str_fixed.get_atoms()]:
            if atom.get_id() == 'CA' or atom.get_id() == 'C' or atom.get_id() == 'N':
                fixed.append(atom)
        for atom in [a for a in str_moving.get_atoms()]:
            if atom.get_id() == 'CA' or atom.get_id() == 'C' or atom.get_id() == 'N':
                moving.append(atom)
        sup = Superimposer()
        sup.set_atoms(fixed, moving)
        sup.apply([a for a in str_moving.get_atoms()])
        return str_moving, sup.rotran

    @staticmethod
    def save_pdb(model, out_path):
        """
        Save residue as pdb
        :param model: structure object
        :param out_path: output pdb file in_dir
        """
        io = PDBIO()
        io.set_structure(model)
        io.save(out_path)

    @staticmethod
    def create_dir(path):
        """
        Create new directory
        :param path: new directory in_dir
        """
        if not os.path.exists(path):
            os.makedirs(path)

    @staticmethod
    def report_elapsed(start):
        """
        Converts seconds to d:h:m:s
        :param start: script starting time
        :return: report
        """
        end = time.time()
        elapsed = round(end - start, 2)
        sec = timedelta(seconds=int(elapsed))
        d = datetime(1, 1, 1) + sec
        formatted = "d:h:m:s   {}:{}:{}:{}".format(d.day - 1, d.hour, d.minute, d.second)
        text = '\nElapsed {} seconds.    {}'.format(elapsed, formatted)
        return text

    @staticmethod
    def read_json(json_file):
        with open(json_file, 'r') as j_file:
            content = json.load(j_file)
        return content


def main():
    parser = argparse.ArgumentParser(description='Create a Mask around the atomic residue')
    parser.add_argument("-j", "--jason", dest="jason_file", required=True, help="Jason file")
    args = parser.parse_args()

    chop = ChopModelMap()
    j_input = chop.read_json(args.jason_file)
    residues = j_input["residues"]
    input_data = j_input["input"]
    parameters = j_input["parameters"]
    out_dir = j_input["out_dir"]
    for record in input_data:
        start = time.time()
        job = '{}_{}_{}_{}'.format(record[0], parameters["chopping_mode"], 'grid', parameters["final_voxel"])
        work_dir = os.path.join(out_dir, job)
        log_path = os.path.join(work_dir, job + '.log')
        chop.set_env(work_dir, record[0], record[1], record[2], log_path)
        chop.set_map_chop_parameters(parameters["cube_radius"], parameters["final_voxel"], parameters["chop_radius"],
                                     parameters["chop_soft_radius"], parameters["inclusion_fraction"])
        chop.chop_model_map(residues, chop_map=parameters["chop_map"], chopping_mode=parameters["chopping_mode"])
        date_time = datetime.now().strftime("%H:%M:%S %Y-%m-%d")
        chop.chop_log.info('Finished: %s', date_time)
        chop.report_elapsed(start)


if __name__ == '__main__':
    main()
