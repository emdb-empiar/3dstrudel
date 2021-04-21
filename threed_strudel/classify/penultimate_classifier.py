"""
penultimateClassifier.py

Groups residues of the same type into Richardson's rotamer classes.
This code uses the Penultimate Rotamer library 'SC Lovell, JM Word,
JS Richardson and DC Richardson (2000) "The Penultimate Rotamer Library"
Proteins: Structure Function and Genetics 40: 389-408.'

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


import sys
import traceback
import math
import os
import copy
import argparse
import json
import time
import numpy as np
from datetime import datetime
from Bio.PDB import calc_dihedral
from Bio.PDB.Superimposer import Superimposer

import threed_strudel.utils.functions as func
import threed_strudel.nomenclature as nomenclature
from threed_strudel.utils import bio_utils
import threed_strudel.lib.strudel.penultimate_rotamer_lib as p_lib


class PenultimateClassifier:
    """
    Class for rotamer classification using a rotamer library in json format
    By default penultimate library (inp/penultimate_lib.json)
    Lovell, S.C., Word, J.M., Richardson, J.S. & Richardson, D.C. "The penultimate rotamer library",
    Proteins: Structure, Function and Genetics, Vol.40, 389-408 (2000)
    """

    def __init__(self, in_dir,
                 out_dir=None,
                 out_prefix='rot_',
                 default_angle_width=30,
                 use_minimal_width=False,
                 warning_level='info',
                 residue_type=None,
                 rotamer_lib_path=None,
                 log_path=None):
        """
        :param in_dir: Input directory containing residues of the same type as separate files
        :param out_dir: Output directory
        :param out_prefix: prefix for output rotamer folders
        :param rotamer_lib_path: penultimate library file (json format)
        :param default_angle_width: rotamer angle width to be used when missing in the library or
               when minimal width flag is True
        :param use_minimal_width: BOOLEAN
        :param warning_level: info, debug
        :param log_path: path for the process_log
        """
        self.in_dir = in_dir
        self.log_path = log_path
        self.corrupted_files_dir = os.path.join(in_dir, 'corrupted_files')
        if not os.path.exists(self.corrupted_files_dir):
            os.mkdir(self.corrupted_files_dir)
        self.model_err_out_path = self.setup_model_read_error_out()
        self.model_err_out = open(self.model_err_out_path, 'w')

        if out_dir is None:
            self.out_dir = in_dir
        else:
            self.out_dir = out_dir
        self.log = self.setup_log(warning_level)

        self.rotamer_lib = self.load_rotamer_lib(rotamer_lib_path)
        self.out_prefix = out_prefix
        # Set True to use the default width as minimal width value
        self.use_minimal_width = use_minimal_width
        # To be used when the width is not defined in the library and when use_minimal_width=True
        self.default_width = default_angle_width
        if residue_type is not None:
            self.residue_type = residue_type.lower()
        else:
            self.residue_type = residue_type
        self.rotamer_classes = None
        self.rotamers = None

        self.all_res_symm_pairs = nomenclature.SYMMETRIC_ATOM_PAIRS
        self.symmetric_chi = nomenclature.SYMMETRIC_CHI
        self.rotamers_data_structure = self.set_rotamers_data_structure()
        # self.run_classification()

    def run_classification(self):
        """

        :return:
        """
        start = time.time()
        date_time = datetime.now().strftime("%H:%M:%S %Y-%m-%d")
        text = '{:_^100}'.format('PenultimateClassifier') + '\n\nStarted: {}'.format(date_time)
        self.log.info('%s\n', text)

        struct_object_list = self.create_model_list()
        if self.residue_type is None:
            self.residue_type = self.get_res_type(struct_object_list[0])

        self.superimpose_all(struct_object_list)
        self.log.info('Classifying %s residues', self.residue_type.upper())

        if self.residue_type.upper() in ('ALA', 'GLY'):
            self.rotamer_classes = [struct_object_list, []]
            self.rotamer_classes[0][0].id = self.rotamer_classes[0][0].id + '_representative'
        else:
            rotamers = self.create_rotamers(self.residue_type)
            self.rotamers = rotamers
            self.rotamer_classes = self.find_rotamers(struct_object_list, rotamers)

            for class_ in self.rotamer_classes:
                if len(class_) > 200:
                    repr_index = self.find_repr_model_fast(class_, rotamers)
                else:
                    repr_index = self.find_repr_model(class_)
                if repr_index is not None:
                    class_[repr_index].id = class_[repr_index].id + '_representative'

        self.save_rotamers(self.rotamer_classes)
        text = func.report_elapsed(start)
        self.log.info(text)

    @staticmethod
    def load_rotamer_lib(rotamer_lib_path=None):
        """
        Reads the rotamer library in json format
        :param rotamer_lib_path: json file
        :return: file content
        """
        if rotamer_lib_path is not None:
            try:
                with open(rotamer_lib_path) as j:
                    data = json.load(j)
                    return data
            except:
                pass
        data = p_lib.PENULTIMATE_2000
        return data

    def setup_log(self, warning_level):
        """
        Setup process_log
        :param warning_level:
        :return:
        """
        class_log = func.setup_logger('process_log', self.log_path, warning_level=warning_level)
        return class_log

    def setup_model_read_error_out(self):
        if self.log_path is not None:
            err_out_file = os.path.splitext(self.log_path)[0] + '.err'
        else:
            err_out_file = os.path.join(self.corrupted_files_dir, 'corrupted_files.err')
        return err_out_file

    def create_model_list(self):
        """
        Reads all .cif files in a directory and creates a list of structure objects
        """
        if self.residue_type is not None:
            r_type = self.residue_type
            type_defined = True
        else:
            type_defined = False
            r_type = None

        files = os.listdir(self.in_dir)
        cif_names = [i for i in files if (i.endswith('.cif') or i.endswith('.pdb') or i.endswith('.ent'))
                     and not i.startswith('.')]

        str_object_list = []
        for element in cif_names:
            try:
                structure = bio_utils.load_structure(os.path.join(self.in_dir, element))
            except:
                # If the residue file could not be read for some reason if will move it
                self.log.info('Failed to process %s file\nMoving it to and related files %s\nCheck %s for details',
                              element, self.corrupted_files_dir, self.model_err_out_path)
                exc_type, exc_value, exc_traceback = sys.exc_info()
                traceback.print_tb(exc_traceback, file=self.model_err_out)
                # os.replace(os.path.join(self.in_dir, element), os.path.join(self.corrupted_files_dir, element))
                prefix = element.split('_')[0]
                related_files = [i for i in files if i.startswith(prefix)]
                if related_files:
                    for file in related_files:
                        os.replace(os.path.join(self.in_dir, file), os.path.join(self.corrupted_files_dir, file))
            else:
                tt = self.get_res_type(structure)
                if tt is not None:
                    if r_type is None:
                        r_type = tt
                    if tt == r_type:
                        str_object_list.append(structure)
                    elif not type_defined:
                        raise Exception('There are multiple residue types in the input and no '
                                        '"residue_type" parameter specified. '
                                        '\nPlease restart with defined "residue_type" parameter')

        if not str_object_list:
            raise Exception(f'There are no atomic residue files in the specified directory: {self.in_dir}')
        return str_object_list

    @staticmethod
    def set_rotamers_data_structure():
        """
        A dictionary to define the rotameric angles for standard amino acids.
        :return: dictionary
        """
        residue_rotamers = nomenclature.ROTAMER_DATA
        for key in residue_rotamers.keys():
            residue_rotamers[key]['id'] = None
            residue_rotamers[key]['name'] = None
        return residue_rotamers

    def create_rotamers(self, residue_3let):
        """
        Creates a list of dictionaries with relevant rotameric parameters for the specified residue type
        :param residue_3let: residue code (3 letter)
        :return: list of dictionaries
        """
        residue_3let = residue_3let.lower()
        rotamers = []

        residue_lib_data = self.rotamer_lib[residue_3let]
        for rotamer in residue_lib_data:
            rot_data_structure = copy.deepcopy(self.rotamers_data_structure[residue_3let])
            for key in rot_data_structure.keys():
                if key.startswith('chi_'):
                    rot_data_structure[key]['angle_val'] = rotamer[key]
                    nr = key.split('_')[-1]
                    try:
                        value = int(rotamer['chi-width_' + nr])
                    except TypeError:
                        value = self.default_width

                    if not self.use_minimal_width:
                        rot_data_structure[key]['angle_width'] = value
                    else:
                        rot_data_structure[key]['angle_width'] = max([value, self.default_width])
                else:
                    rot_data_structure[key] = rotamer[key]
            rotamers.append(rot_data_structure)

        return rotamers

    @staticmethod
    def get_res_type(model):
        """
        Reads the 3 letter residue code from a residue  BIO.PDB object
        :param model:
        :return: 3 letter residue code
        """
        residues = [residue for residue in model.get_residues()]
        if len(residues) == 1:
            return residues[0].get_resname().lower()
        else:
            return None

    @staticmethod
    def _get_dihedral_angle(residue, atom_names):

        if len(atom_names) < 4:
            raise Exception(f"Minimum 4 atom names needed, {len(atom_names)} were given")

        v = [None for _ in range(len(atom_names))]
        atoms = residue.get_atoms()
        for atom in atoms:
            for i, name in enumerate(atom_names):
                if atom.get_name().lower() == name:
                    v[i] = atom.get_vector()

        if any([a is None for a in v]):
            missing = []
            for i, v in enumerate(v):
                if v is None:
                    missing.append(atom_names[i])
            raise Exception(f"Missing atoms in the input residue: {' '.join(missing)}")

        angle = calc_dihedral(v[0], v[1], v[2], v[3])
        angle = math.degrees(angle)
        if len(atom_names) > 4:
            angle1 = calc_dihedral(v[0], v[1], v[2], v[4])
            angle1 = math.degrees(angle1)
            angles = [angle, angle1]
        else:
            angles = [angle]
        return angles

    @staticmethod
    def _get_dihedral_angle1(model, atom_names):
        vector1, vector2, vector3, vector4, vector5 = None, None, None, None, None
        atoms = model.get_atoms()
        for atom in atoms:
            if atom.get_name().lower() == atom_names[0]:
                vector1 = atom.get_vector()
            elif atom.get_name().lower() == atom_names[1]:
                vector2 = atom.get_vector()
            elif atom.get_name().lower() == atom_names[2]:
                vector3 = atom.get_vector()
            elif atom.get_name().lower() == atom_names[3]:
                vector4 = atom.get_vector()
            if len(atom_names) == 5:
                if atom.get_name().lower() == atom_names[4]:
                    vector5 = atom.get_vector()

        if any([a is None for a in [vector1, vector2, vector3, vector4]]):
            missing = []
            for i, v in enumerate([vector1, vector2, vector3, vector4]):
                if v is None:
                    missing.append(atom_names[i])
            raise Exception(f"Missing atoms in the input residue: {' '.join(missing)}")

        angle = calc_dihedral(vector1, vector2, vector3, vector4)
        angle = math.degrees(angle)
        if len(atom_names) == 5:
            angle1 = calc_dihedral(vector1, vector2, vector3, vector5)
            angle1 = math.degrees(angle1)
            angles = [angle, angle1]
        else:
            angles = [angle]
        return angles

    def find_rotamers(self, str_object_list, rotamers):
        """
        Group the input models into Richardson classes
        :return: a list of rotamer lists, the last element contains the non-rotameric models
        """
        classes = []
        for _ in rotamers:
            classes.append([])
        classes.append([])

        for model in str_object_list:
            sorted_ = False
            angles = {}
            for key, value in rotamers[0].items():
                if key not in ('name', 'id'):
                    angles[key] = self._get_dihedral_angle(model, value['atoms'])

            for rotamer in rotamers:
                belongs = True
                rotamer_nr = None
                for key, value in rotamer.items():
                    if key == 'id':
                        rotamer_nr = value
                    elif key == 'name':
                        pass
                    else:
                        min_ = value['angle_val'] - value['angle_width'] / 2
                        max_ = value['angle_val'] + value['angle_width'] / 2
                        if len(angles[key]) > 1:
                            inside1 = self.angle_in_interval(angles[key][0], min_, max_)
                            inside2 = self.angle_in_interval(angles[key][1], min_, max_)
                            inside = inside1 or inside2
                        else:
                            inside = self.angle_in_interval(angles[key][0], min_, max_)
                        if not inside:
                            belongs = False
                if belongs and not sorted_:
                    classes[rotamer_nr - 1].append(model)
                    sorted_ = True
            if not sorted_:
                classes[-1].append(model)
        return classes

    @staticmethod
    def angle_in_interval(n, a, b):
        """
        Checks if the angle n belongs to the interval [a, b]
        :param n: angle (deg)
        :param a: angle (deg)
        :param b: angle (deg)
        :return: bool
        """
        n = n % 360
        a = a % 360
        b = b % 360
        if a < b:
            return a <= n <= b
        return a <= n or n <= b

    def save_rotamers(self, classes):
        """
        Save each rotamer class and the non-rotameric models in separate folders
        :param classes: a list of rotamer lists, the last element contains the non-rotameric models
        """
        nr_rotamers = 0
        nr_non_rotamers = 0
        for index, element in enumerate(classes[0:-1]):
            if self.residue_type.upper() in ('ALA', 'GLY'):
                folder = '{}1_{}'.format(self.out_prefix, 'single-conf')
            else:
                folder = '{}{}_{}'.format(self.out_prefix, index + 1, self.rotamers[index]['name'])
            path = os.path.join(self.out_dir, folder)
            if not os.path.exists(path):
                os.makedirs(path)
            for model in classes[index]:
                model_path = os.path.join(path, model.id + '.cif')
                bio_utils.save_model(model, model_path)
                nr_rotamers += 1
            self.concatenate_residue_models(classes[index], os.path.join(path, 'class_' + str(index + 1) + '.cif'))

        path = os.path.join(self.out_dir, 'non-rotameric')
        if not os.path.exists(path):
            os.makedirs(path)
        for model in classes[-1]:
            model_path = os.path.join(path, model.id + '.cif')
            bio_utils.save_model(model, model_path)
            nr_non_rotamers += 1
        self.concatenate_residue_models(classes[-1], os.path.join(path, self.residue_type + '_non-rotameric.cif'))
        self.log.info("%s %s residues were classified as rotameric\n"
                      "%s %s residues were classified as non rotameric",
                      nr_rotamers, self.residue_type.upper(),
                      nr_non_rotamers, self.residue_type.upper())

    @staticmethod
    def superimpose_n_ca_c(str_fixed, str_moving):
        """
        Superimpose 2 models containing single residues using N, CA, C atom selection
        :param str_fixed: BIO.PDB structure object
        :param str_moving: BIO.PDB structure object
        """
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
        return str_moving

    def superimpose_all(self, str_object_list):
        """
        Superimpose all structures to the first one in the list
        """
        i = 1
        while i < len(str_object_list):
            str_object_list[i] = self.superimpose_n_ca_c(str_object_list[0], str_object_list[i])
            i += 1

    @staticmethod
    def concatenate_residue_models(model_obj_list, out_file):
        """
        Merge a list of structure objects files into a multimodel structure object
        :param model_obj_list: residue directory
        :param out_file: output file name
        """
        if len(model_obj_list) == 1:
            bio_utils.save_model(model_obj_list[0], out_file)
        elif len(model_obj_list) > 1:
            out_structure = copy.deepcopy(model_obj_list[0])
            out_structure[0].serial_num = 1
            for i, structure in enumerate(model_obj_list[1:]):
                structure[0].id = i + 1
                out_structure.add(structure[i + 1])
                out_structure[i + 1].serial_num = i + 2
            for model in out_structure:
                for chain in model:
                    chain.id = 'A'
            bio_utils.save_model(out_structure, out_file)

    def set_default_rotamer_angle_width(self, width):
        """
        Set the default dihedral angle width for rotamer parameters setting
        :param width: dihedral angle width
        """
        self.default_width = width

    def find_repr_model(self, model_list):
        """
        Find the representative models of the rotamer classes
        :param model_list: list of residue objects
        """
        self.log.debug("Searching for the representative residue in each rotamer class")
        if len(model_list) == 0:
            return None
        elif len(model_list) < 3:
            return 0
        rmsd_list = []
        for model1 in model_list:
            rmsd = 0
            for model2 in model_list:
                if model1 != model2:
                    rmsd += self.calc_residue_static_pairwise_rmsd(model1, model2)
            rmsd_list.append(rmsd)
        min_rms = min(rmsd_list)

        return rmsd_list.index(min_rms)

    def find_repr_model_fast(self, model_list, rotamers):
        """
        Find the representative models of the rotamer classes.
        The search is done in two steps:
        First based on the distance between the rotameric angles and the averaged rotameric angles.
        Second, 10% best candidates are evaluated using RMSD
        :param model_list: List of residue objects
        :param rotamers: rotamers dictionary
        :return representative residue index
        """
        self.log.debug("Searching for the representative residue in each rotamer class")
        if len(model_list) == 0:
            return None
        # First stage (find candidates)
        all_angles = {}
        for model in model_list:
            for key, value in rotamers[0].items():
                if key not in ('name', 'id'):
                    angles = self._get_dihedral_angle(model, value['atoms'])
                    angle = min([a % 360 for a in angles])
                    if key in all_angles.keys():
                        all_angles[key].append(angle)
                    else:
                        all_angles[key] = [angle]

        average_angles = {}
        for key in all_angles.keys():
            average_angles[key] = sum(all_angles[key]) / len(all_angles[key])

        keys = all_angles.keys()
        weights = []
        for i in range(len(model_list)):
            weight = 0
            for key in keys:
                weight += math.fabs(average_angles[key] - all_angles[key][i])
            weights.append(weight)
        min_indices = []
        max_weight = max(weights)
        if len(weights) > 10:
            ln = int(len(weights)*0.1)
        else:
            ln = len(weights)
        for i in range(ln):
            index = weights.index(min(weights))
            min_indices.append(index)
            weights[index] = max_weight + 1
        # Second stage accurately evaluate the candidates and find the best one
        rmsds = []
        for i in min_indices:
            ref_model = model_list[i]
            rmsd = 0
            for model in model_list:
                rmsd += self.calc_residue_static_pairwise_rmsd(model, ref_model)
            rmsds.append(rmsd)

        return min_indices[rmsds.index(min(rmsds))]

    def calc_residue_static_pairwise_rmsd(self, model1, model2):
        """
        Calculate pairwise RMSD for residue structure objects without superimposing
        Takes into account residues symmetry
        :param model1: structure object
        :param model2: structure object
        :return: RMSD value
        """
        atoms1 = [atom for atom in model1.get_atoms() if not atom.get_name().upper().startswith('H')]
        atoms2 = [atom for atom in model2.get_atoms() if not atom.get_name().upper().startswith('H')]
        if len(atoms1) != len(atoms2):
            self.log.info("Model %s (%s) and %s (%s) have different number of hevy atoms!!!\n"
                          "The calculated rms might be not accurate",
                          model1.id, len(atoms1), model2.id, len(atoms2))
        names1 = set([atom.get_name() for atom in atoms1])
        names2 = set([atom.get_name() for atom in atoms2])
        names = names1.intersection(names2)
        atoms1 = [atom for atom in atoms1 if atom.get_name() in names]
        atoms2 = [atom for atom in atoms2 if atom.get_name() in names]

        atoms1.sort(key=lambda x: x.get_name())
        atoms2.sort(key=lambda x: x.get_name())

        array1, array2 = np.array([atoms1[0].get_coord()]), np.array([atoms2[0].get_coord()])
        for i in range(1, len(atoms1)):
            array1 = np.concatenate((array1, np.array([atoms1[i].get_coord()])), axis=0)
            array2 = np.concatenate((array2, np.array([atoms2[i].get_coord()])), axis=0)
        rmsd = self.rmsd(array1, array2)

        res_type = self.get_res_type(model1).lower()

        if res_type not in self.all_res_symm_pairs.keys():
            return rmsd

        symmetric_pairs = self.all_res_symm_pairs[res_type]
        swap_indices = [[],[]]
        for i in range(0, len(symmetric_pairs), 2):
            for k, atom in enumerate(atoms2):
                if atom.get_name().lower() == symmetric_pairs[i]:
                    swap_indices[0].append(k)
                elif atom.get_name().lower() == symmetric_pairs[i+1]:
                    swap_indices[1].append(k)
        for i in range(len(swap_indices[0])):
            atoms2[swap_indices[0][i]], atoms2[swap_indices[1][i]] = \
                atoms2[swap_indices[1][i]], atoms2[swap_indices[0][i]]

        array3 = np.array([atoms2[0].get_coord()])
        for i in range(1, len(atoms2)):
            array3 = np.concatenate((array3, np.array([atoms2[i].get_coord()])), axis=0)
        rmsd2 = self.rmsd(array1, array3)
        return min(rmsd, rmsd2)

    @staticmethod
    def rmsd(array1, array2):
        """
        Calculates RMSD between two sets of coordinates without alignment
        :param array1: numpy 2d array
        :param array2: numpy 2d array
        :return: rmsd
        """
        if np.shape(array1) != np.shape(array1):
            raise Exception('The coordinate arrays must have the same dimensions')
        dif = array1 - array2
        return np.sqrt(np.mean(np.sum(dif * dif, axis=1)))

    def report_statistics(self):
        """
        Returns a dictionary with statistics on the number of rotameric and non rotameric residues
        Ex: {'residue': 'PHE', 'rotamers': [4, 4, 5, 7], 'non_rotamers': 10}
        :return: dictionary
        """
        statistics = {'residue': self.residue_type.upper(), 'rotamers': [], 'non_rotamers': None}
        for cl in self.rotamer_classes[:-1]:
            statistics['rotamers'].append(len(cl))
        statistics['non_rotamers'] = len(self.rotamer_classes[-1])
        return statistics


def main():
    parser = argparse.ArgumentParser(description='Classify rotamers based on Penultimate rotamer library')
    parser.add_argument("-i", "--inp_directory", dest="inp_dir", required=True,
                        help="Directory with residues as separate files")
    parser.add_argument("-o", "--out_directory", dest="out_dir", required=True,
                        help="Output directory. If not given than the output will be saved to the input directory")
    parser.add_argument("-p", "--out_prefix", dest="prefix", default='rot_', help="Rotamer folders prefix")
    args = parser.parse_args()
    classify = PenultimateClassifier(args.inp_dir, args.out_dir, args.prefix)
    classify.run_classification()
    # print(classify.report_statistics())
    classify.log.info('%s\n', classify.report_statistics())


if __name__ == '__main__':
    main()
