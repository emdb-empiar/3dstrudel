"""
rotamerClassifier.py

Splits residues of the same type into specified number of rotameric classes
based on RMSD.

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

from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Superimposer import Superimposer
from random import randint
import os
import copy


class RotamerClassifier:
    """
    Class for rotamer classification
    """
    class Atoms:
        res = {'phe': ['CB', 'CG', 'CD1', 'CD2'],
               'val': ['CB', 'CG1', 'CG2']

               }



    def __init__(self, path, rotamer_number):
        self.path = path
        self.map_file_format = 'map'
        self.str_object_list = []
        self.rotamer_number = rotamer_number
        self.create_pdb_list()


    def create_pdb_list(self):
        """
        Create a list of dictionaries with pdb file name as keys

        """
        files = os.listdir(self.path)
        pdb_names = [i for i in files if i.endswith('.pdb')]
        parser = PDBParser(PERMISSIVE=1)
        for element in pdb_names:
            name = element.split('.')[0]
            new_str = parser.get_structure(name, self.path + '/' + element)
            valid_str = True
            atoms = self.Atoms
            for residue in new_str.get_residues():
                type = residue.get_resname().lower()

            for atom in atoms.res[type]:
                if atom not in [a.get_id() for a in new_str.get_atoms()]:
                    valid_str = False
            if valid_str:
                self.str_object_list.append(new_str)
            else:
                print(new_str.id)

    def superimpose_n_ca_c(self, str_fixed, str_moving):
        """
        Superimpose 2 models containing single residues using N, CA, C atom selection
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

    def superimpose_all(self):
        """
        Superimpose all structures to the first one in the list
        """
        i = 1
        while i < len(self.str_object_list):
            self.str_object_list[i] = self.superimpose_n_ca_c(self.str_object_list[0], self.str_object_list[i])
            i += 1

    def calc_pairwise_rmsd(self, model1, model2):
        """
        Calculate pairwise RMSD for structure objects
        :param model1: structure object
        :param model2: structure object
        :return: RMSD value
        """
        summ = 0
        i = 0
        for atom1 in [a for a in model1.get_atoms()]:
            for atom2 in [a for a in model2.get_atoms()]:
                if atom1.get_id() == atom2.get_id():
                    d = atom1 - atom2
                    summ += d * d
                    i += 1
        rmsd = (summ / i) ** (0.5)
        return rmsd

    def save_pdb(self, model, out_path):
        """
        Save residue as pdb
        :param model: BIO.PDB structure object
        :param out_path: out file in_dir
        """
        io = PDBIO()
        io.set_structure(model)
        io.save(out_path)

    def concatenate_pdb(self, model_dir, out_file):
        """
        Merge all pdb files in a directory
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

    def split_rotamers(self):
        """
        Group the models into classes based on pairwise RMSD
        :return: a list of classes
        """
        # Find the initial reference models for each class
        classes = []
        for group in range(self.rotamer_number):
            classes.append([])
        classes[0].append(self.str_object_list[randint(0, len(self.str_object_list))])
        class_index = 1
        while class_index < self.rotamer_number:
            rmsd_list = []

            for structure in self.str_object_list:
                rmsd = 0
                for group in classes:
                    if group != []:
                        rmsd += self.calc_pairwise_rmsd(group[0], structure)
                rmsd_list.append(rmsd)
            max_rmsd = max(rmsd_list)
            classes[class_index].append(self.str_object_list[rmsd_list.index(max_rmsd)])
            class_index += 1
        # Add the remaining models to the corresponding classes
        i = 0
        while i < len(self.str_object_list):
            class_index = 0
            class_rmsd =[]
            while class_index < len(classes):
                rmsd = self.calc_pairwise_rmsd(classes[class_index][0], self.str_object_list[i])
                class_rmsd.append(rmsd)
                class_index += 1
            if 0.0 not in class_rmsd:
                min_rms = min(class_rmsd)
                classes[class_rmsd.index(min_rms)].append(self.str_object_list[i])
            i += 1
        # Find representative models within classes and repeat classification using
        # representative models as references. Execute while representative models change
        new_repr_set = []
        #
        while True:
            # Find representative models within classes
            class_index = 0
            while class_index < len(classes):
                rmsd_list = []
                length = len(classes[class_index])
                k = 0
                while k < length:
                    i = 0
                    rmsd = 0
                    while i < length:
                        rmsd += self.calc_pairwise_rmsd(classes[class_index][k], classes[class_index][i])
                        i += 1
                    rmsd_list.append(rmsd)
                    k += 1
                min_rms = min(rmsd_list)
                representative_model = classes[class_index][rmsd_list.index(min_rms)]
                classes[class_index] = []
                classes[class_index].append(representative_model)
                class_index += 1
            old_repr_set = copy.deepcopy(new_repr_set)
            new_repr_set = [x[0].id for x in classes]
            print(new_repr_set)
            # Repeat classification using representative models as references
            i = 0
            while i < len(self.str_object_list):
                class_index = 0
                class_rmsd =[]
                while class_index < len(classes):
                    rmsd = self.calc_pairwise_rmsd(classes[class_index][0], self.str_object_list[i])
                    class_rmsd.append(rmsd)
                    class_index += 1
                if 0.0 not in class_rmsd:
                    min_rms = min(class_rmsd)
                    classes[class_rmsd.index(min_rms)].append(self.str_object_list[i])
                i += 1
            if old_repr_set == new_repr_set:
                for cl in classes:
                    tmp_name = cl[0].id + '_representative'
                    cl[0].id = tmp_name
                break
        return classes


def main():
    pdb_dir = 'test_new/models_and_maps'
#    pdb_dir = 'x-ray/5o4i/tyr/residue'
    rotamers = 3
    cl = RotamerClassifier(pdb_dir, rotamers)
    cl.superimpose_all()
    classes = cl.split_rotamers()

    for i in range(rotamers):
        path = pdb_dir + '/' + 'rot_' + str(i + 1)
        if not os.path.exists(path):
            os.mkdir(path)
        classes[i].sort(key=lambda x: (int(x.id.split('-')[1]), x.id.split('-')[2]))
        for k in classes[i]:
            cl.save_pdb(k, path + '/' + k.id + '.pdb')
        cl.concatenate_pdb(path, 'class_' + str(i + 1) + '.pdb')


if __name__ == '__main__':
    main()

