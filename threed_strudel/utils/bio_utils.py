"""
bioUtils.py

A collection of functions for atomic models manipulations

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
__date__ = '2019-10-09'


import os
import numpy as np
import string
import copy
from itertools import combinations
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser, FastMMCIFParser
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.mmcifio import MMCIFIO
from scipy.cluster.hierarchy import linkage, fcluster, maxdists
from scipy.spatial import distance as ssd
from . import nomenclature as nm
from threed_strudel.utils import model_map_utils


def load_structure(model_path):
    """
    Creates Biopython structure object from file
    :param model_path: pdb, cif file in_dir
    :return: Biopython structure object
    """
    def copy_label_seq_id(in_path):
        parser = MMCIFParser(QUIET=True)
        parser._mmcif_dict = MMCIF2Dict(in_path)
        if "_atom_site.auth_seq_id" in parser._mmcif_dict:
            for i in parser._mmcif_dict["_atom_site.auth_seq_id"]:
                try:
                    int(i)
                except ValueError:
                    parser._mmcif_dict["_atom_site.auth_seq_id"] = parser._mmcif_dict["_atom_site.label_seq_id"]
        parser._build_structure(in_path.split('/')[-1].split('.')[0])
        return parser._structure_builder.get_structure()

    if model_path.split('.')[-1] == 'pdb' or model_path.split('.')[-1] == 'ent':
        parser = PDBParser(PERMISSIVE=1, QUIET=True)
        structure = parser.get_structure(model_path.split('/')[-1].split('.')[0], model_path)
    elif model_path.split('.')[-1] == 'cif':
        try:
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure(model_path.split('/')[-1].split('.')[0], model_path)
        except ValueError:
            structure = copy_label_seq_id(model_path)
    else:
        raise Exception('Please provide the input residue in pdb or cif format')

    return structure


def load_structure_label_id(model_path):
    """
    Creates Biopython structure object from file
    :param model_path: pdb, cif file in_dir
    :return: Biopython structure object
    """
    # def copy_label_seq_id(in_path):
    parser = MMCIFParser(QUIET=True)
    parser._mmcif_dict = MMCIF2Dict(model_path)
    for field in ['_atom_site.auth_seq_id', '_atom_site.auth_comp_id', '_atom_site.auth_asym_id', '_atom_site.auth_atom_id']:
        parser._mmcif_dict[field] = copy.deepcopy(parser._mmcif_dict[field.replace('auth', 'label')])

    while True:
        for i, rec in enumerate(parser._mmcif_dict['_atom_site.auth_seq_id']):
            try:
                int(rec)
            except ValueError:
                for key, value in parser._mmcif_dict.items():
                    if key.startswith('_atom_site'):
                        del parser._mmcif_dict[key][i]
                break
        else:
            break
    parser._build_structure(model_path.split('/')[-1].split('.')[0])
    structure = parser._structure_builder.get_structure()

    return structure


def fast_load_model(model_path):
    """
    Creates Biopython structure object from file
    :param model_path: pdb, cif file in_dir
    :return: Biopython structure object
    """
    if model_path.split('.')[-1] == 'pdb' or model_path.split('.')[-1] == 'ent':
        parser = PDBParser(PERMISSIVE=1)
        model = parser.get_structure(model_path.split('/')[-1].split('.')[0], model_path)
    elif model_path.split('.')[-1] == 'cif':
        parser = FastMMCIFParser()
        model = parser.get_structure(model_path.split('/')[-1].split('.')[0], model_path)
    else:
        raise Exception('Please provide the input residue in pdb or cif format')
    return model


def save_model(model, out_path, preserve_atom_numbering=False):
    """
    Save residue to file
    :param model: structure object
    :param out_path: output file in_dir
    """
    ext = out_path.split('.')[-1]
    if ext == 'pdb' or ext == 'ent':
        io = PDBIO()
    else:
        io = MMCIFIO()
    io.set_structure(model)
    io.save(out_path, preserve_atom_numbering=preserve_atom_numbering)


def residues2structure(residues):
    """
    Builds a biopython structure object from a residue or a list of residue objects
    :param residues: residue object/ list of residue objects
    :return: structure object
    """

    if type(residues) not in [list, tuple]:
        residues = [residues]
    for res in residues:
        if res.level != 'R':
            raise Exception(f'Residue object expected (R), {res.level} given instead!')

    sb = StructureBuilder()
    sb.init_structure('structure')
    sb.init_seg(' ')
    try:
        res_model_id = set([r.parent.parent.id for r in residues])
    except AttributeError:
        res_model_id = {0}
    res_chain_id = set([r.parent.id for r in residues])
    for m in res_model_id:
        sb.init_model(m)
    for c in res_chain_id:
        sb.init_chain(c)
    for r in residues:
        model_id = r.parent.parent.id
        chain_id = r.parent.id
        sb.structure[model_id][chain_id].add(r.copy())
    for atom in sb.structure.get_atoms():
        atom.disordered_flag = 0

    return sb.structure


def rename_chains(structure):
    """
    Renames the chains using ascii uppercase lowercase and digits in order to minimise the chains names length
    :param structure: Biopython structure object
    :return: Biopython structure object
    """
    ch_ids = []
    for model in structure:
        for chain in model:
            if chain.id not in ch_ids:
                ch_ids.append(chain.id)

    symbols = string.ascii_uppercase + string.ascii_lowercase + string.digits
    base_set = [i for i in symbols]
    ids = base_set
    comb_nr = 2
    while len(ch_ids) > len(ids):
        comb = [''.join(i) for i in combinations(base_set, comb_nr)]
        ids = ids + comb
        comb_nr += 1

    translation = {}
    for i, c in enumerate(ch_ids):
        translation[c] = ids[i]

    for model in structure:
        for chain in model:
            chain.id = translation[chain.id]
    return structure


def get_cif_loop_item_pointer(lines, category, item):
    """
    Finds the position of a data item in a list cif file lines.
    Example: get_cif_loop_item_pointer(lines, '_atom_site', 'pdbx_PDB_model_num')
    :param lines: cif file as list of lines
    :param category: Category name
    :param item: Item name
    :return: Field number (starting from 0),
             Start data line
             End data line
    """
    field = None
    start_line = None
    end_line = None
    for n, line in enumerate(lines):
        if line.rstrip() == 'loop_':
            if lines[n + 1].split('.')[0] == category:
                k = 1
                while True:
                    key_words = lines[n + k].rstrip().split('.')
                    if key_words[0] == category:
                        if key_words[-1] == item:
                            field = k - 1
                        k += 1
                    else:
                        start_line = n + k
                        break
                break

    for n, line in enumerate(lines[start_line:]):
        if line.rstrip() == '#':
            end_line = start_line + n - 1
            break

    return field, start_line, end_line


def make_int_auth_seq_id(model_path, out_path=None):
    """
    Makes auth_seq_id integer numbers
    :param model_path:
    :param out_path:
    :return:
    """
    model_dict = MMCIF2Dict(model_path)
    auth_ids = model_dict['_atom_site.auth_seq_id']
    non_int_entities = {}
    first_id = ''
    entity = ''
    # Find sequences containing non integer auth_seq_id
    for i, a_id in enumerate(auth_ids):
        next_entity = model_dict['_atom_site.label_asym_id'][i]
        if next_entity != entity:
            first_id = a_id
            entity = next_entity
        try:
            int(a_id)

        except ValueError:

            if entity not in non_int_entities.keys():
                non_int_entities[entity] = first_id

    # Rename auth_seq_id. Try to use the first encountered auth_seq_id as starting number

    for _ in range(len(non_int_entities)):
        conversion = {}
        proc_entity = None
        new_id = 1
        for i, a_id in enumerate(auth_ids):
            entity = model_dict['_atom_site.label_asym_id'][i]
            if entity in non_int_entities.keys():
                if proc_entity is None:
                    proc_entity = entity
                    try:
                        new_id = int(non_int_entities[entity])
                    except ValueError:
                        new_id = 1
                    conversion[a_id] = new_id
            if entity == proc_entity:
                if a_id not in conversion.keys():
                    new_id += 1
                    conversion[a_id] = new_id
                auth_ids[i] = str(conversion[a_id])
        del non_int_entities[proc_entity]

    io = MMCIFIO()
    io.set_dict(model_dict)
    if out_path is not None:
        io.save(out_path)
    else:
        io.save(model_path)


def classify_chains(structure, delta=0.1):
    """
    Classifies structure chains based on sequence_dir identity and rms
    :param delta: rms distance between classes
    :param structure: Biopython structure object
    :return: list of chain classes
    """
    # Find pairs of identical chains

    sup = Superimposer()
    edges = []
    if structure.get_level() == 'S':
        model = structure[0]
    elif structure.get_level() == 'M':
        model = structure
    else:
        return None, None

    chain_ids = [chain.get_id() for chain in model]
    if len(chain_ids) < 2:
        return [chain_ids], None
    distance_matrix = np.zeros(shape=(len(chain_ids), len(chain_ids)))

    for chains in combinations(model, 2):
        different = False
        # Exclude hetero residues
        residues1 = chains[0].get_list()
        residues1 = [res for res in residues1 if res.get_id()[0] == " "]
        residues2 = chains[1].get_list()
        residues2 = [res for res in residues2 if res.get_id()[0] == " "]
        if len(residues1) == len(residues2) and len(residues1) != 0:
            for i in range(len(residues1)):
                if residues1[i].get_resname() != residues2[i].get_resname():
                    different = True
                    break
            else:
                # Superimpose chains with identical sequence_dir
                moving = []
                fixed = []
                for res in residues1:
                    for atom in res.get_atoms():
                        moving.append(atom)
                for res in residues2:
                    for atom in res.get_atoms():
                        fixed.append(atom)
                if len(fixed) == len(moving):
                    sup.set_atoms(fixed, moving)
                    edges.append((chains[0].get_id(), chains[1].get_id(), sup.rms))
                    distance_matrix[chain_ids.index(chains[0].get_id()), chain_ids.index(chains[1].get_id())] = sup.rms
                    distance_matrix[chain_ids.index(chains[1].get_id()), chain_ids.index(chains[0].get_id())] = sup.rms
        else:
            different = True

        if different:
            distance_matrix[chain_ids.index(chains[0].get_id()), chain_ids.index(chains[1].get_id())] = delta * 2
            distance_matrix[chain_ids.index(chains[1].get_id()), chain_ids.index(chains[0].get_id())] = delta * 2

    # Cluster using the distance matrix
    # print(distance_matrix)
    Z = linkage(ssd.squareform(distance_matrix), method="complete")
    max_dist = maxdists(Z)
    clusters_ind = fcluster(Z, delta, criterion='distance')

    clusters = []
    for i in range(max(clusters_ind)):
        clusters.append([])
    for ind, c in enumerate(clusters_ind):
        clusters[c - 1].append(chain_ids[ind])



    sorted_clusters = []
    for chain in chain_ids:
        for cluster in clusters:
            if chain in cluster and cluster not in sorted_clusters:
                tmp_cluster = []
                for c in chain_ids:
                    if c in cluster:
                        tmp_cluster.append(c)
                sorted_clusters.append(tmp_cluster)

    return sorted_clusters, list(max_dist[-len(sorted_clusters)+1:])


def classify_chains_(structure, delta=0.1):
    """
    Classifies structure chains based on sequence_dir identity and rms
    :param delta: rms distance between classes
    :param structure: Biopython structure object
    :return: list of chain classes
    """
    # Find pairs of identical chains

    sup = Superimposer()
    # edges = []
    if structure.get_level() == 'S':
        model = structure[0]
    elif structure.get_level() == 'M':
        model = structure
    else:
        return None, None

    chain_ids = [chain.get_id() for chain in model]
    if len(chain_ids) < 2:
        return [chain_ids], None
    distance_matrix = np.zeros(shape=(len(chain_ids), len(chain_ids)))

    twin_classes = []
    has_twins = []

    for chains in combinations(model, 2):
        ch_id1 = chains[0].get_id()
        ch_id2 = chains[1].get_id()
        different = False
        # Exclude hetero residues
        residues1 = chains[0].get_list()
        residues1 = [res for res in residues1 if res.get_id()[0] == " "]
        residues2 = chains[1].get_list()
        residues2 = [res for res in residues2 if res.get_id()[0] == " "]
        if len(residues1) == len(residues2) and len(residues1) != 0:
            for i in range(len(residues1)):
                if residues1[i].get_resname() != residues2[i].get_resname():
                    different = True
                    break
            else:
                duplex = [ch_id1, ch_id2]
                print(duplex)
                print(has_twins)
                if ch_id1 not in has_twins and ch_id2 not in has_twins:
                    twin_classes.append(duplex)
                    has_twins.append(ch_id1)
                    has_twins.append(ch_id2)
                    print(has_twins)
                else:
                    for i, group in enumerate(twin_classes):
                        if ch_id1 in group or ch_id2 in group:
                            group.append(ch_id1)
                            group.append(ch_id2)
                print(twin_classes)

    print(twin_classes)
    for i, group in enumerate(twin_classes):
        twin_classes[i] = list(set(group))
    print(twin_classes)


    # for chains in combinations(residue, 2):
    #     different = False
    #     # Exclude hetero residues
    #     residues1 = chains[0].get_list()
    #     residues1 = [res for res in residues1 if res.get_id()[0] == " "]
    #     residues2 = chains[1].get_list()
    #     residues2 = [res for res in residues2 if res.get_id()[0] == " "]
    #     if len(residues1) == len(residues2) and len(residues1) != 0:
    #         for i in range(len(residues1)):
    #             if residues1[i].get_resname() != residues2[i].get_resname():
    #                 different = True
    #                 break
    #         else:
    #             # Superimpose chains with identical sequence_dir
    #             moving = []
    #             fixed = []
    #             for res in residues1:
    #                 for atom in res.get_atoms():
    #                     moving.append(atom)
    #             for res in residues2:
    #                 for atom in res.get_atoms():
    #                     fixed.append(atom)
    #             if len(fixed) == len(moving):
    #                 sup.set_atoms(fixed, moving)
    #                 # edges.append((chains[0].get_id(), chains[1].get_id(), sup.rms))
    #                 distance_matrix[chain_ids.index(chains[0].get_id()), chain_ids.index(chains[1].get_id())] = sup.rms
    #                 distance_matrix[chain_ids.index(chains[1].get_id()), chain_ids.index(chains[0].get_id())] = sup.rms
    #     else:
    #         different = True
    #
    #     if different:
    #         distance_matrix[chain_ids.index(chains[0].get_id()), chain_ids.index(chains[1].get_id())] = delta * 2
    #         distance_matrix[chain_ids.index(chains[1].get_id()), chain_ids.index(chains[0].get_id())] = delta * 2
    #
    # # Cluster using the distance matrix
    # # print(distance_matrix)
    # Z = linkage(ssd.squareform(distance_matrix), method="complete")
    # max_dist = maxdists(Z)
    # clusters_ind = fcluster(Z, delta, criterion='distance')
    #
    # clusters = []
    # for i in range(max(clusters_ind)):
    #     clusters.append([])
    # for ind, c in enumerate(clusters_ind):
    #     clusters[c - 1].append(chain_ids[ind])
    #
    #
    #
    # sorted_clusters = []
    # for chain in chain_ids:
    #     for cluster in clusters:
    #         if chain in cluster and cluster not in sorted_clusters:
    #             tmp_cluster = []
    #             for c in chain_ids:
    #                 if c in cluster:
    #                     tmp_cluster.append(c)
    #             sorted_clusters.append(tmp_cluster)
    #
    # return sorted_clusters, list(max_dist[-len(sorted_clusters)+1:])

def select_best_chains(structure, map_object, chain_classes, threshold=None):
    """
    Selects best fitted chain in each class
    :param structure: Biopython structure object
    :param map_object: map object created by mapParser
    :param chain_classes: a list of chain classes
    :param threshold: map threshold level
    :return: a list containing best fitted chain from each class
    """

    if threshold is None:
        threshold = model_map_utils.find_threshold(structure, map_object, inclusion_fraction=90)

    if structure.get_level() == 'S':
        model = structure[0]
    elif structure.get_level() == 'M':
        model = structure
    else:
        return None
    filtered = []
    for _class in chain_classes:
        if len(_class) == 1:
            filtered.append(_class[0])
        else:
            best_chain = None
            best_inclusion = -1
            for chain in _class:
                inclusion = model_map_utils.atom_inclusion(model[chain], map_object, threshold)[2]
                if inclusion > best_inclusion:
                    best_chain = chain
                    best_inclusion = inclusion
            filtered.append(best_chain)

    return filtered


def mmcif2pdb_biopython(cif_file, pdb_file=None):
    """
    Converts MMCIF int PDB file when possible
    :param cif_file: MMCIF file
    :param pdb_file: PDB file
    :return: PDB file in_dir
    """
    if pdb_file is None:
        head, tail = os.path.split(cif_file)
        pdb_file = os.path.join(head, tail.split('.')[0] + '.pdb')
    s = load_structure(cif_file)
    l = len([a for a in s.get_atoms()])
    if l > 99999:
        return None
    try:
        save_model(s, pdb_file)
        return pdb_file
    except TypeError:
        pass


def pdb2mmcif(pdb_file, cif_file=None):
    """
    Converts PDB to MMCIF file
    :param cif_file: MMCIF file
    :param pdb_file: PDB file
    :return: MMCIF file in_dir
    """
    if cif_file is None:
        head, tail = os.path.split(pdb_file)
        cif_file = os.path.join(head, tail.split('.')[0] + '.cif')
    s = load_structure(pdb_file)
    save_model(s, cif_file)

    return cif_file


def shift_coord(trans_matrix, structure):
    """
    Apply translation matrix to residue
    :param trans_matrix: translation matrix
    :param structure: biopython structure object
    """
    def shift():
        atom.get_coord()[0] = atom.get_coord()[0] - trans_matrix[0]
        atom.get_coord()[1] = atom.get_coord()[1] - trans_matrix[1]
        atom.get_coord()[2] = atom.get_coord()[2] - trans_matrix[2]
    try:
        for residue in structure.get_residues():
            for atom in residue:
                try:
                    disordered = atom.disordered_get_id_list()
                except AttributeError:
                    shift()
                else:
                    for id_ in disordered:
                        atom.disordered_select(id_)
                        shift()
    except AttributeError:
        for atom in structure:
            try:
                disordered = atom.disordered_get_id_list()
            except AttributeError:
                shift()
            else:
                for id_ in disordered:
                    atom.disordered_select(id_)
                    shift()

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
        if atom.get_name() == "N":
            del residue[atom.get_name()]
    for atom in residue:
        if atom.get_name() == "C":
            del residue[atom.get_name()]
    return residue

def calc_residue_static_pairwise_rmsd(model1, model2):
    """
    Calculate pairwise RMSD for residue structure objects without superimposing
    Takes into account residues symmetry
    :param model1: structure object
    :param model2: structure object
    :return: RMSD value
    """
    atoms1 = [atom for atom in model1.get_atoms() if not atom.get_name().upper().startswith('H')]
    atoms2 = [atom for atom in model2.get_atoms() if not atom.get_name().upper().startswith('H')]
    names1 = [atom.get_name() for atom in atoms1]
    names2 = [atom.get_name() for atom in atoms2]

    if len(atoms1) > len(atoms2):
        atoms1 = [atom for atom in atoms1 if atom.get_name() in names2]
    else:
        atoms2 = [atom for atom in atoms2 if atom.get_name() in names1]
    atoms1.sort(key=lambda x: x.get_name())
    atoms2.sort(key=lambda x: x.get_name())

    array1, array2 = np.array([atoms1[0].get_coord()]), np.array([atoms2[0].get_coord()])
    for i in range(1, len(atoms1)):
        array1 = np.concatenate((array1, np.array([atoms1[i].get_coord()])), axis=0)
        array2 = np.concatenate((array2, np.array([atoms2[i].get_coord()])), axis=0)
    rmsd_ = rmsd(array1, array2)

    res_type = model1.get_resname().lower()
    all_res_symm_pairs = nm.SYMMETRIC_ATOM_PAIRS
    if res_type not in all_res_symm_pairs.keys():
        return rmsd_

    symmetric_pairs = all_res_symm_pairs[res_type]
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
    rmsd2 = rmsd(array1, array3)
    return min(rmsd_, rmsd2)


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
