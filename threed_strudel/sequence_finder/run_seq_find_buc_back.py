from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser, FastMMCIFParser
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.Polypeptide import PPBuilder
import sys
import json
import copy
import os
import argparse
import subprocess

AA_3LET_TO_1LET = {'ALA': 'A',
                   'ARG': 'R',
                   'ASN': 'N',
                   'ASP': 'D',
                   'CYS': 'C',
                   'GLU': 'E',
                   'GLN': 'Q',
                   'GLY': 'G',
                   'HIS': 'H',
                   'ILE': 'I',
                   'LEU': 'L',
                   'LYS': 'K',
                   'MET': 'M',
                   'PHE': 'F',
                   'PRO': 'P',
                   'SER': 'S',
                   'THR': 'T',
                   'TRP': 'W',
                   'TYR': 'Y',
                   'VAL': 'V'}


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


def make_poly_ala(model):
    del_atoms = []
    bb = ['CA', 'C', 'N', 'O', 'CB']
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.get_name().upper() not in bb:
                    del_atoms.append((chain.get_id(), residue.get_id(), atom.get_id()))
    for chain, resi, atom in del_atoms:
        del model[chain][resi][atom]

    for chain in model:
        for residue in chain:
            residue.resname = 'ALA'
    return model


def get_fragments(model, length):
    fragments = []
    all_res = [a for a in model.get_residues()]
    ch = all_res[0].get_parent().get_id()
    l = len(all_res)
    for i in range(l-length):
        tmp_model = copy.deepcopy(model)
        del_ids = []
        for res in all_res:
            if res not in all_res[i:i+length]:
                del_ids.append(res.get_id())
        for _id in del_ids:
            del tmp_model[ch][_id]
        
        fragments.append(tmp_model)
    return fragments


def get_sequence(model):
    code = ''
    for residue in model.get_residues():
        r_c = residue.resname
        code += AA_3LET_TO_1LET[r_c.upper()]
    return code


def save_fragments(fragments_list, frag_path):
    for m in fragments_list:
        seq = get_sequence(m)
        rs = [r for r in m.get_residues()]
        st = rs[0].get_id()[1]
        end = rs[-1].get_id()[1]
        o_p = os.path.join(frag_path, f'{seq}_{st}-{end}.pdb')
        m = make_poly_ala(m)
        save_model(m, o_p)


def run_buc(buc_path, pdb_path, mtzfile, out_path, seq_path):

    command = f'''
{buc_path} -stdin << eof
seqin {seq_path}
mtzin {mtzfile}
colin-fo Fout0,SIGFem
colin-free FreeR_flag
colin-phifom Pout0,FOMem
pdbin {pdb_path}
cycles 1
sequence
eof
'''

    log = open(out_path, 'a')
    log.writelines(command)
    log.flush()
    proc = subprocess.Popen(command,
                            universal_newlines=True, stdout=log, shell=True)
    sys.stdout.flush()
    proc.wait()
    log.flush()
    log.close()


def read_buc_out(log_file, seq):
    with open(log_file, 'r') as f:
        lines = f.readlines()
        place, z = -1, 0
        for line in lines:
            words = line.split()
            if seq in words:
                place = int(words[0])
                z = float(words[1])
                break
    return place, z


def main():
    parser = argparse.ArgumentParser(description='Run Buccaneer')
    parser.add_argument("-buc", "--buccaneer", dest="buc_path", required=False, help="Buccaneer path")
    parser.add_argument("-f", "--fragments_dir", dest="f_dir", required=True, help="PDB fragments directory")
    parser.add_argument("-t", "--task", dest="task", required=True, help="Task to run. F - generate fragments"
                                                                          "R - running buccaneer")
    parser.add_argument("-j", "--json_out", dest="j_out", required=False, help="Json output file")
    parser.add_argument("-o", "--buc_out", dest="buc_out", required=False, help="Buccaneer output")
    parser.add_argument("-p", "--model", dest="m_path", required=True, help="Model path")
    parser.add_argument("-m", "--mtz", dest="mtz_path", required=False, help="Mtz path")
    parser.add_argument("-s", "--seq", dest="seq_path", required=False, help="Sequence path")
    parser.add_argument("-len", "--len", dest="len", required=False, type=int, help="Fragments length")

    args = parser.parse_args()

    for p in [args.f_dir, args.buc_out]:
        if p is not None:
            if not os.path.exists(p):
                os.makedirs(p)

    if args.task == 'F':
        model = load_structure(args.m_path)[0]
        frags = get_fragments(model, args.len)
        save_fragments(frags, args.f_dir)

    if args.task == 'R':
        frag_files = os.listdir(args.f_dir)
        pdb_files = [f for f in frag_files if f.endswith('.pdb') and not f.startswith('.')]
        if len(pdb_files) == 0:
            print('No PDB files found')
            sys.exit(0)
        else:
            results_list = []
            for pdb in pdb_files:
                seq = pdb.split('_')[0]
                lod_file = pdb.split('.')[0] + '.txt'
                log_path = os.path.join(args.buc_out, lod_file)
                pdb_path = os.path.join(args.f_dir, pdb)
                run_buc(args.buc_path, pdb_path, args.mtz_path, log_path, args.seq_path)
                place, z = read_buc_out(log_path, seq)
                results = {'place': place,
                           'z': z,
                           'seq': seq}
                results_list.append(results)
            
                if args.j_out:
                    with open(args.j_out, 'w') as j:
                        json.dump(results_list, j, indent=4)
                else:
                    print(f'{place} {z} {seq}')
                

if __name__ == '__main__':
    main()
