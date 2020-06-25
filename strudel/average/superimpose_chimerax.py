from chimerax.core.commands import Command
from shutil import copy2
import json
import os

rc = Command(session)

align_atoms_dict = {
    'symmetric': {'HIS': [['n', 'c', 'ca', 'cb', 'cg', 'nd1', 'ce1', 'cd2', 'ne2'],
                          ['n', 'c', 'ca', 'cb', 'cg', 'cd2', 'ne2', 'nd1', 'ce1']],
                  'PHE': [['n', 'c', 'ca', 'cb', 'cg', 'cd1', 'ce1', 'cd2', 'ce2', 'cz'],
                          ['n', 'c', 'ca', 'cb', 'cg', 'cd2', 'ce2', 'cd1', 'ce1', 'cz']],
                  'TYR': [['n', 'c', 'ca', 'cb', 'cg', 'cd1', 'ce1', 'cd2', 'ce2', 'cz'],
                          ['n', 'c', 'ca', 'cb', 'cg', 'cd2', 'ce2', 'cd1', 'ce1', 'cz']]
                  },
    'asymmetric': {'ASP': ['n', 'c', 'ca', 'cb', 'cg'],
                   'GLU': ['n', 'c', 'ca', 'cb', 'cg', 'cd'],
                   'ASN': ['n', 'c', 'ca', 'cb', 'cg'],
                   'ARG': ['n', 'c', 'ca', 'cb', 'cg', 'cd'],
                   'LYS': ['n', 'c', 'ca', 'cb', 'cg', 'cd'],
                   'MET': ['n', 'c', 'ca', 'cb', 'cg', 'sd'],
                   'GLN': ['n', 'c', 'ca', 'cb', 'cg', 'cd'],
                   'ILE': ['n', 'c', 'ca', 'cb', 'cg1', 'cd1'],
                   'LEU': ['n', 'c', 'ca', 'cb', 'cg', 'cd1'],
                   'TRP': ['n', 'c', 'ca', 'cb', 'cg', 'cd1'],
                   'PRO': ['n', 'c', 'ca', 'cb', 'cg'],
                   'THR': ['n', 'c', 'ca', 'cb', 'cg2'],
                   'VAL': ['n', 'c', 'ca', 'cb', 'cg1', 'cg2'],
                   'SER': ['n', 'c', 'ca', 'cb', 'og'],
                   'CYS': ['n', 'c', 'ca', 'cb', 'sg'],
                   'GLY': ['n', 'c', 'ca'],
                   'ALA': ['n', 'c', 'ca', 'cb']
                   }
                    }

with open('class_dict.json') as j_file:
    class_dict = json.load(j_file)

global_ref_path = class_dict['global_reference']
grid = class_dict['grid_path']
fit_map = class_dict['fit_map']
inp_dir = class_dict['inp_dir']
out_dir = class_dict['out_dir']
reference_map_path = class_dict['reference_mrc']
reference_model_path = class_dict['reference_model']

# Align the class reference residue and map to the global reference residue
rc.run('open ' + global_ref_path)
rc.run('open ' + reference_map_path)
model = rc.run('open ' + reference_model_path)[0][0]
rc.run('open ' + grid)
rc.run('align #3@ca,n,c to #1@ca,n,c')
rc.run('view position #2 sameAsModels #3')
rc.run('vop  resample #2 onGrid #4')

reference_map_path = os.path.join(out_dir, 'new_ref.mrc')
reference_model_path = os.path.join(out_dir, 'new_ref.cif')
rc.run(f'save {reference_map_path} models #5')
rc.run(f'save {reference_model_path} format mmcif models #3')

res_type = model.residues[0].name

rc.run('close all')

# Align all models and maps to the reference residue [and map]
for pair in class_dict['class_pairs'][:]:
    rc.run('open ' + reference_map_path)
    rc.run('open ' + reference_model_path)
    mrc_file = os.path.join(inp_dir, pair[0])
    pdb_file = os.path.join(inp_dir, pair[1])
    rc.run('open ' + mrc_file)
    rc.run('open ' + pdb_file)
    if res_type.upper() in align_atoms_dict['asymmetric'].keys():
        atoms_list = align_atoms_dict['asymmetric'][res_type.upper()]
        align_atoms = ','.join(atoms_list)
        rc.run(f'align #4@{align_atoms} to #2@{align_atoms}')
    else:
        atoms_list1 = align_atoms_dict['symmetric'][res_type.upper()][0]
        align_atoms1 = ','.join(atoms_list1)
        rms1 = rc.run(f'align #4@{align_atoms1} to #2@{align_atoms1}')[0][2]
        atoms_list2 = align_atoms_dict['symmetric'][res_type.upper()][1]
        align_atoms2 = ','.join(atoms_list2)
        rms2 = rc.run(f'align #4@{align_atoms2} to #2@{align_atoms2}')[0][2]
        if rms1 < rms2:
            rc.run(f'align #4@{align_atoms1} to #2@{align_atoms1}')

    rc.run('view position #3 sameAsModels #4')
    if fit_map:
        rc.run('volume #3 sdLevel 4')
        rc.run('fitmap #3 inMap #1 metric correlation')
    else:
        out_model = os.path.join(out_dir, pair[1])
        rc.run('save {} format mmcif models #4'.format(out_model))
    rc.run('vop resample  #3 onGrid #1')
    out_map = os.path.join(out_dir, pair[0])
    rc.run('save {} models #5'.format(out_map))
    # rc.run('close #3 #4 #5')
    # If just close the maps Chimera keeps the fitmap results in the memory (~1.GB per 10000 fits)
    # Close session deletes the results
    rc.run('close session')

rc.run('close all')
# name = class_dict['reference_model'].split('_')[0] + '_' + class_dict['class'] + '_representative.cif'
# reference_out = os.path.join(inp_dir, class_dict['class'], name)
# copy2(reference_model_path, reference_out)
rc.run('exit')
