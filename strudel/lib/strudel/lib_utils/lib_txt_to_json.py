from strudel.classification.penultimateClassifier import PenultimateClassifier
import copy
import strudel.lib.strudel.nomenclature as nomenclature

lib = {}

for res_name in nomenclature.AA_RESIDUES_3LET_TO_FULL.keys():
    rotamers_data_structure = PenultimateClassifier.set_rotamers_data_structure()
    residue_3let = res_name
    rotamer = rotamers_data_structure[residue_3let]

    dih_angles_index = {}
    with open('penultimate_lib.txt', 'r') as rotamer_file:
        lines = rotamer_file.readlines()
    rr = {}
    rrs = []
    max_chi = 0
    rot_nr = 1
    key_words_length = 0
    for index, line in enumerate(lines):
        words = line.split()
        length = len(words)

        if dih_angles_index != {}:
            if length == key_words_length:
                rr['id'] = rot_nr
                rr['name'] = words[0]
                for key, value in dih_angles_index.items():
                    index1 = dih_angles_index[key]
                    w_key = 'chi-width_' + key.split('_')[1]
                    index2 = dih_angles_index['chi-width_' + key.split('_')[1]]

                    try:
                        rr[key] = int(words[index1])
                    except ValueError:
                        print(res_name, key, words[0])
                        rr[key] = '!!!!!!!wrong!!!!!!'
                    try:
                        rr[w_key] = int(words[index2]) * 2
                    except ValueError:
                        rr[w_key] = None
                tmp = copy.deepcopy(rr)
                rrs.append(tmp)
                rot_nr += 1
            else:
                break
        try:
            if words[0].lower() == nomenclature.AA_RESIDUES_3LET_TO_FULL[residue_3let].lower():
                key_words_length = len(words)
                name_line_words = lines[index - 1].split()
                for word_index, value in enumerate(name_line_words):
                    if 'chi_1' == value:
                        if 'chi_1' not in dih_angles_index.keys():
                            dih_angles_index['chi_1'] = word_index
                            max_chi = 1
                    if 'chi_2' == value:
                        dih_angles_index['chi_2'] = dih_angles_index['chi_1'] + 1
                        max_chi = 2
                    if 'chi_3' == value:
                        dih_angles_index['chi_3'] = dih_angles_index['chi_1'] + 2
                        max_chi = 3
                    if 'chi_4' == value:
                        dih_angles_index['chi_4'] = dih_angles_index['chi_1'] + 3
                        max_chi = 4
                    for i in range(max_chi):
                        key = 'chi-width_' + str(i + 1)
                        dih_angles_index[key] = dih_angles_index['chi_1'] + max_chi + i
        except IndexError:
            pass

    for r in rrs:
        print(r)
    print()

    lib[res_name] = rrs

import json
with open('penultimate_lib2.json', 'w') as j:
    json.dump(lib, j, indent=4, separators=(', ',': '))
