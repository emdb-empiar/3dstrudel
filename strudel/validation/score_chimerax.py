__author__ = 'Andrei Istrate'
__email__ = 'andrei@ebi.ac.uk'
__date__ = '2019-04-10'


import os
import csv
import sys
from multiprocessing import Process, Manager, Lock
from datetime import datetime
import multiprocessing
from shutil import copy
import json
import argparse
import operator
import psutil
# from chimera import runCommand as rc
# from chimera import openModels as om, selection
# from FitMap.fitcmd import fitmap
from chimerax.core.commands import Command
import logging

rc = Command(session)
SD_LEVEL = 3


class DictKeys:
    RES_TYPE = 'residue_type'
    RES_NR = 'residue_nr'
    CHAIN = 'chain'
    TOP_TYPE = 'top_residue_type'
    TOP_CC = 'top_cc'
    RES_MAP = 'residue_map'
    RES_MODEL = 'residue_model'
    SCORES = 'scores'


keys = DictKeys()

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


def run_x(command):
    try:
        return rc.run(command)[0][0]
    except TypeError:
        return None


def read_motif_lib(motif_lib_dir):
    lib_pairs = []
    files = os.listdir(motif_lib_dir)
    pdb_names = [i for i in files if i.endswith('.cif') and not i.startswith('.')]
    for name in pdb_names:
        prefix = name.split('.cif')[0]
        # res = name.split('-')[0].upper()
        map_path = os.path.join(motif_lib_dir, prefix + '.map')
        mrc_path = os.path.join(motif_lib_dir, prefix + '.mrc')
        if os.path.exists(map_path):
            lib_pairs.append([prefix, os.path.join(motif_lib_dir, name), map_path])
        elif os.path.exists(mrc_path):
            lib_pairs.append([prefix, os.path.join(motif_lib_dir, name), mrc_path])
    return lib_pairs


def score_residue(res_map, res_model, loaded_lib, log):
    log.info(f'Scoring {os.path.basename(res_map)}, {os.path.basename(res_model)} ')
    res_score = []
    mod = run_x('open ' + res_model)
    vol = run_x('open ' + res_map)
    for motif, lib_mod, lib_vol in loaded_lib:
        mod_sel = '#{}@c,ca,n'.format(mod.id_string)
        lib_mod_sel = '#{}@c,ca,n'.format(lib_mod.id_string)
        try:
            run_x(f'align {lib_mod_sel} to {mod_sel}')
        except:
            log.error(f'Backbone atoms missing in {os.path.basename(res_model)}. Cannot align! ')
            res_score.append((motif, None))
            continue
        run_x(f'view position #{lib_vol.id_string} sameAsModels #{lib_mod.id_string}')
        run_x(f'volume #{lib_vol.id_string} sdLevel {SD_LEVEL}')
        run_x(f'volume #{vol.id_string} sdLevel {SD_LEVEL}')

        try:
            fit1 = run_x(f'fitmap #{lib_vol.id_string} inMap #{vol.id_string} metric correlation')
            correlation1 = fit1.correlation()
        except Exception as err:
            log.error('Chimerax user error!! %s in map %s at sdLevel =  %s', err, lib_vol.name, SD_LEVEL)
            log.info('Changing threshold level to %s', SD_LEVEL - 0.1)
            try:
                run_x(f'volume #{lib_vol.id_string} sdLevel {SD_LEVEL - 0.1}')
                fit1 = run_x(f'fitmap #{lib_vol.id_string} inMap #{vol.id_string} metric correlation')
                correlation1 = fit1.correlation()
            except Exception as err:
                log.error('Chimerax user error!! %s in map %s at sdLevel =  %s', err, lib_vol.name, SD_LEVEL-0.1)
                log.info('Could not calculate %s and %s correlation', vol.name, lib_vol.name)
                correlation1 = None
        matrix = ','.join(f"{x}" for x in tuple(lib_vol.position.matrix.flat))
        try:
            fit2 = run_x(f'fitmap #{vol.id_string} inMap #{lib_vol.id_string} metric correlation')
            correlation2 = fit2.correlation()
        except Exception as err:
            log.info('User error!! %s %s', err, lib_vol.name)
            try:
                run_x(f'volume #{vol.id_string} sdLevel {SD_LEVEL - 0.1}')
                fit2 = run_x(f'fitmap #{vol.id_string} inMap #{lib_vol.id_string} metric correlation')
                correlation2 = fit2.correlation()
            except Exception as err:
                log.info('User error!! %s %s', err, lib_vol.name)
                correlation2 = None
        if all([correlation1, correlation2]):
            min_corr = min([correlation1, correlation2])
            res_score.append((motif, min_corr, matrix))
        else:
            res_score.append((motif, None, None))
        # This resets the view after fitting
        run_x('view orient; view initial')

    # res_score.sort(reverse=True, key=lambda z: z[1])
    return res_score


def capture_memory_usage():
    p = psutil.Process()
    mem = psutil.Process(p.pid).memory_info()
    return p.pid, mem.rss / 1000000


def slave(pairs_list, lib, score_list, lock, json_out_path, csv_out_path=None, close_rate=50):
    p = psutil.Process()
    log = logging.getLogger(str(p.pid))

    logging.basicConfig(filename='process_logs/' + str(p.pid) + '.log', level=logging.INFO, format='%(levelname)s:  %(message)s')
    def load_lib():
        log.info('Loading motif library')
        l_lib = []
        for mot in lib:
            mod = run_x('open ' + mot[1])
            vol = run_x('open ' + mot[2])
            l_lib.append([mot[0], mod, vol])
        return l_lib

    loaded_lib = load_lib()
    for k, pair in enumerate(pairs_list):
        if k > 0 and k % close_rate == 0:
            run_x('close session')
            loaded_lib = load_lib()

        #     p_id, mem = capture_memory_usage()
        #     with open(f'{p_id}.txt', 'a') as f:
        #         f.write('Closing session\n')
        #         f.write(f'{mem}\n')
        # p_id, mem = capture_memory_usage()
        # with open(f'{p_id}.txt', 'a') as f:
        #     f.write(f'{mem}\n')

        res_scores = score_residue(pair[0], pair[1], loaded_lib, log)

        res_name = os.path.basename(pair[1]).split('.')[0]
        name_lst = res_name.split('-')
        res_type = name_lst[0]
        res_nr = name_lst[1]
        chain = name_lst[2]
        residue_data = {keys.RES_TYPE: res_type,
                        keys.CHAIN: chain,
                        keys.RES_NR: int(res_nr),
                        keys.RES_MAP: os.path.basename(pair[0]),
                        keys.RES_MODEL: os.path.basename(pair[1]),
                        }
        for item in res_scores:
            # motif_name: correlation_m_matrix
            residue_data[item[0]] = f'{item[1]}_m_{item[2]}'

        lock.acquire()
        score_list.append(residue_data)
        if k % 10 == 0:
            scores_list_out = []
            for i in score_list:
                scores_list_out.append(i)
            scores_list_out.sort(key=operator.itemgetter(keys.CHAIN, keys.RES_NR))
            if os.path.exists(json_out_path):
                copy(json_out_path, json_out_path + '_bak')
            save_json(json_out_path, scores_list_out)
            if csv_out_path:
                with open(csv_out_path, mode='w') as csv_file:
                    fieldnames = [key for key in scores_list_out[0].keys()]
                    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                    writer.writeheader()
                    for item in scores_list_out:
                        writer.writerow(item)
        lock.release()


def score_structure(known_correlations, unknown_pairs, lib, json_out_path, csv_out_path=None, np=4):
    log_tmp = 'process_logs'

    if not os.path.exists(log_tmp):
        os.mkdir(log_tmp)

    score_list = Manager().list()
    for entry in known_correlations:
        score_list.append(entry)

    split_pairs = split_list(unknown_pairs, np)
    thread_list = []
    lock = Lock()
    for batch_list in split_pairs:
        t = Process(target=slave, args=(batch_list, lib, score_list, lock, json_out_path, csv_out_path))
        thread_list.append(t)
    for thread in thread_list:
        thread.start()
    for thread in thread_list:
        thread.join()

    scores_list_out = []
    for i in score_list:
        scores_list_out.append(i)
    scores_list_out.sort(key=operator.itemgetter(keys.CHAIN, keys.RES_NR))
    copy(json_out_path, json_out_path + '_bak')
    save_json(json_out_path, scores_list_out)
    if csv_out_path:
        with open(csv_out_path, mode='w') as csv_file:
            fieldnames = [key for key in scores_list_out[0].keys()]
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()
            for item in scores_list_out:
                writer.writerow(item)
    # Combine logs
    combine_files(log_tmp, 'score_' + datetime.now().strftime("%Y-%m-%d") + '.log')
    if os.path.exists(log_tmp):
        os.rmdir(log_tmp)

def combine_files(dir_path, out_file):
    files = os.listdir(dir_path)
    files = [f for f in files if not f.startswith('.')]

    with open(out_file, 'w') as outfile:
        if len(files) > 0:
            for file in files:
                with open(os.path.join(dir_path, file), 'r') as infile:
                    for line in infile:
                        outfile.write(line)



def read_json(json_path):
    with open(json_path) as j:
        data = json.load(j)
    return data


def save_json(file_path, data):
    with open(file_path, 'w') as j:
        json.dump(data, j, indent=4)


def split_list(inp_list, nr):
    """
    Splits evenly a list
    :param inp_list: list
    :param nr: number of parts
    :return: list of "nr" lists
    """
    new_list = []
    nr_el = 1.0 / nr * len(inp_list)
    for i in range(nr):
        start = int(round(i * nr_el))
        end = int(round((i + 1) * nr_el))
        new_list.append(inp_list[start:end])
    return new_list


def check_known_correlations(inp_corr_list, lib):
    incomplete_indices = []
    for index, entry in enumerate(inp_corr_list):
        if len(entry.keys()) - 5 < len(lib) or "None" in entry.values():
            incomplete_indices.append(index)
        incomplete_indices.sort(reverse=True)

    for index in incomplete_indices:
        # del incomplete_indices[index]
        del inp_corr_list[index]
    return inp_corr_list


def exclude_known_pairs(pair_list, inp_corr_list):
    unknown_pairs = []
    for pair in pair_list:
        for entry in inp_corr_list:
            if os.path.basename(pair[0]) in entry.values():
                break
        else:
            unknown_pairs.append(pair)
    return unknown_pairs


def main():
    import time
    start = time.time()
    print(sys.argv)
    st_arg = 1
    for item in sys.argv:
        if str(item).startswith('--'):
            st_arg += 1
    sys.argv = sys.argv[st_arg:]
    print(sys.argv)
    parser = argparse.ArgumentParser(description='Score map')
    parser.add_argument("-p", "--pair_list", dest="pair_list", required=True,
                        help="List of residue map, residue pdb paths as json file")
    parser.add_argument("-l", "--lib", dest="lib", required=True, help="Motif library in_dir")
    parser.add_argument("-np", "--n_processors", dest="np", required=False, help="Number of processors")
    parser.add_argument("-r", "--recompute", dest="recompute", action='store_true',
                        help="Recalculate correlations")
    parser.add_argument("-o", "--out", dest="out", required=True, help="Out json file in_dir")
    args = parser.parse_args()
    # print(args.pair_list)
    pair_list = read_json(args.pair_list)
    lib = read_motif_lib(args.lib)
    if not args.recompute:
        inp_correlations = None
        try:
            inp_correlations = read_json(args.out)
        except (FileNotFoundError, json.decoder.JSONDecodeError) as _:
            try:
                inp_correlations = read_json(args.out + '_bak')
            except (FileNotFoundError, json.decoder.JSONDecodeError) as _:
                pass
        if inp_correlations:
            known_correlations = check_known_correlations(inp_correlations, lib)
            pair_list = exclude_known_pairs(pair_list, known_correlations)
            # for p in pair_list:
            #     print(p)
        else:
            known_correlations = []
    else:
        known_correlations = []

    head, tail = os.path.split(args.out)
    csv_path = os.path.join(head, tail.split('.')[0] + '.csv')
    score_structure(known_correlations, pair_list, lib, args.out,  csv_path, np=int(args.np))
    print(f'Run time = {time.time()-start}')
    run_x('exit')

main()
