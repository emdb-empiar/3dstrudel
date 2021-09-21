import os
import csv
import sys
from multiprocessing import Process, Manager, Lock
from datetime import datetime
import multiprocessing
from shutil import copy, rmtree
import json
import argparse
import operator
import psutil
# from chimera import runCommand as rc
# from chimera import openModels as om, selection
# from FitMap.fitcmd import fitmap
from chimerax.core.commands import Command
from chimerax.map.volume import Volume
from chimerax.atomic.structure import AtomicStructure
import logging


log = logging.getLogger(__name__)
rc = Command(session)
SD_LEVEL = 3


def run_x(command):
    try:
        return rc.run(command)[0][0]
    except TypeError:
        return None


def open_model(path):
    try:
        model_obj = session.open_command.open_data(path)[0][0]
        if type(model_obj) == AtomicStructure or type(model_obj) == Volume:
            session.models.add([model_obj])
            return model_obj
    except AttributeError:
        model_obj = run_x(f'open {path}')
        if type(model_obj) == AtomicStructure or type(model_obj) == Volume:
            return model_obj
        else:
            log.error(f'Could not get model object reference after opening:\n {path}')


def read_motif_lib(motif_lib_dir):
    """
    Reads motif library
    :param motif_lib_dir: library path
    :return: list of map model pair paths
    """
    lib_pairs = []
    files = os.listdir(motif_lib_dir)
    pdb_names = [i for i in files if i.endswith('.cif') and not i.startswith('.')]
    for name in pdb_names:
        prefix = name.split('.cif')[0]
        map_path = os.path.join(motif_lib_dir, prefix + '.map')
        mrc_path = os.path.join(motif_lib_dir, prefix + '.mrc')
        if os.path.exists(map_path):
            lib_pairs.append([prefix, os.path.join(motif_lib_dir, name), map_path])
        elif os.path.exists(mrc_path):
            lib_pairs.append([prefix, os.path.join(motif_lib_dir, name), mrc_path])
    return lib_pairs

def load_lib(lib):
    log.debug('Loading motif library')
    l_lib = []
    for mot in lib:
        # mod = run_x('open ' + mot[1])
        # vol = run_x('open ' + mot[2])
        mod = open_model(mot[1])
        vol = open_model(mot[2])
        l_lib.append([mot[0], mod, vol])
    return l_lib


def main():
    st_arg = 1
    for item in sys.argv:
        if str(item).startswith('--'):
            st_arg += 1
    sys.argv = sys.argv[st_arg:]

    parser = argparse.ArgumentParser(description='Score map')
    parser.add_argument("-sl", "--strudel_lib", dest="lib", required=True, help="Motif library in_dir")
    parser.add_argument("-l", "--log", dest="log", required=True, help="Log file")
    parser.add_argument("-o", "--out", dest="csv", required=True, help="csv_path")

    args = parser.parse_args()

    logging.basicConfig(filename=args.log, level=logging.INFO, format='%(levelname)s:  %(message)s')

    data = []

    lib_pairs = read_motif_lib(args.lib)
    loaded_lib = load_lib(lib_pairs)

    for item in loaded_lib:
        row = {}
        row['target'] = item[0]
        t_mod = item[1]
        t_map = item[2]
        t_mod_sel = '#{}@c,ca,n'.format(t_mod.id_string)

        for mot in loaded_lib:
            if mot[0] == item[0]:
                row[mot[0]] = 1
            else:
                l_mod = mot[1]
                l_map = mot[2]

                l_mod_sel = '#{}@c,ca,n'.format(l_mod.id_string)

                # run_x(f'align {l_mod_sel} to {t_mod_sel}')

                # run_x(f'view position #{l_map.id_string} sameAsModels #{l_mod.id_string}')
                run_x(f'volume #{l_map.id_string} sdLevel {SD_LEVEL}')
                run_x(f'volume #{t_map.id_string} sdLevel {SD_LEVEL}')

                fit1 = run_x(f'fitmap #{l_map.id_string} inMap #{t_map.id_string} metric correlation')
                correlation1 = fit1.correlation()

                fit2 = run_x(f'fitmap #{t_map.id_string} inMap #{l_map.id_string} metric correlation')
                correlation2 = fit2.correlation()

                min_corr = round(min([correlation1, correlation2]), 5)
                if min_corr < 0:
                    log.info(f'cor: {min_corr}, {item[0]}, {mot[0]}')
                row[mot[0]] = min_corr
                run_x('view orient; view initial')

        data.append(row)

    with open(args.csv, mode='w') as csv_file:
        fieldnames = [key for key in data[0].keys()]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for item in data:
            writer.writerow(item)
    run_x('exit')


main()
"""
usage
'ChimeraX_path --nogui  this_script -sl strudel_libpath -l log_path -o csv_path > chimera.log

"""






