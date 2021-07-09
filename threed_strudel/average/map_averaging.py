"""
mapAveraging.py

Superimpose and average maps using atom models as guide.
This code use maps in mrc format
This code uses ChimeraX for superimposing.
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

import logging
import sys
import os
import time
from shutil import copy2
import argparse
from datetime import datetime
import numpy as np
import json
from collections import Counter
import subprocess
import threed_strudel.configure as config
import threed_strudel.utils.functions as func
from threed_strudel.parse.map_parser import MapParser
import threed_strudel.average




class MapAveraging:
    """
    Class for map superimposing using atom models as a guide and map averaging
    """

    def __init__(self, model_dir,
                 rotamer_classes_dir,
                 rot_prefix='rot_',
                 map_class_prefix='class_',
                 box_size=14,
                 motifs_dir=None,
                 scaling='standardise',
                 motif_suffix='',
                 min_reliable_nr=10,
                 log_path=None,
                 warning_level='info',):

        self.scaling = scaling
        self.motif_suffix = motif_suffix
        self.log = self.setup_log(log_path, warning_level)
        self.chimera_path = config.CHIMERA_PATH
        basedir = os.path.dirname(os.path.abspath(__file__))
        self.superimpose_chimera_script = os.path.join(basedir, 'superimpose_chimerax.py')
        self.final_box_size = box_size
        self.global_reference_path = self.get_global_reference_path()
        self.rotamer_classes_dir = os.path.abspath(rotamer_classes_dir)
        self.models_dir = os.path.abspath(model_dir)
        self.motifs_dir = motifs_dir
        self.rot_prefix = rot_prefix
        self.map_class_prefix = map_class_prefix
        self.statistics = []
        self.rotamers = []
        self.classes = []
        self.representative_structures = []
        self.set_env()
        self.map_suffix = None
        self.representative = ''
        self.min_reliable_number = min_reliable_nr
        self.unreliable_motif_folder = 'unreliable'
        self.unreliable_motif_dir = os.path.join(self.motifs_dir, 'unreliable')
        self.final_map = ''
        self.final_model = ''
        if not os.path.exists(self.unreliable_motif_dir):
            os.makedirs(self.unreliable_motif_dir)
        self.run_main()

    @staticmethod
    def setup_log(log_path, warning_level):
        aver_log = func.setup_logger('average_log', log_path, warning_level=warning_level)
        return aver_log

    @staticmethod
    def get_global_reference_path():
        head = os.path.dirname(threed_strudel.average.__file__)
        path = os.path.join(head, 'inp/global_reference.pdb')
        return path

    def run_main(self):
        start = time.time()

        if len(self.classes) > 0:
            for index, class_name in enumerate(self.classes):
                class_dict = self.get_class_info(self.rotamers[index], class_name)

                # Do the averaging using the atomic models as reference
                out_folder = 'model_super'
                out_dir = os.path.join(self.rotamer_classes_dir, class_name, out_folder)
                class_dict['out_dir'] = out_dir
                self.create_dir(out_dir)

                if len(class_dict['class_pairs']) > 0:
                    t = time.time()
                    self.log.info('Superimposing maps using residue as a guide...')
                    chimera_out = self.superimpose(class_dict, fit_map=False)
                    if chimera_out != 0:
                        self.log.error('ChimeraX returned code %s. Check %s for details.',
                                       chimera_out, os.path.join(class_dict['out_dir'], 'chimera.out'))
                    self.log.debug('spent on superimposing %s', func.report_elapsed(t))
                    t = time.time()
                    self.log.info("Averaging the residue superimposed maps")
                    aver_map = self.average(class_dict)
                    self.log.debug('spent on averaging %s', func.report_elapsed(t))
                    # Copy the averaged densities to the motifs directory
                    name = os.path.basename(class_dict['reference_model'])
                    aligned_reference = os.path.join(class_dict['out_dir'], name)
                    copy2(aligned_reference, os.path.join(self.rotamer_classes_dir, class_name))
                    # Additional map to average map fitting ang repeat averaging
                    out_folder = 'map_super'
                    out_dir = os.path.join(self.rotamer_classes_dir, class_name, out_folder)
                    class_dict['out_dir'] = out_dir
                    self.create_dir(out_dir)
                    if aver_map is not None:
                        class_dict['reference_model'] = aligned_reference
                        class_dict['reference_mrc'] = aver_map
                        t = time.time()
                        self.log.info('Superimposing maps using residue as a guide followed by map to map fit...')
                        chimera_out = self.superimpose(class_dict, fit_map=True)
                        if chimera_out != 0:
                            self.log.error('ChimeraX returned code %s. Check %s for details.',
                                           chimera_out, os.path.join(class_dict['out_dir'], 'chimera.out'))
                        self.log.debug('spent on superimposing %s', func.report_elapsed(t))
                        t = time.time()
                        self.log.info("Averaging the map to map superimposed maps")
                        aver_map = self.average(class_dict)
                        self.log.debug('spent on averaging %s', func.report_elapsed(t))

                        # Copy the reference atomic models to the motifs directory
                        prefix = f'{name.split("-")[0]}_rotamer_{class_name.split("_")[-1]}{self.motif_suffix}'
                        map_motif_name = f'{prefix}{self.motif_suffix}.mrc'
                        model_motif_name = f'{prefix}{self.motif_suffix}.cif'

                        if len(class_dict['class_pairs']) >= self.min_reliable_number:
                            dst = self.motifs_dir
                        else:
                            dst = self.unreliable_motif_dir

                        if os.path.exists(aligned_reference) and os.path.exists(aver_map):
                            final_model_path = os.path.join(dst, model_motif_name)
                            final_map_path = os.path.join(dst, map_motif_name)
                            copy2(aligned_reference, final_model_path)
                            copy2(aver_map, final_map_path)
                            self.flag_derived_motif(self.rotamers[index].split('_')[-1], True)
                        else:
                            self.flag_derived_motif(self.rotamers[index].split('_')[-1], False)

        text = func.report_elapsed(start)
        self.log.info('\nElapsed: {}\n{:_^100}'.format(text, ''))

    def flag_derived_motif(self, rotamer_name, flag):
        """
        Flags motif deriving status
        :param rotamer_name: rotamer name
        :param flag: boolean
        """
        for item in self.statistics:
            if item['rotamer'] == rotamer_name:
                item['derived'] = flag

    @staticmethod
    def create_dir(path):
        """
        Create new directory
        :param path: new directory in_dir
        """
        if not os.path.exists(path):
            os.mkdir(path)

    def set_env(self):
        """
        For each folder named rot_prefix* create a folder named map_class_prefix*
        """
        date_time = datetime.now().strftime("%H:%M:%S %Y-%m-%d")
        text = '{:_^100}'.format('MapAveraging') + '\n\nStarted: {}'.format(date_time)
        self.log.info('%s\n\n', text)
        if not config.check_chimerax_executable():
            self.log.error(config.no_chimerax_error())
            sys.exit()
        self.log.info("Searching rotamers in %s directory", self.rotamer_classes_dir)
        rot_folder_lst = [item for item in os.listdir(self.rotamer_classes_dir) if item.startswith(self.rot_prefix)]
        rot_folder_lst.sort()

        for rot_folder in rot_folder_lst:
            self.rotamers.append(rot_folder)
            suffix = rot_folder.split(self.rot_prefix)[-1]
            map_class_folder = self.map_class_prefix + suffix
            self.classes.append(map_class_folder)
            class_dir = os.path.join(self.rotamer_classes_dir, map_class_folder)
            self.create_dir(class_dir)

        if self.motifs_dir is None:
            self.motifs_dir = os.path.join(self.rotamer_classes_dir, 'motifs')
        if not os.path.exists(self.motifs_dir):
            os.makedirs(self.motifs_dir)

        self.log.info("The following rotamer folders were found:")
        for r in self.rotamers:
            self.log.info("%s", r)

    def get_class_info(self, rot_name, class_name):
        """
        Creates a dictionary with the relevant information for the rotamer class
        :param rot_name: input rotamer class folder
        :param class_name: output rotamer class files folder
        :return: dictionary
        """
        class_dict = {'inp_dir': self.models_dir,
                      'class': class_name,
                      'class_pairs': [],
                      'reference_model': 'miss',
                      'reference_mrc': None,
                      }

        all_files = os.listdir(self.models_dir)
        all_model_names = [i for i in all_files if i.endswith('.cif')]
        all_mrc_names = [i for i in all_files if i.endswith('.mrc')]
        self.map_suffix = '_' + all_mrc_names[0].split('_')[-1]

        rotamer_dir = os.path.join(self.rotamer_classes_dir, rot_name)
        files = os.listdir(rotamer_dir)

        rotamer_model_files = [i for i in files if i.endswith('.cif') and not i.startswith('class')]
        for model in rotamer_model_files:
            if model.endswith('representative.cif'):
                tmp = model.split('_representative')[0] + '.cif'
                class_dict['reference_model'] = os.path.join(self.models_dir, tmp)
                class_dict['reference_mrc'] = os.path.join(self.models_dir, f"{tmp.split('_')[0]}{self.map_suffix}")
                pair = (tmp.split('_')[0] + self.map_suffix, tmp)
                class_dict['class_pairs'].append(pair)

            if model in all_model_names:
                base_name = model.split('_')[0]
                mrc = base_name + self.map_suffix
                if mrc in all_mrc_names:
                    pair = (mrc, model)
                    class_dict['class_pairs'].append(pair)
        nr = len(class_dict['class_pairs'])
        self.log.info("\nProcessing %s\n%s models were found for the %s rotamer", rot_name, nr, rot_name)
        self.log.info("The output files will be saved in the %s folder", class_name)
        if len(rotamer_model_files) > 0:
            residue_type = rotamer_model_files[0].split('-')[0]
            motif_stats = {'residue': residue_type,
                           'rotamer': rot_name.split('_')[-1],
                           'fragments_per_entry': dict(self.count_entries(rotamer_model_files)),
                           'derived': False
                           }
            self.statistics.append(motif_stats)

        return class_dict

    @staticmethod
    def count_entries(files_list):
        """
        Counts the number of maps fragments coming from each entry
        :param files_list: list of map or residue files
        :return: counter object
        """
        entries = []
        for file in files_list:
            entry = file.split('_')[0].split('-')[-1]
            entries.append('EMD-' + entry)
        return Counter(entries)

    def create_reference_grid(self, map_example, out_map_path):
        """
        Creates a reference map grid. All densities will be resampled on to this grid after superimposing
        :param map_example: any map from the the list of maps to be averaged
        :param out_map_path: in_dir to the grid
        """
        map_obj = MapParser(map_example)
        c = self.final_box_size
        size = int(c / map_obj.voxel_size[0])
        grid = np.zeros((size, size, size), dtype='float32')
        out_map_obj = MapParser()
        out_map_obj.data = grid
        out_map_obj.cell = (c, c, c)
        out_map_obj.origin = (-c / 2, 0, 0)
        out_map_obj.write_map(out_map_path)

    def superimpose(self, class_dict, fit_map=False):
        """
        Performs map superimposing in Chimera.
        :param class_dict: rotamer class dictionary
        :param fit_map: (Boolean) To do map to map fitting or not
        """
        class_dict['fit_map'] = fit_map
        class_dict['global_reference'] = self.get_global_reference_path()
        out_dir = class_dict['out_dir']
        inp_dir = class_dict['inp_dir']
        # Create a standard grid
        map_sample = os.path.join(inp_dir, class_dict['reference_mrc'])
        grid_file = 'grid_cub_' + str(self.final_box_size) + '.mrc'
        grid_path = os.path.join(out_dir, grid_file)
        self.create_reference_grid(map_sample, grid_path)
        class_dict['grid_path'] = grid_path

        json_path = os.path.join(out_dir, 'class_dict.json')
        with open(json_path, 'w') as j_file:
            json.dump(class_dict, j_file, indent=4)

        copy2(self.superimpose_chimera_script, out_dir)
        py_path = os.path.join(out_dir, os.path.basename(self.superimpose_chimera_script))
        command = self.chimera_path + ' --nogui ' + py_path + ' > chimera.out'
        return subprocess.call(command, cwd=out_dir, shell=True)

    @staticmethod
    def normalize_truncate_neg(xyz, new_max_i):
        """
        Computes image normalisation
        Negative values are set to 0
        :param xyz: numpy array
        :param new_max_i: maximal output voxel intensity
        :return: numpy array
        """
        max_i = float(np.amax(xyz))
        xyz = xyz.clip(0, max_i) * new_max_i / max_i
        return xyz

    @staticmethod
    def normalize(xyz, new_max_i):
        """
        Computes image normalisation
        The negative and positive ranges are computed separately
        The output negative value range is dependant on the positive range
        :param xyz: numpy array
        :param new_max_i: maximal output voxel intensity
        :return: numpy array
        """
        xyz[np.isnan(xyz)] = 0
        min_i = float(np.amin(xyz))
        max_i = float(np.amax(xyz))
        if max_i != 0:
            new_min_i = new_max_i * min_i / max_i
        else:
            return None

        a = xyz.clip(0, max_i) * new_max_i / max_i
        b = (xyz.clip(min_i, 0) - min_i) * new_min_i / min_i + new_min_i

        return a + b

    @staticmethod
    def standardise(xyz, empty):
        """
        Computes image normalisation
        The negative and positive ranges are computed separately
        The output negative value range is dependant on the positive range
        :param xyz: numpy array
        :return: numpy array
        """
        std = np.std(xyz)
        mean = np.mean(xyz)
        xyz[np.isnan(xyz)] = 0
        standardised = (xyz - mean) / std

        return standardised

    def average(self, class_dict):
        """
        Map averaging
        :param class_dict: rotamer class dictionary
        """
        cl_dir = class_dict['out_dir']
        reference_map_path = class_dict['grid_path']

        if self.scaling == 'normalise-all':
            scaling_function = self.normalize
        elif self.scaling == 'normalise-positive':
            scaling_function = self.normalize_truncate_neg
        elif self.scaling == 'standardise':
            scaling_function = self.standardise
        else:
            self.log.info("Error: The normalization flag can be 'normalise-all', "
                          "'normalise-positive' or 'standardise' not %s", self.scaling)
            sys.exit()

        reference_map = MapParser(reference_map_path)
        tmp_grid = np.zeros(reference_map.n_xyz, dtype='float32')
        cella = reference_map.cell
        origin = reference_map.origin
        exceptions = 0
        for pair in class_dict['class_pairs']:
            mrc = pair[0]
            mrc_path = os.path.join(cl_dir, mrc)

            map_obj = MapParser(mrc_path)
            if map_obj.data is None:
                self.log.info("Could not read %s data", mrc_path)
                exceptions += 1
                continue
            map_obj.data = scaling_function(map_obj.data, 1)

            if map_obj.data is not None:
                if np.isnan(np.sum(map_obj.data)):
                    self.log.info("Map %s has NaN values ", mrc_path)
                    map_obj.data[np.isnan(map_obj.data)] = 0
                map_obj.write_map(mrc_path)
                try:
                    tmp_grid = tmp_grid + map_obj.data
                except ValueError:
                    self.log.info("The grid sampling of the \n%s %s\n map is different "
                                  "from the reference map \n%s %s",
                                  mrc_path, map_obj.data.shape, reference_map_path, tmp_grid.shape)
                    exceptions += 1
            else:
                exceptions += 1

        tmp_grid = tmp_grid / (len(class_dict['class_pairs']) - exceptions)
        averaged_map = MapParser()
        averaged_map.cell = cella
        averaged_map.data = np.float32(tmp_grid)
        averaged_map.origin = origin

        head, tail = os.path.split(cl_dir)
        prefix = class_dict['class_pairs'][0][0].split('-')[0] + '_' + class_dict['class']
        out_map = os.path.join(head, prefix + '_aver_' + tail + '.mrc')
        averaged_map.write_map(out_map)

        return out_map


def main():
    parser = argparse.ArgumentParser(description='Create a Mask around the atomic residue')
    parser.add_argument("-d", "--model_dir", dest="m_dir", required=True, help="Models and maps directory")
    parser.add_argument("-r", "--rot_dir", dest="r_dir", required=True, help="Rotamer classes directory")
    parser.add_argument("-p", "--out_rot_prefix", dest="rot_prefix", default='rot_', help="Rotamer folders prefix")
    parser.add_argument("-c", "--out__class_prefix", dest="class_prefix", default='class_',
                        help="Classes folders prefix")
    parser.add_argument("-b", "--box", dest="box", default=14, help="Output motifs map box size")
    parser.add_argument("-m", "--mot_dir", dest="mot", default=None, help="Output motifs directory")
    parser.add_argument("-s", "--scaling", dest="scaling", default='standardise',
                        help="Image scaling normalization: 'normalise-all' - full intensity range, "
                             "'normalise-positive' - negative values are sett to 0 "
                             "'standardise' - standardisation (x-mean)/std")
    parser.add_argument("-l", "--log", dest="log", default=None, help="Log file path")
    parser.add_argument("-nr", "--rel_nr", dest="rel_nr", default=10, help="Minimum number of maps for a reliable motif")
    args = parser.parse_args()

    MapAveraging(args.m_dir, args.r_dir, args.rot_prefix, args.class_prefix, args.box,
                 args.mot, scaling=args.scaling, min_reliable_nr=args.rel_nr)


if __name__ == '__main__':
    main()
