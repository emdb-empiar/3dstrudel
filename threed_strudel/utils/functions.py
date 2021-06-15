"""
functions.py

General purpose functions

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
__date__ = '2019-01-29'

import math
import logging
from datetime import datetime, timedelta
import time
from distutils.spawn import find_executable
import threed_strudel.configure as config
import json
import mrcfile


def global_except_hook(exctype, value, traceback):
    import sys
    try:
        import mpi4py.MPI
        sys.stderr.write("\n*****************************************************\n")
        sys.stderr.write("Uncaught exception was detected on rank {}. \n".format(
            mpi4py.MPI.COMM_WORLD.Get_rank()))
        from traceback import print_exception
        print_exception(exctype, value, traceback)
        sys.stderr.write("*****************************************************\n\n")
        sys.stderr.write("\n")
        sys.stderr.write("Calling MPI_Abort() to shut down MPI processes...\n")
        sys.stderr.flush()
    finally:
        try:
            import mpi4py.MPI
            mpi4py.MPI.COMM_WORLD.Abort(1)
        except Exception as e:
            sys.stderr.write("*****************************************************\n")
            sys.stderr.write("Sorry, we failed to stop MPI, this process will hang.\n")
            sys.stderr.write("*****************************************************\n")
            sys.stderr.flush()
            raise e


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


def split_sorted_list(inp_list, nr):
    new_list = []

    for i in range(nr):
        tmp = []
        for k in range(i, len(inp_list), nr):
            try:
                tmp.append(inp_list[k])
            except IndexError:
                pass
        new_list.append(tmp)
    return new_list


def split_weighted_dict(dictionary, nr):
    """
    From a dictionary ({'ASP': 20, 'GLU': 30, 'PHE': 10, 'ARG': 50, ..)
    Creates a list of lists of keys of the input dictionary such as in each
    group the sum of the corresponding values is minimal.
    :param dictionary:
    :param nr:
    :return:
    """
    out_list = []
    weights = []
    sorted_dict = sorted(dictionary.items(), key=lambda kv: kv[1], reverse=True)

    for key, value in sorted_dict[0:nr]:
        out_list.append([key])
        weights.append(value)

    for key, value in sorted_dict[nr:]:
        min_index = weights.index(min(weights))
        out_list[min_index].append(key)
        weights[min_index] = weights[min_index] + value

    return out_list


def setup_logger(name, log_file=None, warning_level="info"):
    """
    Function setup as many loggers as needed
    :param name: logger name
    :param log_file: process_log file in_dir
    :param warning_level: warning level
    :return: logger
    """
    if warning_level.lower() == 'debug':
        warning_level = logging.DEBUG
    else:
        warning_level = logging.INFO
    formatter = logging.Formatter('%(message)s')
    if log_file:
        handler = logging.FileHandler(log_file)
    else:
        handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(warning_level)
    # Remove old handlers
    for h in logger.handlers:
        logger.removeHandler(h)
    logger.addHandler(handler)

    return logger


def report_elapsed(start):
    """
    Reports elapsed time since start in human readable formatted
    :param start: script starting time
    :return: report
    """
    end = time.time()
    elapsed = round(end - start, 2)
    sec = timedelta(seconds=int(elapsed))
    d = datetime(1, 1, 1) + sec
    formatted = "d:h:m:s {}:{}:{}:{}".format(d.day - 1, d.hour, d.minute, d.second)
    text = '{} seconds   {}'.format(elapsed, formatted)
    return text


def print_time():
    """
    Returns human readable date and time
    :return:
    """
    date_time = datetime.now().strftime("%Y-%m-%d  %H:%M:%S")
    return date_time


def find_good_grid(in_grid):
    def largest_prime_factor(n):
        i = 2
        while i * i <= n:
            if n % i:
                i += 1
            else:
                n //= i
        return n

    out_grid = []
    for dimension in in_grid:
        new_dimension = dimension
        if new_dimension % 2:
            new_dimension += 1
        largest_prime = largest_prime_factor(new_dimension)
        while largest_prime > 19:
            new_dimension += 2
            largest_prime = largest_prime_factor(new_dimension)
        out_grid.append(new_dimension)
    return tuple(out_grid)


def get_map_grid_size(map_path):
    with mrcfile.mmap(map_path, mode='r+', permissive=True) as mrc:
        grid = (mrc.header.mx, mrc.header.my, mrc.header.mz)
    return grid


def calc_fraction_matrix(cell, ang):
    """
    Calculates matrix used to transform Cartesian
    coordinates in the ATOM_SITE category to fractional coordinates
    in the same category
    :param cell: unit cell
    :param ang: unit cell angles in degrees
    :return:
    """
    ang = (ang[0]*math.pi/180, ang[1]*math.pi/180, ang[2]*math.pi/180)

    omega = cell[0]*cell[1]*cell[2] * math.sqrt(1 - math.cos(ang[0])**2 - math.cos(ang[1])**2 - math.cos(ang[2])**2 +
                                                2*math.cos(ang[0])*math.cos(ang[1])*math.cos(ang[0]))
    m11 = 1 / cell[0]
    m12 = - math.cos(ang[2]) / (cell[0]*math.sin(ang[2]))
    m13 = cell[1]*cell[2] * (math.cos(ang[0])*math.cos(ang[2]) - math.cos(ang[1])) / (omega * math.sin(ang[2]))

    m21 = 0
    m22 = 1 / (cell[1]*math.sin(ang[2]))
    m23 = cell[0]*cell[2] * (math.cos(ang[1])*math.cos(ang[2]) - math.cos(ang[0])) / (omega * math.sin(ang[2]))

    m31 = 0
    m32 = 0
    m33 = cell[0]*cell[1]*math.sin(ang[2]) / omega
    matrix = [m11, m12, m13, m21, m22, m23, m31, m32, m33]
    # for i, el in enumerate(matrix):
    #     matrix[i] = round(el, 6)

    return matrix


