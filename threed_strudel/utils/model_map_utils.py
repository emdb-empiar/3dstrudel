"""
modelMapUtils.py

A collection of functions for atomic models and density maps manipulations

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

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def interpolator(map_object):
    """
    Setup scipy regular grid interpolation function
    :return: interpolation function
    """
    nx, ny, nz = map_object.data.shape
    x = range(nx)
    y = range(ny)
    z = range(nz)
    _interpolator = RegularGridInterpolator((x, y, z), map_object.data)
    return _interpolator


def atom_inclusion(structure, map_object, threshold, ignore_atoms=()):
    """
    Counts the number of atoms within and outside map at a given level
    :param threshold:
    :param map_object:
    :param structure:
    :param ignore_atoms:
    :return: [included(nr), not_included(nr), fraction]
    """
    if structure.get_level() == 'S':
        model = structure[0]
    else:
        model = structure

    ignore_atoms = [atom.lower() for atom in ignore_atoms]

    _interpolator = interpolator(map_object)

    included = 0
    not_included = 0
    for atom in model.get_atoms():
        if atom.get_name().lower() not in ignore_atoms:
            map_index = map_object.coord_to_index(atom.coord)
            try:
                val = _interpolator(map_index)
            except ValueError:
                val = threshold

            if val > threshold:
                included += 1
            else:
                not_included += 1
        else:
            included += 1

    fraction = included / (included + not_included)

    return included, not_included, fraction


def find_threshold(structure, map_object, inclusion_fraction, delta=2):
    """
    Calculates the threshold level for which the specified fraction of atoms is inside the map
    :param map_object:
    :param structure:
    :param inclusion_fraction: atom inclusion fraction (%)
    :param delta: atom inclusion fraction precision
    :return: threshold level
    """
    _interpolator = interpolator(map_object)
    upper = np.nanmax(map_object.data)
    lower = np.nanmin(map_object.data)
    values = []
    for atom in structure.get_atoms():
        map_index = map_object.coord_to_index(atom.coord)
        try:
            val = _interpolator(map_index)
        except ValueError:
            val = lower
        values.append(val)
    values = np.array(values)
    total_atoms = len(values)
    current_lvl = 0
    final_lvl = 0.
    delta_lvl = upper - lower
    while delta_lvl > (upper - lower) / 100:
        tmp = np.where(values > current_lvl)
        included = len(values[tmp])
        not_included = total_atoms - included

        current_incl = included / (included + not_included) * 100
        if current_lvl == 0 and current_incl < inclusion_fraction:
            final_lvl = 0.
            break
        elif current_incl < inclusion_fraction:
            upper = current_lvl
            delta_lvl = (upper - lower) / 2
            current_lvl = current_lvl - delta_lvl
        elif current_incl > inclusion_fraction + delta:
            lower = current_lvl
            delta_lvl = (upper - lower) / 2
            current_lvl = current_lvl + delta_lvl
        else:
            final_lvl = round(current_lvl, 8)
            break

    return final_lvl


def coord_to_index(coord, map_obj):
    """
    Converts atomic residue coordinates into map indices
    :param coord: coordinate in angstrom
    :param map_obj: map object
    :return: float indices
    """
    st = map_obj.n_start
    vs = map_obj.voxel_size
    orig = map_obj.origin

    x_index = coord[2] / vs[2] - orig[2] / vs[2] - st[2]
    y_index = coord[1] / vs[1] - orig[1] / vs[1] - st[1]
    z_index = coord[0] / vs[0] - orig[0] / vs[0] - st[0]
    return x_index, y_index, z_index


def coord_to_index_(coord, map_obj):
    """
    Converts atomic residue coordinates into map indices
    :param coord: coordinate in angstrom
    :param map_obj: map object
    :return: float indices
    """
    return coord[2] / map_obj.voxel_size[2] - map_obj.origin[2] / map_obj.voxel_size[2] - map_obj.n_start[2], \
           coord[1] / map_obj.voxel_size[1] - map_obj.origin[1] / map_obj.voxel_size[1] - map_obj.n_start[1], \
           coord[0] / map_obj.voxel_size[0] - map_obj.origin[0] / map_obj.voxel_size[0] - map_obj.n_start[0]


def coord_to_index_int(coord, map_obj):
    """
    Converts atomic residue coordinates into map indices
    :param coord: coordinate in angstrom
    :param map_obj: map object
    :return: float indices
    """
    st = map_obj.n_start
    vs = map_obj.voxel_size
    orig = map_obj.origin

    x_index = int(round(coord[2] / vs[2] - orig[2] / vs[2] - st[2]))
    y_index = int(round(coord[1] / vs[1] - orig[1] / vs[1] - st[1]))
    z_index = int(round(coord[0] / vs[0] - orig[0] / vs[0] - st[0]))
    return x_index, y_index, z_index
