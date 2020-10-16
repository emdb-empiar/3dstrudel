"""
chopMap.py

A collection of methods which allow to chop a map around an atomic model.
The map can be chopped in three different ways:
- using a cub around the atomic residue with hard edges
- using certain radius around atomic residue with hard edges
- using certain radius around atomic residue with soft edges

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

import math
from random import randint
import mrcfile
import numpy as np

from strudel.parse.mapParser import MapParser
from strudel.utils import bioUtils


class ChopMap:
    """
    Class for chopping map around an atomic model
    """

    def __init__(self, in_model=None, in_map=None):
        self.model = in_model
        self.map = in_map
        self.sett_model_map()

    def sett_model_map(self):
        """
        Loads the map and model
        """
        if self.map is not None:
            self.map = MapParser(self.map)
        if self.model is not None:
            self.model = bioUtils.load_structure(self.model)

    @staticmethod
    def chop_cube_list(in_model, in_map, cube_padding, zero_origin=True):
        """
         Chop map using a cubic box around the model
         :param in_model: biopython model object or list of objects
         :param in_map: path to the input map or strudel map object
         :param cube_padding: distance fom any cube edge to the nearest atom
         :param zero_origin: boolean,
         :return: map object or list of map objects, translation matrix or list of matrices depending on input 
         """
        lst = True
        if type(in_map) is str:
            in_map = MapParser(in_map)
        if type(in_model) is not list and type(in_model) is not tuple:
            in_model = [in_model]
            lst = False

        shifts_list = []
        chopped_map_obj_list = []

        for model in in_model:

            atoms_coord = []
            for atom in model.get_atoms():
                atoms_coord.append(atom.coord)

            delta = round(cube_padding / in_map.voxel_size[0])
            # Get the indices in the map grid that correspond to the atom coordinates
            x_indices = []
            y_indices = []
            z_indices = []
            for atom in atoms_coord:
                x_index, y_index, z_index = in_map.coord_to_index(atom)
                if any([x_index, y_index, z_index]) < 0:
                    print("Atom outside map")
                else:
                    x_indices.append(x_index)
                    y_indices.append(y_index)
                    z_indices.append(z_index)
            # Find the voxel located in the middle of the atomic residue
            # and the maximum molecule size in grid points
            deltas = []
            minim = min(x_indices)
            dx = max(x_indices) - minim
            middle_x = int(round(dx / 2 + minim))
            deltas.append(dx)

            minim = min(y_indices)
            dy = max(y_indices) - minim
            middle_y = int(round(dy / 2 + minim))
            deltas.append(dy)

            minim = min(z_indices)
            dz = max(z_indices) - minim
            middle_z = int(round(dz / 2 + minim))
            deltas.append(dz)
            max_d = max(deltas)
            # Calculate the size of the cube
            radius = int(round(max_d / 2 + delta))
            new_dimension = radius * 2
            # Ensure that the new grid size has no prime numbers greater than 19

            # new_dimension = func.find_good_grid([new_dimension, new_dimension, new_dimension])[0]
            # print(new_dimension)
            # hard_radius = int(new_dimension / 2)
            # Create a numpy array to store the chopped voxels
            new_data = np.zeros((new_dimension, new_dimension, new_dimension), dtype='float32')

            # Assign voxel values
            for x in range(new_dimension):
                for y in range(new_dimension):
                    for z in range(new_dimension):
                        try:
                            new_data[x, y, z] = in_map.data[x + middle_x - radius, y + middle_y - radius,
                                                            z + middle_z - radius]
                        except IndexError:
                            pass
            # Calculate the new cell size
            voxel_size = in_map.voxel_size
            new_cell = (round(new_dimension * voxel_size[0], 3),
                        round(new_dimension * voxel_size[1], 3),
                        round(new_dimension * voxel_size[2], 3))

            # Calculate the shifts applied to the chopped map
            if not zero_origin:
                shifts = np.array([0, 0, 0])
            else:
                index_shifts = in_map.index_shifts
                x_shift = (middle_z - radius) * voxel_size[0] + index_shifts[0] * voxel_size[0]
                y_shift = (middle_y - radius) * voxel_size[1] + index_shifts[1] * voxel_size[1]
                z_shift = (middle_x - radius) * voxel_size[2] + index_shifts[2] * voxel_size[2]

                shifts = np.array([x_shift, y_shift, z_shift])

            new_data[np.isnan(new_data)] = 0

            out_map = MapParser('')
            out_map.data = new_data
            out_map.cell = new_cell
            out_map.cellb = in_map.cellb
            out_map.voxel_size = voxel_size
            n_st = in_map.n_start
            if not zero_origin:
                origin = ((middle_z - radius + n_st[0]) * out_map.voxel_size[0] + in_map.origin[0],
                          (middle_y - radius + n_st[1]) * out_map.voxel_size[1] + in_map.origin[1],
                          (middle_x - radius + n_st[2]) * out_map.voxel_size[2] + in_map.origin[2])
                out_map.set_origin(origin)
            chopped_map_obj_list.append(out_map)
            shifts_list.append(shifts)
            # match = self.check_cut_vals(residue, in_map, out_map, shifts)
            # if match:
            # chopped_map_obj_list.append(out_map)
            # shifts_list.append(shifts)
            # else:
            # chopped_map_obj_list.append(None)
            # shifts_list.append(None)
            # print("NO  MATCH")

        if len(shifts_list) > 1 or lst:
            return chopped_map_obj_list, shifts_list
        elif len(shifts_list) == 1:
            return chopped_map_obj_list[0], shifts_list[0]
        else:
            return None, None

    @staticmethod
    def check_cut_vals(model, in_map, out_map, shifts):
        """
        Checks the map values at 5 atoms positions
        :param model:
        :param in_map:
        :param out_map:
        :param shifts:
        :return:
        """
        atoms = [a for a in model.get_atoms()]
        in_coords = []
        out_coords = []

        for i in range(5):
            n = randint(0, len(atoms))
            try:
                in_coords.append(atoms[n].coord)

            except IndexError:
                pass
        for coord in in_coords:
            tmp_coord = list(coord)
            for i in range(3):
                tmp_coord[i] = tmp_coord[i] - shifts[i]
            out_coords.append(tmp_coord)

        in_vals = []
        out_vals = []

        for i in range(len(in_coords)):
            in_indices = in_map.coord_to_index_int(in_coords[i])
            out_indices = out_map.coord_to_index_int(out_coords[i])
            neg_indices = sum(n < 0 for n in in_indices + out_indices)
            if neg_indices > 0:
                in_vals.append(None)
                out_vals.append(None)
            else:
                try:
                    in_vals.append(round(in_map.data[in_indices], 5))
                except IndexError:
                    in_vals.append(None)
                try:
                    out_vals.append(round(out_map.data[out_indices], 5))
                except IndexError:
                    out_vals.append(None)

        if in_vals == out_vals and not all([v is None for v in in_vals]):
            return True
        else:
            print('IN', in_vals)
            print('OU', out_vals)
            return False

    @staticmethod
    def chop_soft_radius(model, in_map, out_map=None, hard_radius=3, soft_radius=2, mask_path=None):
        """
        Chop map using a soft mask with a given radius (hard_radius + soft_radius)
        around the atomic residue. A cosine function is used to create the soft mask.
        :param mask_path: mask output path
        :param model: biopython model object
        :param in_map: path to the input map or strudel map object
        :param out_map: out_map: output map path
        :param hard_radius: hard radius
        :param soft_radius: soft radius, cosine function
        :return strudel map object if out_map not given otherwise None
         """
        # Get atom coordinates
        atoms_coord = []
        for atom in model.get_atoms():
            atoms_coord.append(atom.coord)

        if type(in_map) is str:
            in_map = MapParser(in_map)
        elif type(in_map) is MapParser:
            pass
        else:
            raise TypeError('"in_map" should be a "MapObject" or in_dir to map file')

        voxel_size = in_map.voxel_size
        aver_voxel_size = sum(voxel_size) / 3
        delta1 = hard_radius
        delta2 = hard_radius + soft_radius

        # Create a numpy array for mask
        shape = in_map.data.shape
        mask = np.zeros(shape, dtype='float32')
        for coord in atoms_coord:
            x_index, y_index, z_index = in_map.coord_to_index_int(coord)

            rx = int(round(delta2 / voxel_size[0]))
            ry = int(round(delta2 / voxel_size[1]))
            rz = int(round(delta2 / voxel_size[2]))

            for x in range(x_index - rx, x_index + rx):
                for y in range(y_index - ry, y_index + ry):
                    for z in range(z_index - rz, z_index + rz):
                        # Calculate the distance between the current atom and the current voxel
                        d = aver_voxel_size * math.sqrt((x - x_index) ** 2 + (y - y_index) ** 2 + (z - z_index) ** 2)

                        # Assign mask values based to the distance to the atoms
                        if d <= delta1:
                            try:
                                mask[x, y, z] = 1
                            except IndexError:
                                pass
                        elif delta1 < d < delta2:
                            try:
                                mask[x, y, z] += (math.cos((math.pi / soft_radius) * (d - delta1)) + 1) / 2
                                # if intensity value became > 1 it is set to 1
                                if mask[x, y, z] > 1:
                                    mask[x, y, z] = 1
                            except IndexError:
                                pass

        # Apply the mask to the map data
        final = (mask * in_map.data)

        out_map_obj = MapParser('')
        out_map_obj.copy_header(in_map)
        out_map_obj.data = final

        if mask_path is not None:
            mask_ob = MapParser('')
            mask_ob.copy_header(in_map)
            mask_ob.data = mask
            mask_ob.write_map(mask_path)

        if out_map is not None:
            out_map_obj.write_map(out_map)
        else:
            return out_map_obj

    @staticmethod
    def find_near(atom, model, distance=6):
        close = []
        for atom1 in model.get_atoms():
            if atom1.parent.id != atom.parent.id or atom1.parent.parent.id != atom.parent.parent.id:
                d = atom1 - atom
                if d < distance:
                    close.append(atom1)
        return list(set(close))

    def chop_soft_radius_watershed(self, model, in_map, whole_model, shifts=None, out_map=None,
                                   radius=2, soft_radius=1, mask_path=None):
        """
        Chop map using a soft mask with a given hard_radius around the atomic residue.
        A cosine function is used to create the soft mask. Similar to chop_soft_radius but avoids
        cutting neighboring residues side chains map.
        :param model: biopython atomic residue object
        :param in_map: in_dir to the input map
        :param out_map: out_map: in_dir for the chopped map
        :param radius:
        :param soft_radius:
         """
        # Get atom coordinates
        if shifts is None:
            shifts = np.array([0, 0, 0])

        if type(in_map) is str:
            in_map = MapParser(in_map)
        elif type(in_map) is MapParser:
            pass
        else:
            raise TypeError('"in_map" should be a "MapObject" or in_dir to map file')

        voxel_size = in_map.voxel_size
        aver_voxel_size = sum(voxel_size) / 3

        # Create a numpy array for mask
        shape = in_map.data.shape
        mask = np.zeros(shape, dtype='float32')

        r = int(round((radius + soft_radius) / aver_voxel_size))

        for atom in model.get_atoms():
            xyz = in_map.coord_to_index(atom.coord - shifts)
            xyz_int = in_map.coord_to_index_int(atom.coord - shifts)

            near_atoms = []
            if atom.get_name() not in ['C', 'CA']:
                near_atoms = self.find_near(atom, whole_model, distance=(radius + soft_radius) * 2)

            for x in range(xyz_int[0] - r, xyz_int[0] + r):
                for y in range(xyz_int[1] - r, xyz_int[1] + r):
                    for z in range(xyz_int[2] - r, xyz_int[2] + r):

                        near_ds = [100]
                        for n_atom in near_atoms:
                            n_xyz = in_map.coord_to_index(n_atom.coord - shifts)
                            dn = aver_voxel_size * math.sqrt((x - n_xyz[0]) ** 2 + (y - n_xyz[1]) ** 2
                                                             + (z - n_xyz[2]) ** 2)
                            near_ds.append(dn)
                        dn = min(near_ds)

                        # Calculate the distance between the current atom and the current voxel
                        d = aver_voxel_size * math.sqrt((x - xyz[0]) ** 2 + (y - xyz[1]) ** 2 + (z - xyz[2]) ** 2)
                        if d > dn * 1.3:
                            continue
                        elif dn < radius + soft_radius:
                            delta2 = min((d + dn) * 0.65, radius + soft_radius)
                            delta1 = delta2 - soft_radius
                        else:
                            delta2 = radius + soft_radius
                            delta1 = radius
                        # Assign mask values based to the distance to the atoms
                        if d < delta1:
                            try:
                                mask[x, y, z] = 1
                            except IndexError:
                                pass
                        elif delta1 < d < delta2:
                            try:
                                mask[x, y, z] += (math.cos((math.pi / soft_radius) * (d - delta1)) + 1) / 2
                            except IndexError:
                                pass
        mask[mask > 1] = 1

        final = (mask * in_map.data)

        out_map_obj = MapParser('')
        out_map_obj.copy_header(in_map)
        out_map_obj.data = final

        if mask_path is not None:
            mask_ob = MapParser('')
            mask_ob.copy_header(in_map)
            mask_ob.data = mask
            mask_ob.write_map(mask_path)

        if out_map is not None:
            out_map_obj.write_map(out_map)
        else:
            return out_map_obj

    @staticmethod
    def chop_hard_radius(model, in_map, out_map, radius=3, mask_path=None):
        """
        Chop map using a hard mask with a given hard_radius around the atomic residue.
        :param model: biopython atomic residue object
        :param in_map: in_dir to the input map
        :param out_map: out_map: in_dir for the chopped map
        :param radius:
        """
        # Get atom coordinates
        atoms_coord = []
        for atom in model.get_atoms():
            atoms_coord.append(atom.coord)

        if type(in_map) is str:
            in_map = MapParser(in_map)
        elif type(in_map) is MapParser:
            pass
        else:
            raise TypeError('"in_map" should be a "MapObject" or in_dir to map file')

        voxel_size = in_map.voxel_size
        aver_voxel_size = sum(voxel_size) / 3
        shape = in_map.data.shape
        mask = np.zeros(shape, dtype='float32')

        for coord in atoms_coord:
            x_index, y_index, z_index = in_map.coord_to_index_int(coord)

            rx = int(round(radius / voxel_size[0]))
            ry = int(round(radius / voxel_size[0]))
            rz = int(round(radius / voxel_size[0]))
            for x in range(x_index - rx, x_index + rx):
                for y in range(y_index - ry, y_index + ry):
                    for z in range(z_index - rz, z_index + rz):
                        # Calculate the distance between the current atom and the current voxel
                        d = aver_voxel_size * math.sqrt((x - x_index) ** 2 + (y - y_index) ** 2
                                                        + (z - z_index) ** 2)
                        # Assign mask values based to the distance to the atoms
                        if d < radius:
                            mask[x, y, z] = 1

        # Apply the mask to the map data
        final = (mask * in_map.data)
        # Save the chopped map

        if mask_path is not None:
            mask_ob = MapParser('')
            mask_ob.copy_header(in_map)
            mask_ob.data = mask
            mask_ob.write_map(mask_path)

        out_map_obj = MapParser('')
        out_map_obj.copy_header(in_map)
        out_map_obj.data = final

        if out_map is not None:
            out_map_obj.write_map(out_map)
        else:
            return out_map_obj

    @staticmethod
    def grid_resample_emda(x, num):

        # Forward transform
        X = np.fft.fftn(x)
        X = np.fft.fftshift(X)

        # Placeholder array for output spectrum
        newshape = list(x.shape)
        newshape[0] = num[0]
        newshape[1] = num[1]
        newshape[2] = num[2]
        Y = np.zeros(newshape, X.dtype)

        # upsampling
        if X.shape[0] < newshape[0]:
            dx = abs(newshape[0] - X.shape[0]) // 2
            Y[dx:dx + X.shape[0],
            dx:dx + X.shape[0],
            dx:dx + X.shape[0]] = X
        # downsampling
        if newshape[0] < X.shape[0]:
            dx = abs(newshape[0] - X.shape[0]) // 2
            Y[:, :, :] = X[dx:dx + newshape[0],
                         dx:dx + newshape[0],
                         dx:dx + newshape[0]]

        Y = np.fft.ifftshift(Y)
        return (np.fft.ifftn(Y)).real
