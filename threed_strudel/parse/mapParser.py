"""
mapParser.py

Creates a Map object. This module is built on top of mrcfile library
By default it will transpose the data array so the axis order is x, y, z

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

import os
import mrcfile
import math
import numpy as np


class MapParser:
    def __init__(self, map_path=''):
        self.id = ''
        self.map_path = map_path
        self.data = None
        self.voxel_size = None
        self.n_start = (0, 0, 0)
        self.m = None
        self.cell = None
        self.cellb = (90, 90, 90)
        self.n_xyz = None
        self.origin = (0, 0, 0)
        self.mapc = 1
        self.mapr = 2
        self.maps = 3
        self.info = ''
        self.fraction_matrix = None
        self.index_shifts = (0, 0, 0)
        self.coord_to_index = self._coord_to_index_orthogonal
        self.coord_to_index_int = self._coord_to_index_orthogonal_int
        if os.path.exists(map_path):
            self._load_map()

    def _load_map(self, default_axis=True):
        with mrcfile.open(self.map_path, mode='r+', permissive=True) as mrc:
            self.data = mrc.data
            self.n_start = self._get_n_start(mrc)
            self.m = self._get_m(mrc)
            self.cell = (mrc.header.cella.x, mrc.header.cella.y, mrc.header.cella.z)
            self.cellb = (mrc.header.cellb.alpha, mrc.header.cellb.beta, mrc.header.cellb.gamma)
            self.n_xyz = (mrc.header.nx, mrc.header.ny, mrc.header.nz)
            self.origin = (mrc.header.origin.x, mrc.header.origin.y, mrc.header.origin.z)
            self.id = os.path.basename(self.map_path).split('.')[0]
            self._get_axis_order(mrc)
            if default_axis:
                self._swap_axes()
            self.voxel_size = self._get_voxel_size(mrc)
            self.fraction_matrix = self.calc_fraction_matrix()
            self.index_shifts = self._calc_index_shift()
            if any([z for z in map(lambda k: k != 90, self.cellb)]):
                self.coord_to_index = self._coord_to_index_nonorthogonal
                self.coord_to_index_int = self._coord_to_index_nonorthogonal_int

    def _swap_axes(self):
        crs = (self.mapc, self.mapr, self.maps)
        if crs != (1, 2, 3):
            self.info += 'The input map has non default axis order: '
            order = {1: 'x', 2: 'y', 3: 'z'}
            self.info += 'Changing {} {} {} axis order to x y z'.format(order[crs[0]], order[crs[1]], order[crs[2]])
            # Transpose the data array
            xyz_indices = (crs.index(1), crs.index(2), crs.index(3))
            new_order = [2 - xyz_indices[2 - a] for a in (0, 1, 2)]
            self.data = self.data.transpose(new_order)
            # Rearrange n_start values
            self.n_start = (self.n_start[xyz_indices[0]], self.n_start[xyz_indices[1]], self.n_start[xyz_indices[2]])

            self.mapc = 1
            self.mapr = 2
            self.maps = 3

    @staticmethod
    def _get_voxel_size(mrcfile_object):
        voxel_size = (mrcfile_object.header.cella.x / mrcfile_object.header.mx,
                      mrcfile_object.header.cella.y / mrcfile_object.header.my,
                      mrcfile_object.header.cella.z / mrcfile_object.header.mz)
        return voxel_size

    @staticmethod
    def _get_n_start(mrcfile_object):
        n_start = [int(mrcfile_object.header.nxstart),
                   int(mrcfile_object.header.nystart),
                   int(mrcfile_object.header.nzstart)]
        return n_start

    @staticmethod
    def _get_m(mrcfile_object):
        m = [int(mrcfile_object.header.mx),
             int(mrcfile_object.header.my),
             int(mrcfile_object.header.mz)]
        return m

    def _get_axis_order(self, mrcfile_object):
        self.mapc = int(mrcfile_object.header.mapc)
        self.mapr = int(mrcfile_object.header.mapr)
        self.maps = int(mrcfile_object.header.maps)

    def write_map(self, out_path):
        with mrcfile.new(out_path, overwrite=True) as mr:
            # Assign data and header information
            mr.set_data(self.data)
            mr.update_header_from_data()
            mr.header.cella = self.cell
            mr.header.cellb = self.cellb
            mr.header.nxstart = self.n_start[0]
            mr.header.nystart = self.n_start[1]
            mr.header.nzstart = self.n_start[2]
            try:
                mr.header.mx = self.m[0]
                mr.header.my = self.m[1]
                mr.header.mz = self.m[2]
            except TypeError:
                pass
            mr.header.origin.x = self.origin[0]
            mr.header.origin.y = self.origin[1]
            mr.header.origin.z = self.origin[2]
            mr.header.mapc = self.mapc
            mr.header.mapr = self.mapr
            mr.header.maps = self.maps

    def copy_header(self, map_obj):
        self.voxel_size = map_obj.voxel_size
        self.n_start = map_obj.n_start
        self.cell = map_obj.cell
        self.n_xyz = map_obj.n_xyz
        self.origin = map_obj.origin
        self.cellb = map_obj.cellb

    def _coord_to_index_orthogonal(self, coord):
        """
        Converts atomic residue coordinates into map indices
        :param coord: coordinate in angstrom
        :return: float indices
        """
        x_index = coord[2] / self.voxel_size[2] - self.index_shifts[0]
        y_index = coord[1] / self.voxel_size[1] - self.index_shifts[1]
        z_index = coord[0] / self.voxel_size[0] - self.index_shifts[2]
        return x_index, y_index, z_index

    def _coord_to_index_orthogonal_int(self, coord):
        """
        Converts atomic residue coordinates into map indices
        :param coord: coordinate in angstrom
        :return: int indices
        """
        x, y, z = self._coord_to_index_orthogonal(coord)
        return int(round(x)), int(round(y)), int(round(z))

    def _coord_to_index_nonorthogonal(self, coord):
        coord = np.flip(coord)
        return self.fraction_matrix.dot(coord)*self.n_xyz - self.index_shifts

    def _coord_to_index_nonorthogonal_int(self, coord):
        x, y, z = self._coord_to_index_nonorthogonal(coord)
        return int(round(x)), int(round(y)), int(round(z))

    def calc_fraction_matrix(self):
        """
        Calculates matrix used to transform Cartesian
        coordinates in the ATOM_SITE category to fractional coordinates
        in the same category
        :return: fractional matrix
        """
        cell = self.cell
        ang = self.cellb
        ang = (ang[0] * math.pi / 180, ang[1] * math.pi / 180, ang[2] * math.pi / 180)

        omega = cell[0] * cell[1] * cell[2] * math.sqrt(
            1 - math.cos(ang[0]) ** 2 - math.cos(ang[1]) ** 2 - math.cos(ang[2]) ** 2 +
            2 * math.cos(ang[0]) * math.cos(ang[1]) * math.cos(ang[0]))
        m11 = 1 / cell[0]
        m12 = - math.cos(ang[2]) / (cell[0] * math.sin(ang[2]))
        m13 = cell[1] * cell[2] * (math.cos(ang[0]) * math.cos(ang[2]) - math.cos(ang[1])) / (omega * math.sin(ang[2]))

        m21 = 0
        m22 = 1 / (cell[1] * math.sin(ang[2]))
        m23 = cell[0] * cell[2] * (math.cos(ang[1]) * math.cos(ang[2]) - math.cos(ang[0])) / (omega * math.sin(ang[2]))

        m31 = 0
        m32 = 0
        m33 = cell[0] * cell[1] * math.sin(ang[2]) / omega
        matrix = [m11, m12, m13, m21, m22, m23, m31, m32, m33]

        return np.array(matrix).reshape(3, 3)

    def _calc_index_shift_(self):
        if self.origin == (0, 0, 0):
            x = self.n_start[0]
            y = self.n_start[1]
            z = self.n_start[2]
        elif self.n_start != (0, 0, 0):
            x = self.origin[2] / self.voxel_size[2]
            y = self.origin[1] / self.voxel_size[1]
            z = self.origin[0] / self.voxel_size[0]
        else:
            x, y, z, = 0, 0, 0
        return x, y, z

    def _calc_index_shift(self):
        if self.origin == (0, 0, 0):
            x = self.n_start[0]
            y = self.n_start[1]
            z = self.n_start[2]
        else: #self.n_start != (0, 0, 0):
            x = self.origin[2] / self.voxel_size[2]
            y = self.origin[1] / self.voxel_size[1]
            z = self.origin[0] / self.voxel_size[0]
        # else:
        #     x, y, z, = 0, 0, 0
        return x, y, z

    def set_origin(self, origin):
        self.origin = origin
        self.index_shifts = self._calc_index_shift()

    def grid_resample_emda(self, target_voxel):
        x = self.data

        if x.shape[0] % 2 != 0:
            xshape = list(x.shape)
            xshape[0] = xshape[0] + 1
            xshape[1] = xshape[1] + 1
            xshape[2] = xshape[2] + 1
            temp = np.zeros(xshape, x.dtype)
            temp[:-1, :-1, :-1] = x
            x = temp
            self.cell = (self.cell[0] + self.voxel_size[0], self.cell[1] + self.voxel_size[1], self.cell[2] + self.voxel_size[2])

        scale = self.voxel_size[0] / target_voxel
        old_grid = x.data.shape
        new_grid = [int(round(i * scale)) for i in old_grid]

        # Forward transform
        X = np.fft.fftn(x)
        X = np.fft.fftshift(X)
        # Placeholder array for output spectrum
        newshape = list(x.shape)

        if new_grid[0] % 2 != 0:
            newshape[0] = new_grid[0] + 1
            newshape[1] = new_grid[1] + 1
            newshape[2] = new_grid[2] + 1
        else:
            newshape[0] = new_grid[0]
            newshape[1] = new_grid[1]
            newshape[2] = new_grid[2]

        if x.shape[0] == newshape[0]:
            return x
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

        self.data = np.float32((np.fft.ifftn(Y)).real)
        self.m = self.data.shape
        self.voxel_size = (self.cell[0]/self.m[0], self.cell[1]/self.m[1], self.cell[2]/self.m[2])
        self.index_shifts = self._calc_index_shift()

