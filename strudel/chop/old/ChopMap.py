"""
ChopMap.py

A collection of methods which allow to chop a map around an atomic residue.
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

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
import mrcfile
import numpy as np
import math
import subprocess
import os
from scipy.interpolate import RegularGridInterpolator
from strudel.core.Configure import Configure


class ChopMap:
    """
    Class for chopping map around an atomic residue using soft edge mask.
    """

    def __init__(self):
        config = Configure()
        self.eman2_path = config.eman2_package
        self.model = None
        self.in_map = ''
        self.out_cube = ''
        self.shift = []

    def sett_model_map(self, model, map):
        """
        Creates a biopython structure object
        :param model: Model file in_dir
        :param map: Map file in_dir
        """
        self.in_map = map
        if model.split('.')[-1] == 'pdb' or model.split('.')[-1] == 'ent':
            parser = PDBParser(PERMISSIVE=1)
            self.model = parser.get_structure(model.split('/')[-1].split('.')[0], model)
        elif model.split('.')[-1] == 'cif':
            parser = MMCIFParser()
            self.model = parser.get_structure(model.split('/')[-1].split('.')[0], model)
        else:
            raise Exception('Please provide the input residue in pdb or cif format')

    @staticmethod
    def interpolator(in_map):
        with mrcfile.open(in_map, mode='r+', permissive=True) as mrc:
            nx, ny, nz = mrc.data.shape
            x = range(nx)
            y = range(ny)
            z = range(nz)
            interpolator = RegularGridInterpolator((x, y, z), mrc.data)
        return interpolator

    @staticmethod
    def find_map_parameters(in_map):
        with mrcfile.open(in_map, mode='r+', permissive=True) as mrc:
            x_voxel_size = mrc.header.cella.x / mrc.header.nx
            y_voxel_size = mrc.header.cella.y / mrc.header.ny
            z_voxel_size = mrc.header.cella.z / mrc.header.nz
            nxstart = mrc.header.nxstart
            nystart = mrc.header.nystart
            nzstart = mrc.header.nzstart
            map_max = mrc.data.max()
            map_min = mrc.data.min()

        return (x_voxel_size, y_voxel_size, z_voxel_size), (nxstart, nystart, nzstart), (map_min, map_max)

    @staticmethod
    def atom_inclusion(model, in_map, level):
        """
        Counts the number of atoms within and outside map at a given level
        :param model: Biopython structure object
        :param in_map: map file
        :param level: map level
        :return: [included(nr), not_included(nr)]
        """
        with mrcfile.open(in_map, mode='r+', permissive=True) as mrc:
            inclusion = []
            x_voxel_size = mrc.header.cella.x / mrc.header.nx
            y_voxel_size = mrc.header.cella.y / mrc.header.ny
            z_voxel_size = mrc.header.cella.z / mrc.header.nz
            x = range(mrc.header.nx)
            y = range(mrc.header.ny)
            z = range(mrc.header.nz)
            a = RegularGridInterpolator((x, y, z), mrc.data)

            for atom in model.get_atoms():
                atom_coord = atom.coord
                x_index = atom_coord[2] / x_voxel_size - mrc.header.nxstart
                y_index = atom_coord[1] / y_voxel_size - mrc.header.nystart
                z_index = atom_coord[0] / z_voxel_size - mrc.header.nzstart
                if a([x_index, y_index, z_index]) > level:
                    inclusion.append(1)
                else:
                    inclusion.append(0)

        included = inclusion.count(1)
        not_included = inclusion.count(0)
        return included, not_included

    def find_threshold(self, inclusion, delta=2):
        interpolator = self.interpolator(self.in_map)
        voxel_size, nstart, lvl_range = self.find_map_parameters(self.in_map)
        upper = lvl_range[1]
        lower = 0
        current_lvl = (upper - lower) / 2
        while True:
            included = 0
            not_included = 0
            for atom in self.model.get_atoms():
                atom_coord = atom.coord
                x_index = atom_coord[2] / voxel_size[0] - nstart[0]
                y_index = atom_coord[1] / voxel_size[1] - nstart[1]
                z_index = atom_coord[0] / voxel_size[2] - nstart[2]
                if interpolator([x_index, y_index, z_index]) > current_lvl:
                    included += 1
                else:
                    not_included += 1
            current_incl = included / (included + not_included) * 100
            if current_incl < inclusion:
                upper = current_lvl
                current_lvl = current_lvl - (upper - lower) / 2
            elif current_incl > inclusion + delta:
                lower = current_lvl
                current_lvl = current_lvl + (upper - lower) / 2
            else:
                final_lvl = current_lvl
                break

        return final_lvl

    @staticmethod
    def chop_cube(model, in_map, out_map, cube_size):
        """
        Chop map using a cubic box around the residue
        :param model: biopython atomic residue object
        :param in_map: in_dir to the input map
        :param out_map: in_dir for the chopped map
        :param cube_size: the desired distance fom any cube edge to the nearest atom
        :return: translation matrix
        """
        atoms_coord = []
        for atom in model.get_atoms():
            atoms_coord.append(atom.coord)

        with mrcfile.mmap(in_map, mode='r+', permissive=True) as mrc:
            voxel_size = mrc.header.cella.x / mrc.header.nx
            print(mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart)
            delta = round(cube_size / voxel_size)
            # Get the indices in the map grid that correspond to the atom coordinates
            x_indices = []
            y_indices = []
            z_indices = []
            for atom in atoms_coord:
                x_index = int(round(atom[2] / voxel_size)) - mrc.header.nxstart
                x_indices.append(x_index)
                y_index = int(round(atom[1] / voxel_size)) - mrc.header.nystart
                y_indices.append(y_index)
                z_index = int(round(atom[0] / voxel_size)) - mrc.header.nzstart
                z_indices.append(z_index)
            # Find the voxel located in the middle of the atomic residue
            # and the maximum molecule size in grid points
            deltas = []
            minim = min(x_indices)
            dx = max(x_indices) - minim
            middle_x = int(dx / 2 + minim)
            deltas.append(dx)

            minim = min(y_indices)
            dy = max(y_indices) - minim
            middle_y = int(dy / 2 + minim)
            deltas.append(dy)

            minim = min(z_indices)
            dz = max(z_indices) - minim
            middle_z = int(dz / 2 + minim)
            deltas.append(dz)
            max_d = max(deltas)
            # Calculate the size of the cube
            radius = int(max_d / 2 + delta)
            new_dimension = radius * 2
            # Create a numpy array to store the chopped voxels
            new_grid = np.zeros((new_dimension, new_dimension, new_dimension), dtype='float32')

            # Assign voxel values
            for x in range(new_dimension):
                for y in range(new_dimension):
                    for z in range(new_dimension):
                        try:
                            new_grid[x, y, z] = mrc.data[x + middle_x - radius, y + middle_y - radius,
                                                         z + middle_z - radius]
                        except IndexError:
                            pass
            # Calculate the new cell size
            new_cell = round(new_dimension * voxel_size, 3)
            # Calculate the shifts applied to the chopped map
            x_shift = (middle_z - radius) * voxel_size + mrc.header.nzstart * voxel_size
            y_shift = (middle_y - radius) * voxel_size + mrc.header.nystart * voxel_size
            z_shift = (middle_x - radius) * voxel_size + mrc.header.nxstart * voxel_size
            shifts = [x_shift, y_shift, z_shift]
            # Create a file to store the new map
            with mrcfile.new(out_map, overwrite=True) as mr:
                # Assign data and header information
                mr.set_data(new_grid)
                mr.update_header_from_data()
                mr.header.cella = (new_cell, new_cell, new_cell)
                mr.header.origin = (0, 0, 0)
                mr.header.nxstart = 0
                mr.header.nystart = 0
                mr.header.nzstart = 0

        return shifts

    @staticmethod
    def chop_soft_radius(model, in_map, out_map, radius=3, soft_radius=2):
        """
        Chop map using a soft mask with a given radius around the atomic residue.
        A cosine function is used to create the soft mask.
        :param model: biopython atomic residue object
        :param in_map: in_dir to the input map
        :param out_map: out_map: in_dir for the chopped map
        :param radius:
        :param soft_radius:
         """
        # Get atom coordinates
        atoms_coord = []
        for atom in model.get_atoms():
            atoms_coord.append(atom.coord)

        with mrcfile.open(in_map, mode='r+', permissive=True) as mrc:
            voxel_size = mrc.header.cella.x / mrc.header.nx
            delta1 = radius
            delta2 = radius + soft_radius
            # Create a numpy array for mask
            mask = np.zeros((mrc.header.nx, mrc.header.ny, mrc.header.nz), dtype='float32')
            for atom in atoms_coord:

                # x and z coordinates are flipped in the map and coordinate files
                x_index = int(round(atom[2] / voxel_size)) - mrc.header.nxstart
                y_index = int(round(atom[1] / voxel_size)) - mrc.header.nxstart
                z_index = int(round(atom[0] / voxel_size)) - mrc.header.nxstart

                it = np.nditer(mask, flags=['multi_index'], op_flags=['writeonly'])
                while not it.finished:
                    # Calculate the distance between the current atom and the current voxel
                    d = voxel_size * math.sqrt((it.multi_index[0] - x_index) ** 2 + (it.multi_index[1] - y_index) ** 2
                                               + (it.multi_index[2] - z_index) ** 2)
                    # Assign mask values based to the distance to the atoms
                    if d < delta1:
                        it[0] = 1
                    if delta1 < d < delta2:
                        it[0] += (math.cos((math.pi / soft_radius) * (d - delta1)) + 1) / 2
                        # if intensity value became > 1 it is set to 1
                        if it[0] > 1:
                            it[0] = 1
                    it.iternext()
            # Apply the mask to the map data
            final = (mask * mrc.data)
            # Save the chopped map
            with mrcfile.new(out_map, overwrite=True) as mr:
                mr.set_data(final)
                mr.update_header_from_data()
                mr.header.cella = mrc.header.cella
                mr.header.origin = (0, 0, 0)

    @staticmethod
    def chop_hard_radius(model, in_map, out_map, radius=3):
        """
        Chop map using a hard mask with a given radius around the atomic residue.
        :param model: biopython atomic residue object
        :param in_map: in_dir to the input map
        :param out_map: out_map: in_dir for the chopped map
        :param radius:
        """
        # Get atom coordinates
        atoms_coord = []
        for atom in model.get_atoms():
            atoms_coord.append(atom.coord)

        with mrcfile.open(in_map, mode='r+', permissive=True) as mrc:
            voxel_size = mrc.header.cella.x / mrc.header.nx
            # Create a numpy array for mask
            mask = np.zeros((mrc.header.nx, mrc.header.ny, mrc.header.nz), dtype='float32')
            for atom in atoms_coord:
                # x and z coordinates are flipped in the map and coordinate files
                x_index = int(round(atom[2] / voxel_size)) - mrc.header.nxstart
                y_index = int(round(atom[1] / voxel_size)) - mrc.header.nxstart
                z_index = int(round(atom[0] / voxel_size)) - mrc.header.nxstart

                it = np.nditer(mask, flags=['multi_index'], op_flags=['writeonly'])
                while not it.finished:
                    # Calculate the distance between the current atom and the current voxel
                    d = voxel_size * math.sqrt((it.multi_index[0] - x_index) ** 2 + (it.multi_index[1] - y_index) ** 2
                                               + (it.multi_index[2] - z_index) ** 2)
                    # Assign mask values based to the distance to the atoms
                    if d < radius:
                        it[0] = 1
                    it.iternext()
            # Apply the mask to the map data
            final = (mask * mrc.data)
            # Save the chopped map
            with mrcfile.new(out_map, overwrite=True) as mr:
                mr.set_data(final)
                mr.update_header_from_data()
                mr.header.cella = mrc.header.cella
                mr.header.origin = (0, 0, 0)

    def grid_resample(self, in_map, out_map, desired_voxel):
        """
        Resample the map grid to achieve the desired voxel size. Fourier space interpolation using
        the e2proc3d.py script from EMAN2 package is applied. The image is also normalised.
        :param in_map: in_dir to the input map
        :param out_map: in_dir to the output map
        :param desired_voxel: the desired voxel size
        """
        map_extension = in_map.split('.')[-1]
        if map_extension != 'mrc':
            in_map = in_map.split('.')[0] + '.mrc'
        with mrcfile.open(in_map, mode='r', permissive=True) as mrc:
            voxel_size = mrc.header.cella.x / mrc.header.nx
        multiplier = round(desired_voxel / voxel_size, 3)
        # source_eman = 'export PATH="' + self.eman2_path + 'bin:$PATH" \n'
        # command = 'bash ' + source_eman + 'e2proc3d.py ' + in_map + ' ' + out_map + ' --process=math.fft.resample:n=' + str(multiplier) \
        #          + ' --process=normalize.edgemean'
        self.eman2_path = self.eman2_path + 'bin/'
        command = self.eman2_path + 'e2proc3d.py ' + in_map + ' ' + out_map + ' --process=math.fft.resample:n=' + \
                  str(multiplier) # + ' --process=normalize.edgemean'
        subprocess.call(command, shell=True)
        #os.remove('.eman2log.txt')

    @staticmethod
    def shift_coord(trans_matrix, structure):
        """
        Apply translation matrix to residue
        :param trans_matrix: translation matrix
        :param structure: biopython structure object
        """
        try:
            for residue in structure.get_residues():
                for atom in residue:
                    atom.get_coord()[0] = atom.get_coord()[0] - trans_matrix[0]
                    atom.get_coord()[1] = atom.get_coord()[1] - trans_matrix[1]
                    atom.get_coord()[2] = atom.get_coord()[2] - trans_matrix[2]
        except AttributeError:
            for atom in structure:
                atom.get_coord()[0] = atom.get_coord()[0] - trans_matrix[0]
                atom.get_coord()[1] = atom.get_coord()[1] - trans_matrix[1]
                atom.get_coord()[2] = atom.get_coord()[2] - trans_matrix[2]

    @staticmethod
    def save_pdb(model, out_path):
        """
        Save residue as pdb
        :param model: structure object
        :param out_path: output pdb file in_dir
        :return:
        """
        io = PDBIO()
        io.set_structure(model)
        io.save(out_path)


def main():
    # This is an example how to chop a map
    map = "../test_dir/test_mult/emd_7025.map"
    model_file = "../test_dir/test_mult/arg-115-B-7025_residue.pdb"
    cube_pdb = "../test_dir/test_mult/arg-115-B-7025_residue_side_cube.pdb"
    cube = "../test_dir/test_mult/arg-115-B-7025_residue_side_cube.mrc"
    cube_newgrid = "../test_dir/test_mult/arg-115-B-7025_residue_side_cube05a.mrc"
    soft_map = "../test_dir/test_mult/arg-115-B-7025_residue_side_cube05a_soft.mrc"
    hard_map = "../test_dir/test_mult/arg-115-B-7025_residue_side_cube05a_hard.mrc"

    chop = ChopMap()
    chop.sett_model_map(model_file, map)
    str = chop.model
    matrix = chop.chop_cube(str, chop.in_map, cube, 6)
    chop.shift_coord(matrix, str)
    chop.save_pdb(str, cube_pdb)
    chop.grid_resample(cube, cube_newgrid, 0.25)
    chop.chop_soft_radius(str, cube_newgrid, soft_map, 3, 2)
    chop.chop_hard_radius(str, cube_newgrid, hard_map)


if __name__ == '__main__':
    main()
