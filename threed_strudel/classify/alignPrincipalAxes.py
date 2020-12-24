"""
alignPrincipalAxes.py

Generates rotation-translation matrices (Chimera format) for densities in mrc format
which if loaded in Chimera superimpose the density principal axes of inertia with
the x, y and z axes.

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

import mrcfile
import numpy as np
import os
import json
import argparse


def calc_principal_axis(xyz):
    """
    Calculates the principal axes of inertia of a density. Intensities are
    treated as masses. Theoretical basis: "Classical Mechanics" by J. B. Tatum
    :param xyz: density intensities as a numpy array
    :return: Coordinates of the center of the mass (Mc), principal moments of
    inertia (e_values), principal axes of inertia (e_vectors).
    """
    # Intensities lower than 0.2*maximum are set to 0 to reduce the noise and to
    # exclude negative values.
    max_int = np.amax(xyz)
    it = np.nditer(xyz, flags=['multi_index'], op_flags=['writeonly'])
    while not it.finished:
        if it[0] < max_int*0.2:
            it[0] = 0
        it.iternext()
    # Find coordinates of the centre of the mass, moments and products of inertia
    mass = np.sum(xyz)
    summ_xc = 0
    summ_yc = 0
    summ_zc = 0
    A = 0
    B = 0
    C = 0
    F = 0
    G = 0
    H = 0
    it = np.nditer(xyz, flags=['multi_index'], op_flags=['writeonly'])
    while not it.finished:
        summ_xc += it.multi_index[2] * it[0]
        summ_yc += it.multi_index[1] * it[0]
        summ_zc += it.multi_index[0] * it[0]
        # Find the moments and products of inertia relative to origin (0,0,0)
        A += (it.multi_index[1] ** 2 + it.multi_index[0] ** 2) * it[0]
        B += (it.multi_index[2] ** 2 + it.multi_index[0] ** 2) * it[0]
        C += (it.multi_index[2] ** 2 + it.multi_index[1] ** 2) * it[0]
        F += it.multi_index[1] * it.multi_index[0] * it[0]
        G += it.multi_index[2] * it.multi_index[0] * it[0]
        H += it.multi_index[2] * it.multi_index[1] * it[0]
        it.iternext()
    xc = summ_xc / mass
    yc = summ_yc / mass
    zc = summ_zc / mass
    # Coordinates of the center of the mass
    Mc = np.array([xc, yc, zc])
    # Find the moments and products of inertia relative to origin the center of the mass
    A = A - mass * (yc**2 + zc**2)
    B = B - mass * (xc**2 + zc**2)
    C = C - mass * (xc**2 + yc**2)
    # Products of inertia
    F = F - mass * yc * zc
    G = G - mass * xc * zc
    H = H - mass * xc * yc
    # Find principal moments and axes of inertia
    I = np.array([[A, -H, -G], [-H, B, -F], [-G, -F, C]])
    e_values, e_vectors = np.linalg.eig(I)
    # Sort axes based on principal moments values
    idx = e_values.argsort()[::1]
    e_values = e_values[idx]
    e_vectors = e_vectors[:, idx]
    # Make sure that a right handed set was obtained
    det_e_v = np.linalg.det(e_vectors)
    if det_e_v < 0:
        for i in range(3):
            e_vectors[i,0] = -e_vectors[i,0]
    return Mc, e_values, e_vectors


def calc_rot_trans(Mc, e_vectors, voxel_size, flip=None):
    """
    Calculates rotation-translation matrices in Chimera format. The principal axes
    of inertia are aligned to x, y and z axes of the original box. The center of
    the mass is aligned to the origin of the coordinate system. Since the directions
    of the axes ere arbitrary chosen the axes can be flipped.
    :param Mc: Coordinates of the center of the mass
    :param e_vectors: sorted principal axes of inertia
    :param voxel_size: density voxel size
    :param flip: axes to be flipped (None, x, y or xy)
    :return: rotation-translation matrix
    """
    flip_x = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
    flip_y = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
    if flip is None:
        rot = np.transpose(e_vectors)
        trans = ([0,0,0] - np.dot(Mc, e_vectors)) * voxel_size
        rot_trans = np.hstack((rot, trans.reshape((3,1))))
    elif flip == "x":
        e_vectors = np.dot(e_vectors, flip_x)
        rot = np.transpose(e_vectors)
        trans = ([0, 0, 0] - np.dot(Mc, e_vectors)) * voxel_size
        rot_trans = np.hstack((rot, trans.reshape((3, 1))))
    elif flip == "y":
        e_vectors = np.dot(e_vectors, flip_y)
        rot = np.transpose(e_vectors)
        trans = ([0, 0, 0] - np.dot(Mc, e_vectors)) * voxel_size
        rot_trans = np.hstack((rot, trans.reshape((3, 1))))
    elif flip == "xy":
        e_vectors = np.dot(e_vectors, flip_x)
        e_vectors = np.dot(e_vectors, flip_y)
        rot = np.transpose(e_vectors)
        trans = ([0, 0, 0] - np.dot(Mc, e_vectors)) * voxel_size
        rot_trans = np.hstack((rot, trans.reshape((3, 1))))
    else:
        raise ValueError('Please input "x", "y", or "xy" or skip the flip parameter.')

    rot_trans = np.around(rot_trans, decimals=6)
    return rot_trans


def main():
    parser = argparse.ArgumentParser(description='Generate rotation-translation matrices for '
                                                 'principal component density alignment.')
    parser.add_argument("-f", "--file", dest="mrc_path", required=True,
                        help="Path a directory with densities in mrc format")
    args = parser.parse_args()
    work_dir = args.mrc_path

    files = os.listdir(work_dir)
    mrc_files = [i for i in files if i.endswith('.mrc')]
    length = len(mrc_files)
    moment_list = []
    norm_moment_list = []
    for k, element in enumerate(mrc_files):
        print('Working on', k + 1, 'file out of', length + 1, 'files.      "', element, '"')
        map1 = os.path.join(work_dir, element)
        with mrcfile.open(map1, mode='r+', permissive=True) as mrc:
            xyz = mrc.data
            voxel_size1 = mrc.header.cella.x / mrc.header.nx
        Mc1, e_val1, e_vec1 = calc_principal_axis(xyz)
        e_val1_norm = e_val1 / np.linalg.norm(e_val1)
        moment_list.append([element] + e_val1.tolist())
        norm_moment_list.append([element] + e_val1_norm.tolist())

        rot_trans1 = calc_rot_trans(Mc1, e_vec1, voxel_size1)
        rot_trans2 = calc_rot_trans(Mc1, e_vec1, voxel_size1, 'x')
        rot_trans3 = calc_rot_trans(Mc1, e_vec1, voxel_size1, 'y')
        rot_trans4 = calc_rot_trans(Mc1, e_vec1, voxel_size1, 'xy')
        matrix_list = [rot_trans1.tolist(), rot_trans2.tolist(), rot_trans3.tolist(), rot_trans4.tolist()]
        j_file = os.path.join(work_dir, element.split('.')[0] + '.json')
        # Save matrices as json files
        with open(j_file, 'w')as j:
            json.dump(matrix_list, j, indent=4)
        matrix_file = os.path.join(work_dir, element.split('.')[0] + '.txt')
        # Save matrices as text file loadable by Chimera
        with open(matrix_file, 'w') as fo:
            for i, matrix in enumerate(matrix_list):
                fo.write('Model ' + str(i) + '.0\n')
                for line in matrix:
                    fo.write(str(line[0]) + '  ' + str(line[1]) + '  ' + str(line[2]) + '  ' + str(line[3]) + '\n')

    moment_file = os.path.join(work_dir, 'a_principal_moments.txt')
    with open(moment_file, 'w') as m_f:
        for line in moment_list:
            text = ''
            for k in range(4):
                if k > 0:
                    line[k] = '{:.6f}'.format(round(line[k], 6))
                    text = text + '{0:>{width}}'.format(line[k], width=18)
                else:
                    text = text + '{0:>{width}}'.format(line[k], width=25)
            m_f.write(text + '\n')

    norm_moment_file = os.path.join(work_dir, 'a_principal_moments_norm.txt')
    with open(norm_moment_file, 'w') as m_f:
        for line in norm_moment_list:
            text = ''
            for k in range(4):
                if k > 0:
                    line[k] = '{:.6f}'.format(round(line[k], 6))
                    text = text + '{0:>{width}}'.format(line[k], width=18)
                else:
                    text = text + '{0:>{width}}'.format(line[k], width=25)
            m_f.write(text + '\n')


if __name__ == '__main__':
    main()