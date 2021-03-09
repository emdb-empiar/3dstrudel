
from random import randint

import threed_strudel.utils.bio_utils
# from strudel.chop.chopMap import ChopMap, MapParser
from threed_strudel.parse.map_parser import MapParser
import math
import numpy as np
import time
from numpy import linalg as LA
import os
from scipy.interpolate import InterpolatedUnivariateSpline
from threed_strudel.utils import model_map_utils


class ExtractMap:
    def __init__(self, map_path=None, model_path=None, adapted_map_path=None, work_voxel_size=0.5, tmp_dir='/tmp/segment_density'):
        self.in_map = MapParser(map_path)
        self.cube_map = None
        self.box_size = 10
        self.max_thld = None
        self.min_thld = None
        self.in_model = threed_strudel.utils.bio_utils.load_structure(model_path)
        self.walk_grid = None
        self.operations_counter = 0
        self.max_roll = 5
        self.work_voxel_size = work_voxel_size
        self.interpolator = None
        self.tmp_dir = tmp_dir
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        if adapted_map_path is not None:
            self.adapted_map = MapParser(adapted_map_path)
        else:
            self.adapted_map = None

    @staticmethod
    def cut_cube_around_voxel(map_obj, center, radius):
        radius = int(round(radius / map_obj.voxel_size[0]))
        size = radius * 2 + 1
        new_grid = np.zeros((size, size, size), dtype='float32')
        # Assign voxel values
        for x in range(size):
            for y in range(size):
                for z in range(size):
                    try:
                        new_grid[x, y, z] = map_obj.data[x + center[0] - radius,
                                                         y + center[1] - radius,
                                                         z + center[2] - radius]
                    except IndexError:
                        pass

        out_map = MapParser('')
        out_map.data = new_grid
        out_map.cell = (round(size * map_obj.voxel_size[0], 3),
                        round(size * map_obj.voxel_size[1], 3),
                        round(size * map_obj.voxel_size[2], 3))
        out_map.voxel_size = map_obj.voxel_size
        out_map.n_xyz = new_grid.shape
        n_st = map_obj.n_start
        out_map.n_start = (0, 0, 0)
        out_map.origin = ((center[2] - radius + n_st[0]) * out_map.voxel_size[0] + map_obj.origin[0],
                          (center[1] - radius + n_st[1]) * out_map.voxel_size[1] + map_obj.origin[1],
                          (center[0] - radius + n_st[2]) * out_map.voxel_size[2] + map_obj.origin[2])
        return out_map

    # def map_resample(self, map_obj, file_prefix=''):
    #     """
    #     Resample the grid in the a map object so that the
    #     voxel size become equal to 'self.work_voxel_size'
    #     :param map_obj: map object
    #     :param file_prefix: prefix for the tmp file
    #     :return: map object
    #     """
    #     tmp = os.path.join(self.tmp_dir, file_prefix + 'tmp.mrc')
    #     res_tmp = os.path.join(self.tmp_dir, file_prefix + 'res_tmp.mrc')
    #
    #
    #     map_obj.write_map(tmp)
    #     chop = ChopMap()
    #     chop.grid_resample(tmp, res_tmp, self.work_voxel_size, map_obj.voxel_size[0])
    #     out_map_obj = MapParser(res_tmp)
    #
    #     return out_map_obj

    def update_map_parameters(self, map_obj):
        self.max_thld = np.amax(map_obj.data)
        self.min_thld = np.amin(map_obj.data)

        self.walk_grid = np.zeros(map_obj.data.shape, dtype='float32')

    def get_flank_ca_coord(self, res):
        try:
            prev_ca = self.in_model[0][res.parent.id][res.id[1]-1]['CA'].coord
        except KeyError:
            prev_ca = None
        try:
            next_ca = self.in_model[0][res.parent.id][res.id[1]-1]['CA'].coord
        except KeyError:
            next_ca = None

        return prev_ca, next_ca

    def segment_residue_density(self, model):

        c_coord, ca_coord, n_coord, cb_coord = self._get_bb_coordinates(model)

        ca_index = self.coord_to_index_int(ca_coord, self.in_map)
        ca_map = self.cut_cube_around_voxel(self.in_map, ca_index, self.box_size)
        # print(ca_map.data.shape)
        res_nr = model.id[1]
        prev_ca_coord, next_ca_coord = self.get_flank_ca_coord(model)

        # self.cube_map = self.map_resample(ca_map, str(res_nr))
        ca_map.grid_resample_emda(self.work_voxel_size)
        self.cube_map = ca_map
        # print(self.cube_map.data.shape)

        self.update_map_parameters(self.cube_map)

        if self.adapted_map is not None:
            adapted_ca_map = self.cut_cube_around_voxel(self.adapted_map, ca_index, self.box_size)
            adapted_cube_map = adapted_ca_map.grid_resample_emda(self.work_voxel_size)
            self.interpolator = model_map_utils.interpolator(adapted_cube_map)
        else:
            self.interpolator = model_map_utils.interpolator(self.cube_map)



        start = time.time()
        if cb_coord != []:
            start_coord = cb_coord
        else:
            start_coord = ca_coord
        #start_coord = ca_coord
        print(start_coord)
        side_chain_trace = self.trace_side_chain(start_coord, c_coord, n_coord,
                                                 prev_ca_coord, next_ca_coord,
                                                 self.interpolator, 1.6)

        if len(side_chain_trace) > 0:
            ang = self.path_bb_angle(side_chain_trace, c_coord, n_coord)
            # print(ang)
            if ang < math.pi / 4:
                if cb_coord != []:
                    side_chain_trace = [cb_coord]
                else:
                    side_chain_trace = []




        # print('Trace', time.time() - start)

        #if len(sidechain_trace) > 1:
        #    sidechain_trace = self.exclude_gaps(sidechain_trace)
        residue_trace = [c_coord, ca_coord] + side_chain_trace
        # print(residue_trace)
        residue_trace = self.connect_trace(residue_trace)

        trace_indices = self._coord_lst_to_indices_lst_int(residue_trace, self.cube_map)
        trace_map = self._points_to_map(trace_indices, self.cube_map, 3)

        # print(self.operations_counter)
        start = time.time()
        mask = self.mask_trace(self.cube_map, trace_indices, 1.5, 1)
        # print('Mask', time.time() - start)

        indices = self.coord_to_index_int(c_coord, self.cube_map)
        self.walk_grid[indices] = 1.7
        indices = self.coord_to_index_int(n_coord, self.cube_map)
        self.walk_grid[indices] = 1.7
        indices = self.coord_to_index_int(ca_coord, self.cube_map)
        self.walk_grid[indices] = 1.7

        walk_map = MapParser('')
        walk_map.copy_header(self.cube_map)
        walk_map.data = self.walk_grid

        residue_map = MapParser('')
        residue_map.copy_header(self.cube_map)
        residue_map.data = self.cube_map.data * mask

        return self.cube_map, trace_map, walk_map, residue_map

    @staticmethod
    def unit_vector(v):
        return v / np.linalg.norm(v)

    def path_bb_angle(self, path, c_coord, n_coord):
        aver = sum(path) / len(path)
        v1 = c_coord - n_coord
        v2 = c_coord - aver
        v3 = n_coord - aver

        v1 = self.unit_vector(v1)
        v2 = self.unit_vector(v2)
        v3 = self.unit_vector(v3)

        ang1 = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
        if ang1 > math.pi / 2:
            ang1 = math.pi - ang1
        ang2 = np.arccos(np.clip(np.dot(v1, v3), -1.0, 1.0))
        if ang2 > math.pi / 2:
            ang2 = math.pi - ang2
        return min([ang1, ang2])


    @staticmethod
    def _get_bb_coordinates(model):
        """
        Extracts the C, CA, N and CB coordinates.
        :param model: Biopython residue object
        :return: list of coordinates
        """
        ca = []
        cb = []
        c = []
        n = []
        for atom in model.get_atoms():
            if atom.get_name() == 'CA':
                ca = np.array(atom.coord)
            if atom.get_name() == 'C':
                c = np.array(atom.coord)
            if atom.get_name() == 'N':
                n = np.array(atom.coord)
            if atom.get_name() == 'CB':
                cb = np.array(atom.coord)
        return c, ca, n, cb

    def _coord_lst_to_indices_lst_int(self, coord_lst, map_obj):
        """
        Converts a list of coordinates (in angstroms) to a list of map integer indices
        :param coord_lst: list of coordinates
        :param map_obj: map object
        :return: list of map indices
        """
        indices_lst = []
        for coordinate in coord_lst:
            try:
                indices = self.coord_to_index_int(coordinate, map_obj)
            except IndexError:
                pass
            else:
                indices_lst.append(indices)
        return indices_lst

    @staticmethod
    def _points_to_map(indices_lst, map_template, value):
        """
        Creates a map object of the size of the 'map_template'. Sets the
        map values at 'indices_lst' positions to 'value', all other map values are 0.
        :param indices_lst: list of map indices
        :param map_template: map object
        :return: map object
        """
        out_map = MapParser('')
        grid = np.zeros(map_template.data.shape, dtype='float32')
        i = 0.1
        for indices in indices_lst:
            grid[indices] = value - i
            i += 0.1
        out_map.copy_header(map_template)
        out_map.data = grid
        return out_map

    def mask_trace(self, map_obj, trace, radius, soft_radius):
        mask = np.zeros(map_obj.data.shape, dtype='float32')
        delta2 = radius + soft_radius
        delta1 = radius
        voxel_size = map_obj.voxel_size[0]
        for point in trace:

            rx = int(round(delta2 / voxel_size))
            ry = int(round(delta2 / voxel_size))
            rz = int(round(delta2 / voxel_size))
            for x in range(point[0] - rx, point[0] + rx):
                for y in range(point[1] - ry, point[1] + ry):
                    for z in range(point[2] - rz, point[2] + rz):
                        # Calculate the distance between the current atom and the current voxel
                        d = voxel_size * math.sqrt((x - point[0]) ** 2 + (y - point[1]) ** 2 + (z - point[2]) ** 2)

                        # Assign mask values based to the distance to the atoms
                        try:
                            if d < delta1:
                                mask[x, y, z] = 1
                            elif delta1 < d < delta2:
                                mask[x, y, z] += (math.cos((math.pi / soft_radius) * (d - delta1)) + 1) / 2
                                # if intensity value became > 1 it is set to 1
                                if mask[x, y, z] > 1:
                                    mask[x, y, z] = 1
                        except IndexError:
                            pass
        return mask

    def connect_trace(self, trace, max_dist=0.8):
        """

        :param trace:
        :param max_dist:
        :return:
        """
        connected_trace = trace[:1]
        for index, point in enumerate(trace[1:]):
            point0 = trace[index-1]
            d = self.distance(point0, point)
            if d > max_dist:
                connection_points = int(round(d / max_dist))
                delta = (point - point0) / connection_points
                for i in range(connection_points):
                    new_point = point0 + delta * i
                    connected_trace.append(new_point)
            connected_trace.append(point)
        return connected_trace

    def exclude_gaps(self, trace):
        """

        :param trace:
        :return:
        """
        connected_trace = [trace[0]]
        for index, point in enumerate(trace[1:]):
            point0 = trace[index]
            d = self.distance(point0, point)
            if d > 2:
                deep_min = self.check_minimum(point0, point)
                if deep_min:
                    break
                else:
                    connected_trace.append(point)
            else:
                connected_trace.append(point)
        return connected_trace

    def connect_trace_(self, trace, max_dist=0.8):
        # print(trace)
        connected_trace = [trace[0]]
        for index, point in enumerate(trace[1:]):
            point0 = trace[index-1]
            d = self.distance(point0, point)
            # print('D', d)
            if d > max_dist:
                deep_min = self.check_minimum(point0, point)
                if deep_min:
                    # print('DEEP MIN')
                    break
                connection_points = int(round(d / max_dist))
                delta = (point - point0) / connection_points
                for i in range(connection_points):
                    new_point = point0 + delta * i
                    connected_trace.append(new_point)
                connected_trace.append(point)
            else:
                connected_trace.append(point)
        return connected_trace

    def check_minimum(self, p1, p2, depth=0.2):
        """
        Detects regions where the density of different residues might be connected.
        First takes two points in a map, connects them with a line and calculates the map values along the line. Then
        interpolates the map values using cubic spline and finds the critical points and values on the p1-p2 segment.
        For each minimum checks the flanking maximums and if the average maximum is more than 'depth' higher than
        the minimum returns True else False.
        :param p1: point (np array)
        :param p2: point (np array)
        :param depth: min max relative difference
        :return: True or False
        """
        deep_min = False
        interim_values = []
        v = p2 - p1
        for i in range(20):
            coord = p1 + i * v/20
            indices = self.coord_to_index(coord, self.cube_map)
            val = self.interpolator(indices)
            interim_values.append(val)
        cr_pts, cr_vals = self.critical_points([i for i in range(20)], interim_values)
        min_index = np.argmin(cr_vals)

        if min_index < len(cr_pts) - 2:
            max1 = -1
            max1_ind = None
            max2 = 10000
            max2_ind = None
            for ind, point in enumerate(cr_pts):
                if max1 < point < cr_pts[min_index]:
                    max1 = point
                    max1_ind = ind
                elif cr_pts[min_index] < point < max2:
                    max2 = point
                    max2_ind = ind
            av_max = (cr_vals[max1_ind] + cr_vals[max2_ind]) / 2
            gap = (av_max - cr_vals[min_index]) / av_max
            if gap > depth:
                deep_min = True
        return deep_min

    @staticmethod
    def critical_points(x_axis, y_axis):
        """
        Finds the critical points of a 1D interpolated function
        :param x_axis: x-values, array-like
        :param y_axis: y-values, array-like
        :return: critical_points list, critical_values list
        """
        def quadratic_spline_roots(spl):
            roots = []
            knots = spl.get_knots()
            for a, b in zip(knots[:-1], knots[1:]):
                u, v, w = spl(a), spl((a + b) / 2), spl(b)
                t = np.roots([u + w - 2 * v, w - u, 2 * v])
                t = t[np.isreal(t) & (np.abs(t) <= 1)]
                roots.extend(t * (b - a) / 2 + (b + a) / 2)
            return np.array(roots)

        f = InterpolatedUnivariateSpline(x_axis, y_axis, k=3)
        cr_pts = quadratic_spline_roots(f.derivative())
        cr_pts = np.append(cr_pts, (x_axis[0], x_axis[-1]))
        cr_vals = f(cr_pts)
        return cr_pts, cr_vals

    def find_local_max_map_value(self, index, map_obj, radius):
        # print(map_obj.data[index])
        values = []
        for i in range(index[0]-radius, index[0]+radius):
            for j in range(index[1] - radius, index[1] + radius):
                for k in range(index[2] - radius, index[2] + radius):
                    val = map_obj.data[i, j, k]
                    values.append(val)
        max_val = max(values)
        # print('Mean',max_val)
        return max_val

    def trace_side_chain_(self, start_coord, c_coord, n_coord, prev_ca_coord, next_ca_coord, interpolator,
                         max_radius, thld_step=0.02, thld_delta=0.01):

        path = []
        bb_vector = n_coord - c_coord
        ort_vector = self.orthogonal_vector(bb_vector)
        bb_axis = bb_vector / LA.norm(bb_vector)

        thld = self.max_thld * 0.9
        min_thld = self.max_thld * 0.3

        ca_int = self.interpolator(self.coord_to_index(start_coord, self.cube_map))
        # print('Start I', ca_int)
        start_index = self.coord_to_index_int(start_coord, self.cube_map)
        loc_max = self.find_local_max_map_value(start_index, self.cube_map, 1)


        # print('THLD', thld)
        # print('Min thld', min_thld)

        fused = True
        while fused:
            fused = False
            path = []
            d = 0
            search_center = start_coord
            #thld = self.max_thld * 0.9
            thld = loc_max * 0.9
            radius = 0.4
            while thld > min_thld and d < 10:
                delta = thld_delta * self.max_thld
                thld_range = [thld - delta, thld + delta]
                # print('THLD range', thld_range)
                found = False

                local_max_d = 0
                thld_center = search_center
                while not found:

                    max_dist_coord, d = self._furthest_from_center(c_coord, bb_axis,
                                                                   ort_vector, search_center,
                                                                   prev_ca_coord, next_ca_coord,
                                                                   interpolator, thld_range, radius)
                    # print("d!!", d)

                    if d == 0 and radius < max_radius:
                        radius += 0.1
                    elif d > local_max_d:
                        search_center = max_dist_coord
                        local_max_d = d
                    else:
                        found = True
                        if max_dist_coord != []:
                            path.append(max_dist_coord)

                rolling_dist = self.distance(thld_center, search_center)
                if rolling_dist > self.max_roll:
                    fused = True
                    min_thld = min_thld + thld_step * self.max_thld
                    break

                thld = thld - self.max_thld * thld_step

        return path


    def trace_side_chain(self, start_coord, c_coord, n_coord, prev_ca_coord, next_ca_coord, interpolator,
                         max_radius, thld_step=0.02, thld_delta=0.01):

        path = []
        bb_vector = n_coord - c_coord
        ort_vector = self.orthogonal_vector(bb_vector)
        bb_axis = bb_vector / LA.norm(bb_vector)

        thld = self.max_thld * 0.9
        min_thld = self.max_thld * 0.3

        ca_int = self.interpolator(self.coord_to_index(start_coord, self.cube_map))
        # print('Start I', ca_int)
        start_index = self.coord_to_index_int(start_coord, self.cube_map)
        loc_max = self.find_local_max_map_value(start_index, self.cube_map, 1)


        # print('THLD', thld)
        # print('Min thld', min_thld)

        fused = True
        while fused:
            fused = False
            path = []
            d = 0
            search_center = start_coord
            #thld = self.max_thld * 0.9
            thld = loc_max * 0.9
            radius = 0.4
            while thld > min_thld and d < 10:
                delta = thld_delta * self.max_thld
                thld_range = [thld - delta, thld + delta]
                # print('THLD range', thld_range)
                found = False

                local_max_d = 0
                thld_center = search_center
                while not found:

                    max_dist_coord, d = self._furthest_from_center(c_coord, bb_axis,
                                                                   ort_vector, search_center,
                                                                   prev_ca_coord, next_ca_coord,
                                                                   interpolator, thld_range, radius)
                    # print("d!!", d)

                    # This stops the search if it goes too close to the flanking CA
                    if prev_ca_coord is not None and d:
                        prev_ca_d = self.distance(max_dist_coord, prev_ca_coord)
                    else:
                        prev_ca_d = 3
                    if next_ca_coord is not None and d:
                        next_ca_d = self.distance(max_dist_coord, next_ca_coord)
                    else:
                        next_ca_d = 3
                    if prev_ca_d < 3 and next_ca_d < 3:
                        return path

                    if d == 0 and radius < max_radius:
                        radius += 0.1
                    elif d > local_max_d:
                        search_center = max_dist_coord
                        local_max_d = d
                    else:
                        found = True
                        if max_dist_coord != []:
                            path.append(max_dist_coord)

                rolling_dist = self.distance(thld_center, search_center)
                if rolling_dist > self.max_roll:
                    fused = True
                    min_thld = min_thld + thld_step * self.max_thld
                    break

                thld = thld - self.max_thld * thld_step

        return path


    def _furthest_from_center(self, bb_point, bb_axis, ort_vector, search_center, prev_ca_coord, next_ca_coord,
                              interpolator, thld_range, radius):

        max_dist_coord = []
        max_dist = 0

        x0 = search_center + ort_vector * radius
        x0 = x0 - bb_axis * 0.2 * 4

        for i in range(9):
            x = x0 + bb_axis * 0.2 * i
            for j in range(16):
                self.operations_counter += 1
                alpha = math.pi / 8 * j

                matrix = self.rotation_matrix(bb_axis, alpha)
                x1 = np.dot(matrix, x)
                y1 = np.dot(matrix, search_center)
                t = y1 - search_center
                x1 = x1 - t

                int_indices = self.coord_to_index_int(x1, self.cube_map)
                try:
                    if self.walk_grid[int_indices] == 0:
                        self.walk_grid[int_indices] = 1
                except IndexError:
                    pass

                indices = self.coord_to_index(x1, self.cube_map)
                try:
                    map_value = interpolator(indices)
                except ValueError:
                    continue

                if thld_range[0] < map_value < thld_range[1]:
                    d = self.distance_to_line(x1, bb_point, vector=bb_axis)


                    if d > max_dist:
                        max_dist_coord = x1
                        max_dist = d
                        if self.walk_grid[int_indices] == 0:
                            self.walk_grid[int_indices] = 1.2
            center_indices = self.coord_to_index_int(search_center, self.cube_map)
            self.walk_grid[center_indices] = 1.5
        # print(self.operations_counter)
        return max_dist_coord, max_dist

    @staticmethod
    def distance(p1, p2):
        """
        Calculates the distance between 2 points in 3D
        :param p1: point coordinates
        :param p2: point coordinates
        :return: distance
        """
        # print(p1, p2)
        d = math.sqrt((p1[0] - p2[0]) ** 2
                      + (p1[1] - p2[1]) ** 2
                      + (p1[2] - p2[2]) ** 2)
        return d

    @staticmethod
    def distance_to_line(x0, x1, x2=None, vector=None):
        """
        Calculates the distance between a point and a line. The line can be defined by 2 points or a point and a vector.
        :param x0: point coordinates (np.array)
        :param x1: point (np.array)
        :param x2: point (np.array)
        :param vector: vector (np.array)
        :return: distance
        """
        if vector is None and x2 is not None:
            vector = x2 - x1
        elif x2 is None and vector is None:
            raise Exception("A direction vector and a point or two points are required to define the line")
        d = np.linalg.norm(np.cross(vector, (x1-x0))) / np.linalg.norm(vector)
        return d

    @staticmethod
    def distance_to_line_points(x0, x1, x2):

        d = math.fabs(np.cross(x0-x1, x0-x2)) / math.fabs(x2 - x1)
        return d

    @staticmethod
    def orthogonal_vector_at_point(point, vector):
        """delete"""
        if vector[2] != 0:
            z = (vector[0]*point[0] + vector[1]*point[1]) / vector[2]
            ort_vector = np.array((round(point[0], 2), round(point[1], 2), round(z, 2)))
        elif vector[1] != 0:
            y = (vector[0] * point[0] + vector[2] * point[2]) / vector[1]
            ort_vector = np.array((round(point[0], 2), round(y, 2), round(point[2], 2)))
        else:
            x = (vector[1] * point[1] + vector[2] * point[2]) / vector[0]
            ort_vector = np.array((round(x, 2), round(point[1], 2), round(point[2], 2)))

        return ort_vector

    @staticmethod
    def orthogonal_vector(v):
        """
        Finds a vector orthogonal to the input vector
        :param v: vector
        :return: unit orthogonal vector
        """
        if v[0] != 0:
            y = randint(1, 2)
            z = randint(1, 2)
            x = (v[1] * y + v[2] * z) / -v[0]
            ort_vector = np.array((x, y, z))

        elif v[1] != 0:
            x = randint(1, 2)
            z = randint(1, 2)
            y = (v[2] * z + v[0] * x) / -v[1]
            ort_vector = np.array((x, y, z))
        else:
            x = randint(1, 2)
            y = randint(1, 2)
            z = (v[0] * x + v[1] * y) / -v[2]
            ort_vector = np.array((x, y, z))
        ort_vector = ort_vector / LA.norm(ort_vector)
        return ort_vector

    @staticmethod
    def line_function(point, vector):
        """
        Defines a function F(x) of a line in 3D which
        :param point: a point on the line
        :param vector: line direction vector
        :return: line function
        """
        def f(x):
            y = (x - point[0]) * vector[1] / vector[0] + point[1]
            z = (x - point[0]) * vector[2] / vector[0] + point[2]
            coord = (x, y, z)
            return coord

        return f

    @staticmethod
    def coord_to_index_int(coord, map_obj):
        """
        Converts atomic residue coordinates into integer map indices
        :param coord: coordinate in angstrom
        :param map_obj: map object
        :return: int indices
        """
        map_obj.origin = np.array(map_obj.origin)
        st = map_obj.n_start
        vs = map_obj.voxel_size
        orig = map_obj.origin
        x_index = int(round(coord[2] / vs[2] - orig[2] / vs[2] - st[2]))
        y_index = int(round(coord[1] / vs[1] - orig[1] / vs[1] - st[1]))
        z_index = int(round(coord[0] / vs[0] - orig[0] / vs[0] - st[0]))
        return x_index, y_index, z_index

    @staticmethod
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

    @staticmethod
    def rotation_matrix(axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        axis = axis / math.sqrt(np.dot(axis, axis))
        a = math.cos(theta / 2.0)
        b, c, d = -axis * math.sin(theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                         [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                         [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    @staticmethod
    def get_threshold(map_data, level):
        out_map = np.copy(map_data)
        it = np.nditer(out_map, flags=['multi_index'], op_flags=['writeonly'])
        while not it.finished:
            if it[0] < level:
                it[0] = 0
            it.iternext()
        return out_map

    def trace_length(self, trace):
        length = 0
        if len(trace) > 2:
            a1 = trace[0]
            for a2 in trace:
                d = self.distance(a1, a2)
                length += d
                a1 = a2
        return length

