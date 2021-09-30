import os
import numpy as np
from scipy.interpolate import RegularGridInterpolator

import threed_strudel.utils.functions as func
import threed_strudel.lib as lib_module


def get_molprobity_file_path(residue):
    head = os.path.dirname(lib_module.__file__)
    file_path = os.path.join(head, 'molprobity', 'rota8000-{}.data'.format(residue.lower()))
    if os.path.exists(file_path):
        return file_path
    else:
        raise Exception("could not find molprobity data file for: " + residue)

def generate_molprobity_interpolator(residue):
    '''
    Load a MolProbity data set and return a SciPy RegularGridInterpolator
    object.
    '''
    data_file = get_molprobity_file_path(residue)
    with open(data_file, 'rt', encoding='utf-8') as infile:
        infile.readline()
        # Get number of dimensions
        ndim = int(infile.readline().split()[-1])
        # Throw away the next line - it's just headers
        infile.readline()
        lower_bound = []
        upper_bound = []
        number_of_bins = []

        step_size = []
        first_step = []
        last_step = []
        axis = []

        # Read in the header to get the dimensions and step size for each
        # axis, and initialise the axis arrays
        for i in range(ndim):
            line = infile.readline().split()
            lb = float(line[2])
            lower_bound.append(lb)
            ub = float(line[3])
            upper_bound.append(ub)
            nb = int(line[4])
            number_of_bins.append(nb)

            ss = (ub - lb) / nb
            step_size.append(ss)
            # Values are at the midpoint of each bin
            fs = lb + ss / 2
            first_step.append(fs)
            ls = ub - ss / 2
            last_step.append(ls)
            axis.append(np.linspace(fs, ls, nb))

        data = np.loadtxt(infile)
    full_grid = np.zeros(number_of_bins)
    # Convert each coordinate to an integral number of steps along each
    # axis
    axes = []
    for i in range(ndim):
        axes.append([])

    for i in range(ndim):
        ss = step_size[i]
        fs = first_step[i]
        lb = lower_bound[i]
        axis_vals = data[:, i]
        axes[i] = (((axis_vals - ss / 2 - lb) / ss).astype(int))

    full_grid[tuple(axes)] = data[:, ndim]

    # Replace all zero or negative values with the minimum positive non-zero
    # value, so that we can use logs
    full_grid[full_grid <= 0] = np.min(full_grid[full_grid > 0])
    # Finally, convert the axes to radians
    axis = np.array(axis)
    # axis = axis / 180 * pi

    return RegularGridInterpolator(axis, full_grid, bounds_error=True)


def read_rotamer_names(bounds_file_path, out_json):

    data = {}
    names = []
    with open(bounds_file_path, 'r') as infile:
        lines = infile.readlines()
        for line in lines:
            s_line = line.split('=')
            r_name, rot_name = s_line[0].split()
            if (r_name, rot_name) not in names:
                names.append((r_name, rot_name))
            if r_name not in data.keys():
                data[r_name] = []
                r_id = 1
            bounds = s_line[-1]
            bounds = bounds.strip().strip("\"")

            bounds = list(map(int, bounds.split(", ")))
            print(bounds)
            r_dict = {'id': r_id, 'name': rot_name}
            for i in range(0, len(bounds), 2):
                print(i)
                delta = bounds[i+1] - bounds[i]
                mid = bounds[i] + delta / 2
                r_dict['chi_{}'.format(int((i+2)/2))] = mid
                r_dict['chi-width_{}'.format(int((i + 2) / 2))] = delta
            data[r_name].append(r_dict)
            r_id += 1
        # print(data)
        # print(len(names))
    func.save_json(out_json, data)




# data_dir = '/Users/andrei/Downloads/reference_data-master/Top8000/Top8000_rotamer_pct_contour_grids/'
# c_dir = '/Users/andrei/Downloads/reference_data-master/Top8000/Top8000_rotamer_pct_contour_grids/'
#
# prefix = 'rota8000-arg'
# start = time.time()
# interp = generate_molprobity_interpolator(data_dir, prefix)
# print('built time', time.time()-start)
# start = time.time()
# # print(interp([180.6/180*pi, 77.8/180*pi, 175.2/180*pi, 143.4/180*pi]))
# print(interp([180.6, 77.8, 175.2, 143.4]))
# print('interp time', time.time()-start)