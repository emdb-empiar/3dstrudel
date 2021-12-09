import argparse
import os
import csv
import sys
import matplotlib.pyplot as plt
import itertools
# import seaborn as sns
import collections

import pandas as pd


def statistics(csv_path):
    data = pd.read_csv(csv_path)
    # data = data.sort_values(by=['target'])
    print(data)
    var = data.var()
    print(var[0])
    print(var)
    var_list = []
    for n, val in var.iteritems():
        var_list.append([n, val])
    var_list.sort(key=lambda i: i[0])




    data['variance'] = var
    data = data.sort_values(by=['target'])

    fig, ax = plt.subplots(figsize=(40, 5))
    # ax.bar(data['target'], data['variance'])
    print([i[0] for i in var_list], [i[1] for i in var_list])
    ax.bar([i[0] for i in var_list], [i[1] for i in var_list])
    plt.xticks(rotation=30, ha='right')
    plt.show()

    averaged_var = {}

    for name, val in var_list:
        key = name.split('_')[0]
        if key not in averaged_var.keys():
            averaged_var[key] = [val, 1]
        else:
            averaged_var[key] = [averaged_var[key][0] + val, averaged_var[key][1]+1]

    print(averaged_var)
    for key, val in averaged_var.items():
        averaged_var[key] = val[0]/val[1]
    print(averaged_var)
    fig1, ax1 = plt.subplots(figsize=(10, 5))
    ax1.bar(averaged_var.keys(), averaged_var.values())
    plt.xticks(rotation=30, ha='right')
    plt.show()

    # data.to_csv("/Volumes/data/libs/from_ftp/strudel-libs_ver-2.0_voxel-0.5/var_0.0-2.3_2.csv")


def filter_data1(dataframe):
    dataframe = dataframe.sort_values(by=['target'])
    new_frame = pd.DataFrame()
    length = len(dataframe['target'])

    for index, series in dataframe.iterrows():

        correlations = []
        ref_type = series["target"].split("_")[0]
        for a, b in series.iteritems():
            # print(a, b)
            if a == 'target':
                pass
            else:
                res_type = a.split('_')[0]
                if res_type != ref_type:
                    correlations.append(b)
        if len(correlations) < length:
            for i in range(length-len(correlations)):
                correlations.append(None)
        new_frame[series["target"]] = correlations
    variance = new_frame.mean()
    print(variance)

    averaged_var = {}

    for name, val in variance.iteritems():
        key = name.split('_')[0]
        if key not in averaged_var.keys():
            averaged_var[key] = [val, 1]
        else:
            averaged_var[key] = [averaged_var[key][0] + val, averaged_var[key][1] + 1]
    print(averaged_var)
    for key, val in averaged_var.items():
        averaged_var[key] = val[0] / val[1]
    print(averaged_var)

    fig1, ax1 = plt.subplots(figsize=(10, 5))
    ax1.bar(averaged_var.keys(), averaged_var.values())
    plt.xticks(rotation=30, ha='right')
    ax1.set_ylim(0.6, 1)
    plt.show()
    new_frame.to_csv("/Volumes/data/libs/from_ftp/strudel-libs_ver-2.0_voxel-0.5/filtered_3.5-4.0.csv")
    tmp_dict = {}

def filter_data(dataframe):
    dataframe = dataframe.sort_values(by=['target'])
    new_frame = pd.DataFrame()
    length = len(dataframe['target'])

    for index, series in dataframe.iterrows():

        correlations = []
        ref_type = series["target"].split("_")[0]
        for a, b in series.iteritems():
            # print(a, b)
            if a == 'target':
                pass
            else:
                res_type = a.split('_')[0]
                if res_type != ref_type:
                    correlations.append(b)
        if len(correlations) < length:
            for i in range(length - len(correlations)):
                correlations.append(None)
        new_frame[series["target"]] = correlations
    return new_frame


def extract_max_corr(dataframe):
    dataframe = dataframe.sort_values(by=['target'], ascending=False)
    new_frame = pd.DataFrame()
    # length = len(dataframe['target'])

    for index, series in dataframe.iterrows():
        tmp_dict = {}
        ref_type = series["target"].split("_")[0]
        tmp_dict[ref_type] = 1
        for a, b in series.iteritems():
            # print(a, b)
            if a == 'target':
                pass
            else:
                res_type = a.split('_')[0]
                if res_type != ref_type:
                    if res_type not in tmp_dict.keys():
                        tmp_dict[res_type] = round(b, 5)
                    elif tmp_dict[res_type] < b:
                        tmp_dict[res_type] = round(b, 5)

        ordered_dict = collections.OrderedDict(sorted(tmp_dict.items()))

        new_frame[series["target"]] = ordered_dict.values()
    new_frame.insert(loc=0, column='res', value=ordered_dict.keys())
    new_frame = new_frame.set_index('res')
    return new_frame





def comute_mean(frame):
    mean_series = frame.mean()
    return mean_series

def average_over_rotamers(series):
    averaged_val = {}

    for name, val in series.iteritems():
        key = name.split('_')[0]
        if key not in averaged_val.keys():
            averaged_val[key] = [val, 1]
        else:
            averaged_val[key] = [averaged_val[key][0] + val, averaged_val[key][1] + 1]
    # print(averaged_val)
    for key, val in averaged_val.items():
        averaged_val[key] = val[0] / val[1]
    return averaged_val

def series_to_dict(series):
    dic_data = {}

    for name, val in series.iteritems():
            dic_data[name] = val
    return dic_data



def plot_single(data_dict, resolution_range, plot_path=None):
    fig, ax = plt.subplots(figsize=(80, 10))
    ax.bar(data_dict.keys(), data_dict.values())
    plt.xticks(rotation=45, ha='right')
    ax.set_ylim(0.6, 1)
    plt.title(f'Resolution band: {resolution_range}')
    ax.set_xlabel('Motif')
    ax.set_ylabel('Average correlation')
    plt.subplots_adjust(bottom=0.55)
    if plot_path is None:
        plt.show()
    else:
        fig.savefig(plot_path)


def plot_multiple(data_list, plot_path=None, y_lim=(0, 1)):
    hatches = itertools.cycle(['///', '+++', 'oo', 'XX', 'OO', '.', '--'])

    fig, ax = plt.subplots(figsize=(12, 5))
    width = 0.4
    dd = width*len(data_list)/2 + 1
    x_ticks = [i for i in range(20)]
    x_ticks = [ x * dd for x in x_ticks]
    print(len(x_ticks))
    ax.set_xticks([ x + width * len(data_list)/2 for x in x_ticks])

    for data, res_range in data_list:
        print(res_range)
        print(len(data))
        x_ticks = [x+width for x in x_ticks]
        # ax.bar(data.keys(), data.values(), label=f'{res_range}')
        hatch = next(hatches)
        bar = ax.bar(x_ticks, data.values(), width, label=f'{res_range} $\AA$')
        print(bar)
        for b in bar:
            b.set_hatch(hatch)
    ax.legend(prop={'size': 13})

    ax.set_xticklabels(data_list[0][0].keys(), rotation=0)
    ax.set_xlim(-2*width, x_ticks[-1] + 2*width)
    ax.set_ylim(y_lim)
    x = plt.xlabel('Motif type')
    x.set_fontsize(17)
    # ax.set_ylabel('1-correlation')
    y = plt.ylabel('1-correlation')
    y.set_fontsize(17)

    if plot_path is None:
        plt.show()
    else:
        plt.subplots_adjust(bottom=0.15)
        fig.savefig(plot_path)


def plot_multiple_non_av(data_list, plot_path=None, y_lim=(0, 1)):
    fig, ax = plt.subplots(figsize=(12, 5))
    width = 0.25
    dd = width*len(data_list)/2 + 1
    x_ticks = [i for i in range(len(data_list[0][0]))]
    x_ticks = [ x * dd for x in x_ticks]
    print(len(x_ticks))
    ax.set_xticks([ x + width * len(data_list)/2 for x in x_ticks])

    for data, res_range in data_list:
        print(res_range)
        print(len(data))
        x_ticks = [x+width for x in x_ticks]
        # ax.bar(data.keys(), data.values(), label=f'{res_range}')
        ax.bar(x_ticks, data.values(), width, label=f'{res_range}')
    ax.legend()

    ax.set_xticklabels(data_list[0][0].keys(), rotation=65)
    ax.set_xlim(0, x_ticks[-1] + 5.5)
    ax.set_ylim(y_lim)
    ax.set_xlabel('Motif type')
    ax.set_ylabel('Average correlation')
    if plot_path is None:
        plt.show()
    else:
        plt.subplots_adjust(bottom=0.15)
        fig.savefig(plot_path)

# def filter_max_

def one_minus_val(val_dict):
    for key, val in val_dict.items():
        val_dict[key] = 1 - val
    return val_dict


def plot_all_average_corr(data_path, out_folder=None, action='mean'):
    files = os.listdir(data_path)
    csv_files = [f for f in files if f.endswith('.csv') and not f.startswith('.')]
    csv_files.sort()
    data_list = []
    tmp = csv_files[:3] + [csv_files[-1]]
    # tmp = csv_files
    for file in tmp:
        print(file)
        res_range = file.split('.csv')[0].split('_')[-1]
        path = os.path.join(data_path, file)
        data = pd.read_csv(path)
        filtered = filter_data(data)
        if action == 'mean':
            series = comute_mean(filtered)
        else:
            return
        rot_aver = average_over_rotamers(series)
        print(rot_aver)
        rot_aver = one_minus_val(rot_aver)

        data_list.append([rot_aver, res_range])

    plot_multiple(data_list, y_lim=(0, 0.3), plot_path=os.path.join(out_folder, 'aver_corr_3+1_one_minus_hatch_final.png'))


def plot_all_corr(data_path, out_folder=None, action='mean'):
    files = os.listdir(data_path)
    csv_files = [f for f in files if f.endswith('.csv')]
    csv_files.sort()
    data_list = []

    for file in csv_files:
        res_range = file.split('.csv')[0]
        path = os.path.join(data_path, file)
        data = pd.read_csv(path)
        filtered = filter_data(data)
        if action == 'mean':
            series = comute_mean(filtered)
        else:
            return
        rot_aver = series_to_dict(series)

        data_list.append([rot_aver, res_range])

        plot_single(rot_aver, res_range, plot_path=os.path.join(out_folder,f'{res_range}_all_rot.png'))
        break


def plot_grid_similarity(df, plot_path):
    import numpy as np

    import matplotlib as mpl
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['xtick.labeltop'] = True

    fig, ax = plt.subplots(figsize=(16, 36))

    from matplotlib import colors
    cmap = colors.ListedColormap(['white', 'blue', 'white'])
    bounds = [0, 0.95, 0.999, 1]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    ax.pcolor(df, cmap=cmap, norm=norm)

    y_labels = [l.split('_') for l in df.index]

    form_labels = []
    for l in y_labels:
        if l[0] in ['ala', 'gly']:
            txt = f'{l[0]}'
        else:
            txt = f'{l[0]} {l[-1]}'
        form_labels.append(txt)

    ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=0.5)
    plt.yticks(np.arange(0.99, len(df.index), 1), form_labels, va='top', ha='right')
    print(np.arange(0.5, len(df.columns), 1))
    print(df.columns)
    labels = [l + '  ' for l in df.columns]
    plt.xticks(np.arange(0.99, len(df.columns), 1), labels, ha='right')

    plt.tick_params(axis='x', which='major', labelsize=16)
    plt.tick_params(axis='y', which='major', labelsize=14)
    if plot_path is None:
        plt.show()
    else:
        plt.subplots_adjust(bottom=0.15)
        fig.savefig(plot_path)



def generate_all_drid_similarity(data_path, out_folder=None):
    files = os.listdir(data_path)
    csv_files = [f for f in files if f.endswith('.csv') and not f.startswith('.')]
    csv_files.sort()
    data_list = []
    tmp = csv_files
    # tmp = csv_files
    for file in tmp:
        print(file)
        res_range = file.split('.csv')[0].split('_')[-1]
        path = os.path.join(data_path, file)
        data = pd.read_csv(path)
        filtered = extract_max_corr(data)
        df = filtered.transpose()
        # df = df.sort_index(axis=1)


        plot_grid_similarity(df, os.path.join(out_folder,f'{res_range}_sim_grid.png'))



# path = '/Volumes/data/libs/from_ftp/strudel-libs_ver-2.0_voxel-0.5/2.3-2.5.csv'

# data = pd.read_csv(path)
# filter_data(data)

# statistics(path)
data_path = '/Users/andrei/Documents/Project_data/3D-STRUDEL/libs.3/out'
fig_path = '/Users/andrei/Documents/Project_data/3D-STRUDEL/libs.3/pics_2'
if not os.path.exists(fig_path):
    os.makedirs(fig_path)
plot_all_average_corr(data_path, out_folder=fig_path)
## plot_all_corr(data_path, out_folder=fig_path)

# generate_all_drid_similarity(data_path, out_folder=fig_path)



# tmp_path = '/Volumes/data/libs/from_ftp/strudel-libs_ver-2.0_voxel-0.5/libs_similarity/0.0-2.3.csv'
# tmp_path = '/Volumes/data/libs/from_ftp/strudel-libs_ver-2.0_voxel-0.5/libs_similarity/2.8-3.0.csv'
# data = pd.read_csv(tmp_path)

# new_data = extract_max_corr(data)
# new_data.to_csv('/Volumes/data/libs/from_ftp/strudel-libs_ver-2.0_voxel-0.5/libs_similarity/0.0-2.3_max.csv')
# new_data.to_csv('/Volumes/data/libs/from_ftp/strudel-libs_ver-2.0_voxel-0.5/libs_similarity/2.8-3.0_max.csv')


# no_res = new_data.drop('res', 1)
# np_data = new_data.to_numpy()
# from matplotlib import colors
#
# cmap = colors.ListedColormap(['white', 'blue', 'white'])
# bounds = [0,0.95,0.999, 1]
# norm = colors.BoundaryNorm(bounds, cmap.N)
#
# fig, ax = plt.subplots()
# ax.imshow(np_data, cmap=cmap, norm=norm)
import numpy as np
# # draw gridlines
# ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=0.5)
# ax.set_xticks(np.arange(-.5, 130, 1))
# ax.set_yticks(np.arange(-.5, 20, 1))
#
# plt.show()

# df = new_data.transpose()
# df = df.sort_index(axis=1)
# df.to_csv('/Volumes/data/libs/from_ftp/strudel-libs_ver-2.0_voxel-0.5/libs_similarity/0.0-2.3_max_transp.csv')
#
# print(df.head())
# import matplotlib as mpl
#
# mpl.rcParams['xtick.top'] = True
# mpl.rcParams['xtick.labeltop'] = True
#
# fig, ax = plt.subplots(figsize=(16, 36))
# # plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True
# # plt.pcolor(df)
# from matplotlib import colors
# cmap = colors.ListedColormap(['white', 'blue', 'white'])
# bounds = [0,0.95,0.999, 1]
# norm = colors.BoundaryNorm(bounds, cmap.N)
# # ax.imshow(df, cmap=cmap, norm=norm)
# ax.pcolor(df, cmap=cmap, norm=norm)
# # plt.pcolor(df)
# y_labels = [l.split('_') for l in df.index]
# print(len(y_labels), 'y_lb')
# form_labels = []
# for l in y_labels:
#     if l[0] in ['ala', 'gly']:
#         txt = f'{l[0]}'
#     else:
#         txt = f'{l[0]} {l[-1]}'
#     form_labels.append(txt)
# print(form_labels)
# print(len(form_labels), 'f_lb')
# y_labels = [l[0] + ' ' + l[-1] for l in y_labels]
#
# ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=0.5)
# plt.yticks(np.arange(0.99, len(df.index), 1), form_labels, va='top', ha='right')
# print(np.arange(0.5, len(df.columns), 1))
# print(df.columns)
# labels = [l + '   ' for l in df.columns]
# plt.xticks(np.arange(0.99, len(df.columns), 1), labels, ha='right')
# # def same(x):
# #     return x
# # secaxx = ax.secondary_xaxis('top')
# plt.tick_params(axis='x', which='major', labelsize=16)
# plt.tick_params(axis='y', which='major', labelsize=14)
# plt.show()
#
# # import seaborn as sns
# #
# # sns.heatmap(new_data, annot=True)