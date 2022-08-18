# %%
"""
This module stores some obsolete functionality, including:
    - extracting source data from .dat files from MCNP input
    - Writing that source data to hdf5 files and testing functions
    to check correct implementation
"""

import h5py
import matplotlib.pyplot as plt
import openmc
import numpy as np
import pandas as pd

import os
import re


IDS = ['nPrompt_725g', 'nDelayed', 'gPrompt', 'gDelayed']


def extract_source_data(filename):
    """
    Extracts neutron/gamma source data from MCNP input file,
    Input:
        MCNP input file path
    Returns:
        tuple with four entries:
            string of comments
            list of column names
            np.array of energy bins
            np.array of probabilities for those energy bins (unit varies)
    """
    with open(filename) as f:
        text = f.readlines()

    columns = re.compile(r'^c\s+(\w[\w,/ ]+)\s{2,}(\w[ \w,/]+)', re.M)
    comment = re.compile(r'^#.+')
    data = re.compile(r'^\s+([\d+.E\-+]+)\s+([\d+.E\-+]+)')

    column_list = []
    comment_list = []
    data1 = []
    data2 = []
    for line in text:
        col = re.match(columns, line)
        com = re.match(comment, line)
        dat = re.match(data, line)

        if col:
            # Just tidying some things up
            # Note that we replaced '/' with 'per' because h5py does weird things with '/'
            cols = [col.group(1), col.group(2)]
            cols_clean = [column.strip().replace('/', ' per ')
                          for column in cols]

            column_list.append(cols_clean[0])
            column_list.append(cols_clean[1])
        if com:
            comment_list.append(com.group())
        if dat:
            data1.append(dat.group(1))
            data2.append(dat.group(2))
    return comment_list, column_list, np.array(data1, dtype=np.float64), np.array(data2, dtype='f')


def get_source_file_paths():
    """
        Gives you the filepaths of the .dat MCNP input sources
        and the output paths to write the openmc h5 sources
        Input: 
            No input
        Returns: tuple of two lists.
            First list is input .dat file paths
            Second is output .h5 file paths
    """
    dir = os.getcwd()
    filenames = ['Cf252_' + id + '.dat' for id in IDS]
    input_folder = 'external_data/sources/mcnp_input_source'
    output_folder = 'h5_sources'
    filepaths = [os.path.join(dir, input_folder, filename)
                 for filename in filenames]
    fileouts = []
    for i in range(len(filepaths)):
        fileout = os.path.join(
            dir, output_folder, 'Cf252_' + IDS[i] + '_out.h5')
        fileouts.append(fileout)
    return filepaths, fileouts


def get_source_data_dict(convert=False):
    """
    No input, returns dict of dicts
    corresponding to mcnp input data. Keys are the file ids.
    """
    fileins, fileouts = get_source_file_paths()
    array_dict = {}
    for i, file in enumerate(fileins):
        comments, columns, data1, data2 = extract_source_data(file)
        if convert:
            data1 *= 10**6
            data2 *= 10**6
        if data1[0] < 1.0e-5:
            data1[0] = 1.0e-5
        ########################
        # Experimental stuff
        test_below = 1
        limit = 2.0e7
        data1 = data1[test_below:]
        data1 = data1[data1 <= limit]
        length = len(data1)
        data2 = data2[:length]
        ########################
        array_dict[IDS[i]] = {'x': data1.round(10),
                              'y': data2,
                              'columns': columns}
    return array_dict


def create_hdf5():
    """
        Uses get_file_paths to get input file and output file paths,
        reads in data from .dat input files into .h5 output files
    """
    fileins, fileouts = get_source_file_paths()
    for i in range(len(fileins)):
        with h5py.File(fileouts[i], 'w') as f:
            comments, columns, data1, data2 = extract_source_data(fileins[i])
            col1 = columns[0]
            col2 = columns[1]
            f.create_dataset(col1, data=data1)
            f.create_dataset(col2, data=data2)
            f.attrs['comments'] = comments
            f.attrs['columns'] = columns
    return fileouts


def test_hdf5_output_files(fileouts):
    """
        Tests whether the create_h5 function actually worked properly
        to write the h5 files
    """
    fileins, fileouts = get_source_file_paths()

    # Checking file outputs
    with h5py.File(fileouts[0], 'r') as f:
        dsets = list(f.keys())
        print(dsets)
        dset1 = f['Eup, MeV']
        print(dset1[0:20])
        dset2 = f['Yield, n per bin']
        print(dset2[0:20])
        print(dset1[-20:])
        print(dset2[-20:])
        print(len(dset1))
        print(len(dset2))

    with h5py.File(fileouts[1], 'r') as f:
        dsets = list(f.keys())
        print(dsets)
        dset1 = f['E, MeV']
        print(dset1[0:20])
        dset2 = f['Spectr, n per MeV']
        print(dset2[0:20])
        print(dset1[-20:])
        print(dset2[-20:])
        print(len(dset1))
        print(len(dset2))


def printname(x):
    print(x)


def inspect_openmc_h5():
    f = h5py.File("external_data/Fe56.h5")
    keys = f.keys()
    mat = f['Fe56']
    energy = mat['energy']
    temp = energy['294K']
    print(temp)
    kt = mat['kTs']
    scalar = kt['294K']
    reactions = mat['reactions']
    reaction_002 = reactions['reaction_002']
    reaction_002_294K = reaction_002['294K']['xs']
    print(reaction_002_294K)
    # f.visit(printname)
    reaction = f['Fe56/reactions/reaction_016']
    atts = list(reaction.attrs)
    print(atts)
    for at in atts:
        print(reaction.attrs[at])


def old_way_of_getting_source():
    # Getting the source
    df = pd.read_csv("external_data/sources/new_cf252_source.csv")
    source = openmc.Source()
    x = np.append(df['E_low'].values, df['E_max'].values[-1])
    y = np.append(df['prompt_n_source'].values,
                  df['prompt_n_source'].values[-1])
    source.energy = openmc.stats.Tabular(
        x, y, interpolation='histogram')


def read_partisn_data():
    """
    Reads the partisn data files, and returns a list of dataframes
    for each tally
    Input:
        No input
    Returns:
        Returns list of dataframes with data for each partisn tally
    """
    filepath = "external_data/partisn_model_tallies.txt"
    with open(filepath) as f:
        text = f.read()
    dataset = re.compile(r'(?:^\w.+$\n){2,}', re.MULTILINE)
    columns = re.compile(
        r"^((?:[-\w.()]+ ?)+)(?: {2,}|\t+)((?:[-\w.()]+ ?)+)(?: {2,}|\t+)*$", re.MULTILINE)
    data = re.compile(
        r"^([-\dEe.+]+)\s+([-\dEe.+]+)\s*$", re.MULTILINE)
    datasets = re.findall(dataset, text)
    dfs = []
    for dataset in datasets:
        cols = re.match(columns, dataset)
        col1 = cols.group(1)
        col2 = cols.group(2)
        data1 = []
        data2 = []
        for line in dataset.split('\n'):
            match = re.match(data, line)
            if match:
                data1.append(match.group(1).strip())
                data2.append(match.group(2).strip())
        df = pd.DataFrame({col1: data1, col2: data2})
        dfs.append(df.astype(np.float64))

    df_gamma_partisn, df_neutron_partisn = dfs
    gamma_bins_partisn = np.append(
        df_gamma_partisn.iloc[:, 0].values, 3.0e6)
    neutron_bins_partisn = np.append(
        df_neutron_partisn.iloc[:, 0].values, 20.0e7)
    final_gamma_df = pd.DataFrame({"energy low [eV]": gamma_bins_partisn[:-1],
                                   "energy high [eV]": gamma_bins_partisn[1:],
                                   "Flux": df_gamma_partisn.iloc[:, 1]})
    final_neutron_df = pd.DataFrame({"energy low [eV]": neutron_bins_partisn[:-1],
                                     "energy high [eV]": neutron_bins_partisn[1:],
                                     "Flux": df_neutron_partisn.iloc[:, 1]})
    # an adjustment from 15-15.5 spherical shell to cell 70
    final_neutron_df = read_partisn_neutron_70()
    return final_gamma_df, final_neutron_df


def read_partisn_neutron_70():
    filepath = "external_data/partisn_outputs/partisn-neutron-70.txt"
    df = pd.read_csv(filepath, delimiter='\t')
    new_df = df[["e-low_eV", "e-high_eV", "integral"]]
    new_df = new_df.set_axis(["energy low [eV]", "energy high [eV]", "Flux"],
                             axis=1)
    return new_df


def read_partisn_gamma_52():
    filepath = "external_data/partisn_outputs/partisn-gammas-52.txt"
    df = pd.read_csv(filepath, delimiter='\t')
    new_df = df[["e-low_eV", "e-high_eV", "integral"]]
    new_df = new_df.set_axis(["energy low [eV]", "energy high [eV]", "Flux"],
                             axis=1)
    return new_df


if __name__ == '__main__':
    inspect_openmc_h5()
# %%
