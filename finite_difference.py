# %%
import h5py
import matplotlib.pyplot as plt
import numpy as np
import openmc
import pandas as pd

import importlib
import os
import re


import data_load
import utils
import config

importlib.reload(data_load)


def find_power_folder_max(home_folder):
    directories = os.listdir(home_folder)
    power = re.compile(r'e(\d)')
    powers = []
    for subdir in directories:
        p = re.match(power, subdir)
        if p:
            powers.append(int(p.group(1)))
    max_p = max(powers)
    power_folder = 'e' + str(max_p)
    final_path = os.path.join(home_folder, power_folder)
    return final_path


def get_perturb_table():
    table_file = os.path.join(
        config.MAIN_DIR, 'perturbed_run_data/tables', f'e{config.N}.csv')
    if not os.path.exists(os.path.dirname(table_file)):
        os.makedirs(os.path.dirname(table_file))
    if os.path.exists(table_file):
        table = pd.read_csv(table_file, index_col=[0, 1])
    else:
        table = {"reaction": [],
                 "perturbation": [],
                 "neutron flux sensitivity": [],
                 "gamma flux sensitivity": []}
        table = pd.DataFrame(table)
        table = table.set_index(["reaction", "perturbation"])
        # I know it's sacrilegious to iteratively append to a dataframe
        # But there's so little data going in it that it doesn't matter
    return table


def save(table):
    table_file = os.path.join(
        config.MAIN_DIR, 'perturbed_run_data/tables', f'e{config.N}.csv')
    table.to_csv(table_file)
    return


def compare_perturbation(mt, perturbation, discretization=None,
                         group=None, power=None):
    """Function computes the sensitivity divided by lethargy of the
     neutron flux in the 0.51-2 MeV region following the perturbation
     of a cross section pertaining to a given nuclide, reaction,
     perturbation, discretisation, and discretisation group"""

    # Getting file locations of unperturbed and perturbed runs
    standard_folder = os.path.join(
        config.MAIN_DIR, f'standard_run/output/e{config.N}')
    sens_n_file = 'sens_n.csv'
    sens_g_file = 'sens_g.csv'
    path_sens_n = os.path.join(standard_folder, sens_n_file)
    path_sens_g = os.path.join(standard_folder, sens_g_file)

    if discretization is None:
        id_code = f'mt{mt}-p{perturbation}'
    else:
        id_code = f'mt{mt}-p{perturbation}-d{discretization:03}-g{group+1:03}'
    perturb_home_folder = os.path.join(
        config.MAIN_DIR, 'perturbed_run_data/')
    output_folder = os.path.join(perturb_home_folder, id_code, 'output')
    perturb_folder = os.path.join(output_folder, f'e{config.N}')

    sens_n_file = 'sens_n.csv'
    sens_g_file = 'sens_g.csv'
    path_sens_n_p = os.path.join(perturb_folder, sens_n_file)
    path_sens_g_p = os.path.join(perturb_folder, sens_g_file)

    unperturbed_n = pd.read_csv(path_sens_n)['mean'].values[0]
    unperturbed_g = pd.read_csv(path_sens_g)['mean'].values[0]
    perturbed_n = pd.read_csv(path_sens_n_p)['mean'].values[0]
    perturbed_g = pd.read_csv(path_sens_g_p)['mean'].values[0]

    sens_n = (perturbed_n - unperturbed_n) / unperturbed_n / perturbation
    sens_g = (perturbed_g - unperturbed_g) / unperturbed_g / perturbation
    table = get_perturb_table()
    df = pd.DataFrame({"reaction": [mt],
                       "perturbation": [perturbation],
                       "neutron flux sensitivity": [sens_n],
                       "gamma flux sensitivity": [sens_g]})
    df = df.set_index(["reaction", "perturbation"])
    new_table = pd.concat([table, df])
    new_table = new_table.sort_index()
    new_table = new_table.drop_duplicates()
    save(new_table)
    print(new_table)
    return new_table


if __name__ == '__main__':
    config.N = 7
    config.LIBDIR = "/root/nndc_hdf5"
    config.MAIN_DIR = '/ironbenchmark/Fe-simplified'
    table = compare_perturbation(102, 0.1)

    # %%
