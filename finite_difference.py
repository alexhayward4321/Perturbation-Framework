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
import settings

importlib.reload(data_load)


def find_power_folder(home_folder):
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


def compare_perturbation(nuclide, mt, perturbation, discretization=None,
                         group=None):
    """Function computes the sensitivity divided by lethargy of the
     neutron flux in the 0.51-2 MeV region following the perturbation
     of a cross section pertaining to a given nuclide, reaction,
     perturbation, discretisation, and discretisation group"""

    # Getting file locations of unperturbed and perturbed runs
    standard_folder = find_power_folder(
        os.path.join(settings.MAIN_DIR, 'standard_run/output'))
    sens_n_file = 'sens_n.csv'
    sens_g_file = 'sens_g.csv'
    path_sens_n = os.path.join(standard_folder, sens_n_file)
    path_sens_g = os.path.join(standard_folder, sens_g_file)

    if discretization is None:
        id_code = f'mt{mt}-p{perturbation}'
    else:
        id_code = f'mt{mt}-p{perturbation}-d{discretization:03}-g{group+1:03}'
    perturb_home_folder = os.path.join(
        settings.MAIN_DIR, 'perturbed_run_data/')
    output_folder = os.path.join(perturb_home_folder, id_code, 'output')
    perturb_folder = find_power_folder(output_folder)
    print(perturb_folder)

    sens_n_file = 'sens_n.csv'
    sens_g_file = 'sens_g.csv'
    path_sens_n_p = os.path.join(perturb_folder, sens_n_file)
    path_sens_g_p = os.path.join(perturb_folder, sens_g_file)

    unperturbed_n = pd.read_csv(path_sens_n)['mean'].values[0]
    unperturbed_g = pd.read_csv(path_sens_g)['mean'].values[0]
    perturbed_n = pd.read_csv(path_sens_n_p)['mean'].values[0]
    perturbed_g = pd.read_csv(path_sens_g_p)['mean'].values[0]

    sens_n_leth = (perturbed_n - unperturbed_n) / unperturbed_n
    sens_g_leth = (perturbed_g - unperturbed_g) / unperturbed_g
    print(perturbed_n)
    print(unperturbed_n)
    print(perturbed_g)
    print(unperturbed_g)
    return sens_g_leth, sens_n_leth


if __name__ == '__main__':
    settings.LIBDIR = "/root/nndc_hdf5"
    settings.MAIN_DIR = '/ironbenchmark/Fe-simplified'
    g, n = compare_perturbation('Fe56', 2, 0.01, discretization=None,
                                group=None)


# %%
