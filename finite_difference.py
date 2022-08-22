# %%
import h5py
import matplotlib.pyplot as plt
import numpy as np
import openmc
import pandas as pd

import importlib
import os
import re

import utils
import settings


def group_len():
    test_file = libdir/f"Fe56.h5"
    with h5py.File(test_file) as f:
        bin_struct = f[f"/Fe56/energy/{Temp}K"][:]
        bin_len = len(bin_struct)
    return bin_len


def get_cum_energy_groups(length, discretisation):
    group_size_guide = length // discretization
    group_size = []
    while length > group_size_guide:
        length -= group_size_guide
        group_size.append(group_size_guide)

    i = 0
    while length > 0:
        group_size[i] += 1
        i += 1
        length -= 1

    cumulative = np.cumsum(group_size)
    cumulative = np.insert(cumulative, 0, 0)
    return cumulative


def get_dx(nuclide, mt, perturbation,
           discretization=None, group=None):
    """Function reads in the hdf5 cross section for a given nuclide, reaction,
    perturbation, discretisation, and discretisation group as well as the same
    cross section from the 'standard' (unperturbed) data library"""
    # Implement this by literally just multiplying the standard cross section by
    # 0.01
    data_lib = '/root/nndc_hdf5'
    Temp = 294
    with h5py.File(os.path.join(data_lib, f'{nuclide}.h5')) as f:
        xs = f[f"/{nuclide}/reactions/reaction_{mt:03}/{Temp}K/xs"][:]
        print(type(xs))
        if discretization is None:
            dx = xs * perturbation
        else:
            group_len = group_len()
            cum_index = get_cum_energy_groups(group_len)
            dx = np.zeros(group_len)
            dx += xs[cum_index[group]:cum_index[group+1]] * perturbation
    return dx, xs


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


def finite_difference(nuclide, mt, perturbation, discretization=None,
                      group=None, structure='partisn'):
    """Function computes the derivative following the perturbation of a cross
    section pertaining to a given nuclide, reaction, perturbation, discretisation,
    and discretisation group for a certain energy group structure"""

    dx, xs = get_dx(nuclide, mt, perturbation,
                    discretization, group)

    # Finding folder with perturbation output data
    if discretization is None:
        id_code = f'mt{mt}-p{perturbation}'
    else:
        id_code = f'mt{mt}-p{perturbation}-d{discretization:03}-g{group+1:03}'
    perturb_home_folder = '/ironbenchmark/perturbed_run_data/'
    output_folder = os.path.join(perturb_home_folder, id_code, 'output')
    perturb_folder = find_power_folder(output_folder)
    # Standard folder
    standard_home_folder = '/ironbenchmark/standard_run/output'
    standard_folder = find_power_folder(standard_home_folder)

    if structure == 'partisn':
        partisn_file = 'partisn_n.csv'
        flux_n = pd.read_csv(os.path.join(standard_folder, partisn_file))
        flux_n_p = pd.read_csv(os.path.join(perturb_folder, partisn_file))
    if structure == 'benchmark':
        bench_file = 'bench.csv'
        flux_n = pd.read_csv(os.path.join(standard_folder, bench_file))
        flux_n_p = pd.read_csv(os.path.join(perturb_folder, bench_file))
    if structure == 'mcnp':
        mcnp_file = 'n3.csv'
        flux_n = pd.read_csv(os.path.join(standard_folder, mcnp_file))
        flux_n_p = pd.read_csv(os.path.join(perturb_folder, mcnp_file))

    return dx, xs


if __name__ == '__main__':
    dx, xs = finite_difference(
        'Fe56', 102, 0.01, discretization=None, group=None, structure='partisn')

# %%
