# %%
import numpy as np
import pandas as pd

import os
import re

import config


def read_ng_source(tegrity=False):
    filepath = os.path.join(
        config.MAIN_DIR, "data/external/sources/cf252_newest_ng-source.txt")
    with open(filepath) as f:
        text = f.read()
    dset = re.compile(r'(?:^\d.*$\n)+', re.MULTILINE)
    datasets = re.findall(dset, text)
    order = ['nbin', 'gbin', 'nsource', 'gsource']
    clean_datasets = []
    for dataset in datasets:
        split = dataset.split()
        clean = np.flip(np.array(split, 'f'))
        clean_datasets.append(clean)
    g_bins = clean_datasets[1]
    g_vals = np.append(clean_datasets[3], 0)
    n_bins = clean_datasets[0]
    # CRUCIAL to divide source by energy bins, as that is the form
    # that openmc understands
    n_vals = np.append(clean_datasets[2] / (n_bins[1:] - n_bins[:-1]), 0)
    if tegrity:
        n_vals = np.append(clean_datasets[2], 0)
    return g_bins, g_vals, n_bins, n_vals


def get_raw_mcnp_tally_ebins(text):
    """
    Auxiliary function to get_mcnp_tally_ebins. Returns all tally
    energy bins found in mcnp input file regardless of duplicates.
    Input:
        Text read directly from MCNP input file
    Returns:
        List of lists. Second order lists are energy bins in MeV.
    """
    energy_format = re.compile(
        r'^e\d+\s+.+$\n(?:^\s+.+$\n)+', re.M)
    e_list = re.findall(energy_format, text)
    refine1 = []
    for entry in e_list:
        refine1.append(re.sub(r'\$.+$\n', '', entry, flags=re.M))
    refine2 = []
    for entry in refine1:
        refine2.append(re.sub(r'\s+', ' ', entry, flags=re.M))
    ls = [entry.split() for entry in refine2]
    energy_bins = []
    for entry in ls:
        energy_bin = []
        for i, element in enumerate(entry):
            if i == 0:
                continue
            if 'i' in element:
                pos = element.find('i')
                if len(element) == 2:
                    n = int(element[pos-1])
                else:
                    n = int(element[pos-2] + element[pos-1])
                prev = float(entry[i-1])
                next = float(entry[i+1])
                interval = (next - prev) / (n+1)
                energy_bin += list(np.linspace(prev+interval,
                                   next-interval, n).round(12))
            else:
                energy_bin.append(float(element))
        energy_bins.append(energy_bin)
    return energy_bins


def get_mcnp_tally_ebins_old(filename, convert=False):
    """
        Returns a tuple of np.arrays of the mcnp input card
        tally energy bins for gammas and neutrons.
        Input:
            str: MCNP input file path
        Returns:
            tuple of np.arrays of gamma bins and neutron bins
        Options:
            set convert=True to convert MeV to eV
    """
    with open(filename) as f:
        text = f.read()
    bins = get_raw_mcnp_tally_ebins(text)
    unduped = []
    for ls in bins:
        if ls not in unduped:
            unduped.append(ls)
    bins_next = unduped[:-1]

    if convert:
        true_bins = [np.array(bin)*10**6 for bin in bins_next]
    else:
        true_bins = [np.array(bin) for bin in bins_next]
    gamma_bins, neutron_bins = true_bins
    gamma_bins = np.append(0, gamma_bins)
    neutron_bins = np.append(0, neutron_bins)
    return gamma_bins, neutron_bins


def get_bench_tally_ebins():
    df_bench1 = read_benchmark_data()[0]
    gamma_bins_bench1 = np.append(
        df_bench1['energy low [eV]'].values, df_bench1['energy high [eV]'].values[-1])
    return gamma_bins_bench1


def get_partisn_tally_ebins():
    df_gamma_partisn, df_neutron_partisn = read_partisn_data()
    gamma_bins_partisn = np.append(
        df_gamma_partisn["energy low [eV]"].values,
        df_gamma_partisn["energy high [eV]"].values[-1])
    neutron_bins_partisn = np.append(
        df_neutron_partisn["energy low [eV]"].values,
        df_neutron_partisn["energy high [eV]"].values[-1])
    return gamma_bins_partisn, neutron_bins_partisn


def get_mcnp_tally_ebins():
    dfs = read_mcnp_data()
    gamma_bins = np.append(dfs['gammas_52']["energy low [eV]"].values[0],
                           dfs['gammas_52']["energy high [eV]"])
    neutron_bins = np.append(dfs['neutrons_70']["energy low [eV]"].values[0],
                             dfs['neutrons_70']["energy high [eV]"])
    return gamma_bins, neutron_bins


def read_benchmark_data():
    """
    Reads the benchmark experimental data files, shows flux values
    for spheres of different diameters at particular energy bins.
    Input:
        No input
    Returns:
        Returns list of dataframes with data for each repeat of
        the experiment

    """
    path1 = os.path.join(
        config.MAIN_DIR, "data/external/benchmarks/benchmark1.txt")
    path2 = os.path.join(
        config.MAIN_DIR, "data/external/benchmarks/benchmark2.txt")
    paths = [path1, path2]

    data = re.compile(r'\d')
    columns = re.compile(r'[a-zA-Z]')
    hyphen = re.compile(r'\s+-\s+')

    dfs = []
    for path in paths:
        numbers = []
        with open(path) as f:
            text = f.readlines()
        for line in text:
            if re.match(columns, line):
                cols = line.strip()
                break
        for line in text:
            if re.match(data, line):
                numbers.append(line)

        cols = re.split('\s{2,}', re.sub(hyphen, r'  ', cols))
        numbers_r1 = []
        for line in numbers:
            numbers_r1.append(re.sub(hyphen, r' ', line).split())

        dfs.append(pd.DataFrame(numbers_r1, columns=cols).astype(np.float64))
        for df in dfs:
            df.rename(columns={'Elow (MeV)': 'energy low [eV]',
                               'Eup (MeV)': 'energy high [eV]'},
                      inplace=True)
        df['energy low [eV]'] *= 1.0e6
        df['energy high [eV]'] *= 1.0e6
        df['Fe30'] *= 1.0e-6
    return dfs


def read_partisn_data():
    """
    Reads the partisn data files, and returns a list of dataframes
    for each tally
    Input:
        No input
    Returns:
        Returns list of dataframes with data for each partisn tally
    """
    final_gamma_df = read_partisn_gamma()
    final_neutron_df = read_partisn_neutron()
    return final_gamma_df, final_neutron_df


def read_partisn_gamma():
    filepath = os.path.join(
        config.MAIN_DIR, "data/external/partisn_outputs/partisn_g_final_n.txt")
    df = pd.read_csv(filepath, delimiter='\t')
    df["Energy high [eV]"] = np.append(df["Energy [eV]"].values[1:], 0)
    df = df.drop(len(df)-1)
    new_df = df[["Energy [eV]", "Energy high [eV]", "Flux", "FI/dE [MeV]"]]
    new_df = new_df.set_axis(["energy low [eV]", "energy high [eV]", "Flux", "FI/dU"],
                             axis=1)
    return new_df


def read_partisn_neutron():
    filepath = os.path.join(
        config.MAIN_DIR, "data/external/partisn_outputs/partisn_n_final_n.txt")
    df = pd.read_csv(filepath, delimiter='\t')
    new_df = df.set_axis(["energy low [eV]", "energy high [eV]", "Flux", "FI/dU"],
                         axis=1)
    return new_df


def read_mcnp_data():
    """
    Returns dictionary of dataframes with energy bin, 
    """
    filepath = os.path.join(
        config.MAIN_DIR, "data/external/mcnp_outputs/mcnp_flux_n.txt")
    with open(filepath) as f:
        text = f.read()
    dataset = re.compile(
        r"^.*energy.*$\n(?:^.*\d.*$\n){3,}", re.MULTILINE)
    datasets = re.findall(dataset, text)

    dfs = {}
    labels = ["gammas_52", "neutron_ring", "neutrons_70"]
    for i, dset in enumerate(datasets):
        lines = dset.splitlines()[1:-1]
        row = re.compile(r"\s+([-\d+E.]+)\s+([-\d+E.]+)\s+([-\d+E.]+)")
        ebins = [1.0e-4]
        fluxes = []
        dFi = []
        for line in lines:
            match = re.match(row, line)
            ebins.append(match.group(1))
            fluxes.append(match.group(2))
            dFi.append(match.group(3))
            elow = np.array(ebins[:-1], dtype='f') * 10**6
            ehigh = np.array(ebins[1:], dtype='f') * 10**6
        df = pd.DataFrame({"energy low [eV]": elow,
                           "energy high [eV]": ehigh,
                           "integral": fluxes,
                           "dFi": dFi}, dtype=np.float64)
        df["FI/dU"] = df["integral"] / np.log(df["energy low [eV]"].values /
                                              df["energy high [eV]"].values)
        dfs[labels[i]] = df
    dfs["gammas_52"] = read_mcnp_gammas()
    return dfs


def read_mcnp_gammas():
    filepath = os.path.join(
        config.MAIN_DIR, "data/external/mcnp_outputs/mcnp_g_final_n.txt")
    df = pd.read_csv(filepath, delimiter='\t')
    df = df[["energy", "integral", "Fi/dE [MeV]"]]
    df["Energy high [eV]"] = np.append(df["energy"].values[1:], 0)
    df = df[["energy", "Energy high [eV]", "integral", "Fi/dE [MeV]"]]
    df = df.drop(len(df)-1)
    df = df.set_axis(["energy low [eV]", "energy high [eV]", "integral", "FI/dU"],
                     axis=1)
    df["energy low [eV]"] *= 10**6
    df["energy high [eV]"] *= 10**6
    return df


if __name__ == "__main__":
    a, b, c, d = read_ng_source()
