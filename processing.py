# %%
import matplotlib.pyplot as plt
import numpy as np
import openmc
import pandas as pd

import importlib
import os

import utils

importlib.reload(utils)
N = 8

# utils.settings.OTHERFILE = 'gamma_and_neutron'

# File paths to read data from
filepath_g1 = f'output/e{N}/g1.csv'
filepath_n3 = f'output/e{N}/n3.csv'
filepath_n4 = f'output/e{N}/n4.csv'
filepath_bench1 = f'output/e{N}/bench1.csv'
filepath_partisn_g = f'output/e{N}/partisn_g.csv'
filepath_partisn_n = f'output/e{N}/partisn_n.csv'
# Reading openmc simulation data from csv files
df_g1 = pd.read_csv(filepath_g1)
df_n3 = pd.read_csv(filepath_n3)
df_n4 = pd.read_csv(filepath_n4)
df_bench1 = pd.read_csv(filepath_bench1)
df_partisn_g = pd.read_csv(filepath_partisn_g)
df_partisn_n = pd.read_csv(filepath_partisn_n)
# Reading other data from files
benchmark1 = utils.read_benchmark_data()[0]['Fe30']
benchmark2 = utils.read_benchmark_data()[1]['Fe30']
# old partisn data btw
partisn_g, partisn_n = utils.read_partisn_data()
mcnp_dict = utils.read_mcnp_data()
mcnp_g = mcnp_dict['gammas_52']
mcnp_n = mcnp_dict['neutrons_70']
mcnp_ring_n = mcnp_dict['neutron_ring']

# Normalisation factors for the fluxes (assuming units in cm)
D = 30
R = D/2
g1_norm = 4 * np.pi * (R + 1.55)**2
R_n = 21.0
n3_norm = 4 * np.pi * R_n**2


# %%
# Defining various useful energy bin values
# OpenMC
mid_g_bins_openmc = (df_g1['energy low [eV]'] + df_g1['energy high [eV]']) / 2
mid_n_bins_openmc = (df_n3['energy low [eV]'] + df_n3['energy high [eV]']) / 2
lethargy_g_openmc = np.log(df_g1['energy high [eV]'].values) - \
    np.log(df_g1['energy low [eV]'].values)
diff_g_openmc = df_g1['energy high [eV]'].values - \
    df_g1['energy low [eV]'].values
lethargy_n_openmc = np.log(df_n3['energy high [eV]'].values) - \
    np.log(df_n3['energy low [eV]'].values)
# Experimental benchmark (gammas only)
mid_g_bins_bench = (df_bench1['energy low [eV]'] +
                    df_bench1['energy high [eV]']) / 2
lethargy_g_bench = np.log(df_bench1['energy high [eV]'].values) - \
    np.log(df_bench1['energy low [eV]'].values)
diff_g_bench = df_bench1['energy high [eV]'].values - \
    df_bench1['energy low [eV]'].values
# Partisn
mid_g_bins_partisn = (partisn_g['energy high [eV]'].values +
                      partisn_g['energy low [eV]'].values) / 2
mid_n_bins_partisn = (partisn_n['energy high [eV]'].values +
                      partisn_n['energy low [eV]'].values) / 2
lethargy_g_partisn = np.log(partisn_g['energy high [eV]'].values) - \
    np.log(partisn_g['energy low [eV]'].values)
diff_g_partisn = partisn_g['energy high [eV]'].values - \
    partisn_g['energy low [eV]'].values
lethargy_n_partisn = np.log(partisn_n['energy high [eV]'].values) - \
    np.log(partisn_n['energy low [eV]'].values)

mid_g_bins_mcnp = (mcnp_g["energy low [eV]"].values
                   + mcnp_g["energy high [eV]"].values) / 2
mid_n_bins_mcnp = (mcnp_n["energy low [eV]"].values
                   + mcnp_n["energy high [eV]"].values) / 2
lethargy_g_mcnp = np.log(mcnp_g['energy high [eV]'].values) - \
    np.log(mcnp_g['energy low [eV]'].values)
diff_g_mcnp = mcnp_g['energy high [eV]'].values - \
    mcnp_g['energy low [eV]'].values
lethargy_n_mcnp = np.log(mcnp_n['energy high [eV]'].values) - \
    np.log(mcnp_n['energy low [eV]'].values)
# bins just for histogram
hist_gamma_bins = np.append(df_g1["energy low [eV]"].values,
                            df_g1["energy high [eV]"].values[-1])

# %%
# Direct plotting of openmc tally results
utils.plot_log_axes(mid_g_bins_openmc, df_g1['mean'], 'g1', N)
utils.plot_log_axes(mid_n_bins_openmc, df_n3['mean'], 'n3', N)
utils.plot_log_axes(mid_n_bins_openmc, df_n4['mean'], 'n4', N)
utils.plot_log_axes(mid_g_bins_bench, df_bench1['mean'], 'bench1', N)
# Plotting histogram(s)
utils.plot_histogram(hist_gamma_bins,
                     df_g1['mean'], 'g1', N)

# %%
# Comparing openmc simulation and benchmark data
utils.plot_log_axes([mid_g_bins_bench, mid_g_bins_bench, mid_g_bins_bench],
                    [benchmark1,
                     df_bench1['mean'].values/diff_g_bench,
                     benchmark2],
                    'bench_vs_openmc', N)

openmc_bench_bins_fact = 1/((df_bench1['mean']/diff_g_bench).sum())
bench_fact = 1/benchmark1.sum()
openmc_bench_fact = openmc_bench_bins_fact / bench_fact

utils.plot_log_axes([mid_g_bins_bench, mid_g_bins_bench, mid_g_bins_bench],
                    [benchmark1,
                     df_bench1['mean'].values/diff_g_bench/1000,
                     benchmark2],
                    'bench_vs_openmc_fact', N)

utils.plot_log_axes([mid_g_bins_bench, mid_g_bins_bench, mid_g_bins_bench, mid_g_bins_openmc, mid_g_bins_partisn],
                    [benchmark1,
                     df_bench1['mean'].values/diff_g_bench/1000,
                     benchmark2, df_g1['mean']/diff_g_openmc/1000, df_partisn_g['mean']/diff_g_partisn/1000],
                    'bench_vs_openmc_fact', N)


# %%
# Comparing Partisn and openmc
# Factors
lead = 6
openmc_partisn_fact_g = (partisn_g.iloc[:, 2].values[lead:] /
                         df_partisn_g['mean'].values[lead:]).mean()
kfk_fact_n = 5.54177e+3
kfk_fact_g = 3.44196e+3

utils.plot_log_axes([mid_g_bins_partisn, mid_g_bins_openmc, mid_g_bins_partisn],
                    [partisn_g.iloc[:, 2], df_g1['mean'],
                    df_partisn_g['mean']],
                    "partisn_vs_openmc_g", N,
                    legend=["partisn", "openmc_mcnp_bins",
                            "openmc_partisn_bins"],
                    title="Data direct from output")
utils.plot_log_axes([mid_g_bins_partisn, mid_g_bins_openmc, mid_g_bins_partisn],
                    [partisn_g.iloc[:, 2], df_g1['mean'],
                    df_partisn_g['mean']*openmc_partisn_fact_g],
                    "partisn_vs_openmc_g", N,
                    legend=["partisn", "openmc_mcnp_bins", "openmc_partisn_bins"])
utils.plot_log_axes([mid_g_bins_partisn, mid_g_bins_openmc, mid_g_bins_partisn],
                    [partisn_g.iloc[:, 2]/lethargy_g_partisn,
                     df_g1['mean']*openmc_partisn_fact_g/lethargy_g_openmc,
                    df_partisn_g['mean']*openmc_partisn_fact_g/lethargy_g_partisn],
                    "partisn_vs_openmc_gl", N,
                    legend=["partisn", "openmc_mcnp_bins",
                            "openmc_partisn_bins"],
                    title="Lethargy normalised comparison multiplying openmc results by factor")
utils.plot_log_axes([mid_g_bins_partisn, mid_g_bins_openmc, mid_g_bins_partisn],
                    [partisn_g.iloc[:, 2]/diff_g_partisn,
                     df_g1['mean']*openmc_partisn_fact_g/diff_g_openmc,
                    df_partisn_g['mean']*openmc_partisn_fact_g/diff_g_partisn],
                    "partisn_vs_openmc_nd", N,
                    legend=["partisn", "openmc_mcnp_bins",
                            "openmc_partisn_bins"],
                    title="Energy bin width normalised comparison with factor")
utils.plot_log_axes([mid_n_bins_partisn, mid_n_bins_openmc, mid_n_bins_partisn],
                    [partisn_n.iloc[:, 2], df_n3['mean'],
                    df_partisn_n['mean']],
                    "partisn_vs_openmc_n", N,
                    legend=["partisn", "openmc_mcnp_bins", "openmc_partisn_bins"])
utils.plot_log_axes([mid_n_bins_partisn, mid_n_bins_partisn],
                    [partisn_n.iloc[:, 2],
                    df_partisn_n['mean']*kfk_fact_n],
                    "partisn_vs_openmc_n", N,
                    legend=["partisn", "openmc_mcnp_bins", "openmc_partisn_bins"])


# %%
# Comparing MCNP and openmc
# Factors
tail = 80
openmc_mcnp_fact_g = (mcnp_g.iloc[:, 2].values[:125-tail] /
                      df_g1['mean'].values[:125-tail]).mean()
# openmc_mcnp_fact_g = 8.8
utils.plot_log_axes([mid_g_bins_mcnp, mid_g_bins_mcnp],
                    [mcnp_g['integral'], df_g1['mean']*openmc_mcnp_fact_g],
                    "mcnp_vs_openmc_g", N,
                    legend=["mcnp", "openmc_mcnp_bins"])
utils.plot_log_axes([mid_n_bins_mcnp, mid_n_bins_mcnp],
                    [mcnp_n['integral'], df_n3['mean']*kfk_fact_n],
                    "mcnp_vs_openmc_n", N,
                    legend=["mcnp", "openmc_mcnp_bins"])

# %%
# Comparing multiple all gammas using energy difference (in accordance to benchmark)
importlib.reload(utils)
utils.plot_log_axes([mid_g_bins_partisn, mid_g_bins_bench,
                     mid_g_bins_partisn, mid_g_bins_mcnp],
                    [df_partisn_g['mean']/diff_g_partisn,
                     benchmark1/openmc_bench_fact,
                     partisn_g.iloc[:, 2]/diff_g_partisn /
                     openmc_partisn_fact_g,
                     mcnp_g['Flux']/diff_g_mcnp/openmc_mcnp_fact_g],
                    "overall_comparison_g", N,
                    legend=["openmc", "benchmark", "partisn", "mcnp"])

utils.plot_log_axes([
    mid_g_bins_partisn, mid_g_bins_mcnp],
    [partisn_g.iloc[:, 2]/diff_g_partisn/openmc_partisn_fact_g,
     mcnp_g['Flux']/diff_g_mcnp/openmc_mcnp_fact_g],
    "partsin_openmc_g", N,
    legend=["partisn", "mcnp"])


# Comparing all gammas using energy difference

# %%

# Comparing all neutrons using lethargy

utils.plot_log_axes([mid_n_bins_openmc, mid_n_bins_partisn,
                     mid_n_bins_partisn, mid_n_bins_mcnp],
                    [df_n3['mean']/lethargy_n_openmc*kfk_fact_n,
                     df_partisn_n['mean']/lethargy_n_partisn*kfk_fact_n,
                     partisn_n.iloc[:, 2]/lethargy_n_partisn,
                     mcnp_n['Flux']/lethargy_n_mcnp],
                    "overall_comparison_nl", N,
                    legend=["openmc_mcnp_bins",
                            "openmc_partisn_bins", "partisn", "mcnp"])

utils.plot_log_axes([mid_n_bins_partisn,
                     mid_n_bins_partisn, mid_n_bins_mcnp],
                    [df_partisn_n['mean']/lethargy_n_partisn*kfk_fact_n,
                     partisn_n.iloc[:, 2]/lethargy_n_partisn,
                     mcnp_n['Flux']/lethargy_n_mcnp],
                    "overall_comparison_nl", N,
                    legend=[
                    "openmc_partisn_bins", "partisn", "mcnp"])

utils.plot_log_axes([mid_n_bins_partisn,
                     mid_n_bins_partisn, mid_n_bins_mcnp],
                    [df_partisn_n['mean']/lethargy_n_partisn*kfk_fact_n,
                     partisn_n.iloc[:, 2]/lethargy_n_partisn,
                     mcnp_n['Flux']/lethargy_n_mcnp],
                    "overall_comparison_nl", N,
                    legend=[
                    "openmc_partisn_bins", "partisn", "mcnp"])

# %%

utils.plot_log_axes(
    [partisn_ng_fendl['energy'], partisn_n_fendl['energy'], partisn_n_endfb7['energy']],
    [partisn_ng_fendl['Flux'], partisn_n_fendl['Flux'], partisn_n_endfb7['Flux']],
    'comparing_partisn_data', N,
    legend=['partisn_ng_fendl', 'partisn_n_fendl', 'partisn_n_endfb7'])

utils.plot_log_axes(
    [mid_n_bins_partisn[:-1], mid_n_bins_partisn, mid_n_bins_partisn[:-1]],
    [partisn_ng_fendl['Flux'], partisn_n_fendl['Flux'].values[:-1],
        partisn_n_endfb7['Flux']],
    'comparing_partisn_data', N,
    legend=['partisn_ng_fendl', 'partisn_n_fendl', 'partisn_n_endfb7'])

# print(all(partisn_ng_fendl == partisn_n_endfb7))

# %%
# utils.plot_log_axes(
#     [mid_n_bins_partisn, mid_n_bins_partisn],
#     [ng_source_n, n_source_n['mean']],
#     # 'NG_AND_N_SOURCE_COMPARISON', N,
#     legend=['ng_source', 'n_source'], other=True)

# %%
