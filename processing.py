# %%
import matplotlib.pyplot as plt
import numpy as np
import openmc
import pandas as pd

import importlib
import os

import utils
import settings


def main(funclist):
    N = settings.N

    ###
    # Loading in data
    ###
    # File paths to read data from
    subdir = os.path.join(settings.RUN_ENV, f'output/e{N}')
    filepath_g1 = os.path.join(subdir, 'g1.csv')
    filepath_n3 = os.path.join(subdir, 'n3.csv')
    filepath_n4 = os.path.join(subdir, 'n4.csv')
    filepath_bench = os.path.join(subdir, 'bench.csv')
    filepath_partisn_g = os.path.join(subdir, 'partisn_g.csv')
    filepath_partisn_n = os.path.join(subdir, 'partisn_n.csv')
    # Reading openmc simulation data from csv files
    df_g1 = utils.load_tally(filepath_g1)
    df_n3 = utils.load_tally(filepath_n3)
    df_n4 = utils.load_tally(filepath_n4)
    df_bench = utils.load_tally(filepath_bench)
    df_partisn_g = utils.load_tally(filepath_partisn_g)
    df_partisn_n = utils.load_tally(filepath_partisn_n)
    # Reading other data from files
    benchmark1 = utils.read_benchmark_data()[0]['Fe30']
    benchmark2 = utils.read_benchmark_data()[1]['Fe30']
    # old partisn data btw
    partisn_g, partisn_n = utils.read_partisn_data()
    mcnp_dict = utils.read_mcnp_data()
    mcnp_g = mcnp_dict['gammas_52']
    mcnp_n = mcnp_dict['neutrons_70']
    mcnp_ring_n = mcnp_dict['neutron_ring']

    ###
    # Factors for comparisons between graphs
    ###
    # Factors direct from mcnp input
    kfk_fact_n = 5.54177e+3
    kfk_fact_g = 3.44196e+3
    # Benchmark
    openmc_bench_bins_fact = 1/((df_bench['F/dE']).sum())
    bench_fact = 1/benchmark1.sum()
    openmc_bench_fact = openmc_bench_bins_fact / bench_fact
    # Partisn
    lead = 6
    tail = -3
    openmc_partisn_fact_g = (partisn_g['Flux'].values[lead:tail] /
                             df_partisn_g['mean'].values[lead:tail]).mean()
    # mcnp
    tail = 80
    openmc_mcnp_fact_g = (mcnp_g.iloc[:, 2].values[:125-tail] /
                          df_g1['mean'].values[:125-tail]).mean()


# Comparing openmc simulation and benchmark data

    def inspect_benchmark():
        utils.plot_log_axes(df_bench['mid_bins'],
                            [benchmark1,
                            df_bench['F/dE'].values,
                            benchmark2])
        utils.plot_log_axes(df_bench['mid_bins'],
                            [benchmark1,
                            df_bench['F/dE'].values/1000,
                            benchmark2])
        utils.plot_log_axes([df_bench['mid_bins'], df_bench['mid_bins'],
                            df_bench['mid_bins'], df_g1['mid_bins'], df_partisn_g['mid_bins']],
                            [benchmark1,
                            df_bench['F/dE'].values/1000,
                            benchmark2, df_g1['F/dE']/1000,
                            df_partisn_g['F/dE']/1000],
                            legend=["benchmark1", "openmc_bench_bins", "benchmark2",
                                    "openmc_mcnp_bins", "openmc_partisn_bins"])

    # Comparing Partisn and openmc

    def inspect_partisn():
        utils.plot_log_axes([df_partisn_g['mid_bins'], df_g1['mid_bins'], df_partisn_g['mid_bins']],
                            [partisn_g['Flux'], df_g1['mean'],
                            df_partisn_g['mean']],
                            legend=["partisn", "openmc_mcnp_bins",
                                    "openmc_partisn_bins"],
                            title="Data direct from output")
        utils.plot_log_axes([df_partisn_g['mid_bins'], df_g1['mid_bins'], df_partisn_g['mid_bins']],
                            [partisn_g['Flux'], df_g1['mean'],
                            df_partisn_g['mean']*openmc_partisn_fact_g],
                            legend=["partisn", "openmc_mcnp_bins", "openmc_partisn_bins"])
        utils.plot_log_axes([df_partisn_g['mid_bins'], df_g1['mid_bins'], df_partisn_g['mid_bins']],
                            [partisn_g['Flux']/df_partisn_g['dU'],
                            df_g1['F/dU']*openmc_partisn_fact_g,
                            df_partisn_g['F/dU']*openmc_partisn_fact_g],
                            legend=["partisn", "openmc_mcnp_bins",
                                    "openmc_partisn_bins"],
                            title="Lethargy normalised comparison multiplying openmc results by factor")
        utils.plot_log_axes([df_partisn_g['mid_bins'], df_g1['mid_bins'], df_partisn_g['mid_bins']],
                            [partisn_g['Flux']/df_partisn_g['dE'],
                            df_g1['F/dE']*openmc_partisn_fact_g,
                            df_partisn_g['F/dE']*openmc_partisn_fact_g],
                            legend=["partisn", "openmc_mcnp_bins",
                                    "openmc_partisn_bins"],
                            title="Energy bin width normalised comparison with factor")
        utils.plot_log_axes([df_partisn_n['mid_bins'], df_n3['mid_bins'], df_partisn_n['mid_bins']],
                            [partisn_n.iloc[:, 2], df_n3['mean'],
                            df_partisn_n['mean']],
                            legend=["partisn", "openmc_mcnp_bins", "openmc_partisn_bins"])
        utils.plot_log_axes([df_partisn_n['mid_bins'], df_partisn_n['mid_bins']],
                            [partisn_n.iloc[:, 2],
                            df_partisn_n['mean']],
                            legend=["partisn", "openmc_mcnp_bins", "openmc_partisn_bins"])

    # Comparing MCNP and openmc

    def inspect_mcnp():
        utils.plot_log_axes(df_g1['mid_bins'],
                            [mcnp_g['integral'], df_g1['mean']*openmc_mcnp_fact_g],
                            legend=["mcnp", "openmc_mcnp_bins"])
        utils.plot_log_axes(df_n3['mid_bins'],
                            [mcnp_n['integral'], df_n3['mean']*kfk_fact_n],
                            legend=["mcnp", "openmc_mcnp_bins"])

    # Comparing all gammas using energy difference (in accordance to benchmark)

    def compare_gamma_flux():
        utils.plot_log_axes([df_partisn_g['mid_bins'], df_bench['mid_bins'],
                            df_partisn_g['mid_bins'], df_g1['mid_bins']],
                            [df_partisn_g['F/dE'],
                            benchmark1/openmc_bench_fact,
                            partisn_g.iloc[:, 2]/df_partisn_g['dE'] /
                            openmc_partisn_fact_g,
                            mcnp_g['integral']/df_g1['dE']/openmc_mcnp_fact_g],
                            "overall_comparison_g", N,
                            legend=["openmc", "benchmark", "partisn", "mcnp"],
                            title=f'Comparison of gamma fluxes with 10^{N} particles simulated in openmc')

    # Comparing all neutrons using lethargy

    def compare_neutron_flux():
        utils.plot_log_axes([df_partisn_n['mid_bins'],
                            df_partisn_n['mid_bins'], df_n3['mid_bins']],
                            [df_partisn_n['F/dU'],
                            partisn_n.iloc[:, 2]/df_partisn_n['dU'],
                            mcnp_n['integral']/df_n3['dU']/kfk_fact_n],
                            "overall_comparison_nl", N,
                            legend=[
                            "openmc_partisn_bins", "partisn", "mcnp"],
                            title=f"Comparison of neutron fluxes with 10^{N} particles simulated in openmc")

    def plot_histograms():
        utils.plot_histogram(hist_gamma_bins,
                             df_g1['mean'], 'g1', N)
        hist_gamma_bins = np.append(df_g1["energy low [eV]"].values,
                                    df_g1["energy high [eV]"].values[-1])

    def output_summary():
        utils.plot_log_axes(df_bench['mid_bins'],
                            [benchmark1,
                            df_bench['F/dE'].values/1000,
                            benchmark2],
                            'bench_vs_openmc_fact', N,
                            legend=['benchmark 1',
                                    'openmc_benchmark_bins', 'benchmark 2'],
                            title=f"Openmc and benchmark experiment comparisons, 10^{N} particles")
        compare_gamma_flux()
        compare_neutron_flux()

    for funcstr in funclist:
        func = locals()[funcstr]
        func()


if __name__ == "__main__":
    main(["output_summary"])

# # %%
# inspect_benchmark()
# # %%
# inspect_partisn()
# # %%
# inspect_mcnp()
# # %%
# compare_gamma_flux()
# # %%
# compare_neutron_flux()

# %%
