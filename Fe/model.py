# %%
import matplotlib.pyplot as plt
import numpy as np
import openmc
import pandas as pd

import importlib
import os
import re

import data_load
import settings

importlib.reload(data_load)


def load_model():
    # Specifying settings and source information
    N = settings.N
    run_env = settings.RUN_ENV

    openmc_settings = openmc.Settings()
    openmc_settings.run_mode = "fixed source"
    openmc_settings.particles = int(10**N / 10)
    openmc_settings.batches = 10

    # Getting the source
    g_bins, g_vals, n_bins, n_vals = data_load.read_ng_source()
    g_source = openmc.Source()
    g_source.particle = 'photon'
    g_source.energy = openmc.stats.Tabular(
        g_bins, g_vals, interpolation='histogram')
    n_source = openmc.Source()
    n_source.energy = openmc.stats.Tabular(
        n_bins, n_vals, interpolation='histogram')

    openmc_settings.source = [n_source]
    openmc_settings.photon_transport = True

    openmc_settings.export_to_xml(f"{run_env}/settings.xml")

    ###
    # Specifying tallies
    ###

    # Tally for sensitivity tally (the one that matters most)
    sens_n = openmc.Tally(tally_id=1)
    sens_g = openmc.Tally(tally_id=2)
    sens_n.scores = ['flux']
    sens_g.scores = ['flux']
    gamma_p_filter = openmc.ParticleFilter('photon')
    neutron_p_filter = openmc.ParticleFilter('neutron')
    sens_filter = openmc.EnergyFilter([510000, 2000000])
    sens_n.filters = [openmc.CellFilter(6), neutron_p_filter, sens_filter]
    sens_g.filters = [openmc.CellFilter(6), gamma_p_filter, sens_filter]

    # Complex tallies

    gamma_tally_mcnp_1 = openmc.Tally()
    neutron_tally_mcnp_3 = openmc.Tally()
    neutron_tally_mcnp_4 = openmc.Tally()
    neutron_tally_mcnp_5 = openmc.Tally()

    gamma_bins_mcnp, neutron_bins_mcnp = data_load.get_mcnp_tally_ebins()
    gamma_e_filter_mcnp = openmc.EnergyFilter(gamma_bins_mcnp)
    neutron_e_filter_mcnp = openmc.EnergyFilter(neutron_bins_mcnp)

    gamma_tally_mcnp_1_cell_filter = openmc.CellFilter(52)
    gamma_tally_mcnp_1.scores = ['flux']
    gamma_tally_mcnp_1.filters = [gamma_e_filter_mcnp,
                                  gamma_p_filter,
                                  gamma_tally_mcnp_1_cell_filter]

    neutron_tally_mcnp_3_cell_filter = openmc.CellFilter(70)
    neutron_tally_mcnp_3.scores = ['flux']
    neutron_tally_mcnp_3.filters = [neutron_e_filter_mcnp,
                                    neutron_p_filter,
                                    neutron_tally_mcnp_3_cell_filter]

    neutron_tally_mcnp_4_surface_filter = openmc.SurfaceFilter(89)
    neutron_tally_mcnp_4.scores = ['current']
    neutron_tally_mcnp_4.filters = [neutron_e_filter_mcnp,
                                    neutron_p_filter,
                                    neutron_tally_mcnp_4_surface_filter]

    neutron_tally_mcnp_5_cell_filter = openmc.CellFilter(51)
    neutron_tally_mcnp_5.scores = ['flux']
    neutron_tally_mcnp_5.filters = [neutron_e_filter_mcnp,
                                    neutron_p_filter,
                                    neutron_tally_mcnp_5_cell_filter]

    # Adding tally based on benchmark bins
    gamma_bins_bench = data_load.get_bench_tally_ebins()
    gamma_tally_bench = openmc.Tally()
    gamma_e_filter_bench = openmc.EnergyFilter(gamma_bins_bench)
    gamma_tally_bench.scores = ['flux']
    gamma_tally_bench.filters = [gamma_e_filter_bench,
                                 gamma_p_filter,
                                 gamma_tally_mcnp_1_cell_filter]

    # Adding tallies based on partisn bins
    gamma_bins_partisn, neutron_bins_partisn = \
        data_load.get_partisn_tally_ebins()

    gamma_tally_partisn = openmc.Tally()
    gamma_e_filter_partisn = openmc.EnergyFilter(gamma_bins_partisn)
    gamma_tally_partisn.scores = ['flux']
    gamma_tally_partisn.filters = [
        gamma_e_filter_partisn, gamma_p_filter,
        gamma_tally_mcnp_1_cell_filter]

    neutron_tally_partisn = openmc.Tally()
    neutron_e_filter_partisn = \
        openmc.EnergyFilter(neutron_bins_partisn)
    neutron_tally_partisn.scores = ['flux']
    neutron_tally_partisn.filters = [
        neutron_e_filter_partisn, neutron_p_filter,
        neutron_tally_mcnp_3_cell_filter]

    tallies = openmc.Tallies([
        gamma_tally_mcnp_1,
        neutron_tally_mcnp_3,
        neutron_tally_mcnp_4,
        neutron_tally_mcnp_5,
        gamma_tally_bench,
        gamma_tally_partisn,
        neutron_tally_partisn,
        sens_n, sens_g])
    tallies.export_to_xml(f"{run_env}/tallies.xml")


def process():

    N = settings.N
    run_env = settings.RUN_ENV
    statepoint_path = f'{run_env}/statepoint.10.h5'
    statepoint = openmc.StatePoint(statepoint_path)

    # Constructing filters for identification
    gamma_bins_mcnp, neutron_bins_mcnp = data_load.get_mcnp_tally_ebins()
    gamma_bins_bench = data_load.get_bench_tally_ebins()
    gamma_bins_partisn, neutron_bins_partisn = \
        data_load.get_partisn_tally_ebins()

    id_filter_g1 = openmc.EnergyFilter(gamma_bins_mcnp)
    id_filter_n3 = openmc.EnergyFilter(neutron_bins_mcnp)
    id_filter_n4 = openmc.SurfaceFilter(60)
    id_filter_bench = openmc.EnergyFilter(gamma_bins_bench)
    id_filter_partisn_g = openmc.EnergyFilter(gamma_bins_partisn)
    id_filter_partisn_n = openmc.EnergyFilter(neutron_bins_partisn)
    # Extracting tallies from statepoint file
    tally_g1 = statepoint.get_tally(filters=[id_filter_g1])
    tally_n3 = statepoint.get_tally(
        filters=[id_filter_n3, openmc.CellFilter(70)])
    tally_n4 = statepoint.get_tally(filters=[id_filter_n4])
    tally_bench = statepoint.get_tally(filters=[id_filter_bench])
    tally_partisn_g = statepoint.get_tally(
        filters=[id_filter_partisn_g])
    tally_partisn_n = statepoint.get_tally(
        filters=[id_filter_partisn_n])
    # Getting tally information in dataframes
    df_g1 = tally_g1.get_pandas_dataframe()
    df_n3 = tally_n3.get_pandas_dataframe()
    df_n4 = tally_n4.get_pandas_dataframe()
    df_bench = tally_bench.get_pandas_dataframe()
    df_partisn_g = tally_partisn_g.get_pandas_dataframe()
    df_partisn_n = tally_partisn_n.get_pandas_dataframe()

    # Saving model output for later retrieval
    subdir = os.path.join(run_env, f'output/e{N}')
    filepath_g1 = os.path.join(subdir, 'g1.csv')
    filepath_n3 = os.path.join(subdir, 'n3.csv')
    filepath_n4 = os.path.join(subdir, 'n4.csv')
    filepath_bench = os.path.join(subdir, 'bench.csv')
    filepath_partisn_g = os.path.join(subdir, 'partisn_g.csv')
    filepath_partisn_n = os.path.join(subdir, 'partisn_n.csv')
    if not os.path.exists(os.path.dirname(filepath_g1)):
        os.makedirs(os.path.dirname(filepath_g1))
    df_g1.to_csv(filepath_g1)
    df_n3.to_csv(filepath_n3)
    df_n4.to_csv(filepath_n4)
    df_bench.to_csv(filepath_bench)
    df_partisn_g.to_csv(filepath_partisn_g)
    df_partisn_n.to_csv(filepath_partisn_n)

    # The tallies that matter
    sens_n_tally = statepoint.get_tally(id=1)
    df_sens_n = sens_n_tally.get_pandas_dataframe()
    filepath_sens_n = os.path.join(subdir, 'sens_n.csv')
    df_sens_n.to_csv(filepath_sens_n)
    sens_g_tally = statepoint.get_tally(id=2)
    df_sens_g = sens_g_tally.get_pandas_dataframe()
    filepath_sens_g = os.path.join(subdir, 'sens_g.csv')
    df_sens_g.to_csv(filepath_sens_g)


if __name__ == "__main__":
    settings.N = 7
    settings.MAIN_DIR = '/ironbenchmark/Fe'
    settings.RUN_ENV = '/ironbenchmark/Fe/standard_run'
    load_model()
    process()

# %%
