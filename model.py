# %%
import matplotlib.pyplot as plt
import numpy as np
import openmc
import pandas as pd

import importlib
import os
import re

import utils
import settings

importlib.reload(utils)


def load_model():
    # Specifying settings and source information
    N = settings.N
    run_env = settings.RUN_ENV
    openmc_settings = openmc.Settings()

    openmc_settings.run_mode = "fixed source"

    openmc_settings.particles = int(10**N / 10)
    openmc_settings.batches = 10

    # Getting the source
    g_bins, g_vals, n_bins, n_vals = utils.read_ng_source()
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

    # Specifying tallies

    gamma_tally_mcnp_1 = openmc.Tally()
    neutron_tally_mcnp_3 = openmc.Tally()
    neutron_tally_mcnp_4 = openmc.Tally()
    neutron_tally_mcnp_5 = openmc.Tally()

    gamma_bins_mcnp, neutron_bins_mcnp = utils.get_mcnp_tally_ebins()
    gamma_e_filter_mcnp = openmc.EnergyFilter(gamma_bins_mcnp)
    neutron_e_filter_mcnp = openmc.EnergyFilter(neutron_bins_mcnp)
    gamma_p_filter = openmc.ParticleFilter('photon')
    neutron_p_filter = openmc.ParticleFilter('neutron')

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

    neutron_tally_mcnp_4_surface_filter = openmc.SurfaceFilter(60)
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
    gamma_bins_bench = utils.get_bench_tally_ebins()
    gamma_tally_bench = openmc.Tally()
    gamma_e_filter_bench = openmc.EnergyFilter(gamma_bins_bench)
    gamma_tally_bench.scores = ['flux']
    gamma_tally_bench.filters = [gamma_e_filter_bench,
                                  gamma_p_filter,
                                  gamma_tally_mcnp_1_cell_filter]

    # Adding tallies based on partisn bins
    gamma_bins_partisn, neutron_bins_partisn = \
        utils.get_partisn_tally_ebins()

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

    tallies = openmc.Tallies([gamma_tally_mcnp_1,
                              neutron_tally_mcnp_3,
                              neutron_tally_mcnp_4,
                              neutron_tally_mcnp_5,
                              gamma_tally_bench,
                              gamma_tally_partisn,
                              neutron_tally_partisn])
    tallies.export_to_xml(f"{run_env}/tallies.xml")


def post_process():

    run_env = settings.RUN_ENV
    statepoint_path = f'{run_env}/statepoint.10.h5'
    statepoint = openmc.StatePoint(statepoint_path)

    # Constructing filters for identification
    gamma_bins_mcnp, neutron_bins_mcnp = utils.get_mcnp_tally_ebins()
    gamma_bins_bench = utils.get_bench_tally_ebins()
    gamma_bins_partisn, neutron_bins_partisn = \
        utils.get_partisn_tally_ebins()

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
    N = settings.N
    output_folder = settings.RUN_ENV
    subdir = os.path.join(settings.RUN_ENV, f'output/e{N}')
    print(
        f"This is the directory your openmc simulation output is going to: {subdir}")
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


if __name__ == "__main__":
    load_model()
    post_process()

# %%
