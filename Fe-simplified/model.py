# %%
import matplotlib.pyplot as plt
import numpy as np
import openmc
import pandas as pd

import importlib
import os
import re

import data_load
import config

importlib.reload(data_load)


def materials_geometry():

    Fe = openmc.Material(material_id=1)
    Fe.add_nuclide("Fe56", 0.9138365)
    Fe.add_nuclide("Fe54", 0.05921408)
    Fe.add_nuclide("Fe57", 0.02110447)
    Fe.add_nuclide("Fe58", 0.002808618)
    Fe.set_density("atom/b-cm", 0.0844)

    Air = openmc.Material(material_id=7)
    Air.add_nuclide("N14", 0.07552, 'wo')
    Air.add_nuclide("O16", 0.2314, 'wo')
    Air.add_nuclide("Ar38", 0.0129, 'wo')
    Air.set_density("g/cm3", 0.0012)

    Materials = openmc.Materials([Fe, Air])
    Materials.export_to_xml(f'{config.RUN_ENV}/materials.xml')

    # Creating model geometry
    ball_exterior = openmc.Sphere(r=15, surface_id=60)
    g_detector_inner = openmc.Sphere(r=16.3)
    g_detector_outer = openmc.Sphere(r=16.8)
    n_detector_inner = openmc.Sphere(r=20.5)
    n_detector_outer = openmc.Sphere(r=21.5)
    boundary = openmc.Sphere(r=22.0, boundary_type='vacuum')

    ball = -ball_exterior
    inner_air_layer = +ball_exterior & -g_detector_inner
    g_detector = +g_detector_inner & -g_detector_outer
    middle_air_layer = +g_detector_outer & -n_detector_inner
    n_detector = +n_detector_inner & -n_detector_outer
    outer_air_layer = +n_detector_outer & -boundary

    Fe_ball = openmc.Cell(name="ball", fill=Fe, region=ball, cell_id=6)
    air_layer_inner = openmc.Cell(fill=Air, region=inner_air_layer, cell_id=51)
    detector_g = openmc.Cell(name="detector", fill=Air,
                             region=g_detector, cell_id=52)
    air_layer_middle = openmc.Cell(
        fill=Air, region=middle_air_layer, cell_id=53)
    detector_n = openmc.Cell(fill=Air, region=n_detector, cell_id=70)
    air_layer_outer = openmc.Cell(fill=Air, region=outer_air_layer, cell_id=71)
    Universe = openmc.Universe(cells=[Fe_ball, air_layer_inner, detector_g, air_layer_middle,
                               detector_n, air_layer_outer])
    Geometry = openmc.Geometry(Universe)
    Geometry.export_to_xml(f"{config.RUN_ENV}/geometry.xml")


def settings():
    # Specifying settings and source information
    N = config.N

    openmc_settings = openmc.Settings()
    openmc_settings.run_mode = "fixed source"
    openmc_settings.particles = int(10**N / 10)
    openmc_settings.batches = 10

    # Getting the source
    if os.path.basename(config.RUN_ENV) == 'old_source':
        d = data_load.get_source_data_dict()
        n_bins = d['nPrompt_725g']['x']
        n_vals = d['nPrompt_725g']['y']
        n_source = openmc.Source()
        n_source.energy = openmc.stats.Tabular(
            n_bins, n_vals, interpolation='histogram')
        openmc_settings.source = [n_source]
    elif os.path.basename(config.RUN_ENV) == 'missing_source':
        ...
    else:
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
    openmc_settings.export_to_xml(f"{config.RUN_ENV}/settings.xml")



def tallies():
    # Tally for sensitivity tally (the one that matters most)
    sens_n = openmc.Tally(tally_id=1)
    sens_g = openmc.Tally(tally_id=2)
    sens_n.scores = ['flux']
    sens_g.scores = ['flux']
    gamma_p_filter = openmc.ParticleFilter('photon')
    neutron_p_filter = openmc.ParticleFilter('neutron')
    sens_filter = openmc.EnergyFilter([510000, 2000000])
    sens_n.filters = [openmc.CellFilter(70), neutron_p_filter, sens_filter]
    sens_g.filters = [openmc.CellFilter(52), gamma_p_filter, sens_filter]

    # Complex tallies

    gamma_tally_mcnp_1 = openmc.Tally()
    neutron_tally_mcnp_3 = openmc.Tally()
    neutron_tally_mcnp_4 = openmc.Tally()

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

    neutron_tally_mcnp_4_surface_filter = openmc.SurfaceFilter(60)
    neutron_tally_mcnp_4.scores = ['current']
    neutron_tally_mcnp_4.filters = [neutron_e_filter_mcnp,
                                    neutron_p_filter,
                                    neutron_tally_mcnp_4_surface_filter]

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
        gamma_tally_bench,
        gamma_tally_partisn,
        neutron_tally_partisn,
        sens_n, sens_g])
    tallies.export_to_xml(f"{config.RUN_ENV}/tallies.xml")


def plot_model():
    plot_xy = openmc.Plot()
    plot_xy.width = (50, 50)
    plot_xy.pixels = (2000, 2000)
    plot_xy.colors = {openmc.Cell(6): 'grey',
                      openmc.Cell(51): 'blue',
                      openmc.Cell(52): 'darkseagreen',
                      openmc.Cell(53): 'blue',
                      openmc.Cell(70): 'darkseagreen',
                      openmc.Cell(71): 'blue'}

    plot_xz = openmc.Plot()
    plot_xz.width = (50, 50)
    plot_xz.pixels = (2000, 2000)
    plot_xz.basis = "xz"
    vox_plot = openmc.Plot()
    vox_plot.type = "voxel"
    vox_plot.width = (50, 50, 50)
    vox_plot.pixels = (200, 200, 200)

    # plot_xy.color_by = 'material'
    Plots = openmc.Plots([plot_xy, plot_xz, vox_plot])
    Plots.export_to_xml(f'{config.RUN_ENV}/plots.xml')
    openmc.plot_geometry(cwd=config.RUN_ENV)


def process():

    N = config.N
    statepoint_path = f'{config.RUN_ENV}/statepoint.10.h5'
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
    subdir = os.path.join(config.RUN_ENV, f'output/e{N}')
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
    config.N = 7
    config.RUN_ENV = '/ironbenchmark/Fe-simplified/standard_run'
    config.MAIN_DIR = '/ironbenchmark/Fe-simplified'
    if not os.path.exists(config.RUN_ENV):
        os.makedirs(config.RUN_ENV)
    materials_geometry()
    settings()
    tallies()
    # plot_model()
    openmc.run(cwd=config.RUN_ENV)
    process()

# %%
