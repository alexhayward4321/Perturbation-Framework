# %%
import matplotlib.pyplot as plt
import numpy as np
import openmc
import pandas as pd

import importlib
import os
import re

import model
import processing
import utils
import settings
import modify_materials

# Key examples for run type:
# Overall pattern: sourceInfo_discretisation_discretisationNumber_MTnumberPerturbed
# npd neutron prompt delayed
# gp gamma prompt only
# npdgpd neutron prompt delayed gamma prompt delayed
# We are assuming all perturbations will be of all iron isotopes


def run(N, run_env):
    if os.getcwd() != '/ironbenchmark':
        os.chdir('/ironbenchmark')
    settings.N = N
    settings.RUN_ENV = run_env

    model.load_model()
    openmc.run(cwd=run_env)
    model.post_process()
    processing.main(["output_summary"])


if __name__ == "__main__":
    mt = 102
    perturbation = 0.01
    discretization = None
    standard_run = False
    perturb_folder = '/ironbenchmark/perturbed_run_data/'
    standard_run_folder = '/ironbenchmark/standard_run'

    if standard_run:
        for i in range(6, 8):
            run(i, standard_run_folder)
    elif discretization is None:
        id_code = f'mt{mt}-p{perturbation}'
        run_env = os.path.join(perturb_folder, id_code)
        for i in range(6, 7):
            run(i, run_env)
    else:
        for group in range(discretization):
            id_code = f'mt{mt}-p{perturbation}-d{discretization:03}-g{group+1:03}'
            run_env = os.path.join(perturb_folder, id_code)
            for i in range(6, 8):
                run(i, run_env)


# %%
