# %%
import matplotlib.pyplot as plt
import numpy as np
import openmc
import pandas as pd

import importlib
import os
import sys
import re

import config
import modify_materials
import finite_difference

os.chdir('/ironbenchmark')
importlib.reload(finite_difference)


def execute_perturbation(nuclides, mt, perturbation, discretization=None):
    """
        Runs Ander's modified cross section file perturbation script for specified mt
        and perturbation
    """
    if discretization is None:
        command = f"python3 perturb_xs.py -n {' '.join(nuclides)} -mt {mt} -p {perturbation} -x '{config.XLIB}' \
            -l '{config.LIBDIR}' -d '{config.PERTURB_OUTPUT_DIR}'"
    else:
        command = f"python3 perturb_xs.py -n {' '.join(nuclides)} -mt {mt} -p {perturbation} \
            -di {discretization} -x '{config.XLIB}' -l '{config.LIBDIR}' \
            -d '{config.PERTURB_OUTPUT_DIR}'"
    os.system(command)
    # A bit of a hack to let you change data library instantly without error
    os.environ["OPENMC_CROSS_SECTIONS"] = config.XLIB


def run_single(N, check_repeat):
    config.N = N

    if check_repeat:
        output_path = os.path.join(config.RUN_ENV, f'output/e{N}')
        if os.path.exists(output_path):
            print("Run already performed, loading in results...")
            post_process.main()
            return

    model.settings()
    model.tallies()
    openmc.run(cwd=config.RUN_ENV)
    model.process()
    post_process.main()


def main_run(powers=[6], nuclides=None, mts=None, perturbations=None,
             discretization=None, check_repeat=True):
    """
        Runs multiple simulations depending on varying numbers of particles,
         varying perturbations and discretization within those perturbations
    """

    global model
    global post_process
    global data_load
    import model
    import post_process
    import data_load
    importlib.reload(model)
    importlib.reload(post_process)
    importlib.reload(data_load)

    if perturbations is None:
        for i in powers:
            run_single(i, check_repeat)
    else:
        perturb_folder = os.path.join(config.MAIN_DIR,
                                      'perturbed_run_data')

        if nuclides is None or mts is None or perturbations is None:
            raise Exception("ERROR: either nuclides, mts or perturbations \
                parameter has not been specified. \n All must be \
                    specified to perform a perturbation run.")

        for mt in mts:
            for perturbation in perturbations:
                if discretization is None:
                    execute_perturbation(nuclides, mt, perturbation)
                    modify_materials.main(nuclides, mt, perturbation)
                    id_code = f'mt{mt}-p{perturbation}'
                    config.RUN_ENV = os.path.join(perturb_folder, id_code)
                    for i in powers:
                        run_single(i, check_repeat)
                        finite_difference.compare_perturbation(
                            mt, perturbation)
                else:
                    execute_perturbation(
                        nuclides, mt, perturbation, discretization)
                    modify_materials.main(
                        nuclides, mt, perturbation, discretization)
                    for group in range(discretization):
                        id_code = f'mt{mt}-p{perturbation}d{discretization:03}'
                        group_code = f'g{group+1:03}'
                        config.RUN_ENV = os.path.join(
                            perturb_folder, id_code, group_code)
                        for i in powers:
                            run_single(i, check_repeat)


def load_model(model, run_env=None):
    """
        Function used to change model you are simulating from the default.
        Takes as input the folder name of the model you want to test
    """
    # Where our model is
    config.MAIN_DIR = os.path.join(config.HOME_DIR, model)
    # Which sub-folder in our model to run as our openmc environment
    if run_env is None:
        config.RUN_ENV = os.path.join(config.MAIN_DIR, 'standard_run')
    else:
        config.RUN_ENV = os.path.join(config.MAIN_DIR, run_env)
        if not os.path.exists(config.RUN_ENV):
            os.makedirs(config.RUN_ENV)

    # So modules in sub directories can find important modules defined
    #  in the home directory
    match_items = []
    for item in sys.path:
        if re.match(rf'{config.HOME_DIR}*', item):
            match_items.append(item)
    for item in match_items:
        sys.path.remove(item)
    sys.path.append(config.HOME_DIR)
    sys.path.append(config.MAIN_DIR)


if __name__ == "__main__":
    model = 'H1'
    load_model(model)
    default_nuclides = {'H1': ['H1'],
                        'Fe': ['Fe56'],
                        'Fe-simplified': ['Fe56']}

    main_run(powers=[6], check_repeat=False)
