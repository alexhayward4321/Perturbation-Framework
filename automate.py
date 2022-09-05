# %%
import finite_difference
import modify_materials
import settings
import sys
import re
import matplotlib.pyplot as plt
import numpy as np
import openmc
import pandas as pd

import importlib
import os
os.chdir('/ironbenchmark')

importlib.reload(modify_materials)

# Key examples for run  type:
# Overall pattern: sourceInfo_discretisation_discretisationNumber_MTnumberPerturbed
# npd neutron prompt delayed
# gp gamma prompt only
# npdgpd neutron prompt delayed gamma prompt delayed
# We are assuming all perturbations will be of all iron isotopes


def execute_perturbation(mt, perturbation, nuclides=None, discretization=None):
    """
        Runs Ander's modified cross section file perturbation script for specified mt
        and perturbation
    """
    if nuclides is None:
        if settings.MODEL == 'H1':
            nuclides = ['H1']
        else:
            nuclides = ['Fe56']
    nuclides = ' '.join(nuclides)

    if discretization is None:
        command = f"python3 perturb_xs.py -n {nuclides} -mt {mt} -p {perturbation} -x '{settings.XLIB}' \
            -l '{settings.LIBDIR}' -d '{settings.PERTURB_OUTPUT_DIR}'"
    else:
        command = f"python3 perturb_xs.py -n {nuclides} -mt {mt} -p {perturbation} \
            -di {discretization} -x '{settings.XLIB}' -l '{settings.LIBDIR}' \
            -d '{settings.PERTURB_OUTPUT_DIR}'"
    os.system(command)


def run_single(N, run_env, check_repeat):
    settings.N = N
    settings.RUN_ENV = run_env

    if check_repeat:
        output_path = os.path.join(settings.RUN_ENV, f'output/e{N}')
        if os.path.exists(output_path):
            return

    import model
    import post_process
    model.load_model()
    openmc.run(cwd=run_env)
    model.process()
    post_process.main()


def main_run(powers, mt=None, perturbations=None, discretization=None, check_repeat=True,
             run_env=None):
    """
        Runs multiple simulations depending on varying numbers of particles,
         varying perturbations and discretization within those perturbations
    """
    perturb_folder = os.path.join(settings.MAIN_DIR, 'perturbed_run_data')
    standard_run_folder = os.path.join(settings.MAIN_DIR, 'standard_run')

    if run_env is not None:
        run_env = os.path.join(settings.MAIN_DIR, run_env)
        if not os.path.exists(run_env):
            os.makedirs(run_env)
        for i in powers:
            run_single(i, run_env, check_repeat)
    elif mt is None:
        for i in powers:
            run_single(i, standard_run_folder, check_repeat)
    else:
        for perturbation in perturbations:
            if discretization is None:
                execute_perturbation(mt, perturbation)
                modify_materials.main(perturbation=perturbation, mt=mt)
                id_code = f'mt{mt}-p{perturbation}'
                run_env = os.path.join(perturb_folder, id_code)
                for i in powers:
                    run_single(i, run_env, check_repeat)
            else:
                execute_perturbation(mt, perturbation, discretization)
                modify_materials.main(
                    perturbation=perturbation, discretization=discretization)
                for group in range(discretization):
                    id_code = f'mt{mt}-p{perturbation}d{discretization:03}'
                    group_code = f'g{group+1:03}'
                    run_env = os.path.join(perturb_folder, id_code, group_code)
                    for i in powers:
                        run_single(i, run_env, check_repeat)


def load_config(model, n=3):
    """
        Loads various important settings needed by module functions
    """
    os.system(
        'export OPENMC_CROSS_SECTIONS=/root/neutron_perturbed/cross_sections_perturbed.xml')
    if os.uname().sysname == 'Linux':
        # Full path to original unperturbed nuclear data library
        settings.LIBDIR = "/root/nndc_hdf5"
        # File where new cross section file for perturbed data is stored
        settings.XLIB = '/root/neutron_perturbed/cross_sections_perturbed.xml'
        # File where perturbed cross section hdf5 files are stored
        settings.PERTURB_OUTPUT_DIR = '/root/neutron_perturbed'
        # Full path to main folder
        settings.HOME_DIR = '/ironbenchmark'
    elif os.uname().sysname == 'Darwin':
        settings.LIBDIR = "/Users/user1/Documents/Summer Internship 2022/Nuclear Data/endfb71_hdf5"
        settings.XLIB = '/Users/user1/Documents/Summer Internship 2022/Nuclear Data/neutron_perturbed/cross_sections_perturbed.xml'
        settings.PERTURB_OUTPUT_DIR = '/Users/user1/Documents/Summer Internship 2022/Nuclear Data/neutron_perturbed'
        settings.HOME_DIR = '/Users/user1/Documents/Summer Internship 2022/Python code/Iron'
    # which model we are simulating, Fe or H1
    settings.MODEL = model
    settings.MAIN_DIR = f'/ironbenchmark/{model}'
    if model == 'H1':
        settings.n = n

    # So modules in sub directories can find important modules defined in the home directory
    for item in sys.path:
        if re.match(rf'{settings.HOME_DIR}*', item):
            sys.path.remove(item)
    sys.path.append(settings.HOME_DIR)
    if settings.MAIN_DIR not in sys.path:
        sys.path.append(settings.MAIN_DIR)
    if os.getcwd() != settings.HOME_DIR:
        os.chdir(settings.HOME_DIR)


if __name__ == "__main__":
    # Tells you which model folder has all of the information you want and loads settings
    # for that
    load_config('Fe-simplified')
    main_run(powers=[7], mt=2, perturbations=[0.1], check_repeat=False)
    

    # %env OPENMC_CROSS_SECTIONS /root/neutron_perturbed/cross_sections_perturbed.xml


# %%
