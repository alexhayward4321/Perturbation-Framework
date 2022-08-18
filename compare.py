# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import importlib
import os
import re

import model
import processing
import utils
import settings


# Key examples for run type:
# Overall pattern: sourceInfo_discretisation_discretisationNumber_MTnumberPerturbed
# npd neutron prompt delayed
# gp gamma prompt only
# npdgpd neutron prompt delayed gamma prompt delayed
# We are assuming all perturbations will be of all iron isotopes


def run(N, output_folder):
    if os.getcwd() != '/ironbenchmark':
        os.chdir('/ironbenchmark')
    settings.N = N
    settings.RUN_TYPE = output_folder
    model.load_model()
    model.run()
    model.post_process()
    processing.main(["output_summary"])


if __name__ == "__main__":
    for i in range(6, 9):
        run(i, "npd_d000_000_MT0")
    



# %%
