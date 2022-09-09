# %%
import automate
import settings

import importlib
importlib.reload(settings)


def run():
    # Telling the system what model to run
    model = 'H1'
    automate.load_model(model, run_env=None)
    # If you are lazy and so when you switch models you want some
    #  default parameters just add stuff to / modify the dictionary
    #  below and access it as necessary like the example for nuclides
    default_nuclides = {'H1': ['H1'],
                        'Fe': ['Fe56'],
                        'Fe-simplified': ['Fe56']}

    # Running the model you want to run
    # automate.main_run(powers=[6], check_repeat=False)
    automate.main_run(powers=[6], nuclides=default_nuclides[model],
                      mts=[2], perturbations=[0.1], check_repeat=False)
    # automate.main_run(powers=[7], mts=[2, 102, 4],
    #          perturbations=[0.01, 0.01, 0.3, 1.0], check_repeat=False)


if __name__ == '__main__':
    run()
