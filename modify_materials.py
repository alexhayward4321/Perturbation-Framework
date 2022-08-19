# %%
import openmc

import collections as c
import copy
import os


disc = 10
g = 1
perturb = 0.01
MT = 102


def main(nuclides, mt, perturbation, discretization):
    for i in range(discretization):
        generate_materials_xml(nuclides, mt, perturbation, discretization, i)


def generate_materials_xml(nuclide_list, mt, perturbation, discretization, group):
    mat_folder = '/ironbenchmark/perturbed_run_data/'
    mat_file = f'materials-mt{MT}-p{perturbation}-d{discretization:03}-g{group+1:03}.xml'
    output_file = os.path.join(mat_folder, mat_file)

    materials = openmc.Materials.from_xml('openmc_model_data/materials.xml')
    Iron_copy = copy.deepcopy(materials[0])
    nuclides = Iron_copy.nuclides
    Iron_new = openmc.Material(material_id=Iron_copy.id)
    for nuclide in nuclides:
        if nuclide.name in nuclide_list:
            new_nuclide_name = nuclide.name + \
                f'-mt{mt}-p{perturbation}-d{discretization:03}-g{group+1:03}'
            Iron_new.add_nuclide(
                new_nuclide_name, nuclide.percent, nuclide.percent_type)

        else:
            Iron_new.add_nuclide(nuclide.name, nuclide.percent,
                                 nuclide.percent_type)

    Iron_new.set_density(Iron_copy.density_units, Iron_copy.density)
    Iron_new.temperature = Iron_copy.temperature
    materials[0] = Iron_new
    materials.export_to_xml(output_file)

# %%
