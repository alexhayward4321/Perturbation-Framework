# %%
import openmc

import collections as c
import copy
import os
import shutil


disc = 10
g = 1
perturb = 0.01
MT = 102


def main(nuclides, mt, perturbation, discretization=None):
    if discretization is None:
        generate_materials_xml(nuclides, mt, perturbation)
    else:
        for i in range(discretization):
            generate_materials_xml(
                nuclides, mt, perturbation, discretization, i)


def generate_materials_xml(nuclide_list, mt, perturbation,
                           discretization=None, group=None):
    # Creating openmc simulation run folder and specifying geometry and material xml file paths
    geom_file = '/ironbenchmark/standard_run/geometry.xml'
    mat_file = '/ironbenchmark/standard_run/materials.xml'
    perturb_folder = '/ironbenchmark/perturbed_run_data/'
    if discretization is None:
        id_code = f'mt{MT}-p{perturbation}'
    else:
        id_code = f'mt{MT}-p{perturbation}-d{discretization:03}-g{group+1:03}'
    folder_path = os.path.join(perturb_folder, id_code)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    output_file = os.path.join(folder_path, "materials.xml")
    shutil.copyfile(geom_file, os.path.join(folder_path, 'geometry.xml'))

    materials = openmc.Materials.from_xml(mat_file)
    Iron_copy = copy.deepcopy(materials[0])
    nuclides = Iron_copy.nuclides
    Iron_new = openmc.Material(material_id=Iron_copy.id)
    for nuclide in nuclides:
        if nuclide.name in nuclide_list:
            new_nuclide_name = nuclide.name + '-' +\
                id_code
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
