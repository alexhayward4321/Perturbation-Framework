# %%
import openmc

import collections as c
import copy
import os
import shutil

import config


disc = 10
g = 1
perturb = 0.01
MT = 102


def main(nuclides, mt, perturbation, discretization=None):
    if discretization is None:
        create_folder_env(nuclides, mt, perturbation)
    else:
        for i in range(discretization):
            create_folder_env(
                nuclides, mt, perturbation, discretization, i)


def create_folder_env(nuclides, mt, perturbation,
                      discretization=None, group=None):
    # Creating openmc simulation run folder and specifying geometry and material xml file paths
    geom_file = os.path.join(config.MAIN_DIR, 'standard_run/geometry.xml')
    mat_file = os.path.join(config.MAIN_DIR, 'standard_run/materials.xml')
    perturb_folder = os.path.join(config.MAIN_DIR, 'perturbed_run_data/')
    if discretization is None:
        id_code = f'mt{mt}-p{perturbation}'
        folder_path = os.path.join(perturb_folder, id_code)
    else:
        id_code = f'mt{mt}-p{perturbation}d{discretization:03}'
        group_code = f'g{group+1:03}'
        folder_path = os.path.join(perturb_folder, id_code, group_code)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    output_file = os.path.join(folder_path, "materials.xml")
    shutil.copyfile(geom_file, os.path.join(folder_path, 'geometry.xml'))

    # creating new perturbed materials file from the standard one
    materials = openmc.Materials.from_xml(mat_file)
    new_materials = openmc.Materials()
    for material in materials:
        new_mat = openmc.Material(material_id=material.id)
        for nuclide in material.nuclides:
            if nuclides is None:
                new_mat.add_nuclide(
                    nuclide.name, nuclide.percent, nuclide.percent_type)
            elif nuclide.name in nuclides:
                new_nuclide_name = nuclide.name + '-' + id_code
                new_mat.add_nuclide(
                    new_nuclide_name, nuclide.percent, nuclide.percent_type)
            else:
                new_mat.add_nuclide(
                    nuclide.name, nuclide.percent, nuclide.percent_type)
        new_mat.set_density(material.density_units, material.density)
        new_mat.temperature = material.temperature
        new_materials.append(new_mat)
    new_materials.export_to_xml(output_file)


if __name__ == '__main__':
    main(nuclides=None, mt=102, perturbation=0.01, discretization=None)
    # %%
