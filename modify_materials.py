# %%
import openmc
import collections as c
import copy

disc = 10
g = 1
perturb = 0.01
MT = 102


def generate_materials_xml(mt, perturbation, discretization, group):
    materials = openmc.Materials.from_xml('openmc_model_data/materials.xml')
    Iron_copy = copy.deepcopy(materials[0])
    nuclides = Iron_copy.nuclides
    nuclide_change_list = ['Fe54', 'Fe56', 'Fe57', 'Fe58']
    Iron_new = openmc.Material(material_id=Iron_copy.id)
    new_nuclides = []
    for nuclide in nuclides:
        if nuclide.name in nuclide_change_list:
            new_nuclide_name = nuclide.name + \
                f'-mt{MT}-p{perturbation}-d{discretization:03}-g{group:03}'
            Iron_new.add_nuclide(
                new_nuclide_name, nuclide.percent, nuclide.percent_type)

        else:
            Iron_new.add_nuclide(nuclide.name, nuclide.percent,
                                nuclide.percent_type)


    Iron_new.set_density(Iron_copy.density_units, Iron_copy.density)
    Iron_new.temperature = Iron_copy.temperature
    print(materials)
    materials[0] = Iron_new
    print(materials)

    materials.export_to_xml('tests/perturbed_materials.xml')











# %%
