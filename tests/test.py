# %%
import numpy as np
import openmc

openmc.Materials.cross_sections = "/root/Lib_Fe_all_0.01/cross_sections_perturbed.xml"

mat = openmc.Material()
mat.add_nuclide('Fe56-102', 1.0)
mat.set_density('g/cm3', 7.97)

materials = openmc.Materials([mat])
materials.export_to_xml()

# %%


length = 404
ones = np.ones(length)
perturbation = 0.01
discretization = 10

group_size_guide = length // discretization
group_size = []
while length > group_size_guide:
    length -= group_size_guide
    group_size.append(group_size_guide)

i = 0
while length > 0:
    group_size[i] += 1
    i += 1
    length -= 1

print(group_size)
print(np.array(group_size).sum())
cumulative = np.cumsum(group_size)
cumulative = np.insert(cumulative, 0, 0)
print(cumulative)


perturbed_groups = []
for i in range(discretization):
    perturbed_group = ones.copy()
    perturbed_group[cumulative[i]:cumulative[i+1]] *= perturbation + 1
    perturbed_groups.append(perturbed_group)

print(perturbed_groups)


# %%
