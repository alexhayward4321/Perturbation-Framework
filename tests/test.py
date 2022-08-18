# %%
import openmc

openmc.Materials.cross_sections = "/root/Lib_Fe_all_0.01/cross_sections_perturbed.xml"

mat = openmc.Material()
mat.add_nuclide('Fe56-102', 1.0)
mat.set_density('g/cm3', 7.97)

materials = openmc.Materials([mat])
materials.export_to_xml()

# %%
