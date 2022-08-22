#!/usr/bin/env python3

import argparse
import h5py
import numpy as np
import glob
import os
import sys
import shutil
import matplotlib.pyplot as plt
from pathlib import Path

import openmc.data
import modify_materials

import importlib
importlib.reload(modify_materials)

description = """
This script generates perturbed cross sections for local sensitivity analysis. Script
generates a cross_sections_perturbed.xml file with the standard library plus the perturbed evaluations.
"""


class CustomFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    pass


parser = argparse.ArgumentParser(
    description=description, formatter_class=CustomFormatter)

parser.add_argument("-n", "--nuclides", nargs="+",
                    default=["Fe56"], help="The nuclide(s) to be sampled")
parser.add_argument("-d", "--destination", default=None,
                    help="Directory to create new library in")
parser.add_argument("-mt", "--cross_sections", default=["2"], nargs="+",
                    help="MT numbers to perturb (default = 2)")
parser.add_argument("-l", "--libdir", default=None,
                    help="Directory of endf library to sample eg. nndc-b7.1-nndc folder")
parser.add_argument("-x", "--xlib", default=None,
                    help="cross_section.xml library to add random evaluations to. Default is OPENMC_CROSS_SECTIONS")
parser.add_argument("-p", "--perturbation", default=0.01,
                    help="perturbation of the xs (default = 0.01)")
parser.add_argument("-t", "--temp", default="294",
                    help="Only format previously sampled files to HDF5")
parser.add_argument("-di", "--discretization", default=None,
                    help="Number of groups to discretize energy bins into")


args = parser.parse_args()
script_dir = Path.cwd()

output_dir = args.destination
if output_dir == None:
    # output_dir = script_dir / "neutron_perturbed"
    output_dir = Path('/root/neutron_perturbed').resolve()
else:
    output_dir = Path(output_dir).resolve()

output_dir.mkdir(exist_ok=True)

libdir = args.libdir
if libdir == None:
    # libdir = script_dir / "neutron" # Just redefining this to be easier
    libdir = Path("/root/nndc_hdf5").resolve()
else:
    libdir = Path(libdir).resolve()

xlib = args.xlib
if xlib == None:
    # xlib = os.getenv("OPENMC_CROSS_SECTIONS") # Don't have combine_libraries.py
    xlib = Path('/root/neutron_perturbed/cross_sections_perturbed.xml').resolve()
else:
    xlib = Path(xlib).resolve()


nuclides = args.nuclides
mts = [int(mt) for mt in args.cross_sections]
perturbation = float(args.perturbation)
Temp = int(args.temp)
discretization = args.discretization
###
#   Path deffs
###

###
#
###


###
#   Parameter deffs
###
REACTION = openmc.data.reaction.REACTION_NAME

lib = openmc.data.DataLibrary()
lib = lib.from_xml(xlib)

def valid_discretization_test():
    test_file = libdir/f"Fe56.h5"
    f = h5py.File(test_file)                
    bin_struct = f[f"/Fe56/energy/{Temp}K"][:]
    bin_len = len(bin_struct)
    f.close()
    if discretization <= 0 or discretization > bin_len:
        print(f'{discretization} is not a valid discretization number')
        sys.exit(None)        
    return bin_len


def get_cum_energy_groups(length):
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

    cumulative = np.cumsum(group_size)
    cumulative = np.insert(cumulative, 0, 0)
    return cumulative



if discretization:
    # Tests if valid discretization number has been input and returns total length of
    # energy groups
    discretization = int(discretization)
    energy_group_length = valid_discretization_test()
    cum_index = get_cum_energy_groups(energy_group_length)

    for nuc in nuclides:
        for MT in mts:
            if not MT in REACTION:
                print('{} is not a valid MT number'.format(MT))
                sys.exit(None)

            reaction = REACTION[MT]

            for i in range(discretization):

                print(
                    f""""Perturbing {nuc} {reaction} (MT = {MT}) at {Temp} by {perturbation}
                    with discretization {discretization} in energy group {i+1}""")

                ###
                #   Perturb nuclide
                ###
                filename = libdir/f"{nuc}.h5"
                filename_new = output_dir/f"{nuc}-mt{MT}-p{perturbation}-d{discretization:03}-g{i+1:03}.h5"
                if os.path.exists(filename_new):
                    print("This perturbation already exists so no need to repeat")
                    continue

                shutil.copyfile(filename, filename_new)         # Make a copy of file.

                f = h5py.File(filename_new, 'r+')                # open new file
                energy = f[f"/{nuc}/energy/{Temp}K"][:]
                cross_section = f[f"/{nuc}/reactions/reaction_{MT:03}/{Temp}K/xs"][:]
                cross_section_perturb = cross_section.copy()
                cross_section_perturb[cum_index[i]:cum_index[i+1]] *= perturbation + 1

                # Overwrite the chosen nuclide
                f[f"/{nuc}/reactions/reaction_{MT:03}/{Temp}K/xs"][:] = cross_section_perturb

                f.close()

                ###
                #   Overwrite the name of the new file, and add to library
                ###
                data = openmc.data.IncidentNeutron.from_hdf5(filename_new)
                data.name = f"{nuc}-mt{MT}-p{perturbation}-d{discretization:03}-g{i+1:03}"
                data.export_to_hdf5(filename_new, "w")

                lib.register_file(filename_new)
else:
    for nuc in nuclides:
        for MT in mts:
            if not MT in REACTION:
                print('{} is not a valid MT number'.format(MT))
                sys.exit(None)

            reaction = REACTION[MT]

            print(
                f"Perturbing {nuc} {reaction} (MT = {MT}) at {Temp} by {perturbation}")

            ###
            #   Perturb nuclide
            ###
            filename = libdir/f"{nuc}.h5"
            filename_new = output_dir/f"{nuc}-mt{MT}-p{perturbation}.h5"
            if os.path.exists(filename_new):
                    print("This perturbation already exists so no need to repeat")
                    continue

            shutil.copyfile(filename, filename_new)         # Make a copy of file.

            f = h5py.File(filename_new, 'r+')                # open new file
            energy = f[f"/{nuc}/energy/{Temp}K"][:]
            cross_section = f[f"/{nuc}/reactions/reaction_{MT:03}/{Temp}K/xs"][:]

            cross_section_perturb = cross_section + cross_section * perturbation

            # Overwrite the chosen nuclide
            f[f"/{nuc}/reactions/reaction_{MT:03}/{Temp}K/xs"][:] = cross_section_perturb

            f.close()

            ###
            #   Overwrite the name of the new file, and add to library
            ###
            data = openmc.data.IncidentNeutron.from_hdf5(filename_new)
            data.name = f"{nuc}-mt{MT}-p{perturbation}"
            data.export_to_hdf5(filename_new, "w")

            lib.register_file(filename_new)


pre = output_dir / "cross_sections_perturbed_pre.xml"
post = output_dir / "cross_sections_perturbed.xml"

# lib.export_to_xml(pre)
# if post.exists():
#     command = f"combine_libraries.py -l {pre} {post} -o {post}"
#     os.system(command)
# else:
#     lib.export_to_xml(post)

lib.export_to_xml(post)
# Uses function from alternative module to automatically generate materials.xml files
# NOTE THIS BREAKS FOR MULTIPLE MT NUMBERS FOR NOW, THIS IS A KNOWN BUG
# Or does it? Check later
modify_materials.main(nuclides, mts, perturbation, discretization)

# pre.unlink()


###
#   Could add plotting option here
###