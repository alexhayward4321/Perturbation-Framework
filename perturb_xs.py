#!/usr/bin/env python3
# %%
import re
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
import settings

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


def inelastic(nuc):
    # Modifying mts for MT = 4
    filename = libdir/f"{nuc}.h5"
    f = h5py.File(filename, 'r')
    groups = f[f"/{nuc}/reactions"]
    groups = groups.keys()
    p = re.compile('\w+_(\d{3})')
    inelastic_mt = []
    for reaction in groups:
        m = p.match(reaction)
        if m:
            inelastic_mt.append(m.group(1))
    inelastic_mt = [int(r) for r in inelastic_mt]
    # for sorting and removing duplicates
    inelastic_mt = list(set(inelastic_mt))
    final = []
    for mt in inelastic_mt:
        if mt > 50 and mt < 92:
            final.append(mt)
    final = [int(r) for r in final]
    f.close()
    return final


###
#   Parameter deffs
###
REACTION = openmc.data.reaction.REACTION_NAME


lib = openmc.data.DataLibrary()
if not os.path.exists(xlib):
    lib = lib.from_xml(libdir/'cross_sections.xml')
else:
    lib = lib.from_xml(xlib)


for nuc in nuclides:

    for MT in mts:
        if not MT in REACTION and MT != 4:
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

        # Make a copy of file.
        shutil.copyfile(filename, filename_new)

        f = h5py.File(filename_new, 'r+')                # open new file
        energy = f[f"/{nuc}/energy/{Temp}K"][:]

        if MT == 4:
            inelastic_mts = inelastic(nuc)
            for imt in inelastic_mts:
                cross_section = f[f"/{nuc}/reactions/reaction_{imt:03}/{Temp}K/xs"][:]
                cross_section_perturb = cross_section + cross_section * perturbation
                f[f"/{nuc}/reactions/reaction_{imt:03}/{Temp}K/xs"][:] = cross_section_perturb
        else:
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

post = output_dir / "cross_sections_perturbed.xml"
lib.export_to_xml(post)


# %%
