# %%
import re
import os
import sys

# ----------------------------
# Crucial configuration variables

# Home folder under which the whole repository sits
HOME_DIR = '/ironbenchmark'
# Standard cross section library to use when creating perturbed
#  cross sections
LIBDIR = '/root/jeff33_hdf5'


# ------------------------
# Defaults you may wish to change (opportunities to change elsewhere in
#  model)

#  MAIN_DIR refers to the model environment where you store a base
#  openmc model with the requisite openmc .xml files inside the
#  'standard_run' folder. Perturbation folders are automatically copied
#  and modified from that standard run folder and put elsewhere inside
#  your MAIN_DIR folder.
MAIN_DIR = os.path.join(HOME_DIR, 'Fe-simplified')

# --------------------
# Defaults that you amost likely will not want / need to change

# Default place to send perturbed cross section files
PERTURB_OUTPUT_DIR = LIBDIR + '_perturbed'
# cross_section.xml file for OPENMC to use with the location of all
#  unperturbed and perturbed files
XLIB = os.path.join(PERTURB_OUTPUT_DIR, 'cross_sections_perturbed.xml')
# By default the code will run your unperturbed model
RUN_ENV = os.path.join(MAIN_DIR, 'standard_run')

# Determines how xml files will be generated each time when you do a run
# with a run_env that is neither the standard_run folder or a perturbation
# of your standard_run_folder. Three options for each xml file:
# - File directly copied from your standard run folder (every run)
# - File generated using function(s) in your model.py file
# - Leave the files as they are (allows for by-hand modification of
#   existing xml files)


# --------------------------------
# Automatically setting sys.path variables for module location,
#  OPENMC_CROSS_SECTIONS variable and changing your directory to
#  the home directory
if os.path.exists(XLIB):
    os.environ["OPENMC_CROSS_SECTIONS"] = XLIB
else:
    os.environ["OPENMC_CROSS_SECTIONS"] = os.path.join(
        LIBDIR, 'cross_sections.xml')


if os.getcwd() != HOME_DIR:
    os.chdir(HOME_DIR)
# Cleaning the sys.path pallette
for item in sys.path:
    if re.match(rf'{HOME_DIR}*', item):
        sys.path.remove(item)
sys.path.append(HOME_DIR)
sys.path.append(MAIN_DIR)


# ----------------------------------
# Add any further desired global configuration variables below
