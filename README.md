# Openmc cross section perturbation and sensitivity analysis toolkit

This repository houses a python code used to compute the sensitivity of gamma and neutron fluxes to perturbations of cross section data. It automates the perturbations of cross sections of a set of nuclides, reactions, and perturbation strengths. Further, if you specify a 'base' model in a folder, it automatically generates all of the materials, geometry, settings and tallies files within properly labelled folders for each perturbation, runs openmc simulations for them for a specified number of particles, then stores the results of the simulations as well as processed graphs in a standardised and easily accessible file structure and format. 

The main power of this simulation running structure is indeed the fact that it automates perturbations and openmc simulation calculations for any general model then stores them in an accessible format.

## Disclaimer

It is still in a very early stage of development, so there are doubtless numerous bugs, but if determined a person should be able to pick up how the code works and modify it in order to improve it and make it more generally applicable. Note also that the primary developer of the code was very inexperienced in building large interconnected python projects, and this explains the conspicuous absence of defined objects and classes, as well as some experimental ways of implementing features that were based on spontaneous problem solving while potentially ignorant of standard (and perhaps much simpler) ways of doing things.

If you have any questions please feel free to email me at alexhayward4321@gmail.com.

## Getting Started

First clone the repository using the command:
```
git clone https://github.com/alexhayward4321/Iron <folder_name>
```
inside the directory you would like to perform your data analysis (if you do not have git installed, download it). 

Then open the configuration file config.py located in your newly downloaded repository. The crucial configuration variables you must set are the `HOME_DIR` and `LIBDIR` variables. Read the instructions in the file to find out what they mean. Full instructions on all variables are given in the file, thus no section has been dedicated in this README.md for it, though reading the rest of this document will likely illuminate some strange-seeming decisions or incomplete explanations.

Go to run.py file and specify a model in the eponymous variable according to the main folders defined under the repository directory. Go to the automate.main_run() function and modify the parameters there to configure the run how you want it. Then simply run the file. Congratulations! You have made your first automated monte carlo run.

## How it works

### model folders

Underneath the main repository folder you create folders with files within of a specific format in order to specify your models. Let us work through this with the example model Fe-simplified. Within Fe-simplified you will find several other folders and some python files. Pay attention in particular to the following five files/folders:
- model.py
- post_process.py
- data_load.py
- standard_run
- perturbed_run_data

These are the five folders that are automatically built into the automation of the perturbations, and within them you may need to define functions in a certain way to get them to work how you want them.

#### standard_run 
Firstly, the standard_run folder is where your 'base' model which you would like to run perturbations for is stored. This is where information for the perturbed runs will be copied from and into the perturbed_run_data folder. This is also where a run of your model will automatically go (if running from run.py) if you specify no perturbations, and no separate run_env variable in the load_config() function in run.py. 

#### model.py

The model.py file is where all of the openmc xml file generating functions live. These are, by convention named:
- materials_geometry() loads materials.xml and geometry.xml file for inside a given config.RUN_ENV
- settings() loads config.xml file
- tallies() loads tallies.xml file
- [Optional] 

Do not change this convention for the naming of the settings() and tallies() functions, as those are built into the automation (though technically removing the tallies() function from the automation might be a sensible idea in future). The materials_geometry() function however is not built in to the automation. All of these functions redirect their output xml file to the config.RUN_ENV configuration variable (make sure to do this with new models).

Within this model.py file is also a process() function, this processes the tallies and sends them to a folder called output (within config.RUN_ENV) in csv format, and separates the runs by how many particles were simulated. These e[N] folders mean 10^[N] particles were simulated. It's important to store properly the output after very long monte carlo runs so that you can identify it later, since you would typically want to avoid re-running them. 

You may also want to define other functions, like the plot_model() function I have written. These can be arbitrarily named and won't affect the automation, it's just that to run them you will need to run it from within the model.py folder itself, and make sure your configuration variables at the bottom are properly defined to direct your output where you want it.

#### post_process.py

The post_process.py takes the output from the csv files stored in the output folder and performs processing on them. It has a main() function defined that takes as input a list of the functions that you have defined within the main() function to run as input. automate.py calls the main function without argument, so you can easily change what functions run by default just modifying the func_list default in the post_process.py file itself. In my processing files I have defined plotting functions, but in theory you could define anything you wish and save your output to anywhere you want. I however defined the `utils.plot_log_axes()` function which plots two variables on log-log axes - very convenient for plotting flux spectra. It also has parameters you can specify that can automatically redirect your graphs into your current `config.RUN_ENV` environment in the graphs folder, again in a subdirectory that tells you the number of particles in the simulation. You can also specify filename of course.

#### data_load.py

This is the least interconnected file with the rest of the automation. This is typically only used for loading data output files to compare with your openmc output in the post_processing.py file, or external information needed to specify the openmc model (e.g. tally bin structures, source information). The only thing you may want to do is specify the file location of the data you are extracting relative to config.MAIN_DIR, so that if you decide to share your work on a remote repository and clone it for another PC or your laptop, or for someone else to take a look at you won't have to undo a lot of mess with your file naming (could alternatively specify relative file locations).

#### perturbed_run_data

This is where all the data from your runs with perturbed cross sections go. Types of perturbation are specified according to a code. The code is not difficult to understand. If you look in the Fe_simplified/perturbed_run/data folder the code will read "mt#p#" where the number following mt is the mt number of the reaction, and the number following p is the fractional perturbation of the total cross section of that reaction. Admittedly, perhaps the code should include nuclide as well, depending on the application, but for the applications I was working in the nuclides you perturbed never changed.


### File Descriptions

Below is a list of the files in the repository and a brief description of how they relate to each other and what they do:

- run.py : this is your main run file, and will act as a sort of 'black box' when it it comes to running all of your simulations. 
- automate.py: A module that draws functions from most other modules all together to automate the file and folder construction, cross section perturbations, OpenMC simulation runs, and post-processing of output / plotting of graphs.
- finite_difference.py: runs the senstivity calculations following any perturbation runs and tabulates the results 
- modify_materials.py: copies the materials.xml file from the standard_run folder then modifies it to insert the perturbed materials. Then it places this materials.xml file as well as the geometry.xml file from the standard_run folder in a new folder with a name that allows easy identification of the perturbation performed.
- perturb_xs.py written originally by Ander Gray, and subsequently modified, perturbs the cross sections that are stored in a way that is accessible to openmc.
- utils.py: a few utilities to do with plotting and loading data of general applicability to all models (add more if you want)

#### Within a model

- model.py: Loads the openmc model xml files for the simulations. Also processes data from statepoint file output from openmc model for storage.    
- post_process.py: post-processes output data that was stored from model.processing() function. Used primarily for plotting in the current models defined.
- data_load.py: loads any data you want for comparison with the model you are simulated. Might be PARTISN, MCNP or other output you want to compare. Also might load source (a particle source) information.


### Running a simulation

Once a simulation has been set up, whether by you or someone else, and you have been careful to obey the needed conventions and define the right functions, running automated perturbation runs is relatively simple. Head to the run.py file. the run() function contains two key module functions, the first, automate.load_model() loads configuration variables for the particular model you want to run an in a specific run environment. Second comes the automate.main_run() function. These functions are important enough that a fuller description of their functionality and parameters is given below

**automate.load_model(model, run_env=None)**

Loads the configuration variables needed for modules to know where to direct simulation outputs, graphs, and locate comparison data. Also changes sys.path variable to be able to reload necessary model.py, data_load.py and post_process.py modules when switching between models. The run_env parameter is a crucial understanding point to expand your flexibility in modifying your simulations for experimentation of model conditions. 

The run_env parameter sets the config.RUN_ENV variable. This tells the code where to find the xml files to actually run the simulation, or where to generate them if necessary. Let us say, for example, you would like to see what happens when you make an unperturbed run with a different source to the one you specified, just to check the difference. You could add an if clause to your settings() function to load a certain source instead of the normal one after it checks that the config.RUN_ENV variable, then run the run.py file as usual. The code will automatically create the new run environment folder for you, copy the materials.xml, geometry.xml, settings.xml and tallies.xml files from your standard_run in there, then run your model.settings() and model.tallies() functions to overwrite the ones that have just been copied to add your new source information. 

There are some caveats to this functionality. Since the automate.py file only calls the settings() and tallies() commands while making sequential runs, if you want to change the materials or geometry, you will need to run the materials and geometry generating functions within the model.py file in python separately inside the model you want to simulate and point its output to the config.RUN_ENV environment. As an example, just change what comes after `config.RUN_ENV =` at the bottom of the model.py file within Fe-simplified.

Note, however, that the copying of the xml files from the standard_run folder only occurs when the program cannot detect the filename of the run_env you have specified. So you can change the materials.xml file once and it will stay the same however many times you run run.py

**automate.main_run(powers=[6], nuclides=None, mts=None, perturbations=None,
             discretization=None, check_repeat=True)**

This function performs the actual simulations, output generation, cross section perturbation and etc. When specifying a perturbation run, all of nuclides, mts and perturbations must not be None. Either an error will be raised, or if perturbations is None then a standard_run will be performed

Parameters: 

**powers**: a list of the the powers of 10 of the number of particles you would like to simulate in your simulations.

**nuclides**: list of nuclides to perturb. Nuclides are strings in same format as openmc E.g. ["Fe56", "Fe57", "Fe54].

**mts**: list of mt number of the reactions you want to perturb. ***Note*** the program will not perturb all of these mts and do a single run with them for a given perturbation, but will do successive runs with one mt number perturbed at a time.

**perturbations**: A list of perturbations indicating how much you would like the total cross section of a particular reaction to all the nuclides to be perturbed. E.g. a perturbation of 0.01 will multiply the total cross section of a given nuclide and reaction by 1.01.

**discretization**: an unimplemented half-completed potential feature. See proposed extensions section.

**check_repeat**: checks if a run has already been performed (i.e. does a folder with the corresponding run_env exist), and if so skips all the simulations and perturbations and just loads the graph output.

#### Some illustrative examples:
```
automate.main_run(powers=[6, 7, 8], nuclides=['Fe56', 'Fe57']],
    mts=[2], perturbations=[0.1], check_repeat=False)
```
Will run 3 simulations with three different powers. Each simulation both nuclides Fe56 and Fe57 will each be perturbed by 0.1 for the elastic reaction was has an mt of 2 (i.e. cross sections of mt=2 will be multiplied by 1.1)
```
automate.main_run(powers=[6, 7], nuclides=['Fe56', 'Fe57','Fe58'], 
    mts=[2, 4], perturbations=[0.1, 0.3, 0.01], 
    check_repeat=False)
```
Will have 12 runs, 6 for each power. The first run will involve all nuclides having their elastic (mt=2) cross sections perturbed by 0.1 and the run will have 10^6 particles. The second will be the same excpet the perturbation will now be 0.3, and so on...






