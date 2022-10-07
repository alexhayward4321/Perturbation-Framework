# Openmc Cross Section Perturbation and Sensitivity Analysis Toolkit

## Description

This repository houses a python code used to compute the sensitivity of gamma and neutron fluxes to perturbations in nuclear cross section data for a user-specified 'base' physical model in [OpenMC](https://docs.openmc.org/en/stable/). 

## Why It's Useful

The power of this toolkit is you can import into it an OpenMC model you have built with the python API, then automate the hdf5 nuclear data file modification; OpenMC input file preparation; running of simulations; output post-processing; and clearly labelled storage of output files / graphs for easy access. All these processes occur for your physical model according to:
- How many particles you want to simulate
- What nuclides you want to perturb
- For what reactions you would like to perturb
- How strongly you would like to perturb them

Not only is this automated for a single run, but you can specify multiple perturbation strengths for multiple reaction types and multiple numbers of particles that will perform several simulations from a single run of a command. 

## Disclaimer

This is still in a very early stage of development and written by a relatively inexperienced developer. Expect some bugs, unimplemented features and convoluted existing feature implementation. There should nevertheless be sufficient information for a determined a person to pick up how the code works and improve it. 

If you have any questions please feel free to email me at alexhayward4321@gmail.com.

## Getting Started

First clone the repository using the command:
```
git clone https://github.com/alexhayward4321/Iron <folder_name>
```
inside the directory you would like to perform your data analysis (if you do not have git installed, install it [here](https://git-scm.com/downloads)). 

Then open the configuration file config.py located in your newly downloaded repository. The crucial configuration variables you must set are the `HOME_DIR` and `LIBDIR` variables. Full instructions on all variables are given in the file, thus no section has been dedicated in this README for it, though reading the rest of this document will likely illuminate some incomplete explanations.

Go to the run.py file and specify a model in the eponymous variable according to the main folders defined under the repository directory. Go to the `automate.main_run()` function and modify the parameters there to configure the run how you want it. Then simply run the file. Congratulations! You have made your first automated monte carlo run.

## How it works

### Model folders

Underneath the main folder where you installed this repository (`config.MAIN_DIR`) you create folders for each model you would like to run simulations for with files within of a specific format to facilitate the automation. Let us work through this with the example model Fe-simplified. Fe-simplified is an openmc model of an iron sphere with a fissile source in the centre, and several tallies specified around it. Within Fe-simplified you will find several other folders and some python files. Pay attention in particular to the following five files/folders:
- model.py
- post_process.py
- data_load.py
- standard_run
- perturbed_run_data

These are the five folders that are automatically built into the automation of the perturbations, and within the python files you may need to define functions in a certain way to get them to work how you want them.

#### **standard_run**
Firstly, the standard_run folder is where your 'base' model which you would like to run perturbations for is stored. This is where information for the perturbed runs will be copied from and into the perturbed_run_data folder. This is also where a run of your model will automatically go (if running from run.py) if you specify no perturbations, and no separate run_env variable in the `automate.load_model()` function in run.py. 

#### **model.py**

The model.py file is where all of the openmc xml file generating functions live. These are, by convention named:
- materials_geometry() loads materials.xml and geometry.xml file for inside a given `config.RUN_ENV`
- settings() loads config.xml file
- tallies() loads tallies.xml file

Do not change the convention for the naming of the `settings()` and `tallies()` functions specifically, as those are built into the automation (though technically removing the tallies() function from the automation might be a sensible idea in future as it is not strictly necessary). The `materials_geometry()` function however is not built in to the automation, so you may choose to name this anything you wish. All of these functions redirect their output xml file to the `config.RUN_ENV` configuration variable.

Within this model.py file is also a `process()` function, this processes the tallies fromt the statepoint file after a simulation and sends them to a folder called output (within `config.RUN_ENV`) in csv format, and separates the runs by how many particles were simulated. These e[N] folders mean 10^[N] particles were simulated. It's important to store properly the output after very long monte carlo runs so that you can identify it later, since you would typically want to avoid re-running them. Note that statepoint files are not currently stored in the output folders.

You may also want to define other functions, like the `plot_model()` function I have written. These can be arbitrarily named and won't affect the automation, it's just that to run them you will need to run the model.py folder itself, and make sure your configuration variables below `if __name__=='__main__':` are properly defined to direct your output where you want it.

#### post_process.py

The post_process.py script takes the output from the csv files stored in the output folder and performs processing operations on them. It has a `main()` function, within which are defined several nested functions defined (not mandatory, but true for Fe-simplified and I like this workflow). `main()` takes as input a list of these nested functions to run when called. automate.py calls the main function without argument, so you can easily change what functions run by default just modifying the func_list default in the post_process.py file itself. In my processing files I have defined plotting functions, but in theory you could define anything you wish and save your output to anywhere you want. I however defined the `utils.plot_log_axes()` function which plots two variables on log-log axes - very convenient for plotting flux spectra. It also has parameters you can specify to save your graphs relative to your current `config.RUN_ENV` environment in the graphs folder, again in a subdirectory that tells you the number of particles in the simulation.

#### data_load.py

This is the least interconnected file with respect to the rest of the automation. Its functions are called only within the model.py and post_process.py scripts within the model, and the functionality it can be used for is flexible. This has so far been used for loading data output files to compare with your openmc output in the post_processing.py file, or for loading external information needed to specify the openmc model (e.g. tally bin structures, source information). The only way you may want to connect the information in this file to the automating structure for the perturbations is to specify the file location of the data you are extracting relative to `config.MAIN_DIR`, so that if you decide to share your work on a remote repository and clone it for another PC or your laptop, or for someone else to take a look at you won't have to undo a lot of mess with your file naming (could alternatively specify relative file locations).

#### perturbed_run_data

This is where all the data from your runs with perturbed cross sections goes. Runs of differing perturbation are stored in folders according to a code. The code is simple. If you look in the Fe_simplified/perturbed_run/data folder the code will read "mt#p#" where the number following mt is the mt number of the reaction, and the number following p is the fractional perturbation of the total cross section of that reaction. Admittedly, perhaps the code should include nuclides perturbed as well, but for the applications I was working in the nuclides being perturbed never changed thus this change was not made. If necessary you could just check the materials.xml file to check which nuclides had been perturbed.


### File Descriptions

Below is a list of the files in the repository and a brief description of how they relate to each other and what they do:

- run.py : this is your main run file, and is intended to act as a sort of 'black box' when it it comes to running all of your simulations - particularly when it comes to the perturbations. 
- automate.py: A module that draws functions from most other modules all together to automate the file and folder construction, cross section perturbations, OpenMC simulation runs, and post-processing of output / plotting of graphs.
- finite_difference.py: runs the senstivity calculations following any perturbation runs and tabulates the results 
- modify_materials.py: copies the materials.xml file from the standard_run folder then modifies it to insert the perturbed materials. Then it places this materials.xml file as well as the geometry.xml file from the standard_run folder in a new folder. This folder is named according to a code that allows easy identification of the perturbation performed.
- perturb_xs.py written originally by Ander Gray, and subsequently modified, perturbs the cross sections that are stored in a way that is accessible to openmc.
- utils.py: a few utilities to do with plotting and loading data of general applicability to all models (add more if you want)

#### Within a model

- model.py: Loads the openmc xml files needed for the simulations. Also processes data from statepoint file output from the simulations for storage.    
- post_process.py: post-processes output data that was stored from `model.processing()` function. Used primarily for plotting in the current models defined.
- data_load.py: loads any data you want for comparison with the model you are simulated. Might be PARTISN, MCNP or other output you want to compare. Also might load source (a particle source) information. The functionality of this file is flexible.


### Running a simulation

Once a simulation has been set up, whether by you or someone else, and you have been careful to obey the needed conventions and define the right functions, running automated perturbation runs is relatively simple. Head to the run.py file. the `run()` function contains two key module functions, the first, `automate.load_model()` loads configuration variables for the particular model you want to run an in a specific run environment. Second comes the `automate.main_run()` function. These functions are important enough that a fuller description of their functionality and parameters is given below.

`automate.load_model(model, run_env=None)`

Loads the configuration variables needed for modules to know where to direct simulation outputs, graphs, and to locate comparison data. Also changes `sys.path` variable to be able to reload necessary model.py, data_load.py and post_process.py modules when switching between models. The `run_env` parameter is a crucial understanding point to expand your flexibility in modifying your simulations for experimentation of model conditions. 

The run_env parameter sets the `config.RUN_ENV` variable. This tells the code where to find the xml files to actually run the simulation, or where to generate them if necessary. Let us say, for example, you would like to see what happens when you make an unperturbed run with a different source to the one you specified, just to check the difference. You could add an if clause to your `model.settings()` function to load a certain source instead of the normal one after it checks that the `config.RUN_ENV` variable is pointing to the appropriate folder, then run the run.py file as usual. The code will automatically create the new run environment folder for you, copy the materials.xml, geometry.xml, settings.xml and tallies.xml files from your standard_run in there, then run your `model.settings()` and `model.tallies()` functions to overwrite the ones that have just been copied to add your new source information. 

There are some caveats to this functionality. Since the automate.py file only calls the `model.settings()` and `model.tallies()` commands while making sequential runs, if you want to change the materials or geometry, you will need to run the materials and geometry generating functions within the model.py file in python separately inside the model you want to simulate and point its output to the `config.RUN_ENV` environment. As an example, just change what comes after `config.RUN_ENV =` at the bottom of the model.py file within Fe-simplified. 

Note, however, that the copying of the xml files from the standard_run folder only occurs when the program cannot detect the filename of the run_env you have specified. So you can change the materials.xml file once and it will stay the same however many times you run run.py (I found this a useful feature while experimenting with my simulations).

`automate.main_run(powers=[6], nuclides=None, mts=None,   
    perturbations=None, discretization=None, check_repeat=True)`

This function performs the actual simulations, output generation, cross section perturbation and etc. When specifying a perturbation run, all of nuclides, mts and perturbations must not be None. Either an error will be raised, or if perturbations is None then a run using the xml files in the standard_run file will be performed (unless another run_env is specified in the `automate.load_model()` function).

Parameters: 

**powers**: a list of the the powers of 10 of the number of particles you would like to simulate in your simulations.

**nuclides**: list of nuclides to perturb. Nuclides are strings in same format as openmc E.g. ["Fe56", "Fe57", "Fe54].

**mts**: list of mt number of the reactions you want to perturb. ***Note*** the program will not perturb all of these mts and do a single run with them for a given perturbation, but will do successive runs with one mt number perturbed at a time.

**perturbations**: A list of perturbations indicating how much you would like the total cross section of a particular reaction to all the nuclides to be perturbed. E.g. a perturbation of 0.01 will multiply the total cross section of a given nuclide and reaction by 1.01.

**discretization**: an unimplemented half-completed potential feature. See proposed extensions section.

**check_repeat**: checks if a run has already been performed (i.e. does a folder with the corresponding run_env exist), and if so skips all the simulations and perturbations and just loads the graph output.

**NOTE**, running perturbations at the same time as specifying a new run_env is not supported. A normal perturbation run while taking xml files from the standard_run file will be performed. If you really want perturbations with that model configuration, just create a new model folder (at the level of the H1 and Fe-simplified folders), copy and paste from your base model and modify what is necessary.

#### Some illustrative examples:
<br>
```
automate.main_run(powers=[6, 7, 8], nuclides=['Fe56', 'Fe57']],
    mts=[2], perturbations=[0.1], check_repeat=False)
```
Will run 3 simulations with three different powers. Each simulation both nuclides Fe56 and Fe57 will each be perturbed by 0.1 for elastic nuclear reactions which have an mt of 2 (i.e. cross sections of mt=2 will be multiplied by 1.1)

```
    automate.main_run(powers=[6, 7], 
        nuclides=['Fe56', 'Fe57','Fe58'], 
        mts=[2, 4], perturbations=[0.1, 0.3, 0.01], 
        check_repeat=False)
```
Will have 12 runs, 6 for each power. The first run will involve all nuclides having their elastic (mt=2) cross sections perturbed by 0.1 and the run will have 10^6 particles. The second will be the same except the perturbation will now be 0.3, and so on...

#### Known Problems and extensions

1. In finite_difference.py file, irritatingly runs of the same reaction and perturbation with only slight differences (e.g. 10th decimal place) in sensitivies are registered as different runs. Refine the criteria for whether or not a row should be duplicated. Note, sometimes it is useful to have data from a historical run with the same perturbation but a slightly different simulation to track the effects of any changes, so we don't necessarily want to automatically delete runs with the same perturbation
1. Finish implementing the discretization feature for the cross section data. The idea would be that somehow you would input a energy group structure into the perturbation file for perturbations to occur at each group individually along the cross section. This is what is needed to verify Ivan Kodeli's code fully. There are partially started sections of code to implement this functionality already (should in theory have been separated in a different branch, I know, this is another extension if you choose).
1. In theory no need for the command line interface for the perturb_xs.py module. Clean it up to make it fit better with the rest of the code and more understandable.
1. Adapt the code so there is no dependence at all on the python API. It only directly modifies copies of xml files in the standard_run folder when running perturbations. Python files can still live in the folder to make changes to the xml files.


### Further information

For further information or questions please email me at alexhayward4321@gmail.com, and if you would like to know some of the specifics of the motivations behind the models already implemented ask for me my report on the subject.









