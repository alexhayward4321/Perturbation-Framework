# Openmc cross section perturbation and sensitivity analysis toolkit

This repository houses a python code used to compute the sensitivity of gamma and neutron fluxes to perturbations of cross section data. It automates the perturbations of cross sections of a set of nuclides, reactions, and perturbation strengths. Further, if you specify a base model in a folder, it automatically generates all of the materials, geometry, settings and tallies files within properly labelled folders for
each perturbation, runs openmc simulations for them for a specified number of particles, then stores the results of the simulation as well as processed graphs in a standardised and easily accessible file structure and format.

## Disclaimer
---
It is still in a very early stage of development, so there are doubtless numerous bugs, but if determined a person should be able to pick up how the code works and modify it in order to improve it and make it more generally applicable. Note also that the primary developer of the code was very inexperienced in building large interconnected python projects, and this explains the conspicuous absence of defined objects and classes, as well as some experimental ways of implementing features that were based on spontaneous problem solving while potentially ignorant of standard (and perhaps much simpler) ways of doing things.

If you have any questions please feel free to email me at alexhayward4321@gmail.com.

## How to use it
---
First clone the repository using the command:
```
git clone https://github.com/alexhayward4321/Iron <folder_name>
```
inside the directory you would like to perform your data analysis (if you do not have git installed, download it). 

Then open the configuration file settings.py located in your newly downloaded repository.

