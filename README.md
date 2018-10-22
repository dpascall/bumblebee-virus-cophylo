# bumblebee-virus-cophylo
The code for the analysis perfomed in "Host evolutionary history predicts virus prevalence across bumblebee species". Some code is currently optimised for running on a cluster other for running locally. Feel free to contact me on twitter (@davidandthebees) or via email (david.pascall@glasgow.ac.uk) for help in repeating these analyses or usage of this code to analyse new datasets.

R code for data handling is a slightly modified version of code by Jarrod Hadfield for his paper "A tale of two phylogenies". The tongue length analysis uses code extracted from the source code of the package "Prevalence" outside of that package (prevmodel.txt). The Stan model code was originally built off the base code by Diogo Melo in their Stan animal model here (https://github.com/diogro/stanAnimal). Thanks to multiple members on the Stan Discourse for help in implementing the correct marginalisation to account for phylogenetic uncertainty. All errors are my own.

Due to the base off Diogo Melo's MIT licensed code, this code is MIT licensed, and any code build off it must also be MIT licensed (as I understand it, 100% not a lawyer here, please don't sue me!).
